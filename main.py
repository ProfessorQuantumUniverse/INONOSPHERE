import georinex as gr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.spatial import cKDTree
import pymap3d as pm
from tqdm import tqdm
import warnings
import os

warnings.filterwarnings('ignore')

EARTH_RADIUS_KM = 6371.0
IONO_HEIGHT_KM = 350.0
GM = 3.986005e14          
OMEGA_E = 7.2921151467e-5 

def get_tec_multiplier(sys_id):
    if sys_id == 'G': return 9.52
    elif sys_id == 'E': return 7.76
    elif sys_id == 'C': return 8.99 
    return 9.52

def mapping_function(elevation_deg):
    elevation_rad = np.radians(elevation_deg)
    sin_z_prime = (EARTH_RADIUS_KM / (EARTH_RADIUS_KM + IONO_HEIGHT_KM)) * np.cos(elevation_rad)
    return np.sqrt(1.0 - sin_z_prime**2)

def calculate_ipp(lat_sta, lon_sta, azimuth_deg, elevation_deg):
    az_rad = np.radians(azimuth_deg)
    el_rad = np.radians(elevation_deg)
    
    val1 = (EARTH_RADIUS_KM / (EARTH_RADIUS_KM + IONO_HEIGHT_KM)) * np.cos(el_rad)
    psi_rad = (np.pi / 2.0) - el_rad - np.arcsin(np.clip(val1, -1.0, 1.0))
    
    lat_sta_rad = np.radians(lat_sta)
    lon_sta_rad = np.radians(lon_sta)
    
    val2 = np.sin(lat_sta_rad) * np.cos(psi_rad) + np.cos(lat_sta_rad) * np.sin(psi_rad) * np.cos(az_rad)
    lat_ipp_rad = np.arcsin(np.clip(val2, -1.0, 1.0))
    
    lon_ipp_rad = lon_sta_rad + np.arcsin(np.clip((np.sin(psi_rad) * np.sin(az_rad)) / np.cos(lat_ipp_rad), -1.0, 1.0))
    
    return np.degrees(lat_ipp_rad), np.degrees(lon_ipp_rad)

def get_satellite_position_kepler(nav_sv, t_obs):
    try:
        nav_epoch = nav_sv.sel(time=t_obs, method='nearest')
        t_nav = nav_epoch['time'].values
        if abs((t_obs - t_nav) / np.timedelta64(1, 'h')) > 4.0:
            return None, None, None

        def safe_get(*keys):
            for k in keys:
                if k in nav_epoch.data_vars or k in nav_epoch.coords:
                    val = nav_epoch[k].values
                    if isinstance(val, np.ndarray): return val.flatten()[0] 
                    return val
            raise KeyError(f"Missing {keys}")

        A = safe_get('sqrtA')**2
        n0 = np.sqrt(GM / A**3)
        n = n0 + safe_get('DeltaN')
        tk = (t_obs - t_nav) / np.timedelta64(1, 's')
        
        M = safe_get('M0') + n * tk
        E = M
        e = safe_get('Eccentricity', 'e')
        for _ in range(7): E = M + e * np.sin(E)
            
        v = np.arctan2(np.sqrt(1 - e**2) * np.sin(E), np.cos(E) - e)
        phi = v + safe_get('omega')
        
        du = safe_get('Cus') * np.sin(2*phi) + safe_get('Cuc') * np.cos(2*phi)
        dr = safe_get('Crs') * np.sin(2*phi) + safe_get('Crc') * np.cos(2*phi)
        di = safe_get('Cis') * np.sin(2*phi) + safe_get('Cic') * np.cos(2*phi)
        
        u = phi + du
        r = A * (1 - e * np.cos(E)) + dr
        i = safe_get('Io', 'i0') + safe_get('IDOT', 'idot') * tk + di
        
        x_prime = r * np.cos(u)
        y_prime = r * np.sin(u)
        
        Omega = safe_get('Omega0') + (safe_get('OmegaDot', 'omegadot') - OMEGA_E) * tk - OMEGA_E * safe_get('Toe', 'toe')
        X = x_prime * np.cos(Omega) - y_prime * np.cos(i) * np.sin(Omega)
        Y = x_prime * np.sin(Omega) + y_prime * np.cos(i) * np.cos(Omega)
        Z = y_prime * np.sin(i)
        
        return X, Y, Z
    except Exception:
        return None, None, None

def get_absolute_base_vtec(ionex_file, lat_sta, lon_sta, target_time):
    """
    Eigenbau IONEX-Parser: Sucht die passende Epoche zur Messzeit und 
    extrahiert den exakten VTEC-Wert für den Standort der GNSS-Station.
    """
    if not ionex_file or not os.path.exists(ionex_file):
        print(f"⚠️ IONEX Datei nicht gefunden: {ionex_file}. Nutze Fallback (20.0).")
        return 20.0

    try:
        exponent = -1 
        best_epoch = None
        min_time_diff = float('inf')
        
        # 1. Finde die passendste Uhrzeit (Epoche) in der Datei
        with open(ionex_file, 'r') as f:
            for line in f:
                if "EXPONENT" in line:
                    exponent = int(float(line.split()[0]))
                elif "EPOCH OF CURRENT MAP" in line:
                    parts = line.split()
                    y, m, d, h, mn, s = map(int, parts[:6])
                    if y < 100: y += 2000
                    epoch_time = pd.Timestamp(year=y, month=m, day=d, hour=h, minute=mn, second=s)
                    
                    diff = abs((epoch_time - target_time).total_seconds())
                    if diff < min_time_diff:
                        min_time_diff = diff
                        best_epoch = epoch_time

        if best_epoch is None:
            raise ValueError("Keine Epochen in IONEX Datei gefunden.")

        # 2. Lade nur die TEC-Werte für exakt diese Epoche
        parsing_target_map = False
        current_lat = None
        lons, tec_values = [],[]
        tec_map = {}
        
        with open(ionex_file, 'r') as f:
            for line in f:
                if "EPOCH OF CURRENT MAP" in line:
                    parts = line.split()
                    y, m, d, h, mn, s = map(int, parts[:6])
                    if y < 100: y += 2000
                    epoch_time = pd.Timestamp(year=y, month=m, day=d, hour=h, minute=mn, second=s)
                    parsing_target_map = (epoch_time == best_epoch)
                    
                elif parsing_target_map and "LAT/LON1/LON2/DLON/H" in line:
                    if current_lat is not None and tec_values:
                        tec_map[current_lat] = dict(zip(lons, tec_values))
                    
                    s_line = line[:40].replace('-', ' -')
                    parts = s_line.split()
                    current_lat = float(parts[0])
                    lon1, lon2, dlon = float(parts[1]), float(parts[2]), float(parts[3])
                    lons = np.arange(lon1, lon2 + dlon/2.0, dlon)
                    tec_values =[]
                    
                elif parsing_target_map and "END OF TEC MAP" in line:
                    if current_lat is not None and tec_values:
                        tec_map[current_lat] = dict(zip(lons, tec_values))
                    break
                    
                elif parsing_target_map and current_lat is not None:
                    if not any(x in line for x in["START", "END", "LAT/LON", "EPOCH"]):
                        for c in range(0, min(80, len(line.rstrip('\n'))), 5):
                            chunk = line[c:c+5]
                            if chunk.strip():
                                tec_values.append(int(chunk))
                                
        # 3. Nächsten Gitterpunkt zur Station finden
        closest_lat = min(tec_map.keys(), key=lambda k: abs(k - lat_sta))
        closest_lon = min(tec_map[closest_lat].keys(), key=lambda k: abs(k - lon_sta))
        
        raw_tec = tec_map[closest_lat][closest_lon]
        actual_tec = raw_tec * (10 ** exponent)
        
        print("\n" + "🌟"*22)
        print(" 🌍 IONEX KALIBRIERUNG ERFOLGREICH!")
        print("🌟"*22)
        print(f"Gewählte Karte: {best_epoch.strftime('%H:%M')} UTC")
        print(f"Gitter-Punkt:   {closest_lat}°N, {closest_lon}°E")
        print(f"Anker-TEC Wert: {actual_tec:.1f} TECU")
        print("🌟"*22 + "\n")
        
        return actual_tec
        
    except Exception as e:
        print(f"⚠️ IONEX Parser fehlgeschlagen: {e} -> Nutze Fallback (20.0)")
        return 20.0

def process_rinex_advanced(obs_file, nav_files, ionex_file):
    obs = gr.load(obs_file)
    if isinstance(nav_files, str): nav_files = [nav_files]
    nav_datasets = [gr.load(f) for f in nav_files]
    
    if 'position' not in obs.attrs:
        lat_sta, lon_sta, alt_sta = 49.87, 8.62, 100.0
    else:
        x, y, z = obs.attrs['position']
        lat_sta, lon_sta, alt_sta = pm.ecef2geodetic(x, y, z)

    # --- Die magische IONEX Kalibrierung ---
    time_min = pd.to_datetime(obs.time.values.min())
    time_max = pd.to_datetime(obs.time.values.max())
    target_time = time_min + (time_max - time_min) / 2
    
    BASE_VTEC_LEVEL = get_absolute_base_vtec(ionex_file, lat_sta, lon_sta, target_time)
    # ---------------------------------------

    l1_bands =['C1C', 'C1P', 'C1W', 'C1X', 'C1S', 'L1C', 'C2I']
    l2_bands =['C2W', 'C2C', 'C2L', 'C5Q', 'C7Q', 'C8Q', 'L2C', 'L5Q', 'C7I']

    all_lats, all_lons, all_vtecs = [], [],[]
    valid_svs =[]
    nav_map = {} 
    
    stats = {'total': 0, 'no_nav': 0, 'no_dual_freq': 0, 'low_elevation': 0, 'success': 0}
    
    for sv in obs.sv.values:
        stats['total'] += 1
        sys_id = str(sv)[0]
        if sys_id in ['G', 'E', 'C']: 
            found_nav = False
            for nav_ds in nav_datasets:
                if sv in nav_ds.sv.values:
                    valid_svs.append(sv)
                    nav_map[sv] = nav_ds
                    found_nav = True
                    break
            if not found_nav: stats['no_nav'] += 1

    for sv in tqdm(valid_svs, desc="🛰️ Analysiere & Filtere Satelliten"):
        sv_id = str(sv)
        sv_data = obs.sel(sv=sv)
        nav_sv = nav_map[sv].sel(sv=sv) 
        
        best_p1, best_p2 = None, None
        max_valid = 0
        
        for b1 in l1_bands:
            for b2 in l2_bands:
                if b1 in sv_data.data_vars and b2 in sv_data.data_vars:
                    valid_count = int((~np.isnan(sv_data[b1] - sv_data[b2])).sum().item())
                    if valid_count > max_valid:
                        max_valid = valid_count
                        best_p1, best_p2 = b1, b2

        if max_valid < 5: 
            stats['no_dual_freq'] += 1
            continue 
            
        stec_raw = get_tec_multiplier(sv_id[0]) * (sv_data[best_p2] - sv_data[best_p1])
        valid_mask = np.isfinite(stec_raw.values)
        
        stec_vals = stec_raw.values[valid_mask]
        time_vals = stec_raw.time.values[valid_mask]

        track_lats_raw, track_lons_raw, track_vtecs_raw = [], [],[]
        had_low_elevation = False
        
        for t_obs, stec in zip(time_vals, stec_vals):
            X, Y, Z = get_satellite_position_kepler(nav_sv, t_obs)
            if X is None: continue
            az, el, _ = pm.ecef2aer(X, Y, Z, lat_sta, lon_sta, alt_sta)
            
            if el < 25.0:
                had_low_elevation = True
            else:
                lat_ipp, lon_ipp = calculate_ipp(lat_sta, lon_sta, az, el)
                vtec = stec * mapping_function(el)
                
                if np.isfinite(vtec) and np.isfinite(lat_ipp) and np.isfinite(lon_ipp):
                    track_lats_raw.append(lat_ipp)
                    track_lons_raw.append(lon_ipp)
                    track_vtecs_raw.append(vtec)
                
        if len(track_vtecs_raw) > 5:
            x_idx = np.arange(len(track_vtecs_raw))
            poly_coeffs = np.polyfit(x_idx, track_vtecs_raw, 1) 
            poly_func = np.poly1d(poly_coeffs)
            vtecs_smoothed = poly_func(x_idx)
            
            vtecs_leveled = vtecs_smoothed - np.mean(vtecs_smoothed) + BASE_VTEC_LEVEL
            
            all_lats.extend(track_lats_raw)
            all_lons.extend(track_lons_raw)
            all_vtecs.extend(vtecs_leveled)
            stats['success'] += 1
        elif had_low_elevation:
            stats['low_elevation'] += 1

    all_vtecs = np.array(all_vtecs)
    if len(all_vtecs) > 0:
        all_vtecs = np.maximum(all_vtecs, 0.1)

    time_str = f"{time_min.strftime('%H:%M')} - {time_max.strftime('%H:%M')} UTC"

    print("\n" + "="*35)
    print(" 📊 SATELLITEN-DIAGNOSE TABELLE")
    print("="*35)
    print(f"Total im OBS-File gesichtet:     {stats['total']}")
    print(f"Erfolgreich geplottet:           {stats['success']}")
    print("-" * 35)
    print("Aussortiert weil:")
    print(f"- Keine Navigationsdaten (.26N/L): {stats['no_nav']}")
    print(f"- Fehlendes L2 Signal (Multipath): {stats['no_dual_freq']}")
    print(f"- Zu nah am Horizont (<25° Elev.): {stats['low_elevation']}")
    print("="*35 + "\n")

    return np.array(all_lats), np.array(all_lons), all_vtecs, lat_sta, lon_sta, time_str


def plot_professional_ionosphere_map(lats, lons, vtec, lat_sta, lon_sta, time_str):
    print("🌍 Erzeuge VTEC-Karte (2D Spatial Plane Model)...")
    
    valid = np.isfinite(lons) & np.isfinite(lats) & np.isfinite(vtec)
    lons = lons[valid]
    lats = lats[valid]
    vtec = vtec[valid]
    
    if len(vtec) < 3:
        print("❌ ABBRUCH: Zu wenige gültige Messpunkte für ein 2D-Flächenmodell.")
        return

    lons_centered = lons - lon_sta
    lats_centered = lats - lat_sta
    
    A = np.c_[lons_centered, lats_centered, np.ones(len(lons_centered))]
    C, _, _, _ = np.linalg.lstsq(A, vtec, rcond=None)

    grid_lon, grid_lat = np.mgrid[lon_sta-6:lon_sta+6:400j, lat_sta-5:lat_sta+5:400j]
    grid_vtec = C[0] * (grid_lon - lon_sta) + C[1] * (grid_lat - lat_sta) + C[2]

    tree = cKDTree(np.c_[lons, lats])
    dist, _ = tree.query(np.c_[grid_lon.ravel(), grid_lat.ravel()])
    dist = dist.reshape(grid_lon.shape)
    
    # ERHÖHT auf 2.0 für eine schönere, zusammenhängende Wolke
    MAX_DISTANCE_DEG = 2.0 
    grid_vtec[dist > MAX_DISTANCE_DEG] = np.nan

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    ax.set_extent([lon_sta - 5, lon_sta + 5, lat_sta - 4, lat_sta + 4], crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='#333333')
    ax.add_feature(cfeature.BORDERS, linestyle='--', alpha=0.5)
    ax.add_feature(cfeature.LAND, facecolor='#f9f9f9')
    ax.add_feature(cfeature.OCEAN, facecolor='#e6f2ff')
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.4, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    
    vmin = np.nanmin(grid_vtec)
    vmax = np.nanmax(grid_vtec)
    if np.isnan(vmin) or np.isnan(vmax) or vmin == vmax: 
        vmax = vmin + 1.0 
        
    levels = np.linspace(vmin, vmax, 50) 
    
    contour = ax.contourf(grid_lon, grid_lat, grid_vtec, levels=levels, cmap='plasma', alpha=0.85, transform=ccrs.PlateCarree())
    ax.contour(grid_lon, grid_lat, grid_vtec, levels=10, colors='white', alpha=0.3, linewidths=0.5, transform=ccrs.PlateCarree())
    
    ax.scatter(lons, lats, color='white', s=15, transform=ccrs.PlateCarree(), alpha=0.9)
    ax.scatter(lons, lats, color='black', s=3, transform=ccrs.PlateCarree(), alpha=1.0, label='IPP Messkorridor')
    ax.plot(lon_sta, lat_sta, marker='^', color='red', markersize=14, markeredgecolor='white', markeredgewidth=1.5, transform=ccrs.PlateCarree(), label="GNSS Station")
    
    cbar = plt.colorbar(contour, ax=ax, shrink=0.75, pad=0.04)
    # NEU: Formatierung der Farbskala (nur 1 Nachkommastelle)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    cbar.set_label('Kalibrierter VTEC [TECU]', fontsize=13, weight='bold')
    
    plt.title(f'Lokales Ionosphärenmodell (LIM) ({time_str})\nAbsolut kalibriert (IONEX) | Filter: 25° Elev.', pad=15, fontsize=14, weight='bold')
    plt.legend(loc='lower right', framealpha=0.95, fontsize=11)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # HIER DEINE 3 DATEIPFADE EINTRAGEN:
    RINEX_OBS_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26O"
    RINEX_NAV_FILES =[
        r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26N", 
        r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26L"
    ]
    IONEX_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\c1pg1060.26i"
    
    # Der IONEX File Pfad wird nun korrekt in die Hauptfunktion übergeben!
    lats, lons, vtecs, sta_lat, sta_lon, time_str = process_rinex_advanced(RINEX_OBS_FILE, RINEX_NAV_FILES, IONEX_FILE)
    plot_professional_ionosphere_map(lats, lons, vtecs, sta_lat, sta_lon, time_str)