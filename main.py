import georinex as gr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import Rbf
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
import pymap3d as pm
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

EARTH_RADIUS_KM = 6371.0
IONO_HEIGHT_KM = 350.0
GM = 3.986005e14          
OMEGA_E = 7.2921151467e-5 

def get_tec_multiplier(sys_id):
    if sys_id == 'G': return 9.52
    elif sys_id == 'E': return 7.76
    return 9.52

def mapping_function(elevation_deg):
    elevation_rad = np.radians(elevation_deg)
    sin_z_prime = (EARTH_RADIUS_KM / (EARTH_RADIUS_KM + IONO_HEIGHT_KM)) * np.cos(elevation_rad)
    return np.sqrt(1.0 - sin_z_prime**2)

def calculate_ipp(lat_sta, lon_sta, azimuth_deg, elevation_deg):
    az_rad = np.radians(azimuth_deg)
    el_rad = np.radians(elevation_deg)
    psi_rad = (np.pi / 2.0) - el_rad - np.arcsin((EARTH_RADIUS_KM / (EARTH_RADIUS_KM + IONO_HEIGHT_KM)) * np.cos(el_rad))
    lat_sta_rad = np.radians(lat_sta)
    lon_sta_rad = np.radians(lon_sta)
    lat_ipp_rad = np.arcsin(np.sin(lat_sta_rad) * np.cos(psi_rad) + np.cos(lat_sta_rad) * np.sin(psi_rad) * np.cos(az_rad))
    lon_ipp_rad = lon_sta_rad + np.arcsin((np.sin(psi_rad) * np.sin(az_rad)) / np.cos(lat_ipp_rad))
    return np.degrees(lat_ipp_rad), np.degrees(lon_ipp_rad)

def get_satellite_position_kepler(nav_sv, t_obs):
    try:
        nav_epoch = nav_sv.sel(time=t_obs, method='nearest')
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
        t_nav = safe_get('time')
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
    except Exception as e:
        # --- FIX D: Bessere Fehlerbehandlung ---
        # Falls unerwartet viele Bahndaten fehlen, kann man dies hier ausgeben lassen
        # print(f"Orbit Error für {nav_sv.sv.values} um {t_obs}: {e}")
        return None, None, None

def process_rinex_advanced(obs_file, nav_files):
    obs = gr.load(obs_file)
    if isinstance(nav_files, str): nav_files = [nav_files]
    nav_datasets = [gr.load(f) for f in nav_files]
    
    if 'position' not in obs.attrs:
        lat_sta, lon_sta, alt_sta = 49.87, 8.62, 100.0
    else:
        x, y, z = obs.attrs['position']
        lat_sta, lon_sta, alt_sta = pm.ecef2geodetic(x, y, z)

    l1_bands =['C1C', 'C1P', 'C1W', 'C1X', 'C1S', 'L1C']
    l2_bands =['C2W', 'C2C', 'C2L', 'C2I', 'C5Q', 'C7Q', 'C8Q', 'L2C', 'L5Q']

    all_lats, all_lons, all_vtecs =[], [], []
    valid_svs =[]
    nav_map = {} 
    
    for sv in obs.sv.values:
        sys_id = str(sv)[0]
        if sys_id in ['G', 'E']: 
            for nav_ds in nav_datasets:
                if sv in nav_ds.sv.values:
                    valid_svs.append(sv)
                    nav_map[sv] = nav_ds
                    break

    for sv in tqdm(valid_svs, desc="🛰️ Verarbeite GPS & Galileo"):
        sv_id = str(sv)
        sv_data = obs.sel(sv=sv)
        nav_sv = nav_map[sv].sel(sv=sv) 
        
        # --- FIX A: Dynamische & robuste Band-Selektion ---
        best_p1, best_p2 = None, None
        max_valid = 0
        
        for b1 in l1_bands:
            for b2 in l2_bands:
                if b1 in sv_data.data_vars and b2 in sv_data.data_vars:
                    # Zählt, wie viele Werte für dieses Bandpaar NICHT NaN sind
                    valid_mask = ~np.isnan(sv_data[b1] - sv_data[b2])
                    valid_count = int(valid_mask.sum().item())
                    
                    if valid_count > max_valid:
                        max_valid = valid_count
                        best_p1, best_p2 = b1, b2

        if max_valid < 5: 
            continue # Überspringen, wenn der Satellit weniger als 5 nutzbare Epochen hat
            
        p1_var, p2_var = best_p1, best_p2
        # --------------------------------------------------
            
        stec_raw = get_tec_multiplier(sv_id[0]) * (sv_data[p2_var] - sv_data[p1_var])
        valid_mask = ~np.isnan(stec_raw.values)
        stec_vals = stec_raw.values[valid_mask]
        time_vals = stec_raw.time.values[valid_mask]

        track_lats, track_lons, track_vtecs = [], [],[]
        for t_obs, stec in zip(time_vals, stec_vals):
            X, Y, Z = get_satellite_position_kepler(nav_sv, t_obs)
            if X is None: continue
            az, el, _ = pm.ecef2aer(X, Y, Z, lat_sta, lon_sta, alt_sta)
            
            if el > 15.0: # Elevationsmaske
                lat_ipp, lon_ipp = calculate_ipp(lat_sta, lon_sta, az, el)
                vtec = stec * mapping_function(el)
                track_lats.append(lat_ipp)
                track_lons.append(lon_ipp)
                track_vtecs.append(vtec)
                
        if len(track_vtecs) > 0:
            track_vtecs = np.array(track_vtecs)
            
            # --- FIX C: Per-Satellite Leveling (Pseudo-DCB Korrektur) ---
            # Statt pauschal am Ende alles zu korrigieren, leveln wir JEDEN 
            # Satelliten einzeln auf einen physikalisch realistischen Minimalwert.
            # Das verhindert die extremen Sprünge (z.B. 10 TECU vs 60 TECU).
            track_vtecs = track_vtecs - np.min(track_vtecs) + 5.0 
            # -------------------------------------------------------------
            
            all_lats.extend(track_lats)
            all_lons.extend(track_lons)
            all_vtecs.extend(track_vtecs)

    # --- Zeitfenster für den Plot-Titel extrahieren ---
    time_min = pd.to_datetime(obs.time.values.min())
    time_max = pd.to_datetime(obs.time.values.max())
    time_str = f"{time_min.strftime('%H:%M')} - {time_max.strftime('%H:%M')} UTC"

    return np.array(all_lats), np.array(all_lons), np.array(all_vtecs), lat_sta, lon_sta, time_str

def plot_professional_ionosphere_map(lats, lons, vtec, lat_sta, lon_sta, time_str):
    print("\n🌍 Erzeuge VTEC-Karte (RBF Smoothed & Masked)...")
    if len(vtec) == 0:
        print("❌ ABBRUCH: Keine Daten.")
        return

    grid_lon, grid_lat = np.mgrid[lon_sta-10:lon_sta+10:300j, lat_sta-8:lat_sta+8:300j]
    
    # RBF Interpolation
    rbf = Rbf(lons, lats, vtec, function='linear', epsilon=2.0)
    grid_vtec = rbf(grid_lon, grid_lat)

    # Masking mit KDTree (Begrenzung des Einflussbereichs)
    tree = cKDTree(np.c_[lons, lats])
    dist, _ = tree.query(np.c_[grid_lon.ravel(), grid_lat.ravel()])
    dist = dist.reshape(grid_lon.shape)
    
    MAX_DISTANCE_DEG = 2.5 
    grid_vtec[dist > MAX_DISTANCE_DEG] = np.nan

    grid_vtec_smooth = grid_vtec.copy()
    mask = np.isnan(grid_vtec)
    grid_vtec_smooth[mask] = np.nanmean(grid_vtec) 
    grid_vtec_smooth = gaussian_filter(grid_vtec_smooth, sigma=1.0)
    grid_vtec_smooth[mask] = np.nan 

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    ax.set_extent([lon_sta - 8, lon_sta + 8, lat_sta - 6, lat_sta + 6], crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor='#333333')
    ax.add_feature(cfeature.BORDERS, linestyle='--', alpha=0.6)
    ax.add_feature(cfeature.LAND, facecolor='#f5f5f5')
    ax.add_feature(cfeature.OCEAN, facecolor='#e0f0ff')
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    
    # Levels dynamisch an echte Daten anpassen (Verhindert zu starke Sättigung)
    vmax = max(15, np.nanmax(grid_vtec_smooth))
    levels = np.linspace(0, vmax, 40)
    
    contour = ax.contourf(grid_lon, grid_lat, grid_vtec_smooth, levels=levels, cmap='plasma', alpha=0.85, transform=ccrs.PlateCarree())
    
    ax.scatter(lons, lats, color='white', s=15, transform=ccrs.PlateCarree(), alpha=0.8)
    ax.scatter(lons, lats, color='black', s=5, transform=ccrs.PlateCarree(), alpha=1.0, label='IPP Messpunkte')
    ax.plot(lon_sta, lat_sta, marker='^', color='red', markersize=14, markeredgecolor='white', markeredgewidth=1.5, transform=ccrs.PlateCarree(), label="GNSS Station")
    
    cbar = plt.colorbar(contour, ax=ax, shrink=0.75, pad=0.04)
    cbar.set_label('Relativer VTEC [TECU]', fontsize=13, weight='bold') # Titel leicht angepasst
    cbar.ax.tick_params(labelsize=11)
    
    # Titel um Uhrzeit ergänzt!
    plt.title(f'Lokale VTEC Ionosphären-Map ({time_str})\nStation: {lon_sta:.2f}E, {lat_sta:.2f}N | Höhe: 350km', pad=15, fontsize=15, weight='bold')
    plt.legend(loc='lower right', framealpha=0.95, fontsize=11)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # ====================================================================
    RINEX_OBS_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26O"
    RINEX_NAV_FILES =[r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26N", r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 5\V3RJ1060.26L"]  
    # ====================================================================
    
    lats, lons, vtecs, sta_lat, sta_lon, time_str = process_rinex_advanced(RINEX_OBS_FILE, RINEX_NAV_FILES)
    plot_professional_ionosphere_map(lats, lons, vtecs, sta_lat, sta_lon, time_str)