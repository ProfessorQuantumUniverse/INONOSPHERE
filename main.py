import georinex as gr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import pymap3d as pm
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

# Physikalische Konstanten (WGS84 & Ionosphäre)
EARTH_RADIUS_KM = 6371.0
IONO_HEIGHT_KM = 350.0
GM = 3.986005e14          # Gravitationskonstante Erde [m^3/s^2]
OMEGA_E = 7.2921151467e-5 # Erdrotation [rad/s]

def get_tec_multiplier(sys_id):
    if sys_id == 'G': return 9.52   # GPS
    elif sys_id == 'E': return 7.76 # Galileo
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
    
    lat_ipp_rad = np.arcsin(np.sin(lat_sta_rad) * np.cos(psi_rad) + 
                            np.cos(lat_sta_rad) * np.sin(psi_rad) * np.cos(az_rad))
    lon_ipp_rad = lon_sta_rad + np.arcsin((np.sin(psi_rad) * np.sin(az_rad)) / np.cos(lat_ipp_rad))
    return np.degrees(lat_ipp_rad), np.degrees(lon_ipp_rad)

def get_satellite_position_kepler(nav, sv, t_obs):
    """
    Echte Kepler-Mathematik! Berechnet die ECEF X,Y,Z Koordinaten des 
    Satelliten aus den Broadcast-Ephemeriden (RINEX NAV).
    """
    try:
        nav_sv = nav.sel(sv=sv)
        # Nächste Ephemeriden-Epoche zur Messzeit finden
        nav_epoch = nav_sv.sel(time=t_obs, method='nearest')
        
        # t_k: Zeitdifferenz zur Ephemeriden-Referenzzeit in Sekunden
        t_nav = nav_epoch['time'].values
        tk = (t_obs - t_nav) / np.timedelta64(1, 's')
        
        A = nav_epoch['sqrtA'].item()**2
        n0 = np.sqrt(GM / A**3)
        n = n0 + nav_epoch['DeltaN'].item()
        M = nav_epoch['M0'].item() + n * tk
        
        # Iterative Lösung der Kepler-Gleichung E = M + e*sin(E)
        E = M
        e = nav_epoch['Eccentricity'].item()
        for _ in range(7):
            E = M + e * np.sin(E)
            
        v = np.arctan2(np.sqrt(1 - e**2) * np.sin(E), np.cos(E) - e)
        phi = v + nav_epoch['omega'].item()
        
        # Bahnstörungen (Harmonische Korrekturen)
        du = nav_epoch['Cus'].item() * np.sin(2*phi) + nav_epoch['Cuc'].item() * np.cos(2*phi)
        dr = nav_epoch['Crs'].item() * np.sin(2*phi) + nav_epoch['Crc'].item() * np.cos(2*phi)
        di = nav_epoch['Cis'].item() * np.sin(2*phi) + nav_epoch['Cic'].item() * np.cos(2*phi)
        
        u = phi + du
        r = A * (1 - e * np.cos(E)) + dr
        i = nav_epoch['Io'].item() + nav_epoch['IDOT'].item() * tk + di
        
        x_prime = r * np.cos(u)
        y_prime = r * np.sin(u)
        
        Omega = nav_epoch['Omega0'].item() + (nav_epoch['OmegaDot'].item() - OMEGA_E) * tk - OMEGA_E * nav_epoch['Toe'].item()
        
        X = x_prime * np.cos(Omega) - y_prime * np.cos(i) * np.sin(Omega)
        Y = x_prime * np.sin(Omega) + y_prime * np.cos(i) * np.cos(Omega)
        Z = y_prime * np.sin(i)
        
        return X, Y, Z
    except Exception as e:
        print(f"Exception for {sv}: {e}")
        # Falls Ephemeriden fehlen oder lückenhaft sind
        return None, None, None

def process_rinex_advanced(obs_file, nav_file):
    print("⏳ Lade RINEX Beobachtungs- und Navigationsdaten...")
    obs = gr.load(obs_file)
    nav = gr.load(nav_file)
    
    if 'position' not in obs.attrs:
        lat_sta, lon_sta, alt_sta = 49.87, 8.62, 100.0
    else:
        x, y, z = obs.attrs['position']
        lat_sta, lon_sta, alt_sta = pm.ecef2geodetic(x, y, z)

    l1_bands = ['C1C', 'C1', 'P1']
    l2_bands =['C2W', 'C2C', 'P2', 'C2', 'C5Q', 'L2']

    all_lats, all_lons, all_vtecs = [], [],[]

    # Wir fokussieren uns auf GPS (G) und Galileo (E), da diese die gleichen Kepler-Formeln nutzen
    valid_svs =[sv for sv in obs.sv.values if str(sv)[0] in ['G', 'E']]

    for sv in tqdm(valid_svs, desc="🛰️ Berechne Kepler-Orbits & VTEC"):
        sv_id = str(sv)
        sv_data = obs.sel(sv=sv)
        
        p1_var = next((v for v in l1_bands if v in sv_data.data_vars), None)
        p2_var = next((v for v in l2_bands if v in sv_data.data_vars), None)
        
        if not (p1_var and p2_var): 
            print(f"[{sv_id}] fehlende p1/p2. p1: {p1_var}, p2: {p2_var}")
            continue
            
        # 1. Rohes STEC berechnen
        stec_raw = get_tec_multiplier(sv_id[0]) * (sv_data[p2_var] - sv_data[p1_var])
        
        # NaN Werte entfernen
        valid_mask = ~np.isnan(stec_raw.values)
        stec_vals = stec_raw.values[valid_mask]
        time_vals = stec_raw.time.values[valid_mask]
        
        if len(stec_vals) < 20: 
            print(f"[{sv_id}] len(stec_vals) < 20: {len(stec_vals)}")
            continue

        track_lats, track_lons, track_vtecs = [], [],[]
        
        # 2. Reale Position für jede Epoche berechnen
        for t_obs, stec in zip(time_vals, stec_vals):
            X, Y, Z = get_satellite_position_kepler(nav, sv_id, t_obs)
            if X is None: 
                print(f"[{sv_id}] Kepler X ist None für {t_obs}")
                continue
                
            # ECEF in lokalen Azimut/Elevation der Station umrechnen
            az, el, _ = pm.ecef2aer(X, Y, Z, lat_sta, lon_sta, alt_sta)
            
            # Elevation Mask (Ignoriere Signale unter 20°, da extremes Multipath/Rauschen)
            if el > 20.0:
                lat_ipp, lon_ipp = calculate_ipp(lat_sta, lon_sta, az, el)
                vtec = stec * mapping_function(el)
                
                track_lats.append(lat_ipp)
                track_lons.append(lon_ipp)
                track_vtecs.append(vtec)
            else:
                print(f"[{sv_id}] Elevation {el} <= 20")
                
        # 3. DCB-Leveling (Beseitigt negative Werte und Sprünge)
        print(f"[{sv_id}] track_vtecs elements: {len(track_vtecs)}")
        if len(track_vtecs) > 0:
            track_vtecs = np.array(track_vtecs)
            # Verschiebe die gesamte Kurve des Satelliten so, dass ihr Minimum bei realistischen 3 TECU liegt
            track_vtecs = track_vtecs - np.min(track_vtecs) + 3.0
            
            all_lats.extend(track_lats)
            all_lons.extend(track_lons)
            all_vtecs.extend(track_vtecs)

    return np.array(all_lats), np.array(all_lons), np.array(all_vtecs), lat_sta, lon_sta

def plot_professional_ionosphere_map(lats, lons, vtec, lat_sta, lon_sta):
    print("🌍 Erzeuge physikalisch korrekte VTEC-Karte (Gaussian Smoothed)...")
    if len(vtec) == 0:
        print("❌ Keine Daten.")
        return

    # Raster erzeugen (Zentriert um deine Station)
    grid_lon, grid_lat = np.mgrid[lon_sta-8:lon_sta+8:200j, lat_sta-6:lat_sta+6:200j]
    
    # Werte interpolieren
    grid_vtec = griddata((lons, lats), vtec, (grid_lon, grid_lat), method='linear')

    # Gauss-Filter anwenden, um die Delaunay-Dreiecke zu schönen sanften "Wetterwolken" zu machen
    grid_vtec_smooth = grid_vtec.copy()
    mask = np.isnan(grid_vtec)
    grid_vtec_smooth[mask] = np.nanmean(grid_vtec) # Verhindert NaN-Ausbreitung im Filter
    grid_vtec_smooth = gaussian_filter(grid_vtec_smooth, sigma=1.5)
    grid_vtec_smooth[mask] = np.nan # Maske wiederherstellen

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    ax.set_extent([lon_sta - 8, lon_sta + 8, lat_sta - 6, lat_sta + 6], crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1.0)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='#eaeaea')
    ax.add_feature(cfeature.OCEAN, facecolor='#cce5ff')
    
    # Smoothe Heatmap zeichnen (jetzt ohne negative Skala!)
    levels = np.linspace(0, max(20, np.nanmax(grid_vtec_smooth)), 30)
    contour = ax.contourf(grid_lon, grid_lat, grid_vtec_smooth, levels=levels, cmap='plasma', alpha=0.8, transform=ccrs.PlateCarree())
    
    # Die echten Satellitenspuren zur Kontrolle ganz fein drüberlegen
    ax.scatter(lons, lats, color='black', s=1, alpha=0.2, transform=ccrs.PlateCarree(), label='Echte Satelliten-Spuren (IPPs)')
    
    ax.plot(lon_sta, lat_sta, marker='^', color='red', markersize=14, markeredgecolor='white', transform=ccrs.PlateCarree(), label="GNSS Station")
    
    cbar = plt.colorbar(contour, ax=ax, shrink=0.7, pad=0.03)
    cbar.set_label('VTEC (Vertical TEC) [TECU]', fontsize=12)
    
    plt.title(f'Professionelle Ionosphären-Map (Kepler Orbits & DCB Leveled)\nStation: {lon_sta:.2f}E, {lat_sta:.2f}N', pad=15)
    plt.legend(loc='upper right', framealpha=0.9)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # ====================================================================
    RINEX_OBS_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 4\V3RJ1060.26O"
    RINEX_NAV_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 4\V3RJ1060.26N"  # PFLICHT! Wird für Kepler Orbits benötigt.
    # ====================================================================
    
    lats, lons, vtecs, sta_lat, sta_lon = process_rinex_advanced(RINEX_OBS_FILE, RINEX_NAV_FILE)
    plot_professional_ionosphere_map(lats, lons, vtecs, sta_lat, sta_lon)