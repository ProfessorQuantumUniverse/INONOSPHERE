import georinex as gr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

def ecef_to_latlon(x, y, z):
    """Konvertiert ECEF-Koordinaten in Breitengrad und Längengrad."""
    a = 6378137.0
    e2 = 0.00669437999014
    b = np.sqrt(a**2 * (1 - e2))
    ep = np.sqrt((a**2 - b**2) / b**2)
    p = np.sqrt(x**2 + y**2)
    th = np.arctan2(a * z, b * p)
    lon = np.arctan2(y, x)
    lat = np.arctan2((z + ep**2 * b * np.sin(th)**3), (p - e2 * a * np.cos(th)**3))
    return np.degrees(lat), np.degrees(lon)

def get_tec_multiplier(sys_id):
    """Gibt den Umrechnungsfaktor für P2-P1 zu STEC basierend auf den Frequenzen zurück."""
    # Formel: M = (f1^2 * f2^2) / (40.3 * 10^16 * (f1^2 - f2^2))
    if sys_id == 'G': # GPS L1 & L2 (Standard)
        return 9.52
    elif sys_id == 'E': # Galileo E1 & E5a (L1/L5 equivalent)
        return 7.76
    elif sys_id == 'R': # GLONASS G1 & G2 (angenäherte Center-Frequenzen)
        return 9.25
    else:
        return 9.52 # Fallback

def calculate_stec(obs_file, nav_file=None):
    """Liest die RINEX Datei und berechnet STEC für GPS, Galileo und GLONASS."""
    print("⏳ Lade RINEX Dateien in den Arbeitsspeicher (das kann kurz dauern)...")
    obs = gr.load(obs_file)
    
    if nav_file:
        print(f"⏳ Lade Navigationsdaten ({nav_file})...")
        nav = gr.load(nav_file)
        # Hier könnten künftig mit einer Library wie 'gnsscalcs' die exakten 
        # Kepler-Orbits berechnet werden.
    
    stec_results = {}
    
    # Listen der bevorzugten Observablen für L1 und L2/L5
    l1_bands = ['C1C', 'C1', 'P1', 'L1']
    l2_bands =['C2W', 'C2C', 'P2', 'C2', 'C5Q', 'C5', 'L2', 'L5']

    sv_list = obs.sv.values
    
    # Progress Bar für die Verarbeitung der Satelliten
    for sv in tqdm(sv_list, desc="🛰️ Verarbeite Satelliten"):
        sv_id = str(sv)
        sys_id = sv_id[0] # 'G' für GPS, 'E' für Galileo, 'R' für GLONASS
        
        if sys_id not in ['G', 'E', 'R']:
            continue
            
        sv_data = obs.sel(sv=sv)
        
        # Finde die besten verfügbaren Frequenzen für diesen Satelliten
        p1_var = next((v for v in l1_bands if v in sv_data.data_vars), None)
        p2_var = next((v for v in l2_bands if v in sv_data.data_vars), None)
        
        if p1_var and p2_var:
            P1 = sv_data[p1_var]
            P2 = sv_data[p2_var]
            
            # STEC berechnen
            multiplier = get_tec_multiplier(sys_id)
            stec = multiplier * (P2 - P1)
            
            # NaN-Werte entfernen
            stec_clean = stec.dropna(dim='time')
            
            if len(stec_clean) > 10: # Nur Satelliten mit ausreichend Daten
                stec_results[sv_id] = {
                    'time': stec_clean.time.values,
                    'stec': stec_clean.values,
                    'mean_stec': float(stec_clean.mean().values)
                }
                
    return stec_results, obs

def plot_temporal_heatmap(stec_dict):
    """Erstellt eine Matrix-Heatmap (Satelliten vs. Zeit) wie in RXTools."""
    print("📊 Erstelle Zeit-Satelliten Heatmap...")
    
    svs = list(stec_dict.keys())
    svs.sort() # Satelliten sortieren (E.., G.., R..)
    
    if not svs:
        print("Keine ausreichenden TEC Daten für die Heatmap gefunden.")
        return

    # Sammle alle existierenden Zeitstempel
    all_times =[]
    for sv in svs:
        all_times.extend(stec_dict[sv]['time'])
    all_times = np.unique(all_times)
    all_times.sort()
    
    # Matrix mit NaNs initialisieren
    heatmap_data = np.full((len(svs), len(all_times)), np.nan)
    
    for i, sv in enumerate(tqdm(svs, desc="🎨 Generiere Heatmap Matrix")):
        times = stec_dict[sv]['time']
        vals = stec_dict[sv]['stec']
        # Finde die Indizes der Zeitstempel
        indices = np.searchsorted(all_times, times)
        heatmap_data[i, indices] = vals

    plt.figure(figsize=(12, 6))
    plt.imshow(heatmap_data, aspect='auto', cmap='plasma', interpolation='nearest',
               extent=[0, len(all_times), len(svs), 0])
    
    plt.colorbar(label='STEC [TECU]')
    plt.yticks(np.arange(len(svs)) + 0.5, svs, fontsize=8)
    plt.xlabel('Zeitverlauf (Epochen)')
    plt.ylabel('Satelliten (PRN)')
    plt.title('Ionosphären-Aktivität (Heatmap: Zeit vs. Satellit)')
    plt.tight_layout()
    plt.show()

def plot_spatial_heatmap(stec_dict, obs):
    """Erstellt eine geografische Heatmap über einer Europakarte."""
    print("🌍 Erstelle geografische Spatial-Heatmap...")
    
    # Stationskoordinaten ermitteln
    if 'position' not in obs.attrs:
        lat_sta, lon_sta = 50.0, 10.0 # Zentrum Deutschland als Fallback
    else:
        x, y, z = obs.attrs['position']
        lat_sta, lon_sta = ecef_to_latlon(x, y, z)

    # Wir erzeugen simulierte Durchstoßpunkte (IPPs) um die Station herum
    # (Für echte IPPs müssten die Navigationsdaten (.26n) mittels Kepler-Gleichung aufgelöst werden)
    lats, lons, tec_vals = [], [],[]
    
    np.random.seed(42) # Für reproduzierbare Pseudoverteilung
    for sv, data in stec_dict.items():
        # Simuliere einen Himmels-Einstichpunkt ca. 2-6 Grad entfernt von der Station
        # Abhängig vom Satelliten-Namen
        angle = hash(sv) % 360
        distance = 1.5 + (hash(sv[::-1]) % 40) / 10.0 
        
        ipp_lat = lat_sta + distance * np.cos(np.radians(angle))
        ipp_lon = lon_sta + distance * np.sin(np.radians(angle))
        
        lats.append(ipp_lat)
        lons.append(ipp_lon)
        tec_vals.append(data['mean_stec'])

    lats = np.array(lats)
    lons = np.array(lons)
    tec_vals = np.array(tec_vals)

    # 1. Erstelle ein Grid (Gitter) für die Heatmap-Interpolation
    grid_lon, grid_lat = np.mgrid[lon_sta-8:lon_sta+8:100j, lat_sta-6:lat_sta+6:100j]
    
    # 2. Interpoliere die Werte der einzelnen Satelliten zu einer weichen Heatmap
    grid_tec = griddata((lons, lats), tec_vals, (grid_lon, grid_lat), method='cubic', fill_value=np.nan)

    # 3. Karte zeichnen
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    ax.set_extent([lon_sta - 8, lon_sta + 8, lat_sta - 6, lat_sta + 6], crs=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, color='#e0e0e0')
    ax.add_feature(cfeature.OCEAN, color='#b0d0e0')
    
    # Heatmap (Contourf) über die Karte legen
    contour = ax.contourf(grid_lon, grid_lat, grid_tec, levels=20, cmap='jet', alpha=0.6, transform=ccrs.PlateCarree())
    
    # Satelliten-Punkte und Station eintragen
    ax.scatter(lons, lats, color='black', s=10, marker='x', transform=ccrs.PlateCarree(), label='Simulierte IPP-Durchstoßpunkte')
    ax.plot(lon_sta, lat_sta, marker='^', color='red', markersize=12, transform=ccrs.PlateCarree(), label="GNSS Station")
    
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', shrink=0.7)
    cbar.set_label('Mittlerer STEC [TECU]')
    
    plt.title(f'Geografische Ionosphären-Heatmap (VTEC-Näherung)\nStation: {lon_sta:.2f}E, {lat_sta:.2f}N')
    plt.legend(loc='lower left')
    plt.show()

if __name__ == '__main__':
    # ====================================================================
    # DEINE DATEIEN HIER EINTRAGEN
    RINEX_OBS_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 4\V3RJ1060.26O"
    RINEX_NAV_FILE = r"C:\Users\loren\OneDrive\Documenten\GNSS\LOG 4\V3RJ1060.26N" # Optional. Wenn nicht vorhanden: None
    # ====================================================================
    
    try:
        # 1. Daten berechnen (mit Multi-Constellation Logik & Progress Bars)
        stec_data, rinex_obs = calculate_stec(RINEX_OBS_FILE, RINEX_NAV_FILE)
        
        if len(stec_data) > 0:
            print(f"✅ TEC für {len(stec_data)} Satelliten erfolgreich berechnet.")
            
            # 2. Zeitliche Matrix-Heatmap (Zeit vs. Satellit)
            plot_temporal_heatmap(stec_data)
            
            # 3. Geografische Heatmap (über Europa-Karte)
            plot_spatial_heatmap(stec_data, rinex_obs)
        else:
            print("❌ Konnte keine ausreichenden Daten für eine Darstellung extrahieren.")
            
    except FileNotFoundError as e:
        print(f"❌ Fehler: Datei nicht gefunden. ({e})")
    except Exception as e:
        print(f"❌ Ein Fehler ist aufgetreten: {e}")