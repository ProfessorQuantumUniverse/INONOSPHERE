# Was man mit deinen aufgezeichneten RINEX-Daten noch machen kann (neben VTEC)

## Kurz zur Ausgangslage
Im Repo gibt es bereits ein Python-Skript (`/home/runner/work/INONOSPHERE/INONOSPHERE/main.py`) zur VTEC-/Ionosphärenauswertung mit `georinex`, IPP-Berechnung und Kartenplot.

> Hinweis: Ich konnte nur **1 Bild** auslesen. Die Ideen unten decken dessen Inhalte vollständig ab und ergänzen sie um weitere sinnvolle Auswertungen.

---

## A) Fragen aus dem Bild als konkrete Auswertungsprojekte

## 1) Wie genau ist professionelles GNSS vs. Smartphone?
**WIE?**
- Gleiche Strecke/Standorte mit professionellem Empfänger und Smartphone aufzeichnen (möglichst gleichzeitig).
- Beide Datensätze mit gleichem Verfahren lösen (z. B. SPP und PPP), dann gegen Referenzposition vergleichen.
- Metriken: 2D/3D-RMSE, 95%-Fehler, Medianfehler, Ausreißerquote.
- Ergebnisse nach Umgebung aufteilen (offen, urban, Bäume, Indoor-nah).

**Existierende Skripte/Programme?**
- Ja: **RTKLIB**, **gLAB**, **GNSS-SDR** (für Auswertung/Processing).
- Für „fertigen 1-Klick-Vergleich professionell vs. Smartphone“ gibt es i. d. R. **kein einzelnes Standardskript** → meist eigene Auswertepipeline (Python/Notebook).

---

## 2) Welche Orte haben die schlechteste Genauigkeit – und warum?
**WIE?**
- Pro Standort Positionsfehler berechnen.
- Parallel Qualitätsindikatoren extrahieren: Satellitenanzahl, DOP (GDOP/PDOP/HDOP), Elevationsverteilung, C/N0/SNR, Cycle Slips, Multipath-Indikatoren.
- Standort-Ranking erstellen (best/worst).
- Ursachenanalyse: schlechte Geometrie, Abschattung, Mehrwegeffekte, Signalunterbrechungen.

**Existierende Skripte/Programme?**
- Teilweise: **RTKLIB (rtkplot/qc)**, **teqc** (historisch), **georinex + eigene Python-Auswertung**.
- Vollständige automatische „Warum schlechtester Standort?“-Diagnose: meist **kein fertiges End-to-End-Tool**, eher Kombination mehrerer Tools.

---

## 3) Wie verbessern mehr Satelliten/Frequenzen die Genauigkeit?
**WIE?**
- Szenarien rechnen: 
  - nur GPS,
  - GPS+Galileo,
  - Single-Frequency vs. Dual-Frequency,
  - begrenzte Satellitenzahl (Subsampling).
- Für jedes Szenario gleiche Trajektorie/Station neu lösen.
- Fehlerstatistik und Konvergenzzeit vergleichen.

**Existierende Skripte/Programme?**
- Ja: **RTKLIB**, **gLAB**, **CSRS-PPP** (servicebasiert, für PPP-Vergleich), teils über Konfigurationsprofile.
- Automatisierte Batch-„Was-wäre-wenn“-Studie häufig **eigene Skripte**.

---

## B) Themenblöcke aus dem Bild

## 4) SPP (Single Point Positioning)
**WIE?**
- Pseudorange-basierte Lösung aus OBS+NAV rechnen.
- Laufend Positionfehler, DOP und Satellitenanzahl protokollieren.
- Genauigkeit je Umgebung/Frequenzkonfiguration vergleichen.

**Existierende Skripte/Programme?**
- Ja: **RTKLIB**, **gLAB**, **GNSS-SDR**.
- In Python gibt es Beispiele, aber selten ein universelles „fertiges Wissenschafts-Skript“ für alle Konstellationen.

---

## 5) PPP (Precise Point Positioning)
**WIE?**
- Statische Daten im PPP-Modus verarbeiten.
- Final/Rapid Orbit- und Clock-Produkte (IGS/GSSC) verwenden.
- Lösung mit **CSRS-PPP** oder lokal mit **gLAB** rechnen.
- PPP mit RTK vergleichen (Genauigkeit, Konvergenzzeit, Infrastrukturbedarf).
- Abschlusstabelle: Methode | Genauigkeit | Konvergenz | Aufwand | Einschränkungen.

**Existierende Skripte/Programme?**
- Ja: **CSRS-PPP** (Webservice), **gLAB**, **RTKLIB PPP**.
- Vollständiger Vergleichsbericht (inkl. Tabellen/Plots) bleibt meist **eigene Analysearbeit**.

---

## 6) RTK
**WIE?**
- Rover-RINEX mit Basisstation-RINEX (z. B. ESOC/GSSC) kombinieren.
- Fixed/Float-Rate auswerten, Baseline-Länge berücksichtigen.
- Genauigkeit gegen Referenz vergleichen.

**Existierende Skripte/Programme?**
- Ja: **RTKLIB (rtkrcv/rnx2rtkp)**, **goGPS**.
- Für deine spezifische Messkampagne ist meist Konfigurations-/Skriptanpassung nötig.

---

## 7) PVT (Position, Velocity, Time) aus kinematischen Daten
**WIE?**
- Zeitreihen für Position und Geschwindigkeit erzeugen.
- Geschwindigkeitsprofile (z. B. Stop-and-Go, Beschleunigung, Kurven) plotten.
- Vergleich der Geschwindigkeitsqualität je Umgebung (offen vs. urban).

**Existierende Skripte/Programme?**
- Ja: **RTKLIB**, **GNSS-SDR**, teils **gLAB**-Workflows.
- Für saubere Bewegungssegmentierung + Qualitätsmetriken oft **eigene Python-Skripte**.

---

## 8) Genauigkeit, Multipath und Satellitengeometrie
**WIE?**
- Aus RINEX C/N0 (SNR), Cycle Slips, Elevation/Azimut, DOP extrahieren.
- Plots erstellen: 
  - SNR vs. Elevation,
  - Fehler vs. DOP,
  - Cycle-Slip-Rate pro Standort,
  - Skyplot / Sichtbarkeitsfenster.
- Multipath-Indikatoren berechnen (z. B. code-carrier Kombis) und mit Fehler korrelieren.

**Existierende Skripte/Programme?**
- Teilweise: **georinex** (Einlesen), **RTKLIB-Tools**, verschiedene Forschungsnotebooks auf GitHub.
- Einheitliches „komplett fertiges Multipath-Labor“ für alle RINEX-Varianten ist selten → häufig **eigene Pipeline**.

---

## C) Zusätzliche sinnvolle Ideen (eigene Ergänzungen)

## 9) ROT/ROTI und ionosphärische Aktivität (aus deinen TEC-Zeitreihen)
**WIE?**
- Aus STEC/VTEC-Zeitreihen ROT und ROTI berechnen.
- Zeitfenster mit starker ionosphärischer Unruhe markieren.
- Mit GNSS-Fehlern (Position, Cycle Slips) korrelieren.

**Existierende Skripte/Programme?**
- Forschungsbeispiele vorhanden, aber oft projektspezifisch.
- Für deine Daten wahrscheinlich **eigene Implementierung** am sinnvollsten.

---

## 10) DCB-Bias-Abschätzung (Empfänger/Satellit)
**WIE?**
- Dual-Frequency-Kombinationen für Bias-Schätzung nutzen.
- Tägliche Stabilität der Biases prüfen.
- Einfluss auf TEC- und Positionsgenauigkeit quantifizieren.

**Existierende Skripte/Programme?**
- Teils in Forschungssoftware vorhanden.
- Für robuste, reproduzierbare DCB-Auswertung im eigenen Setup oft **kein direkt fertiges Universaltool**.

---

## 11) Qualitätsmonitoring/Anomalie-Erkennung deiner Station
**WIE?**
- Tägliche KPI-Pipeline: Datenvollständigkeit, Slip-Rate, SNR-Statistik, Satellitenverfügbarkeit, Ausreißer.
- Schwellenwerte/Alarmregeln definieren.
- Trends über Wochen/Monate auswerten.

**Existierende Skripte/Programme?**
- Teilweise mit QC-Tools möglich.
- Vollständiges, auf deine Station angepasstes Monitoring ist meist **Eigenbau**.

---

## 12) Troposphäre (ZTD) und Wetterbezug
**WIE?**
- Aus PPP-Lösungen ZTD-Parameter extrahieren.
- Mit lokalen Wetterdaten (Druck/Feuchte/Regen) vergleichen.
- Prüfen, ob Wetterereignisse die GNSS-Qualität beeinflussen.

**Existierende Skripte/Programme?**
- PPP-Software liefert teils ZTD.
- Datensynchronisation + Analyse meist **eigene Skriptarbeit**.

---

## Praktische Empfehlung für dein Repo (nächste Schritte)
1. Bestehendes VTEC-Skript als Modul aufteilen (Einlesen, Geometrie, Plotting).
2. Neues Auswerte-Notebook/Script für **SPP/PPP/RTK-Vergleichstabelle** anlegen.
3. Gemeinsames KPI-Schema definieren (RMSE, 95%, DOP, SNR, Slip-Rate, Konvergenzzeit).
4. Automatische Reports als Markdown/CSV/Plots erzeugen.

Damit hättest du aus denselben RINEX-Daten nicht nur Ionosphäre, sondern auch Positionierungs-Performance, Multipath-Qualität und Umgebungsabhängigkeit wissenschaftlich vergleichbar dokumentiert.
