# RINEX-Analyseideen (zusätzlich zu VTEC)

Ausgangsstruktur deiner Daten:
- `.26O` = Observation (Code/Phase/SNR, Zeitreihe)
- `.26N` = GPS Navigation
- `.26L` = Galileo Navigation
- `.26P` = Mixed Navigation (mehrere GNSS)

Das Repo enthält bereits VTEC-Logik (`main.py`).  
Die folgenden Ideen sind **grundsätzlich sinnvoll**. Der Hauptunterschied ist der benötigte Aufwand und die Datenvoraussetzungen.

## Kurzfazit: Machen die Ideen Sinn?

Ja – alle 11 Ideen sind fachlich sinnvoll. Priorisierung:

- **Sofort sinnvoll (mit vorhandenen RINEX-Daten):** 1, 3, 4, 8, 9, 11
- **Sinnvoll mit zusätzlichem Aufwand/Daten:** 5 (präzise Produkte), 6 (Base/NTRIP), 7 (kinematische Referenz)
- **Sinnvoll mit externer Zusatzquelle:** 2 (Smartphone-Daten), 10 (Space-Weather-Indizes optional)

Wichtig für belastbare Ergebnisse:
- Einheitliche Zeitsysteme (GPS/GAL), saubere Referenzkoordinaten, konsistente Qualitätsfilter
- Klare Trennung zwischen statischen und kinematischen Sessions
- Statistische Kennzahlen immer mit Stichprobengröße (`N`) berichten

---

## 1) SPP (Single Point Positioning) – Basisgenauigkeit
**Sinnhaftigkeit:** Sehr hoch. Beste Start-Baseline vor PPP/RTK.

**Ziel:** Wie genau ist eine einfache Pseudorange-Lösung?

**Erweiterte Anleitung**
1. `.26O` + passende Navigation (`.26N/.26L/.26P`) laden.
2. Vorfilter definieren (Elevation-Maske, SNR-Minimum, Ausreißergrenze).
3. Nur Code-Messungen nutzen (keine Carrier-Phase).
4. Pro Epoche Position + Empfängeruhr schätzen.
5. In ENU gegen Referenzposition auswerten.
6. Kennzahlen je Modus berichten: RMS, Median, 95%-Perzentil, max. Fehler.
7. Getrennt rechnen: GPS-only, Galileo-only, Multi-GNSS.

**Ergebnisartefakte**
- Zeitreihe Ost/Nord/Höhe
- Tabelle der Fehlerkennzahlen je Konfiguration

**Existierendes Script/Programm**
- Ja: RTKLIB (`rnx2rtkp`), gLAB
- Python-Bausteine: `georinex`, `gnss-lib-py` (teilweise)

---

## 2) Professionelles GNSS vs. Smartphone
**Sinnhaftigkeit:** Hoch, wenn identische Messbedingungen eingehalten werden.

**Ziel:** Genauigkeitsvergleich nach Gerätetyp.

**Erweiterte Anleitung**
1. Gleiche Strecke/Ort/Zeit messen (möglichst simultan).
2. Beide Datensätze auf denselben Analysepfad bringen (gleiche Filter, gleiche Metriken).
3. SPP-Benchmark für beide rechnen, optional PPP ergänzen.
4. Vergleichsmetriken: Horizontal/Vertical RMS, CEP50/95, Datenverfügbarkeit.
5. Qualitätsvergleich: SNR, Cycle-Slip-Rate, Tracking-Ausfälle.
6. Ergebnis nach Umgebung splitten (offen, urban, teilverdeckt).

**Achtung**
- Smartphone-Antenne und Duty-Cycle können Ergebnisse stark verzerren.
- Nur fair vergleichen, wenn Zeitfenster und Satellitenkonfiguration ähnlich sind.

**Existierendes Script/Programm**
- Bausteine vorhanden (RTKLIB, gLAB, Beispielnotebooks), meist kein einheitliches End-to-End-Tool

---

## 3) Welche Orte sind am schlechtesten – und warum?
**Sinnhaftigkeit:** Sehr hoch, liefert direkt umsetzbare Erkenntnisse.

**Ziel:** Standortabhängige Fehlerquellen quantifizieren.

**Erweiterte Anleitung**
1. Mehrere Standortklassen definieren (offen, urban canyon, Baumdach etc.).
2. Pro Standort identische Auswerteparameter verwenden.
3. Fehlerkennzahlen pro Standort berechnen.
4. Qualitätsmerkmale ableiten: SNR-Verteilung, DOP, Satellitenanzahl, Slip-Rate.
5. Korrelationen berechnen (z. B. Fehler vs. DOP/SNR).
6. Ranking erstellen und je Standort die Hauptursache dokumentieren.

**Ergebnisartefakte**
- Ranking-Tabelle „schlechtester bis bester Ort“
- Kurzbegründung pro Standort (Geometrie, Abschattung, Multipath)

---

## 4) Einfluss von Satellitenanzahl und Frequenzen
**Sinnhaftigkeit:** Hoch, methodisch gut zur Sensitivitätsanalyse.

**Ziel:** Wie stark verbessern mehr Satelliten/Signale die Genauigkeit?

**Erweiterte Anleitung**
1. Experimente mit Teilmengen aufsetzen (4/6/8/10+ Satelliten).
2. Single- vs Dual-Frequency getrennt auswerten.
3. GPS-only, Galileo-only, Multi-GNSS separat rechnen.
4. Pro Fall DOP + Fehlerstatistiken berichten.
5. Optional Robustheit prüfen (z. B. künstliches Entfernen niedriger Elevationen).

**Achtung**
- Sehr kleine Satellitenmengen können numerisch instabil werden.
- Vergleich nur mit konsistenten Filtern interpretieren.

---

## 5) PPP (Precise Point Positioning) inkl. Vergleich zu RTK
**Sinnhaftigkeit:** Sehr hoch für Präzision, aber deutlich komplexer.

**Ziel:** Präzisionsgewinn, Konvergenzzeit, Grenzen.

**Erweiterte Anleitung**
1. Statische RINEX-Sessions mit ausreichender Dauer auswählen.
2. Präzise Orbit/Clock/ggf. Bias-Produkte (IGS) laden.
3. PPP-Lauf durchführen (gLAB, RTKLIB-PPP oder CSRS-PPP).
4. Konvergenzphase klar definieren (z. B. bis Fehler < Schwellwert stabil).
5. Stationäre Phase separat statistisch auswerten.
6. Mit SPP und ggf. RTK in einer Vergleichsmatrix gegenüberstellen.

**Achtung**
- PPP-Konvergenz kann je Umgebung und Datenqualität stark schwanken.
- Exakte Ergebnisreproduzierbarkeit benötigt dokumentierte Produktversionen.

---

## 6) RTK-Auswertung (falls Basisstation/Referenzdaten vorhanden)
**Sinnhaftigkeit:** Sehr hoch, wenn Basisdaten sauber synchronisiert sind.

**Ziel:** Zentimeterbereich prüfen.

**Erweiterte Anleitung**
1. Rover- und Base-Daten zeitlich synchronisieren.
2. Referenzkoordinate der Base verifizieren.
3. RTK-Lösung rechnen (Float/Fix), Qualitätskriterien definieren.
4. Kennzahlen ausgeben: Fix-Rate, Time-to-Fix, RMS horizontal/vertikal.
5. Fehler über Zeit und bei Fix↔Float-Übergängen analysieren.

**Achtung**
- Ohne saubere Basisreferenz sind cm-Aussagen nicht belastbar.
- NTRIP-Latenz/Signalabbrüche können Fix-Rate stark drücken.

---

## 7) PVT-Zeitreihe (Position, Velocity, Time) aus kinematischen Daten
**Sinnhaftigkeit:** Sinnvoll bei Bewegungsdaten; mit statischen Daten nur begrenzt.

**Ziel:** Bewegungsprofil und Geschwindigkeitsgüte.

**Erweiterte Anleitung**
1. Kinematische Session wählen, epochweise Position bestimmen.
2. Geschwindigkeit entweder aus Doppler oder differenzierter Position schätzen.
3. Zeitreihen glätten/filtern (ohne Dynamik zu verfälschen).
4. Mit Referenztrajektorie (IMU/Odometer/Survey-Track) vergleichen.
5. Fehler nach Bewegungsphase auswerten (Anfahren, Kurven, Stop-and-go).

**Achtung**
- Reine Differenzen verstärken Rauschen; Doppler ist oft stabiler für Geschwindigkeit.

---

## 8) Multipath- und Qualitätsanalyse
**Sinnhaftigkeit:** Sehr hoch, erklärt viele Ausreißer direkt.

**Ziel:** Einfluss von Mehrwegeffekten auf Genauigkeit.

**Erweiterte Anleitung**
1. Aus `.26O` SNR/CN0 und relevante Beobachtungstypen extrahieren.
2. Skyplot und SNR-vs-Elevation erstellen.
3. Cycle-Slips und Datenlücken markieren.
4. Multipath-Indikatoren berechnen (z. B. code-minus-carrier, MP1/MP2).
5. Fehlerstatistik für „alle Epochen“ vs „bereinigte Epochen“ vergleichen.

**Ergebnisartefakte**
- Multipath-Hotspots nach Azimut/Elevation
- Quantifizierter Einfluss auf Positionsfehler

---

## 9) Zusätzliche Idee: DOP-/Geometrie-Heatmaps
**Sinnhaftigkeit:** Hoch, schnell interpretierbar und sehr nützlich im Reporting.

**Ziel:** Sichtbar machen, wann Geometrie statt Signalqualität limitiert.

**Erweiterte Anleitung**
1. GDOP/PDOP/HDOP/VDOP pro Epoche berechnen.
2. Mit Satellitenanzahl und Fehlern zusammenführen.
3. Heatmap und Zeitfenster mit schlechter Geometrie markieren.
4. Prüfen, ob Fehlerspitzen zeitlich mit DOP-Spitzen zusammenfallen.

---

## 10) Zusätzliche Idee: Ionosphäre-Dynamik (ROTI/Storm-Monitoring)
**Sinnhaftigkeit:** Hoch, gute Ergänzung zu VTEC.

**Ziel:** Neben VTEC kurzfristige ionosphärische Unruhe erfassen.

**Erweiterte Anleitung**
1. Dual-Frequency-Phasen für Slant-TEC vorbereiten.
2. ROT (zeitliche TEC-Änderung) pro Satellit berechnen.
3. ROTI über gleitende Fenster bilden.
4. Ruhige und gestörte Perioden vergleichen.
5. Optional: Korrelation mit Kp/Dst/Solar-Wind-Indizes.

**Achtung**
- Strenge Outlier-/Cycle-Slip-Behandlung nötig, sonst wird ROTI verfälscht.

---

## 11) Zusätzliche Idee: Qualitäts-Dashboard pro Datei (.26O/.26N/.26L/.26P)
**Sinnhaftigkeit:** Sehr hoch, ideal für schnelle Datensatz-Freigabe.

**Ziel:** Sofort sehen, ob ein Datensatz wissenschaftlich „gut genug“ ist.

**Erweiterte Anleitung**
1. Format- und Parser-Checks (Header, Beobachtungstypen, Epochendichte).
2. Kennzahlen berechnen: Lückenquote, Slip-Rate, mittlere SNR, Dual-Frequency-Anteil.
3. Schwellenwerte definieren und Scorecard (A/B/C + Ampel) erzeugen.
4. Pro Datei Kurzdiagnose mit Hauptproblem ausgeben.
5. Dashboard als CSV/JSON + kompakte Grafik exportieren.

**Nutzen**
- Spart Zeit vor aufwendigen PPP/RTK-Läufen.
- Schafft reproduzierbare Qualitätskriterien.

---

## Empfohlene Reihenfolge für dieses Repo
1. **SPP-Baselines** (schnell, robust, direkt vergleichbar)
2. **Multipath/SNR/Cycle-Slip-Diagnostik** (Ursachen sichtbar machen)
3. **DOP-Heatmaps + Standortvergleich** (Geometrie vs Umgebung trennen)
4. **PPP-Vergleich** (Präzisionsgewinn und Konvergenz bewerten)
5. **Optional RTK** (wenn Base/NTRIP sauber vorliegt)
6. **Danach PVT/ROTI/Dashboard als Ausbau**

So entwickelst du das Repo von „nur VTEC“ zu einer vollständigen GNSS-Qualitäts- und Genauigkeitsanalyse mit klarer Priorisierung.
