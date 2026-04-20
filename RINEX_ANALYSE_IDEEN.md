# RINEX-Analyseideen (zusätzlich zu VTEC)

Ausgangsstruktur Deiner Daten:
- `.26O` = Observation (Code/Phase/SNR, Zeitreihe)
- `.26N` = GPS Navigation
- `.26L` = Galileo Navigation
- `.26P` = Mixed Navigation (mehrere GNSS)

Dein aktuelles Repo hat bereits ein VTEC-Skript (`main.py`).  
Hier sind weitere sinnvolle Auswertungen (inkl. **WIE?**) basierend auf den genannten Themen aus den Bildern (u. a. SPP, PPP, RTK, PVT, Multipath, Genauigkeitsvergleich) plus zusätzlichen Vorschlägen:

---

## 1) SPP (Single Point Positioning) – Basisgenauigkeit
**Ziel:** Wie genau ist eine einfache Pseudorange-Lösung?

**WIE?**
1. Lade `.26O` + passende Navigation (`.26N/.26L/.26P`).
2. Nutze nur Code-Messungen (Pseudorange), keine Carrier-Phase.
3. Schätze Position + Empfängeruhr pro Epoche.
4. Vergleiche gegen Referenzkoordinate (z. B. bekannte Antennenposition) und berechne RMS/95%-Fehler.
5. Wiederhole getrennt für GPS-only, Galileo-only, Multi-GNSS.

**Existierendes Script/Programm?**
- **Ja:** RTKLIB (`rnx2rtkp` SPP), gLAB.
- **Python-Bibliotheken vorhanden:** `gnss-lib-py` (teils), `georinex` fürs Parsing.
- **Hinweis:** Komplett fertiges, universelles SPP-Python-Repo für *deine* Dateistruktur ist meist nicht 1:1 plug-and-play.

---

## 2) Professionelles GNSS vs. Smartphone
**Ziel:** Genauigkeitsvergleich nach Gerätetyp.

**WIE?**
1. Zwei Datensätze aufnehmen (gleicher Ort/Zeit): geodätischer Empfänger vs Smartphone-RINEX.
2. Für beide identische Auswertepipeline nutzen (SPP und/oder PPP).
3. Fehlerstatistiken vergleichen (Horizontal/Vertical RMS, CEP50/95).
4. Unterschiede in SNR, Mehrwege, Tracking-Unterbrechungen auswerten.

**Existierendes Script/Programm?**
- **Teilweise:** Viele Papers/Notebooks, aber selten als einheitliches End-to-End-Tool.
- **Ja (Bausteine):** RTKLIB, gLAB, Google Smartphone Decimeter Datensätze/Beispiele.

---

## 3) Welche Orte sind am schlechtesten – und warum?
**Ziel:** Standortabhängige Fehlerquellen quantifizieren.

**WIE?**
1. Daten an mehreren Umgebungen sammeln (freies Feld, urban canyon, Bäume, Innenstadtnähe).
2. Pro Standort SPP/PPP-Fehler berechnen.
3. Zusätzliche Qualitätsmerkmale ableiten: SNR-Verteilung, Cycle Slips, Satellitenanzahl, DOP.
4. Korrelation Fehler ↔ Qualitätsmerkmale auswerten.
5. Ergebnis als Ranking + Ursachenanalyse dokumentieren.

**Existierendes Script/Programm?**
- **Kein einzelnes Standard-Tool**, das dir automatisch „schlechtester Ort + Ursache“ liefert.
- **Machbar mit:** Python (`pandas`, `georinex`, `matplotlib`) + RTKLIB/gLAB Outputs.

---

## 4) Einfluss von Satellitenanzahl und Frequenzen
**Ziel:** Wie stark verbessern mehr Satelliten/Signale die Genauigkeit?

**WIE?**
1. Auswertung mit künstlichen Teilmengen: nur 4, 6, 8, 10+ Satelliten.
2. Single-Frequency vs Dual-Frequency vergleichen.
3. GPS-only vs Galileo-only vs Multi-GNSS vergleichen.
4. Pro Fall Positionsfehler und DOP-Werte ausgeben.

**Existierendes Script/Programm?**
- **Teilweise:** RTKLIB/gLAB können Konstellationen/Frequenzen selektieren.
- **Kein universelles One-click-Repo** für alle Vergleichsmatrizen.

---

## 5) PPP (Precise Point Positioning) inkl. Vergleich zu RTK
**Ziel:** Präzisionsgewinn, Konvergenzzeit, Grenzen.

**WIE?**
1. Statische `.26O` + Nav-Dateien vorbereiten.
2. Präzise Orbit/Clock Produkte (Final/Rapid) laden (z. B. IGS/GSSC).
3. PPP rechnen (gLAB oder CSRS-PPP Webservice).
4. Kennzahlen: Konvergenzzeit, stationäre Genauigkeit, Ausreißer.
5. RTK-Resultate danebenstellen (falls Basisdaten vorhanden).
6. Vergleichstabelle erstellen (SPP vs PPP vs RTK).

**Existierendes Script/Programm?**
- **Ja:** gLAB, RTKLIB PPP, CSRS-PPP.
- **Hinweis:** CSRS-PPP ist Webdienst (Account nötig).

---

## 6) RTK-Auswertung (falls Basisstation/Referenzdaten vorhanden)
**Ziel:** Zentimeterbereich prüfen.

**WIE?**
1. Rover `.26O` + Base-RINEX synchronisieren.
2. Gemeinsame Epochen, gemeinsame Satelliten, Qualitätsfilter setzen.
3. RTK-Fix/Float-Lösung berechnen.
4. Fix-Rate, Time-to-Fix, Positionsfehler auswerten.
5. Falls vorhanden: Vergleich mit ESOC/GSSC-Referenzantenne.

**Existierendes Script/Programm?**
- **Ja:** RTKLIB (Standard), auch BKG Ntrip Client + RTK-Workflows.
- **Kein einzelnes einfaches GitHub-Script** deckt alle Kampagnenfälle automatisch ab.

---

## 7) PVT-Zeitreihe (Position, Velocity, Time) aus kinematischen Daten
**Ziel:** Bewegungsprofil und Geschwindigkeitsgüte.

**WIE?**
1. Epoche-für-Epoche Position schätzen (SPP/PPP/RTK).
2. Geschwindigkeit aus Positionsdifferenz oder Doppler ableiten.
3. Zeitverlauf für Position/Geschwindigkeit plotten.
4. Gegen Referenz (z. B. bekannte Strecke/IMU/odometer) validieren.
5. Fehler nach Umgebung vergleichen.

**Existierendes Script/Programm?**
- **Ja (Bausteine):** RTKLIB Outputs + Python-Postprocessing.
- **Selten als fertiges Komplettskript** exakt für deine Datenpipeline.

---

## 8) Multipath- und Qualitätsanalyse
**Ziel:** Einfluss von Mehrwegeffekten auf Genauigkeit.

**WIE?**
1. Nutze `.26O` Felder (SNR/CN0, Code-Phase-Kombinationen).
2. Erzeuge Skyplot + SNR-vs-Elevation + Zeitplots.
3. Detektiere Cycle Slips und Phasenunterbrechungen.
4. Bilde Multipath-Indikatoren (MP1/MP2, code-minus-carrier).
5. Vergleiche Fehlerstatistik mit/ohne problematische Epochen.

**Existierendes Script/Programm?**
- **Ja (teilweise):** teqc (historisch), RTKLIB-Tools, georinex + eigene Python-Skripte.
- **Kein universelles modernes Kompletttool** für alle Multipath-KPIs „out of the box“.

---

## 9) Zusätzliche eigene Idee: DOP-/Geometrie-Heatmaps
**Ziel:** Sichtbar machen, wann Geometrie statt Signalqualität limitiert.

**WIE?**
1. Berechne GDOP/PDOP/HDOP/VDOP pro Epoche.
2. Kombiniere mit Positionsfehler und Satellitenzahl.
3. Plotte Tages-Heatmap (Zeit vs DOP vs Fehler).

**Existierendes Script/Programm?**
- **Ja:** RTKLIB liefert DOP; Auswertung meist per eigenem Python-Skript.

---

## 10) Zusätzliche eigene Idee: Ionosphäre-Dynamik (ROTI/Storm-Monitoring)
**Ziel:** Neben VTEC auch kurzfristige ionosphärische Unruhe erfassen.

**WIE?**
1. Aus dual-frequency Phase Slant TEC ableiten.
2. ROT/ROTI pro Satellit und Zeitfenster berechnen.
3. Ruhe- vs Störungsperioden vergleichen.
4. Optional mit Space-Weather-Indizes (Kp, Dst) korrelieren.

**Existierendes Script/Programm?**
- **Teilweise:** Forschungsrepos vorhanden, aber oft nicht direkt produktionsreif.
- **Meist eigener Analysecode nötig** (Python/Matlab).

---

## 11) Zusätzliche eigene Idee: Qualitäts-Dashboard pro Datei (.26O/.26N/.26L/.26P)
**Ziel:** Sofort sehen, ob ein Datensatz wissenschaftlich „gut genug“ ist.

**WIE?**
1. Parser-Checks: fehlende Epochen, Beobachtungstypen, Satellitenverfügbarkeit.
2. Kennzahlen: Datenlücken, Slip-Rate, mittlere SNR, nutzbare Dual-Frequency-Anteile.
3. Ausgabe als Scorecard (A/B/C) + Ampel.

**Existierendes Script/Programm?**
- **Kein verbreitetes Standardtool** mit genau dieser kompakten Scorecard.
- **Sehr gut umsetzbar** als eigenes Python-Modul in diesem Repo.

---

## Empfohlene Reihenfolge für dein Repo
1. **SPP-Baselines** (schnell, robust, sofort vergleichbar)  
2. **Multipath/SNR/Cycle-Slip Diagnostik** (Ursachen verstehen)  
3. **PPP (gLAB/CSRS-PPP) + Vergleichstabelle**  
4. **Optional RTK**, wenn Basisdaten/NTRIP vorliegen  
5. **Erweiterung auf PVT/ROTI/Dashboard**

Damit bekommst Du von „nur VTEC“ zu einer vollständigen GNSS-Qualitäts- und Genauigkeitsanalyse.
