# Begleitdokument: Hyperelastizität & Arruda-Boyce Modell

Dieses Dokument erklärt die theoretischen Hintergründe der Scilab-Implementierung des Arruda-Boyce Materialmodells für hyperelastische Werkstoffe. Es dient als Brücke zwischen der akademischen Herleitung und dem praktischen Code, mit Fokus auf die einzelnen Rechenschritte, physikalischen Konzepte und deren Umsetzung.

---

## Table of Contents

1. [Aufgabe 1: Die Materialroutine (`F_to_sig`)](#1-aufgabe-1-die-materialroutine-f_to_sig)
   - 1.1 Der deviatorische Teil (Taylor-Reihe)
   - 1.2 Der volumetrische Teil
   - 1.3 Die 1. Piola-Kirchhoff Spannung
   - 1.4 Zusammensetzen der Cauchy-Spannung

2. [Aufgabe 2: Spannungs-Dehnungs-Kurven](#2-aufgabe-2-spannungs-dehnungs-kurven)
   - 2.1 Der Datenfluss im Programm
   - 2.2 Warum die Rücktransformation zu PK1?
   - 2.3 Initiale Werte für den Plot
   - 2.4 Übersicht der Lastfälle

3. [Aufgabe 3: Ploterstellung in Scilab](#3-aufgabe-3-programmierung-ploterstellung-in-scilab)
   - 3.1 Fenstersteuerung (`scf` & `clf`)
   - 3.2 Fenster unterteilen (`subplot`)
   - 3.3 Zeichnen und Beschriften
   - 3.4 Fenstergröße anpassen

4. [Aufgabe 4: Data Fitting (`lsqrsolve`)](#4-aufgabe-4-data-fitting-lsqrsolve)
   - 4.1 Der Solver-Befehl `lsqrsolve`
   - 4.2 Die Fehlerfunktion (Residuen)
   - 4.3 Multimodales Fitting

5. [Aufgabe 5: Tangentialsteifigkeit ($E_T$)](#5-aufgabe-5-tangentialsteifigkeit-e_t)
   - 5.1 Die 2. Piola-Kirchhoff Spannung
   - 5.2 Der Steifigkeitstensor

6. [Aufgabe 6: Glyphen-Visualisierung ($E_{nnnn}$)](#6-aufgabe-6-glyphen-visualisierung-e_nnnn)
   - 6.1 Programmierung der Glyphen
   - 6.2 Interpretation der Glyphen-Formen
   - 6.3 Colorbar-Steuerung
   - 6.4 Kamera und Perspektive

7. [Aufgabe 7: Material- vs. Räumliche Steifigkeit](#7-aufgabe-7-analyse-material--vs-räumliche-steifigkeit)
   - 7.1 Referenz- vs. Momentankonfiguration
   - 7.2 Berücksichtigung von Querschnitt und Volumen
   - 7.3 Analytisch vs. Numerisch

8. [Nomenklatur & Abkürzungen](#nomenklatur--abkürzungen)

*Hinweis zur Notation: Die verwendeten kontinuumsmechanischen Herleitungen und Tensordefinitionen folgen primär der Nomenklatur von G. A. Holzapfel (Nonlinear Solid Mechanics, 2000) sowie dem Skript "Solid mechanics of continua".*

---

## 1. Aufgabe 1: Die Materialroutine (`F_to_sig`)

Das Arruda-Boyce-Modell (Skript Gl. 2.242) beschreibt das elastische Potenzial $W$ als Summe aus einem deviatorischen (Gestaltänderung) und einem volumetrischen (Volumenänderung) Teil [--> Holzapfel, Gl. 6.85]:

$$W(\boldsymbol{C}) = W_{dev}(\bar{I}_1) + W_{vol}(J)$$

Dieses Modell ist ideal für Elastomere (Gummi), da es die charakteristische S-förmige Spannungs-Dehnungs-Kurve mit Steifigkeitszunahme bei großen Dehnungen (Locking-Effekt) korrekt reproduziert.

### 1.1 Der deviatorische Teil (Taylor-Reihe)

Die inverse Langevin-Funktion, die exakt das Verhalten von Polymerketten beschreibt, ist analytisch schwer zu handhaben. Daher nutzen wir eine Taylor-Entwicklung für $W_{dev}$ um den isotropen Referenzzustand [--> Skript, Gl. 2.242; vgl. Holzapfel, Gl. 6.136]:

$$W_{dev} = \mu \left[ \frac{1}{2}(\bar{I}_1 - 3) + \frac{1}{20\beta^2}(\bar{I}_1^2 - 9) + \frac{11}{1050\beta^4}(\bar{I}_1^3 - 27) + \dots \right]$$

**Physikalische Interpretation:**
- Der erste Term (linear in $\bar{I}_1 - 3$) repräsentiert das Neo-Hooksche Verhalten bei kleinen Dehnungen
- Die höheren Terme (quadratisch, kubisch, etc.) modellieren die zunehmende Steifigkeit bei großen Streckungen
- Der Parameter $\beta$ (Sperrparameter) steuert, bei welcher Streckung der Locking-Effekt dominant wird

**Rechenschritte im Code (`Material_Model.sci`):**

1. **Invariante berechnen:** [--> Holzapfel, Gl. 6.109]
   $$\bar{I}_1 = J^{-2/3} \text{Sp}(\boldsymbol{C}) = \text{Sp}(\bar{\boldsymbol{B}})$$
   Die isochore (volumenerhaltende) erste Invariante beschreibt die reine Gestaltänderung ohne Kompression.

2. **Erste Ableitung bilden ($W_1$):** 
   Um zur Spannung zu gelangen, benötigen wir die Ableitung des Potentials nach der Invariante:
   $$W_1 = \frac{\partial W_{dev}}{\partial \bar{I}_1} = \mu \left( \frac{1}{2} + \frac{1}{10\beta^2}\bar{I}_1 + \frac{33}{1050\beta^4}\bar{I}_1^2 + \dots \right)$$

3. **Taylor-Koeffizienten berechnen:** 
   Die Funktion `calculate_taylor_expansion_coefficients` berechnet die Vorfaktoren $a_1, a_2, a_3, a_4$ aus den Termen der Reihe. Diese werden später für schnelle Auswertung tabelliert.

### 1.2 Der volumetrische Teil

Für die Kompressibilität des Materials verwenden wir einen quadratischen Ansatz, der eine logarithmische Energiefunktion approximiert:

$$W_{vol} = \frac{K}{2}(J-1)^2$$

Die zugehörige hydrostatische Druckreaktion ist:
$$p_J = \frac{\partial W_{vol}}{\partial J} = K(J-1)$$

**Interpretation:**
- Für $J = 1$ (Inkompressibilität) ist $p_J = 0$
- Für $J > 1$ (Expansion) wirkt ein positiver innerer Druck $p_J > 0$
- Der Modul $K$ kontrolliert die Steifigkeit gegen Volumenänderung; typisch wird $K \gg \mu$ gewählt, um nahezu inkompressibles Verhalten zu erzwingen

### 1.3 Der Zwischenschritt: Die 1. Piola-Kirchhoff Spannung ($\boldsymbol{\Pi}$)

Bevor wir die Cauchy-Spannung im aktuellen (verformten) Konfiguration berechnen, ist es oft hilfreich, den Weg über die Nennspannung (1. PK-Spannung) zu gehen. Diese ist definiert als die Ableitung des Potentials nach dem Deformationsgradienten [--> Holzapfel, Gl. 6.1]:

$$\boldsymbol{\Pi} = \frac{\partial W}{\partial \boldsymbol{F}} = \frac{\partial W_{dev}}{\partial \bar{I}_1} \frac{\partial \bar{I}_1}{\partial \boldsymbol{F}} + \frac{\partial W_{vol}}{\partial J} \frac{\partial J}{\partial \boldsymbol{F}}$$

**Wichtige Teil-Ableitungen** (aus Kettenregel):

- **Volumen:** 
  $$\frac{\partial J}{\partial \boldsymbol{F}} = J \boldsymbol{F}^{-T}$$

- **Invariante:** 
  $$\frac{\partial \bar{I}_1}{\partial \boldsymbol{F}} = 2 J^{-2/3} \left( \boldsymbol{F} - \frac{1}{3} I_1 \boldsymbol{F}^{-T} \right)$$

Eingesetzt ergibt sich die **1. PK-Spannung:**

$$\boldsymbol{\Pi} = 2 W_1 J^{-2/3} \left( \boldsymbol{F} - \frac{1}{3} I_1 \boldsymbol{F}^{-T} \right) + p_J J \boldsymbol{F}^{-T}$$

Diese Darstellung ist numerisch stabil und trennt die deviatorischen und volumetrischen Beiträge auf natürliche Weise.

### 1.4 Zusammensetzen der Cauchy-Spannung $\boldsymbol{\Sigma}$

Die wahre Spannung $\boldsymbol{\Sigma}$ (Cauchy-Spannung) erhalten wir durch den **"Push-forward"** der 1. PK-Spannung von der Referenzkonfiguration in die aktuelle Konfiguration:

$$\boldsymbol{\Sigma} = \frac{1}{J} \boldsymbol{F} \boldsymbol{\Pi}$$

Der Faktor $1/J$ berücksichtigt dabei, dass die Spannung nun auf die aktuelle (verformte) Fläche bezogen wird, nicht mehr auf die ursprüngliche Fläche.

Unter Verwendung von $\boldsymbol{B} = \boldsymbol{F}\boldsymbol{F}^T$ (linker Cauchy-Green Tensor) und $\boldsymbol{I} = \boldsymbol{F}^{-T}\boldsymbol{F}^T$ vereinfacht sich dies zur Cauchy-Spannung für das Arruda-Boyce Modell [--> Skript, Gl. 2.243; vgl. Holzapfel, Gl. 6.145]:

$$\boldsymbol{\Sigma} = \frac{2}{J^{5/3}} W_1 \bar{\boldsymbol{B}} + p_J \boldsymbol{I}$$

Im Code wird dies wie folgt umgesetzt:

```scilab
// Deviatorischer Anteil (formänderungsbedingt)
deviatoric_stress = (mu / J^(5/3)) * langevin_expansion * (B - (1/3)*trace(B)*eye(3,3));

// Volumetrischer Anteil (druckbedingt)
volumetric_stress = K * (J - 1) * eye(3,3);

// Gesamtspannung
sigma = deviatoric_stress + volumetric_stress;
```

---

## 2. Aufgabe 2: Spannungs-Dehnungs-Kurven

In dieser Aufgabe erzeugen wir die theoretischen Kurven, die später mit den Messdaten verglichen werden. Obwohl die Materialroutine Cauchy-Spannungen ($\sigma$) liefert, erfordert der Vergleich mit experimentellen Daten eine Umrechnung in technische (Nominal) Spannungen.

### 2.1 Der Datenfluss im Programm

Die Simulation eines Lastfalls (z.B. `simulate_uniaxial_tension`) folgt immer diesem Schema:

1. **Input:** 
   - Ein Vektor von Streckungen `lambda_range` (z.B. von 1.0 bis 7.0)
   - Initiale Materialparameter ($\mu, \beta, K$)

2. **Schleife:** Für jeden Wert $\lambda_i$:
   - Aufstellen des Deformationsgradienten $\boldsymbol{F}(\lambda_i)$ (siehe [Abschnitt 2.4](#24-übersicht-der-lastfälle-deformationsgradienten))
   - Berechnen der Cauchy-Spannung $\boldsymbol{\sigma}_{raw}$ über die Routine aus [Abschnitt 1](#1-aufgabe-1-die-materialroutine-f_to_sig)
   - **Druck-Korrektur:** Da das Modell kompressibel ist, aber der Versuch "frei" (ohne Seitendruck) stattfindet, ziehen wir den Seitendruck ab: 
     $$\Sigma_{true} = \Sigma_{raw} - \Sigma_{22} I$$
   - **Transformation zu PK1:** Umrechnung der wahren Spannung in die Nennspannung (siehe [Abschnitt 2.2](#22-warum-die-rücktransformation-zu-pk1))
   - **Speichern:** Extraktion der Hauptspannung $P_{11}$ in einen Ergebnisvektor

3. **Output:** Ein Vektor technischer Spannungen, der direkt mit Messdaten vergleichbar ist

### 2.2 Warum die Rücktransformation zu PK1?

Experimentelle Daten (wie in `Hyperelastic.txt`) werden fast immer als **Nennspannung** (Nominal Stress / Engineering Stress) aufgezeichnet. Das bedeutet:

$$\text{Nennspannung} = \frac{\text{Aktuelle Kraft}}{\text{Ursprünglicher Querschnitt } A_0}$$

Unsere Materialroutine berechnet jedoch die **wahre Spannung** (Cauchy-Spannung):

$$\sigma = \frac{\text{Aktuelle Kraft}}{\text{Aktueller Querschnitt } A(t)}$$

**Das Problem:** Bei Gummi ändern sich die Querschnitte extrem stark durch Querkontraktion. Eine Zugprobe mit Ausgangsquerschnitt $A_0$ wird beim Dehnen um $\lambda = 7$ in der Querrichtung um den Faktor $\approx 1/\sqrt{\lambda} \approx 0.378$ gestaucht. Der aktuelle Querschnitt ist also deutlich kleiner, und die Cauchy-Spannung kann um ein Mehrfaches höher liegen als die Nennspannung.

**Die Lösung:** Wir transformieren die Cauchy-Spannung über die **1. Piola-Kirchhoff Spannung** (Nennspannung) zurück [--> Holzapfel, Gl. 3.8 für die Piola-Transformation]:

$$\boldsymbol{P} = J \boldsymbol{\sigma} \boldsymbol{F}^{-T}$$

Nur die Komponente $P_{11}$ dieses Tensors entspricht der im Labor gemessenen Kraft pro Ausgangsfläche und ist damit direkt mit experimentellen Daten vergleichbar.

### 2.3 Initiale Werte für den Plot

Um einen ersten Plot zu erstellen (bevor das Fitting läuft), werden Startwerte benötigt. Diese sollten physikalisch sinnvoll sein:

- **Schubmodul $\mu$ [MPa]:** 
  Kann grob aus der Anfangssteigung der Uniaxial-Kurve geschätzt werden. Im elastischen Limit gilt $\sigma \approx 3\mu \cdot (\lambda - 1/\lambda^2)$ bei kleinen Dehnungen, also $E_0 \approx 3\mu$ bei reiner Zugbelastung [--> Skript, Gl. 2.243; vgl. Holzapfel, S. 238, neo-Hookean model]. Typischer Startwert: $\mu_0 \approx 0.3$–$1.0$ MPa.

- **Sperrstreckung $\beta$ (dimensionslos):** 
  Bestimmt den "Aufstieg" der Kurve bei hohen Dehnungen (Locking-Effekt). Ein typischer Startwert ist lt. Skript: $\beta_0 = 10$ (beliebte Wahl). Höhere Werte verschieben das Locking zu größeren Streckungen.

- **Kompressionsmodul $K$ [MPa]:** 
  Wird oft sehr hoch angesetzt (z.B. $K_0 = 1000 \cdot \mu$), um nahezu inkompressibles Verhalten zu erzwingen. Dies ist typisch für Elastomere, die sich kaum volumetrisch ändern.

Ein guter Initial Guess beschleunigt die Konvergenz des Fittings massiv (siehe [Abschnitt 4](#4-aufgabe-4-data-fitting-lsqrsolve)).

### 2.4 Übersicht der Lastfälle (Deformationsgradienten)

| Lastfall | Matrix $\boldsymbol{F}$ | Anwendung | Physik |
| :--- | :--- | :--- | :--- |
| **Uniaxial** | $\text{diag}(\lambda, 1/\sqrt{\lambda}, 1/\sqrt{\lambda})$ | Standard Zugversuch | Einfacher Zug mit freier Querkontraktion [--> Holzapfel, Gl. 2.129] |
| **Biaxial** | $\text{diag}(\lambda, \lambda, \lambda^{-2})$ | Aufblasen einer Membran | Ebene Dehnung in zwei Richtungen [--> Holzapfel, Abs. 2.6 nach Gl. 2.130] |
| **Pure Shear** | $\text{diag}(\lambda, 1, \lambda^{-1})$ | Breiter Streifen | Scherung ohne Volumenänderung [--> Holzapfel, Gl. 2.131] |

Diese drei Lastfälle sind mathematisch unabhängig und bilden zusammen den vollen Parameterraum des Modells ab. Ihre Kombination im Fitting (siehe [Abschnitt 4.3](#43-multimodales-kombiniertes-fitting)) führt zu robusten Materialkennwerten.

---

## 3. Aufgabe 3: Programmierung, Ploterstellung in Scilab

Um die Ergebnisse der Simulation und des Fittings zu visualisieren, nutzen wir die Grafik-Engine von Scilab. Hier ist eine Erklärung der wichtigsten Befehle anhand praktischer Beispiele aus dem Code.

### 3.1 Fenstersteuerung (`scf` & `clf`)

- **`scf(n)` (Set Current Figure):** 
  Dieser Befehl öffnet oder aktiviert ein Grafikfenster mit der Nummer `n`. 
  - *Beispiel:* `scf(0);` stellt sicher, dass alle folgenden Zeichenbefehle in Fenster Nr. 0 landen
  - Ohne diesen Befehl würde Scilab einfach das zuletzt aktive Fenster überschreiben
  - Nützlich, um mehrere Fenster parallel zu verwalten

- **`clf()` (Clear Figure):** 
  Löscht den Inhalt des aktuellen Fensters. Das ist wichtig, wenn man das Skript mehrfach ausführt, damit alte Kurven verschwinden und nicht überzeichnet werden.

### 3.2 Fenster unterteilen (`subplot`)

- **`subplot(m, n, p)`:** 
  Teilt das Fenster in ein Raster aus `m` Zeilen und `n` Spalten auf und aktiviert das `p`-te Teilfenster
  - *Beispiel:* `subplot(1, 3, 1);` erstellt drei Fenster nebeneinander und wählt das linke aus
  - Zählung läuft zeilenweise: `subplot(2,2,1)` ist oben-links, `subplot(2,2,4)` ist unten-rechts
  - Useful für Vergleichsplots (z.B. Initial Guess vs. Fitted vs. Experiment in einer Zeile)

### 3.3 Zeichnen und Beschriften (`plot`, `xtitle`, `legend`)

```scilab
plot((lambda-1)*100, Pi_u_init, "r-", "LineWidth", 1.5); 
xtitle("Initial Guess (Zug/Schub)");
xlabel("$\large \text{techn. Dehnung,}\ \varepsilon_{11}\ [\%]$");
ylabel("Nennspannung $P_{11}$ [MPa]");
legend(["Uniaxial", "Biaxial", "Planar Shear"], 4); 
xgrid();
```

**Erklärung der Befehle:**

- **`plot(x, y, "Stil", ...)`:** 
  - Zeichnet die Daten `y` gegen `x`
  - "r-" steht für rote durchgehende Linie; andere Optionen: "b--" (blau gestrichelt), "g:" (grün gepunktet)
  - "LineWidth", 1.5 setzt die Liniendicke (Standard: 1.0)

- **`xtitle("Titel")`:** 
  Setzt den Haupttitel des (Teil-)Plots

- **`xlabel` / `ylabel`:** 
  - Beschriften die Achsen
  - Unterstützen LaTeX-Code zwischen "$...$" für mathematische Symbole
  - "\large" vergrößert den Text; "\varepsilon" erzeugt das Epsilon-Symbol

- **`legend([...], Position)`:** 
  - Erstellt die Legende
  - Die Zahl gibt die Position im Fenster an: "1" = oben-rechts, "2" = oben-links, "3" = unten-links, "4" = unten-rechts

- **`xgrid()`:** 
  Fügt ein Koordinatengitter zur besseren Lesbarkeit hinzu; `ygrid()` für y-Achse, `grid()` für beide

### 3.4 Fenstergröße anpassen

```scilab
gcf().figure_size = [1300, 590];  // Breite x Höhe in Pixel
```

- **`gcf()` (Get Current Figure):** 
  Greift auf die Eigenschaften des aktuellen Fensters zu
- **`figure_size`:** 
  Setzt Breite und Höhe in Pixeln, damit die Plots in deinen Bericht oder Präsentation passen
- Typische Werte: `[1000, 600]` für breite Layouts, `[800, 800]` für quadratische Vergleiche

---

## 4. Aufgabe 4: Data Fitting (`lsqrsolve`)

Hierbei werden die Parameter $\mu, \beta, K$ so optimiert, dass die Summe der Fehlerquadrate zwischen Modell und Messung minimal wird. Dies ist die Materialcharakterisierung: Wir passen das Modell an experimentelle Daten an.

### 4.1 Der Solver-Befehl: `lsqrsolve`

In Scilab nutzen wir den **Levenberg-Marquardt-Algorithmus** über den Befehl `lsqrsolve`. Die Syntax im Code:

```scilab
[optimized_params, v, info] = lsqrsolve(initial_params, calculate_residuals, num_data_points);
```

**Parameter:**

- **`initial_params`**: 
  Ein Vektor mit den Startwerten für $[\mu, \beta, K]$ (siehe [Abschnitt 2.3](#23-initiale-werte-für-den-plot))
  - Ein guter "Initial Guess" beschleunigt die Konvergenz massiv und vermeidet lokale Minima
  - Beispiel: `initial_params = [0.5; 10; 5000];`

- **`calculate_residuals`**: 
  Der Name der Funktion (unserer "Fehlerfunktion"), die der Solver immer wieder mit neuen Test-Parametern aufruft
  - Diese Funktion muss die Abweichung zwischen Modell und Experiment quantifizieren

- **`num_data_points`**: 
  Die Anzahl der Messpunkte
  - Der Solver muss wissen, wie viele Fehlerwerte er minimieren soll
  - Beispiel: Bei 50 Messwerten je Lastfall und 3 Lastfällen: `num_data_points = 150`

**Rückgabewerte:**

- **`optimized_params`:** 
  Der optimierte Parametervektor $[\mu^*, \beta^*, K^*]$
- **`v`:** 
  Residuen-Vektor (sollte gegen Null konvergieren)
- **`info`:** 
  Konvergenzinformation (0 = erfolgreich, >0 = Warnungen/Fehler)

### 4.2 Die Fehlerfunktion (Residuen)

Der Kern des Fittings ist die Berechnung der **Residuen** (Abweichungen). Innerhalb der Funktion wird für jeden Messpunkt die Differenz zwischen Modellvorhersage und Experiment gebildet:

```scilab
function [residual] = calculate_residuals(p, m)
    mu_p = p(1); 
    beta_p = p(2); 
    K_p = p(3);
    
    // 1. Modellspannung für den aktuellen Parametersatz berechnen
    [model_stresses] = simulate_uniaxial_tension(lambdas, mu_p, beta_p, K_p);
    
    // 2. Differenz zum Experiment bilden
    residual = model_stresses - experimental_data(:,1);
endfunction
```

**Wie der Solver funktioniert:**

1. Der Solver versucht verschiedene Parameterkombinationen
2. Für jede Kombination ruft er `calculate_residuals(p, m)` auf
3. Der Residuenvektor wird berechnet: `residual = model - measurement`
4. Scilab minimiert $\sum_i (\text{residual}_i)^2$ mittels Levenberg-Marquardt
5. Iteration bis Konvergenz oder maximale Iterationen erreicht

Ein guter Residuenvektor sollte:
- Klein sein (Modell passt gut)
- Keine systematischen Trends aufweisen (z.B. nicht immer zu hoch oder zu niedrig)

### 4.3 Multimodales (kombiniertes) Fitting

Ein einzelner Versuch (z.B. nur Uniaxial) reicht oft nicht aus, um alle Parameter eindeutig zu bestimmen. Deshalb nutzen wir das **kombinierte Fitting** mit allen drei Lastfällen (`perform_combined_parameter_fitting`).

Dabei werden die Residuen aller drei Lastfälle einfach untereinander in einen langen Vektor geschrieben ("gestapelt"):

```scilab
function [combined_residual] = perform_combined_parameter_fitting(p, m)
    mu_p = p(1); beta_p = p(2); K_p = p(3);
    
    // Simuliere alle drei Lastfälle mit neuen Parametern
    [sig_u] = simulate_uniaxial_tension(lambdas_u, mu_p, beta_p, K_p);
    [sig_b] = simulate_biaxial_tension(lambdas_b, mu_p, beta_p, K_p);
    [sig_s] = simulate_pure_shear(lambdas_s, mu_p, beta_p, K_p);
    
    // Gestapelte Residuen: Alle drei Lastfälle in einem Vektor
    combined_residual = [sig_u - uni_data(:,1); 
                         sig_b - biax_data(:,1); 
                         sig_s - shear_data(:,1)];
endfunction
```

**Vorteil des kombinierten Fittings:**

Der Solver wird gezwungen, einen Parametersatz zu finden, der für **alle** Belastungsarten gleichzeitig gut funktioniert. Dies macht die Materialkarte deutlich robuster für komplexe 3D-Belastungen im realen Einsatz (z.B. in FEM-Simulationen).

**Mathematischer Hintergrund:**

- Uniaxiale Daten allein können auf 2–3 Parameter passen (degenerate solutions möglich)
- Durch Hinzufügen von Biaxial- und Scherdaten entsteht ein überbestimmtes System
- Das überbestimmte System erzwingt eine *physikalisch realistische* Lösung
- Das Fit-Residuum ist größer, aber die Materialparameter sind verallgemeinerbarer

---

## 5. Aufgabe 5: Tangentialsteifigkeit ($E_T$)

Die Tangentialsteifigkeit ist die **4. Stufe Ableitung der Deformationsenergie**. Sie beschreibt, wie sich die Spannung bei infinitesimaler zusätzlicher Verformung ändert. Dies ist essentiell für numerische Stabilität in FEM-Simulationen.

### 5.1 Die 2. Piola-Kirchhoff Spannung $\boldsymbol{P}$

Zuerst berechnen wir analog zu [Abschnitt 1.3](#13-der-zwischenschritt-die-1-piola-kirchhoff-spannung-boldsymbolπ) [--> Skript, Gl. 2.197; Holzapfel, Gl. 6.13]: 

$$\boldsymbol{P} = 2 \frac{\partial W}{\partial \boldsymbol{C}}$$

Dies führt auf die Zerlegung in deviatorischen und volumetrischen Teil [--> Holzapfel, Gl. 6.88]:

$$\boldsymbol{P} = \underbrace{\frac{2 W_1}{J^{2/3}}  \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)}_{\boldsymbol{P}_{dev}} + \underbrace{K(J-1) J \boldsymbol{C}^{-1}}_{\boldsymbol{P}_{vol}}$$

Die 2. PK-Spannung ist die "Reaktionskraft" des Materials in der **Referenzkonfiguration** (unverformt). Der volumetrische hydrostatische Anteil folgt dabei [--> Holzapfel, Gl. 6.91].

### 5.2 Der Steifigkeitstensor $\mathbb{C}_T = 2 \frac{\partial \boldsymbol{P}}{\partial \boldsymbol{C}}$

Der Steifigkeitstensor 4. Stufe ist die Änderung der Spannung bei Änderung der Metrik (des Deformationsmasses) [--> Skript, Gl. 2.198; Holzapfel, Gl. 6.157]. In `Tangentialsteifigkeit.sci` wird dieser analytisch zusammengesetzt. Wir nutzen die Produktregel:

$$\mathbb{C}_T = \frac{\partial^2 W}{\partial \boldsymbol{C} \partial \boldsymbol{C}}$$

**A) Volumetrischer Teil (`ET_vol`):**

Aus der Ableitung von $P_{vol} = K(J-1) J \boldsymbol{C}^{-1}$ folgt mittels Kettenregel [--> Holzapfel, Gl. 6.166]:

$$\mathbb{C}_{T, vol} = 2 \left[ \frac{\partial (p_J J)}{\partial \boldsymbol{C}} \otimes \boldsymbol{C}^{-1} + (p_J J) \frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}} \right]$$

Im Code wird dies über zwei Terme realisiert:

- **`term_K1`**: $2 K (J - 0.5) J \boldsymbol{C}^{-1} \otimes \boldsymbol{C}^{-1}$ 
  (Resultat aus der Kettenregel für $p_J = K(J-1)$ und $J$)

- **`term_K2`**: $2 p_J J \mathbb{I}_{C^{-1}}$ 
  (Geometrischer Anteil der inversen Metrik)

Dabei ist $\boldsymbol{I}_{C^{-1}}$ der **symmetrische Identitätstensor der inversen Metrik** (`fourth_order_inv_symm`), dessen Ableitung $\frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}} = -\mathbb{I}_{C^{-1}}$ folgt [--> Holzapfel, Gl. 6.164]:

$$\mathbb{I}_{ijkl} = \frac{1}{2} \left( \boldsymbol{C}^{-1}_{ik} \boldsymbol{C}^{-1}_{jl} + \boldsymbol{C}^{-1}_{il} \boldsymbol{C}^{-1}_{jk} \right)$$

Dieser Tensor ist zentral für die korrekte Darstellung geometrischer Nichtlinearitäten.

**B) Deviatorischer Teil (`ET_dev`):**

$$\mathbb{C}_{T, dev} = 4 W_{11} \left( \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} \otimes \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} \right) + 4 W_1 \left( \frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C}^2} \right)$$

Hier sind die **Zwischenschritte** für die 2. Ableitung der Invariante entscheidend (`d2I1bar_dC2`):

$$\frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C}^2} = J^{-2/3} \left[ \text{Coupling} + \text{Product} + \text{Symmetry} \right]$$

Mit den drei Komponenten:

- **`term_vol_coupling`**: $-\frac{1}{3} \left( \boldsymbol{C}^{-1} \otimes \boldsymbol{I} + \boldsymbol{I} \otimes \boldsymbol{C}^{-1} \right)$
  (Kopplung zwischen Volumen und Gestalt)

- **`term_inv_prod`**: $\frac{1}{9} I_1 \boldsymbol{C}^{-1} \otimes \boldsymbol{C}^{-1}$
  (Selbst-Wechselwirkung der inversen Metrik)

- **`term_inv_symm`**: $\frac{1}{3} I_1 \mathbb{I}_{C^{-1}}$
  (Symmetrischer Anteil)

Diese mathematisch exakte Zerlegung stellt sicher, dass der FEM-Solver (oder die 3D-Glyphe in [Abschnitt 6](#6-aufgabe-6-glyphen-visualisierung-e_nnnn)) die korrekte richtungsabhängige Steifigkeit anzeigt. Ohne die 2. Ableitung der Invariante (geometrische Nichtlinearitäten) würde die Steifigkeit bei großen Dehnungen unphysikalisch.

---

## 6. Aufgabe 6: Glyphen-Visualisierung ($E_{nnnn}$)

Um die Steifigkeit im Raum zu verstehen, berechnen wir für jede Raumrichtung $\boldsymbol{n}$ (Einheitsvektor) den sog. **Richtungsmodul** (Directional Stiffness):

$$E_{nnnn} = \sum_{i,j,k,l} n_i n_j \mathbb{E}_{ijkl} n_k n_l$$

Diese Größe gibt an, wie steif das Material in Richtung $\boldsymbol{n}$ gegen Streckung widersteht.

### 6.1 Programmierung der Glyphen (`visualize_stiffness_glyph_3d`)

Die Glyphe visualisiert den Richtungsmodul $E_{nnnn}$ als 3D-Körper. Der Radius der Glyphe in eine bestimmte Richtung entspricht exakt der Steifigkeit in dieser Richtung.

**Der Algorithmus im Code:**

1. **Sphärisches Gitter:** 
   Wir definieren zwei Winkel-Vektoren `phi_grid` (0 bis $2\pi$, Azimut) und `theta_grid` (0 bis $\pi$, Polar)

2. **Richtungsvektor:** 
   Für jede Winkelkombination berechnen wir den Einheitsvektor in sphärischen Koordinaten:
   ```scilab
   nx = sin(theta) * cos(phi); 
   ny = sin(theta) * sin(phi); 
   nz = cos(theta);
   ```

3. **Steifigkeit:** 
   Wir rufen `calculate_directional_stiffness` auf, welche die 4-fach Kontraktion durchführt:
   ```scilab
   radial_stiffness = n(1)*n(1)*ET(1,1,1,1) + n(1)*n(2)*ET(1,1,1,2) + ... // alle 81 Terme
   ```
   In Praxis wird dies mit verschachtelten Schleifen realisiert.

4. **Skalierung:** 
   Die kartesischen Koordinaten für den Plot werden durch Multiplikation des Einheitsvektors mit der berechneten Steifigkeit erzeugt:
   ```scilab
   X = radial_stiffness * nx; // Radius = Steifigkeit
   Y = radial_stiffness * ny;
   Z = radial_stiffness * nz;
   ```

5. **Darstellung:** 
   `surf(X, Y, Z)` zeichnet die Hülle. `gca().isoview = "on"` ist hierbei **essenziell**, damit die Glyphe nicht verzerrt dargestellt wird und die physikalischen Proportionen erhalten bleiben.

### 6.2 Interpretation der Glyphen-Formen

Die Form der Glyphe gibt sofort Auskunft über den Materialzustand:

- **Kugel (Undeformiert, $\lambda=1$):** 
  Das Material ist **isotrop**. In jede Richtung ist der Widerstand gegen Dehnung gleich groß (isotrope Steifigkeit). Der Radius entspricht dem Anfangs-E-Modul.

- **Verzerrte Glyphe (Deformiert, $\lambda > 1$):** 
  Unter unaxialem Zug richten sich die Polymerketten in Zugrichtung aus (strain-induced anisotropy). Die Glyphe nimmt eine elliptische form an:
  
  - **Längsachse:** Der große Radius zeigt die massive Versteifung in Zugrichtung (Locking-Effekt). Das Material wird hier regelrecht hart.
  - **Taille/Äquator:** In Querrichtung ist die Steifigkeit geringer, da hier kaum Kettenausrichtung stattfindet. Das Material bleibt relativ weich.

### 6.3 Die Colorbar-Steuerung (UI-Tricks)

Da die Steifigkeit oft in sehr großen Wertebereichen liegt (oder bei Isotropie fast konstant ist), haben wir die Colorbar manuell "gezähmt":

- **Dynamische Beschriftung (`msprintf`):** 
  Wir nutzen `msprintf`, um Min- und Max-Werte der Steifigkeit direkt in den Titel der Colorbar zu schreiben:
  ```scilab
  c = colorbar();
  c.title.text = msprintf("%s\nMin: %.2f\nMax: %.2f", title, min_val, max_val);
  ```
  Dies macht die Skalierung transparent und dokumentiert die numerischen Werte direkt im Plot.

- **Präzise Positionierung (`axes_bounds`):** 
  Anstatt die Colorbar standardmäßig irgendwo zu platzieren, berechnen wir ihre Position relativ zum aktuellen Subplot:
  ```scilab
  c.position = [gca().axes_bounds(1) + gca().axes_bounds(3) + 0.02, ...];
  ```
  Dies verhindert, dass Colorbars sich überlagern oder Achsenbeschriftungen verdecken.

- **Konstanz-Check:** 
  Wenn das Material fast isotrop ist ($\text{Min} \approx \text{Max}$), erzwingen wir eine kleine künstliche Spanne für die Farbskala, damit Scilab die Colormap nicht mit nur einer Farbe füllt:
  ```scilab
  if (max_val - min_val) < eps then
      max_val = min_val + small_perturbation;
  end
  ```

### 6.4 Kamera und Perspektive

- **`view_angles`**: 
  Wir übergeben der Funktion einen Vektor `[alpha, theta]`, um beide Glyphen (undeformiert und deformiert) aus exakt dem gleichen Winkel zu betrachten. Dies ermöglicht den visuellen Vergleich ohne Verwechslungen durch unterschiedliche Perspektiven.

- **`isoview = "on"`**: 
  Verhindert, dass Scilab die Achsen unterschiedlich skaliert. Nur so bleibt eine Kugel eine Kugel und eine Ellipse eine Ellipse.

---

## 7. Aufgabe 7: Analyse: Material- vs. Räumliche Steifigkeit

In Aufgabe 7 haben wir den aufwendigen Weg gewählt, um die Aufgabenstellung exakt und vollständig abzudecken. Hierbei wird zwischen zwei Sichtweisen der Steifigkeit unterschieden:

### 7.1 Referenz- vs. Momentankonfiguration

Der Professor gibt in der Angabe einen wichtigen Hinweis: *"Oder analog im Deformierten System: C -> B; bzw. P -> SIG"*. Unser Code setzt dies über den **Push-Forward** um:

1. **Materieller Tensor ($\mathbb{C}_T$):** 
   Wir leiten analytisch die Änderung der 2. PK-Spannung (**$\boldsymbol{P}$**) nach dem rechten Cauchy-Green-Tensor ($\boldsymbol{C}$) ab (siehe [Abschnitt 5.2](#52-der-steifigkeitstensor-mathbb-c_t--2-fracpartial-boldsymbbolppartial-boldsymbolc)). Dies ist die **materielle Formulierung**, gültig in der **Referenzkonfiguration** (unverformter Zustand).

2. **Räumlicher Tensor ($\mathbb{c}_T$):** 
   Anstatt die direkte Ableitung $\frac{\partial \Sigma}{\partial B}$ per Hand aufzustellen, transformieren wir $\mathbb{C}_T$ mittels Push-Forward in die **Momentankonfiguration** (aktueller deformter Zustand) [--> Holzapfel, Gl. 6.159]:

   $$c_{ijkl} = \frac{1}{J} F_{iI} F_{jJ} F_{kK} F_{lL} \mathbb{C}_{IJKL}$$
   
   Dieser Tensor $\mathbb{c}_T$ beschreibt die Steifigkeit im deformierten System (Momentankonfiguration) und ist die **räumliche Formulierung**.

**Physikalische Interpretation:**

- Die **materielle Steifigkeit** $\mathbb{C}_T$ beschreibt die Steifigkeit aus der Sicht eines Beobachters, der mit dem Material mitbewegt (Lagrange-Beschreibung)
- Die **räumliche Steifigkeit** $\mathbb{c}_T$ beschreibt die Steifigkeit im Labor-Koordinatensystem (Euler-Beschreibung)
- Für FEM-Formulierungen wird meist die räumliche Steifigkeit verwendet, da die aktuelle Konfiguration relevant ist

### 7.2 Berücksichtigung von Querschnitt und Volumen

Ein kritischer Punkt bei Hyperelastizität ist der Bezug der Kraft:

- **Nennspannung (PK1):** 
  Bezieht die Kraft auf den **ursprünglichen Querschnitt** ($A_0$)
  $$P = \frac{F}{A_0}$$

- **Wahre Spannung (Cauchy):** 
  Bezieht die Kraft auf den **aktuellen, eingeschnürten Querschnitt** ($A_{\text{aktuell}}$)
  $$\sigma = \frac{F}{A_{\text{aktuell}}} = \frac{F}{A_0 / \lambda_{\text{quer}}} = J \cdot P / \lambda$$

Der fundamentale Zusammenhang zwischen der wahren Hauptspannung (Cauchy) und der 1. PK-Hauptspannung folgt dabei [--> Holzapfel, Gl. 6.48]:

$$P_a = J \lambda_a^{-1} \sigma_a$$

Durch den Push-Forward (insbesondere den Faktor $1/J$) wird diese physikalische Realität automatisch korrekt abgebildet. Wir berechnen die Steifigkeit also exakt so, wie sie in einem realen Bauteil unter Last wirkt:

$$\mathbb{c}_T \propto \frac{1}{J} \mathbb{C}_T$$

Der Faktor $1/J$ berücksichtigt, dass der aktuelle Querschnitt kleiner ist und daher die Spannung (und damit auch die Steifigkeit) höher wirkt.

### 7.3 Analytisch vs. Numerisch

Während man die Steifigkeit auch numerisch anhand der Dehnung approximieren könnte (Finite Differenzen: $\Delta P / \Delta C \approx \frac{P(\boldsymbol{C} + \delta) - P(\boldsymbol{C})}{\delta}$), nutzt unser Code die **analytische Ableitung**. Dies bietet drei Vorteile:

1. **Mathematisch exakt:** Keine Approximationsfehler durch Wahl von $\delta$
2. **Schneller:** Analytische Formeln sind schneller auszuwerten als numerische Differenzen
3. **Stabiler:** Keine numerische Instabilität durch Division durch sehr kleine $\delta$-Werte

Für FEM-Implementierungen ist die analytische Steifigkeit essentiell, da Newton-Raphson Solver sie in jedem Iterationsschritt benötigen.

---

## Nomenklatur & Abkürzungen

### Kinematik

| Symbol | Code-Variable | Bedeutung |
| :--- | :--- | :--- |
| $\boldsymbol{F}$ | `F` | Deformationsgradient (Abbildung Referenz $\to$ Momentan) |
| $\boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F}$ | `C` | Rechter Cauchy-Green Tensor |
| $\boldsymbol{B} = \boldsymbol{F}\boldsymbol{F}^T$ | `B` | Linker Cauchy-Green Tensor |
| $\bar{\boldsymbol{B}} = J^{-2/3}\boldsymbol{B}$ | `B_bar` | Isochorer linker Cauchy-Green Tensor |
| $J = \det(\boldsymbol{F})$ | `J` | Volumenverhältnis / Determinante von $\boldsymbol{F}$ |
| $\lambda$ | `lam`, `lambda` | Streckung (aktuelle Länge / Ausgangslänge) |
| $I_1 = \text{Sp}(\boldsymbol{C})$ | `I1` | Erste Invariante von $\boldsymbol{C}$ |
| $\bar{I}_1 = J^{-2/3} I_1$ | `I1bar` | Erste isochore Invariante |
| $\boldsymbol{n}$ | `n_dir`, `n` | Richtungs-Einheitsvektor im Raum |

### Materialparameter

| Symbol | Code-Variable | Bedeutung | Typischer Bereich |
| :--- | :--- | :--- | :--- |
| $\mu$ | `mu`, `MY` | Initialer Schubmodul [MPa] | 0.1–1.0 MPa (Gummi) |
| $\beta$ | `beta`, `BETA` | Sperrparameter (Inverse Langevin Streckgrenze) | 5–20 (dimensionslos) |
| $K$ | `K`, `KMODUL` | Kompressionsmodul [MPa] | 1000–5000 MPa (oder $K = 1000\mu$) |
| $a_1, a_2, a_3, a_4$ | `a` | Taylor-Entwicklungskoeffizienten | — |

### Energie & Spannung

| Symbol | Code-Variable | Bedeutung | Formel / Kontext |
| :--- | :--- | :--- | :--- |
| $W$ | — | Deformationsenergiedichte (Potential) | $W = W_{\text{dev}} + W_{vol}$ |
| $W_1$ | `W1` | 1. Ableitung des Potentials | $W_1 = \frac{\partial W}{\partial \bar{I}_1}$ |
| $W_{11}$ | `W11` | 2. Ableitung des Potentials | $W_{11} = \frac{\partial^2 W}{\partial \bar{I}_1^2}$ |
| $p_J$ | `pressure_J` | Hydrostatischer Druck | $p_J = K(J-1)$ |
| $\boldsymbol{\Sigma}$ (oder $\boldsymbol{\sigma}$) | `sigma`, `sig` | **Cauchy-Spannung** (wahre Spannung) | Kraft / aktuelle Fläche |
| $\boldsymbol{\Pi}$ | `Pi`, `P_nom` | **1. Piola-Kirchhoff Spannung** | Kraft / $A_0$ |
| $\boldsymbol{P}$ | `P_material` | **2. Piola-Kirchhoff Spannung** | Materielle Konfiguration |

### Steifigkeit

| Symbol | Code-Variable | Bedeutung | Dimension |
| :--- | :--- | :--- | :--- |
| $\mathbb{C}_T$ | `ET_tensor` | Tangentiale Steifigkeit (Referenzkonfiguration) | [4. Stufe Tensor] |
| $\mathbb{c}_T$ | `c_spatial` | Räumliche tangentiale Steifigkeit (Momentankonfiguration) | [4. Stufe Tensor] |
| $\mathbb{I}_{C^{-1}}$ | `fourth_order_inv_symm` | Sym. Identitätstensor 4. Stufe der inversen Metrik | [4. Stufe Tensor] |
| $E_{nnnn}$ | `radial_stiffness` | Richtungs-E-Modul (Steifigkeit in Richtung $\boldsymbol{n}$) | [MPa] |
| $E_{T,\text{true}}$ | `ET_true_vec` | Wahrer Tangentenmodul | $d\sigma / d\epsilon_{\text{true}}$ |
| $E_{T,\text{nom}}$ | `ET_tech_vec` | Nominaler Tangentenmodul | $dP / d\lambda$ |
| $\boldsymbol{M}$ | `Voigt_matrix` | Steifigkeitstensor in 6×6 Voigt-Notation | [6×6 Matrix] |

### Abkürzungen

| Abkürzung | Bedeutung | Kontext |
| :--- | :--- | :--- |
| **PK1** | 1. Piola-Kirchhoff | Nennspannung, technisch messbar |
| **PK2** | 2. Piola-Kirchhoff | Materielle Formulierung, FEM |
| **FEM** | Finite-Elemente-Methode | Numerische Simulation |
| **Isotrop** | Richtungsunabhängig | Material ohne Vorzugsrichtung |
| **Anisotrop** | Richtungsabhängig | Hier: Deformationsinduziert durch Kettenausrichtung |
| **Voigt** | Voigt-Notation | Abbildung 4. Stufe Tensor auf 6×6 Matrix |
| **LM** / **L-M** | Levenberg-Marquardt | Optimierungsalgorithmus in `lsqrsolve` |
