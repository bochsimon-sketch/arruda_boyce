# Begleitdokument: Hyperelastizität & Arruda-Boyce Modell
Dieses Dokument erläutert die theoretischen Grundlagen der Scilab-Implementierung des Arruda-Boyce-Materialmodells für hyperelastische Werkstoffe. Es dient als Brücke zwischen der akademischen Herleitung und der numerischen Umsetzung, mit Fokus auf die einzelnen Rechenschritte, die physikalischen Konzepte und deren algorithmische Realisierung in den beiliegenden Skripten.

*Hinweis zur Notation: Die verwendeten kontinuumsmechanischen Herleitungen und Tensordefinitionen folgen primär dem Skript "Solid mechanics of continua" (Reiter, 2025). Das Lehrbuch "Nonlinear Solid Mechanics" (Holzapfel, 2000) dient als ergänzende Quelle für tiefergehende Vergleiche.*

--------------------------------------------------------------------------------

## Inhaltsverzeichnis
1. [Aufgabe 1: Die Materialroutine (`calculate_Cauchy_stress`)](#1-aufgabe-1-die-materialroutine-calculate_Cauchy_stress)
    *  1.1 Der isochore Anteil (Taylor-Reihe)
    *  1.2 Der volumetrische Anteil
    *  1.3 Die 1. Piola-Kirchhoff-Spannung ($\boldsymbol{\Pi}$)
    *  1.4 Push-Forward zur Cauchy-Spannung ($\boldsymbol{\sigma}$)
2. [Aufgabe 2: Lastfall-Simulation & Querdehnung](#2-aufgabe-2-lastfall-simulation-querdehnung)
    *  2.1 Iterative Bestimmung der transversalen Streckung
    *  2.2 Der kompressible Deformationsgradient
    *  2.3 Transformation in messbare Nennspannungen
    *  2.4 Übersicht der Lastfälle
3. [Aufgabe 3: Visualisierung in Scilab](#3-aufgabe-3-visualisierung-in-scilab)
    *  3.1 Fenstersteuerung und Grafik-Engine
    *  3.2 Strukturierte Darstellung via Subplots
    *  3.3 Formatierung und LaTeX-Unterstützung
4. [Aufgabe 4: Parameteridentifikation (`lsqrsolve`)](#4-aufgabe-4-parameteridentifikation-lsqrsolve)
    *  4.1 Der Levenberg-Marquardt-Algorithmus
    *  4.2 Definition der Residuenfunktion
    *  4.3 Multimodales (kombiniertes) Fitting
5. [Aufgabe 5: Tangentialsteifigkeit ($\mathbb{C}_T$) & Voigt-Notation](#5-aufgabe-5-tangentialsteifigkeit-voigt-notation)
    *  5.1 Der materielle Steifigkeitstensor
    *  5.2 Effizienzsteigerung durch 6x6-Reduktion
6. [Glyphen-Visualisierung ($E_{nn}$)](#6-glyphen-visualisierung)
    *  6.1 Vektorisierte Richtungsmodul-Berechnung
    *  6.2 Interpretation der richtungsabhängigen Steifigkeit
    *  6.3 Skalierung und Farbraumsteuerung
7. [Transformation der Steifigkeitstensorik](#7-tranformation-der-steifigkeitstensorik)
    *  7.1 Push-Forward via Voigt-Transformationsmatrix ($\mathbf{Q}$)
    *  7.2 Räumliche vs. Materielle Formulierung
8. [Nomenklatur & Abkürzungen](#8-nomenklatur-abkürzungen)

--------------------------------------------------------------------------------

## 1. Aufgabe 1: Die Materialroutine (`calculate_Cauchy_stress`)
Das Arruda-Boyce-Modell (vgl. Arruda & Boyce, 1993) beschreibt das elastische Potenzial $W$ als Summe aus einem isochoren (gestaltändernden) und einem volumetrischen (volumenändernden) Anteil [--> Skript, Gl. 2.193; vgl. Holzapfel, Gl. 6.85]:

$$W(C) = W_{dev}(\bar{I}_1) + W_{vol}(J)$$

Dieses Modell basiert auf der statistischen Mechanik von Polymerketten und ist ideal für Elastomere (Gummi), da es die charakteristische entropische Versteifung bei großen Dehnungen (Locking-Effekt) präzise reproduziert.

### 1.1 Der isochore Anteil (Taylor-Reihe)
Die exakte mathematische Beschreibung der Kettenstatistik erfolgt über die inverse Langevin-Funktion $\mathcal{L}^{-1}$. Da diese nicht analytisch geschlossen integrierbar ist, nutzt die Implementierung eine Taylor-Entwicklung für das deviatorische Potenzial $W_{dev}$ um den isotropen Referenzzustand [--> Skript, Gl. 2.242; vgl. Holzapfel, Gl. 6.136]:

$$W_{dev} = \mu \left[ \frac{1}{2}(\bar{I}_1 - 3) + \frac{1}{20\beta^2}(\bar{I}_1^2 - 9) + \frac{11}{1050\beta^4}(\bar{I}_1^3 - 27) + \dots \right]$$

**Physikalische Interpretation & Code-Umsetzung (`Material_Model.sci`):**
*   Der erste Term (linear in $\bar{I}_1 - 3$) repräsentiert das Neo-Hookesche Verhalten bei kleinen Dehnungen.
*   Die höheren Terme modellieren die zunehmende Steifigkeit bei großen Streckungen.
*   Der Parameter $\beta = \sqrt{N}$ (Sperrparameter) steuert, bei welcher Streckung der Locking-Effekt dominant wird.
*   Die Funktion `calculate_taylor_expansion_coefficients` berechnet genau diese Vorfaktoren für eine effiziente algorithmische Auswertung.

Dabei ist die modifizierte isochore erste Invariante definiert als [--> Skript, Gl. 2.194]:

$$\bar{I}_1 = J^{-2/3} \text{Sp}(\boldsymbol{C}) = \text{Sp}(\bar{\boldsymbol{B}})$$

Um später die Spannung zu berechnen, benötigen wir die erste Ableitung des isochoren Potentials nach dieser Invariante:

$$W_1 = \frac{\partial W_{dev}}{\partial \bar{I}_1} = \mu \left( \frac{1}{2} + \frac{1}{10\beta^2}\bar{I}_1 + \frac{33}{1050\beta^4}\bar{I}_1^2 + \dots \right)$$

```scilab
// Code-Ausschnitt: Berechnung der modifizierten Invariante I1_bar
invariant_1 = trace(F*F'); // I1 = tr(C)
invariant_3 = det(F);      // J = det(F)
invariant_1_dev = invariant_3^(-2/3)*invariant_1; // I1_bar
```

### 1.2 Der volumetrische Anteil
Um die Kompressibilität des Werkstoffs zu berücksichtigen (Penalty-Ansatz), wird ein quadratischer Ansatz für die volumetrische Energiedichte verwendet:

$$W_{vol} = \frac{K}{2}(J-1)^2$$

Die zugehörige hydrostatische Druckreaktion (im Code `calculate_PK1_vol_from_F`) ist: 

$$p_J = \frac{\partial W_{vol}}{\partial J} = K(J-1)$$

Der Kompressionsmodul $K$ wirkt hierbei als Strafparameter, der bei Elastomeren typischerweise hoch gewählt wird ($K \gg \mu$), um ein nahezu inkompressibles Verhalten zu erzwingen.

### 1.3 Der Zwischenschritt: Die 1. Piola-Kirchhoff-Spannung ($\boldsymbol{\Pi}$)
Bevor die wahre Spannung im aktuellen Zustand berechnet wird, ist es algorithmisch robuster, den Weg über die Nennspannung (1. PK-Spannung) zu gehen. Diese resultiert aus der Ableitung des Gesamtpotentials nach dem Deformationsgradienten $\boldsymbol{F}$ [--> vgl. Skript, Gl. 2.199]:

$$\boldsymbol{\Pi} = \frac{\partial W}{\partial \boldsymbol{F}} = \frac{\partial W_{dev}}{\partial \bar{I}_1} \frac{\partial \bar{I}_1}{\partial \boldsymbol{F}} + \frac{\partial W_{vol}}{\partial J} \frac{\partial J}{\partial \boldsymbol{F}}$$

Unter Verwendung der Kettenregel ergeben sich folgende grundlegende Teil-Ableitungen:
*   **Volumen** [--> Skript, Gl. 2.212]: $\frac{\partial J}{\partial \boldsymbol{F}} = J \boldsymbol{F}^{-T}$
*   **Invariante** [--> Skript, Gl. 2.213]: $\frac{\partial \bar{I}_1}{\partial \boldsymbol{F}} = 2 J^{-2/3} \left( \boldsymbol{F} - \frac{1}{3} I_1 \boldsymbol{F}^{-T} \right)$

Eingesetzt ergibt sich exakt der im Code (`calculate_PK1_stress`) implementierte Ausdruck, der deviatorische und volumetrische Beiträge trennt:

$$\boldsymbol{\Pi} = \underbrace{2 W_1 J^{-2/3} \left( \boldsymbol{F} - \frac{1}{3} I_1 \boldsymbol{F}^{-T} \right)}_{\boldsymbol{\Pi}_{dev}} + \underbrace{K(J-1) J \boldsymbol{F}^{-T}}_{\boldsymbol{\Pi}_{vol}}$$

```scilab
// Berechnung des deviatorischen Anteils der 1. Piola-Kirchhoff Spannung (PK1).
// Resultiert aus der Formänderung des molekularen Netzwerks unter konstantem Volumen.
function PK1_dev = calculate_PK1_dev_from_F(F, invariants, my, Beta)
    I1 = invariants(1); I1_dev = invariants(2); J = invariants(3);
    
    // Taylor-Approximation der Netzwerkversteifung (dW/dI1_dev)
    a = calculate_taylor_expansion_coefficients(Beta);
    d_hyperel_potential_dev = my*(a(1) + a(2)*2.0*I1_dev^1 + a(3)*3.0*I1_dev^2 ..
                                       + a(4)*4.0*I1_dev^3 + a(5)*5.0*I1_dev^4);
    
    // Ableitung der modifizierten Invariante nach F (Kettenregel nach Holzapfel)
    d_invariant_1_dev = 2.0*J^(-2/3)*(F-(1/3)*I1*inv(F)'); 

    // PK1_dev = dW_bar/dI1_dev * dI1_dev/dF
    PK1_dev = d_hyperel_potential_dev * d_invariant_1_dev;
endfunction

// Berechnung des volumetrischen Anteils der PK1-Spannung.
// Modelliert die hydrostatische Reaktion basierend auf dem Kompressionsmodul K.
function PK1_vol = calculate_PK1_vol_from_F(F, J, K)
    // PK1_vol = dW_vol/dJ * dJ/dF = K*(J-1) * J*F^-T
    PK1_vol = K*(J-1.0) * J*inv(F)';
endfunction

// Aggregation der Spannungsbeiträge zur 1. Piola-Kirchhoff Spannung (Nennspannung).
function PK1 = calculate_PK1_stress(F, mat_params)
    my = mat_params(1); Beta = mat_params(2); K = mat_params(3);
    invariants = calculate_invariants_of_F(F);
    
    PK1_dev = calculate_PK1_dev_from_F(F, invariants, my, Beta);
    PK1_vol = calculate_PK1_vol_from_F(F, invariants(3), K);
    PK1     = PK1_dev + PK1_vol;
endfunction
```

### 1.4 Push-Forward zur Cauchy-Spannung ($\boldsymbol{\Sigma}$)
Die wahre Spannung (Cauchy-Spannung) erhalten wir durch den **"Push-Forward"** der 1. PK-Spannung von der Referenzkonfiguration in die aktuelle Konfiguration. Im Skript ist dies über die inverse Piola-Transformation definiert [--> Skript, Gl. 1.149]:

$$\boldsymbol{\Sigma} = \frac{1}{J} \boldsymbol{F} \boldsymbol{\Pi}$$

Der Faktor $1/J$ berücksichtigt dabei, dass die Spannung nun auf die aktuelle (verformte) Fläche bezogen wird, nicht mehr auf den ursprünglichen Querschnitt.

Setzt man die Herleitung für das Arruda-Boyce-Modell konsequent fort (mit dem linken Cauchy-Green-Tensor $\boldsymbol{B} = \boldsymbol{F}\boldsymbol{F}^T$), ergibt sich die kompakte Form der Cauchy-Spannung [--> Skript, Gl. 2.243]:

$$\boldsymbol{\Sigma} = \frac{2}{J^{5/3}} W_1 \bar{\boldsymbol{B}} + p_J \boldsymbol{I}$$

Die Routine `calculate_Cauchy_stress` führt exakt diese Piola-Rücktransformation algorithmisch durch.

```scilab
// Code-Ausschnitt: Push-Forward zur wahren Spannung
function Sigma = calculate_Cauchy_stress(F, mat_params)
    J       = det(F);
    PK1     = calculate_PK1_stress(F, mat_params);
    Sigma   = (1/J) * F * PK1;
endfunction
```

## 2. Aufgabe 2: Lastfall-Simulation & Querdehnung
Ein zentraler Aspekt bei der Simulation experimenteller Lastfälle ist die Berücksichtigung der realen Kinematik. Da Elastomere in unserem Modell kompressibel angesetzt werden ($K < \infty$), muss die transversale Streckung $\lambda_t$ aktiv bestimmt werden, um die physikalischen Randbedingungen der Prüfaufbauten exakt zu erfüllen.

### 2.1 Iterative Bestimmung der transversalen Streckung
In freien Zugversuchen (z. B. uniaxial) sind die Seitenflächen der Probe spannungsfrei, es gilt also für die Normalspannungen in Querrichtung $\Sigma_{22} = \Sigma_{33} = 0$ bzw. $\Pi_{22} = \Pi_{33} = 0$. 
In der Simulation wird dies durch die Funktion `solve_transverse_stretch` realisiert. Anstatt ein perfekt inkompressibles Verhalten ($J=1$) vorauszusetzen, wird für jede vorgegebene Hauptstreckung $\lambda$ die zugehörige Querdehnung $\lambda_t$ iterativ (über ein Newton-Verfahren) so bestimmt, dass das Kraftgleichgewicht in Transversalrichtung exakt erfüllt ist:

$$\Pi_{33}(\lambda, \lambda_t) = 0$$

Dieser Algorithmus stellt sicher, dass der Deformationsgradient $\boldsymbol{F}$ konsistent mit den Materialparametern ($\mu, \beta, K$) berechnet wird.

```scilab
// Code-Ausschnitt: Startwert-Schätzung für die iterative Querdehnungs-Berechnung
function lam_t = solve_transverse_stretch(lam, test_type, mat_params)
    // Startwert aus dem inkompressiblen Grenzfall (J=1)
    select test_type
        case 'uni' then   lam_t = 1/lam^(0.5);
        case 'biax' then  lam_t = 1/(lam^2);
        case 'shear' then lam_t = 1/lam;
    end
    // ... gefolgt von der Newton-Raphson Iteration für Pi_33 = 0
endfunction
```

### 2.2 Der kompressible Deformationsgradient
Durch die explizite numerische Berechnung von $\lambda_t$ erhalten wir einen Deformationsgradienten, der die Volumenänderung $J = \det(\boldsymbol{F})$ korrekt wiedergibt. Für den einachsigen Zug ergibt sich damit:

$$\boldsymbol{F} = \text{diag}(\lambda, \lambda_t, \lambda_t) \quad \text{mit} \quad \lambda_t \neq \frac{1}{\sqrt{\lambda}}$$

Erst im Grenzfall eines unendlich hohen Kompressionsmoduls ($K \to \infty$) nähert sich dieses Ergebnis der klassischen Inkompressibilitätsannahme $\lambda_t = \lambda^{-1/2}$ an.

```scilab
// Code-Ausschnitt: Aufbau des kompressiblen Deformationsgradienten
function F = define_deformation_gradient_for_loadcase(lam, lam_t, test_type)
    select test_type
        case 'uni' then   F = [lam, 0, 0;   0, lam_t, 0;    0, 0, lam_t];
        case 'biax' then  F = [lam, 0, 0;   0,   lam, 0;    0, 0, lam_t];
        case 'shear' then F = [lam, 0, 0;   0,     1, 0;    0, 0, lam_t];
        case 'comp' then  F = [lam, 0, 0;   0,   lam, 0;    0, 0,   lam];
    end
endfunction
```

### 2.3 Transformation in messbare Nennspannungen
Experimentelle Daten (wie in `Hyperelastic.txt`) werden fast immer als **Nennspannung** (Engineering Stress) aufgezeichnet. Das bedeutet:

$$\text{Nennspannung} = \frac{\text{Aktuelle Kraft}}{\text{Ursprünglicher Querschnitt } A_0}$$

**Das Problem:** Unsere Materialroutine berechnet basierend auf der Deformation die **wahre Spannung** (Cauchy-Spannung $\boldsymbol{\Sigma}$). Bei hochelastischen Materialien wie Gummi ändert sich der Querschnitt durch die Querkontraktion beim Ziehen extrem. Eine Zugprobe wird in Querrichtung massiv gestaucht, wodurch der aktuelle Querschnitt deutlich kleiner ist und die Cauchy-Spannung um ein Mehrfaches über der experimentellen Nennspannung liegt.

**Die Lösung:** Wir transformieren die berechnete Cauchy-Spannung über die **Piola-Transformation** in die 1. Piola-Kirchhoff-Spannung (Nennspannung) zurück [--> Skript, Gl. 1.148; vgl. Holzapfel, Gl. 3.8]:
$$\boldsymbol{\Pi} = J \boldsymbol{\Sigma} \boldsymbol{F}^{-T}$$
Nur die Hauptkomponente $\Pi_{11}$ dieses Tensors entspricht der im Labor gemessenen Kraft pro Ausgangsfläche und ist damit direkt mit den experimentellen Daten abgleichbar.

*Tipp für die Startwerte (Initial Guess):* Um die unkalibrierten Modellantworten erstmals plotten zu können, wird im Code ein Startwert für $\mu$ benötigt. Dieser lässt sich grob aus der Anfangssteigung $E_0$ des uniaxialen Zugversuchs schätzen. Im elastischen Limit geht das Modell in ein Neo-Hookesches Material über, für das $\Sigma \approx 3\mu \cdot (\lambda - 1/\lambda^2)$ und somit $E_0 \approx 3\mu$ gilt [--> vgl. Holzapfel, S. 238].

### 2.4 Übersicht der Lastfälle
Die Routine `define_deformation_gradient_for_loadcase` bildet die klassischen Belastungszustände ab. Die Struktur der Deformationsgradienten stützt sich auf die fundamentalen kinematischen Annahmen [--> vgl. Holzapfel, S. 226 f.]:

| Lastfall | Matrix $\boldsymbol{F}$ | Anwendung | Physik / Randbedingung |
| :--- | :--- | :--- | :--- |
| **Uniaxial** | $\text{diag}(\lambda, \lambda_t, \lambda_t)$ | Standard Zugversuch | Einfacher Zug mit freier Querkontraktion, $\Sigma_{22} = \Sigma_{33} = 0$. |
| **Biaxial** | $\text{diag}(\lambda, \lambda, \lambda_t)$ | Aufblasen einer Membran | Ebene (gleiche) Dehnung in zwei Richtungen, transversale Schrumpfung. |
| **Planar Shear** | $\text{diag}(\lambda, 1, \lambda_t)$ | Breiter Streifen | Reine Scherung (Pure Shear). Eine Querrichtung ist blockiert ($\lambda_2=1$). |
| **Compression** | $\text{diag}(\lambda, \lambda, \lambda)$ | Volumetrischer Test | Reine hydrostatische Kompression (volumetrischer Druck). |

--------------------------------------------------------------------------------

## 3. Aufgabe 3: Visualisierung in Scilab
Die grafische Aufbereitung der Simulationsergebnisse erfolgt über die integrierte Grafik-Engine von Scilab. Zur Erstellung publikationsreifer Abbildungen werden spezifische Befehle zur Steuerung der Grafikfenster und der Achsen-Eigenschaften eingesetzt.

### 3.1 Fenstersteuerung (scf & clf)
*   **`scf(n)` (Set Current Figure):** Aktiviert das Grafikfenster mit der Kennung `n`. Falls das Fenster nicht existiert, wird es neu initialisiert. Dies erlaubt die simultane Verwaltung mehrerer Analyseansichten (z. B. Initial Guess in `scf(0)`, validierte Ergebnisse in `scf(1)`).
*   **`clf()` (Clear Figure):** Löscht den aktuellen Inhalt des aktiven Fensters. Dieser Befehl ist essentiell für die iterative Skriptentwicklung, um das Überzeichnen von Kurven aus vorherigen Simulationsläufen zu vermeiden.

### 3.2 Strukturierte Darstellung via Subplots
Die Funktion `subplot(m, n, p)` unterteilt die Grafikfläche in eine Matrix aus `m` Zeilen und `n` Spalten. Der Parameter `p` adressiert das jeweilige Teilfenster (z. B. `subplot(1, 2, 1)` für die linke und `subplot(1, 2, 2)` für die rechte Spalte). Dies ermöglicht den direkten visuellen Vergleich unterschiedlicher Lastfälle oder den Abgleich zwischen Modellvorhersage und experimentellen Daten in einer einzigen kompakten Abbildung.

### 3.3 Formatierung und LaTeX-Unterstützung
Zur präzisen Beschriftung werden Funktionen wie `title`, `xlabel` und `ylabel` verwendet. Scilab unterstützt hierbei die native Einbindung von LaTeX-Strings (eingeschlossen in `$ ... $`), was die korrekte Darstellung mathematischer Symbole wie $\varepsilon_{11}$ oder $\Pi_{11}$ ermöglicht. Durch `gcf().figure_size` wird die Ausgabegröße der Abbildungen in Pixeln exakt definiert, um eine konsistente Dokumentationsqualität für Berichte zu gewährleisten.

```scilab
// Code-Ausschnitt: Ploterstellung mit LaTeX-Beschriftung und Fenstersteuerung
scf(0); clf(); // Aktiviere und leere Fenster 0
gcf().figure_size = ; // Definiere feste Ausgabegröße

subplot(1,2,1);
plot((lambda_coarse-1)*100, Pi_u_init, "r-", "LineWidth", 1.5);
plot((lambda_coarse-1)*100, Pi_b_init, "b-", "LineWidth", 1.5);
plot((lambda_coarse-1)*100, Pi_s_init, "g-", "LineWidth", 1.5);

// LaTeX-Beschriftung der Achsen
xlabel(" $\large \text{technische Dehnung,}\ \Large \varepsilon_{11}\ \normalsize\left[\%\right]$ ");
ylabel(" $\large \text{Nennspannung,}\  \Pi_{11}\ \normalsize\left[MPa\right]$ ");
xgrid(1,1,7); // Hilfsgitter zur besseren Lesbarkeit einblenden
```

---

## 4. Aufgabe 4: Parameteridentifikation (`lsqrsolve`)
Die Bestimmung der optimalen Materialparameter ($\mu, \beta, K$) erfolgt durch die Minimierung der Fehlerquadrate zwischen der theoretischen Modellvorhersage und den experimentellen Messwerten. Dies ist der Kern der Materialcharakterisierung.

### 4.1 Der Levenberg-Marquardt-Algorithmus
In Scilab nutzen wir für diese nichtlineare Optimierung die Funktion `lsqrsolve`, welche den Levenberg-Marquardt-Algorithmus implementiert. Dieser Algorithmus kombiniert die Robustheit des Gradientenverfahrens (für weit entfernte Startwerte) mit der Konvergenzgeschwindigkeit des Gauß-Newton-Verfahrens (in der Nähe des Minimums). Ein physikalisch sinnvoller Startvektor (Initial Guess, siehe Kap. 2.3) ist dennoch essenziell, um lokale Minima zu vermeiden.

### 4.2 Definition der Residuenfunktion
Die Qualität des Fittings wird durch Residuenfunktionen bestimmt. Diese berechnen für einen gegebenen, iterierten Parametersatz den Differenzvektor zwischen den simulierten Spannungen und den experimentellen Daten:
$$r_i = \Pi_{11, \text{Modell}}(\lambda_i) - \Pi_{11, \text{Experiment}, i}$$
Der Solver variiert die Materialparameter iterativ so lange, bis die Norm des Residuenvektors ($\sum r_i^2$) ein definiertes Minimum erreicht.

### 4.3 Robuste, mehrstufige Fitting-Strategie
Ein einzelner Versuch (z.B. nur Uniaxial) reicht bei hyperelastischen Modellen oft nicht aus, um alle Parameter eindeutig und physikalisch sinnvoll zu bestimmen. Daher implementiert der Code in `Main.sci` und `Curve_Fit.sci` eine robuste 2-stufige Kalibrierungsstrategie:

*   **Schritt A (Volumetrisches Fitting):** Zunächst wird der Kompressionsmodul $K$ isoliert an den Daten des volumetrischen Druckversuchs (`Compression_Test.txt`) kalibriert (`compression_fitting`). Da $K$ bei Elastomeren meist um Potenzen größer ist als der Schubmodul $\mu$, würde ein simultanes Fitten aller Parameter zu numerischen Instabilitäten führen.
*   **Schritt B (Kombiniertes deviatorisches Fitting):** Mit dem nun fixierten $K$-Modul werden die deviatorischen Parameter ($\mu, \beta$) simultan über die drei restlichen Lastfälle (Uniaxial, Biaxial, Planar Shear) optimiert (`combined_deviatoric_fitting`). Die Residuenvektoren aller drei Lastfälle werden hierbei zu einem globalen Vektor aggregiert. Dies zwingt den Algorithmus, eine physikalisch generalisierbare "Materialkarte" zu finden, die den gesamten Deformationsraum bestmöglich abbildet.

```scilab
// Code-Ausschnitt: 2-Stufige Fitting-Strategie (aus Main.sci)

// Schritt A: Volumetrisches Fitting zur isolierten Bestimmung von K
p_comp_only = compression_fitting(p_init, Comp_data);
K_fitted    = p_comp_only(3);

// Schritt B: Kombiniertes Fitting (mu & beta) bei fixiertem K-Modul
p_comb_init   = [my_init, Beta_init];
p_comb_fitted = combined_deviatoric_fitting(p_comb_init, Uni_data, Biax_data, Shear_data, K_fitted);
```

--------------------------------------------------------------------------------

## 5. Aufgabe 5: Tangentialsteifigkeit ($\mathbb{C}_T$) & Voigt-Notation
Die Tangentialsteifigkeit (der materielle Steifigkeitstensor 4. Stufe) ist die zweite Ableitung der Deformationsenergiedichte $W$ nach dem Deformationsmaß. Sie beschreibt das inkrementelle Materialverhalten und ist im Rahmen der Finite-Elemente-Methode (FEM) essenziell für die numerische Stabilität und die quadratische Konvergenz des globalen Newton-Raphson-Solvers.

### 5.1 Die 2. Piola-Kirchhoff-Spannung ($\boldsymbol{P}$)
Der Ausgangspunkt für den materiellen Steifigkeitstensor ist die 2. Piola-Kirchhoff-Spannung $\boldsymbol{P}$. Diese "Reaktionskraft" des Materials in der Referenzkonfiguration ergibt sich aus der ersten Ableitung des Potentials nach dem rechten Cauchy-Green-Tensor $\boldsymbol{C}$ [--> Skript, Gl. 2.197; vgl. Holzapfel, Gl. 6.13]:

$$\boldsymbol{P} = 2 \frac{\partial W}{\partial \boldsymbol{C}}$$

Aufgrund der additiven Zerlegung des Arruda-Boyce-Potentials führt dies direkt auf eine Zerlegung in einen deviatorischen und einen volumetrischen Spannungsteil [--> vgl. Holzapfel, Gl. 6.88 u. 6.91]:

$$\boldsymbol{P} = \underbrace{\frac{2 W_1}{J^{2/3}} \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)}_{\boldsymbol{P}_{dev}} + \underbrace{K(J-1) J \boldsymbol{C}^{-1}}_{\boldsymbol{P}_{vol}}$$

### 5.2 Der materielle Steifigkeitstensor ($\mathbb{C}_T$)
Der Steifigkeitstensor 4. Stufe folgt aus der erneuten Ableitung der Spannung nach $\boldsymbol{C}$ [--> Skript, Gl. 2.198; vgl. Holzapfel, Gl. 6.157]:

$$\mathbb{C}_T = 2 \frac{\partial P}{\partial C} = 4 \frac{\partial^2 W}{\partial C \partial C}$$

**A) Volumetrischer Teil:**
Aus der Ableitung von $P_{vol} = K(J-1) J \boldsymbol{C}^{-1}$ folgt mittels Produkt- und Kettenregel der volumetrische Steifigkeitstensor [--> vgl. Holzapfel, Gl. 6.166]:

$$\mathbb{C}_{T, vol} = 2 \left[ \frac{\partial (p_J J)}{\partial C} \otimes C^{-1} + (p_J J) \frac{\partial C^{-1}}{\partial C} \right]$$

Dabei ist die Ableitung der inversen Metrik negativ und führt auf den symmetrischen Identitätstensor 4. Stufe ($I_{C^{-1}}$) [--> vgl. Holzapfel, Gl. 6.164]: $\frac{\partial C^{-1}}{\partial C} = -I_{C^{-1}}$. Im Code (`Tangentialsteifigkeit.sci`) spiegelt sich dies in den Termen für $K$ wider, welche exakt diese geometrische Nichtlinearität der inversen Metrik abbilden.

**B) Deviatorischer (isochorer) Teil:**
Die analytische Ableitung des isochoren Anteils erfordert die Anwendung der Kettenregel auf die modifizierte Invariante $\bar{I}_1$:

$$\mathbb{C}_{T, dev} = 4 W_{11} \left( \frac{\partial \bar{I}_1}{\partial C} \otimes \frac{\partial \bar{I}_1}{\partial C} \right) + 4 W_1 \left( \frac{\partial^2 \bar{I}_1}{\partial C^2} \right)$$

Die korrekte algorithmische Implementierung der zweiten Ableitung $\frac{\partial^2 I_1}{\partial C^2}$ ist kompliziert und besteht aus der Kopplung zwischen Volumen und Gestaltänderung, der Selbst-Wechselwirkung der inversen Metrik sowie dem symmetrischen Anteil $\mathbb{I}_{C^{-1}}$. Diese mathematisch exakte Zerlegung stellt sicher, dass die richtungsabhängige Steifigkeit bei finiten Dehnungen physikalisch konsistent bleibt.

### 5.3 Effizienzsteigerung durch 6x6-Voigt-Notation
Ein Tensor 4. Stufe besitzt in 3D grundsätzlich $3^4 = 81$ Komponenten. Da sowohl der Spannungs- als auch der Verzerrungstensor symmetrisch sind, weist auch der Steifigkeitstensor $\mathbb{C}_T$ weitreichende Symmetrien auf (Minor- und Major-Symmetries).
Um die enorme Rechenlast bei Tensoroperationen 4. Stufe in Scilab zu umgehen, nutzt die Funktion `calculate_tangent_stiffness_tensor` die **Voigt-Notation**. Dabei wird der Tensor verlustfrei auf eine $6 \times 6$-Matrix abgebildet. Dies beschleunigt die nachfolgenden Push-Forward-Transformationen und die Visualisierung der 3D-Glyphen massiv.

```scilab
// Code-Ausschnitt: Mapping des 4. Stufe Tensors auf eine 6x6 Matrix
function [Voigt_matrix] = convert_tensor_to_voigt_6x6(ET_tensor)
    Voigt_matrix = zeros(6,6);
    // Voigt-Mapping: 11->1, 22->2, 33->3, 23->4, 13->5, 12->6
    voigt_map = [1,1; 2,2; 3,3; 2,3; 1,3; 1,2]; 
    
    for I=1:6,
        for J=1:6
            Voigt_matrix(I,J) = ET_tensor(voigt_map(I,1), voigt_map(I,2), voigt_map(J,1), voigt_map(J,2));
        end
    end
endfunction
```

## 6. Aufgabe 6: Glyphen-Visualisierung ($E_{nn}$)
Zur tiefgehenden Analyse der richtungsabhängigen Materialsteifigkeit wird der Richtungsmodul $E_{nnnn}$ berechnet. Dieser skalare Wert beschreibt den Widerstand des Materials gegen eine inkrementelle Streckung in einer beliebigen Raumrichtung $\boldsymbol{n}$.

### 6.1 Vektorisierte Richtungsmodul-Berechnung
Die aktuelle Implementierung (`visualize_stiffness_glyph_3d`) nutzt eine hochgradig effiziente, vektorisierte Berechnungsmethode. Für jede Raumrichtung $\boldsymbol{n}$ wird ein Voigt-Projektionsvektor $\mathbf{v}_n$ konstruiert:

$$\mathbf{v}_n = [n_1^2, n_2^2, n_3^2, 2n_2n_3, 2n_1n_3, 2n_1n_2]^T$$

Der Richtungsmodul ergibt sich dann schlicht als quadratische Form der $6 \times 6$-Voigt-Matrix:

$$E_{nnnn} = \mathbf{v}_n^T \cdot \mathbf{M}_{6 \times 6} \cdot \mathbf{v}_n$$

Dieser Ansatz erlaubt es, die Steifigkeit in Scilab für ein gesamtes sphärisches Gitter simultan zu berechnen, ohne auf rechenintensive Schleifen zurückgreifen zu müssen.

```scilab
// Code-Ausschnitt: Berechnung des Richtungsmoduls via Voigt-Matrix
function E_direction = calculate_directional_stiffness_voigt(Voigt_matrix, n_vector)
    n = n_vector ./ norm(n_vector);
    vn = [n(1)^2; n(2)^2; n(3)^2; 2*n(2)*n(3); 2*n(1)*n(3); 2*n(1)*n(2)];
    E_direction = vn' * Voigt_matrix * vn;
endfunction
```

### 6.2 Interpretation der richtungsabhängigen Steifigkeit
Die Form der resultierenden 3D-Glyphe gibt sofort Auskunft über den physikalischen Materialzustand:
*   **Kugel (Undeformiert, $\lambda=1$):** Das Material ist initial isotrop. In jede Richtung ist der Widerstand gegen Dehnung gleich groß. Der Radius der Kugel entspricht dem initialen Tangentenmodul.
*   **Verzerrtes Ellipsoid (Deformiert, $\lambda > 1$):** Unter unaxialem Zug richten sich die Polymerketten in Zugrichtung aus (Dehnungsinduzierte Anisotropie). 
    *   *Längsachse:* Der massiv anwachsende Radius visualisiert die entropische Versteifung in Zugrichtung (Locking-Effekt). Das Material wird hier regelrecht hart.
    *   *Taille / Äquator:* In Querrichtung ist die Steifigkeit deutlich geringer, da hier die Kettenausrichtung fehlt und das Material quer zur Zugrichtung weich bleibt.

### 6.3 Scilab Grafik-Tools: Farbraum, Kamera und Perspektive
Um die 3D-Glyphen publikationsreif und unverzerrt darzustellen, werden spezifische Post-Processing-Tools von Scilab verwendet:
*   **`surf(X,Y,Z, color)`:** Zeichnet die 3D-Oberfläche. Durch die Übergabe der radialen Steifigkeit als 4. Parameter (`r_stiffness`) wird die Farbe an die lokale Steifigkeit gekoppelt.
*   **`gcf().color_map = jet(100)` & `colorbar`:** Definiert einen Farbverlauf von Blau (weich) nach Rot (steif) in 100 Abstufungen und blendet eine Farbskala zur quantitativen Ablesung ein.
*   **`gca().isoview = "on"`:** Dies ist **essenziell**! Es verhindert, dass Scilab die Achsen automatisch unterschiedlich skaliert, was die physikalischen Proportionen der Glyphe (Kugel vs. Ellipse) visuell zerstören würde.
*   **`gca().rotation_angles = my_view`:** Setzt den Kamerawinkel (Azimut und Polar) für beide Subplots exakt gleich, um einen unverfälschten Vorher-Nachher-Vergleich zu garantieren.

--------------------------------------------------------------------------------

## 7. Aufgabe 7: Transformation der Steifigkeitstensorik
In der Kontinuumsmechanik muss streng zwischen der Referenzkonfiguration (materiell) und der Momentankonfiguration (räumlich) unterschieden werden.

### 7.1 Push-Forward via Voigt-Transformationsmatrix ($\mathbf{Q}$)
Der räumliche Steifigkeitstensor $\mathbb{c}_T$ beschreibt die Steifigkeit im deformierten Labor-Koordinatensystem. Er wird über die Push-Forward-Operation des materiellen Tensors $\mathbb{C}_T$ berechnet [--> vgl. Holzapfel, Gl. 6.159]:

$$c_{ijkl} = \frac{1}{J} F_{iI} F_{jJ} F_{kK} F_{lL} \mathbb{C}_{IJKL}$$

Um diese rechenintensive 4.-Stufe-Operation in Schleifen zu umgehen, nutzt der Code eine $6 \times 6$ Voigt-Transformationsmatrix $\mathbf{Q}$. Diese bildet die Tensor-Transformation exakt auf den Voigt-Raum ab:

```scilab
// Code-Ausschnitt: Hocheffizienter Push-Forward der Steifigkeit im Voigt-Raum
function Voigt_c = push_forward_stiffness_voigt(Voigt_C, F)
    J = det(F);
    Q = get_voigt_transformation_matrix(F);
    Voigt_c = (1/J) * Q * Voigt_C * Q'; // M_räumlich = (1/J) * Q * M_materiell * Q^T
endfunction
```

### 7.2 Räumliche vs. Materielle Formulierung (Der 1D-Fall)
Da das Materialmodell kompressibel angesetzt ist ($K < \infty$) und die tatsächliche Querkontraktion numerisch exakt über das transversale Kraftgleichgewicht bestimmt wird, verbietet sich eine vereinfachte analytische 1D-Auswertung der Steifigkeit (wie sie bei perfekter Inkompressibilität möglich wäre).
Stattdessen extrahiert die Routine `calculate_uniaxial_stiffness_curves` die Steifigkeitsentwicklung direkt aus den physikalischen Tensorkomponenten:
*   **Materielle Steifigkeit ($\mathbb{C}_{11}$):** Bezogen auf die Referenzkonfiguration. Zeigt das fundamentale Versteifen des Materialnetzwerks durch Kettenstreckung.
*   **Räumliche Steifigkeit ($\mathbb{c}_{11}$):** Bezogen auf die Momentankonfiguration. Hier ist der "geometrische Erweichungseffekt" durch die schrumpfende Querschnittsfläche und die Volumenänderung ($1/J$) bereits vollständig und exakt einkondensiert.

```scilab
// Code-Ausschnitt: Extraktion der Tensorkomponenten über den Dehnungspfad
C_comps.C11(i) = Voigt_C_mat(1,1); // Referenzkonfiguration (materiell)
C_comps.C22(i) = Voigt_C_mat(2,2); 

c_comps.c11(i) = Voigt_c_mat(1,1); // Momentankonfiguration (räumlich)
c_comps.c22(i) = Voigt_c_mat(2,2); 
```

--------------------------------------------------------------------------------

## 8. Nomenklatur & Abkürzungen

### Kinematik & Tensorik
| Symbol | Code-Variable | Bedeutung |
| :--- | :--- | :--- |
| $\boldsymbol{F}$ | `F` | Deformationsgradient (Abbildung Referenz $\to$ Momentan) |
| $\boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F}$ | `C` | Rechter Cauchy-Green Tensor |
| $\boldsymbol{B} = \boldsymbol{F}\boldsymbol{F}^T$ | `B` | Linker Cauchy-Green Tensor |
| $\bar{\boldsymbol{B}} = J^{-2/3}\boldsymbol{B}$ | `B_bar` | Isochorer linker Cauchy-Green Tensor |
| $J = \det(\boldsymbol{F})$ | `J` | Volumenverhältnis / Determinante von $\boldsymbol{F}$ |
| $\lambda, \lambda_t$ | `lam`, `lam_t` | Längsstreckung / transversale Streckung |
| $I_1 = \text{Sp}(\boldsymbol{C})$ | `I1` | Erste Invariante von $\boldsymbol{C}$ |
| $\bar{I}_1 = J^{-2/3} I_1$ | `I1_dev` | Erste modifizierte isochore Invariante |
| $\boldsymbol{n}$ | `n_vector`, `n` | Richtungs-Einheitsvektor im Raum |

### Materialparameter & Energie
| Symbol | Code-Variable | Bedeutung |
| :--- | :--- | :--- |
| $\mu$ | `my`, `mu` | Initialer Schubmodul im Referenzzustand |
| $\beta$ | `Beta` | Sperrparameter (Inverse Langevin Kettensteifigkeit) |
| $K$ | `K` | Kompressionsmodul (Bulk-Modulus) |
| $W$ | — | Deformationsenergiedichte (Potential $W = W_{\text{dev}} + W_{vol}$) |
| $W_1$ | `W1` | Erste Ableitung des isochoren Potentials nach $\bar{I}_1$ |
| $p_J$ | `p_J` | Hydrostatischer Reaktionsdruck aus Kompressibilität |

### Spannungen & Steifigkeit
| Symbol | Code-Variable | Bedeutung |
| :--- | :--- | :--- |
| $\boldsymbol{\sigma}$ | `Sigma`, `sig` | Cauchy-Spannung (Wahre Spannung) |
| $\boldsymbol{\Pi}$ | `PK1`, `Pi` | 1. Piola-Kirchhoff-Spannung (Nennspannung) |
| $\boldsymbol{P}$ | `P_material` | 2. Piola-Kirchhoff-Spannung |
| $\mathbb{C}_T$ | `ET_tensor` | Materielle tangentiale Steifigkeit (Referenzzustand) |
| $\mathbb{c}_T$ | `c_spatial` | Räumliche tangentiale Steifigkeit (Momentanzustand) |
| $\mathbf{M}_{6 \times 6}$ | `Voigt_C`, `Voigt_c`| Steifigkeitstensor in $6 \times 6$ Voigt-Notation |
| $\mathbf{Q}$ | `Q` | $6 \times 6$ Voigt-Transformationsmatrix |
| $E_{nnnn}$ | `radial_stiffness`| Vektorisierter Richtungsmodul |

### Abkürzungen
| Abkürzung | Bedeutung | Kontext |
| :--- | :--- | :--- |
| **PK1 / PK2** | 1. bzw. 2. Piola-Kirchhoff | Nennspannung vs. Materielle Formulierung |
| **FEM** | Finite-Elemente-Methode | Numerische Integrationsumgebung |
| **Voigt** | Voigt-Notation | Effiziente Reduktion von Tensoren 4. Stufe |
| **Locking** | Locking-Effekt | Steiler Anstieg der Kraft durch Kettenstreckung |
| **Isochor** | Volumenerhaltend | $J=1$, reine Gestaltänderung (Deviator) |
| **LM-Solver** | Levenberg-Marquardt | Optimierungsalgorithmus (`lsqrsolve`) |
