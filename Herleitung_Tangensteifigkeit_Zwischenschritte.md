# Herleitung der 1. und 2. Ableitung für das Arruda-Boyce Materialmodell

## 1. Ausgangssituation: Das elastische Potential $W$

Das Arruda-Boyce-Modell (Skript Gl. 2.242) teilt sich in einen deviatorischen (volumenerhaltenden) und einen volumetrischen Teil auf:

$$
W(\boldsymbol{C}) = W_{dev}(\bar{I}_1) + W_{vol}(J)
$$

Die Definitionen der Invarianten basieren auf dem rechten Cauchy-Green-Tensor $\boldsymbol{C}$:
* **Volumenänderung:** $J = \sqrt{\det(\boldsymbol{C})}$
* **Modifizierte 1. Invariante:** $\bar{I}_1 = J^{-2/3} I_1 = J^{-2/3} \text{Sp}(\boldsymbol{C})$

Das Potential für Arruda-Boyce lautet:

$$
W_{dev} = \mu \left( \frac{1}{2}(\bar{I}_1 - 3) + \frac{1}{20\beta^2}(\bar{I}_1^2 - 9) + \frac{11}{1050\beta^4}(\bar{I}_1^3 - 27) + \dots \right)
$$

$$
W_{vol} = \frac{K}{2}(J-1)^2
$$

---

## 2. Hilfsmittel: Tensor-Ableitungen der Invarianten
Um die Kettenregel anzuwenden, brauchen wir die Ableitungen der Invarianten nach $\boldsymbol{C}$. 

$$\frac{\partial J}{\partial \boldsymbol{C}} = \frac{1}{2}(\det\boldsymbol{C})^{-1/2} \cdot \left( J^2 \boldsymbol{C}^{-1} \right) = \frac{1}{2J} J^2 \boldsymbol{C}^{-1} = \frac{1}{2} J \boldsymbol{C}^{-1}$$

$$
\frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} = \frac{\partial (J^{-2/3})}{\partial \boldsymbol{C}} I_1 + J^{-2/3} \frac{\partial I_1}{\partial \boldsymbol{C}} 
$$

Herleitung für $$\frac{\partial \bar{I}_1}{\partial \boldsymbol{C}}$$:

$$\bar{I}_1 = J^{-2/3} I_1$$
$$\frac{\partial I_1}{\partial \boldsymbol{C}} = \boldsymbol{I}$$

Kettenregel anwenden: da $({J(\boldsymbol{C})})^{-2/3}$ 

Innere mal äußere Ableitung:

$$\frac{\partial (J^{-2/3})}{\partial \boldsymbol{C}} = -\frac{2}{3} J^{-5/3} \cdot \frac{\partial J}{\partial \boldsymbol{C}}$$

ergibt:

$$
\frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} = \frac{\partial (J^{-2/3})}{\partial \boldsymbol{C}} I_1 + J^{-2/3} \frac{\partial I_1}{\partial \boldsymbol{C}} = J^{-2/3} \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)
$$

## 3. die 2. Piola-Kirchhoff Spannung $\boldsymbol{P}$

Nach Skript (Gl. 2.197) ist die 2. PK-Spannung definiert als:

$$
\boldsymbol{P} = 2 \frac{\partial W}{\partial \boldsymbol{C}}
$$

Aufgeteilt in deviatorische und volumetrische Anteile des Potentials:

$$
\boldsymbol{P} = 2 \left( \frac{\partial W_{dev}}{\partial \bar{I}_1} \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} + \frac{\partial W_{vol}}{\partial J} \frac{\partial J}{\partial \boldsymbol{C}} \right)
$$

Wir definieren nun die Ableitungen der Arruda-Boyce-Energiefunktion:

$$
W_1 := \frac{\partial W_{dev}}{\partial \bar{I}_1} = \mu \left( \frac{1}{2} + \frac{2}{20\beta^2}\bar{I}_1 + \frac{33}{1050\beta^4}\bar{I}_1^2 + \dots \right)
$$

$$
p_J := \frac{\partial W_{vol}}{\partial J} = K(J-1)
$$

Eingesetzt in die Spannungsgleichung ergibt das:

$$
\boldsymbol{P} = \frac{2 \cdot W_1}{J^{2/3}}  \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right) + p_J \cdot J \boldsymbol{C}^{-1}
$$

---

# Ableitung der PK2-Spannung nach dem Cauchy-Green-Tensor $\boldsymbol{C}$

Unsere Ausgangsgleichung lautet:

$$
\boldsymbol{P} = \boldsymbol{P}_{dev} + \boldsymbol{P}_{vol}
$$

Mit den beiden Anteilen:

$$
\boldsymbol{P}_{dev} = \frac{2 W_1}{J^{2/3}} \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)
$$

$$
\boldsymbol{P}_{vol} = p_J J \boldsymbol{C}^{-1}
$$

Wir suchen die Ableitung $\mathbb{C}_T = 2 \frac{\partial \boldsymbol{P}}{\partial \boldsymbol{B}}$.

---

## 1. Der volumetrische Teil (Produktregel)

Wir müssen den Term $\boldsymbol{P}_{vol} = (p_J \cdot J) \boldsymbol{C}^{-1}$ nach $\boldsymbol{C}$ ableiten. Da hier ein Skalar $(p_J \cdot J)$ mit einem Tensor $\boldsymbol{C}^{-1}$ multipliziert wird, greift die Produktregel: $(f \cdot \boldsymbol{T})' = f' \otimes \boldsymbol{T} + f \cdot \boldsymbol{T}'$.

**Schritt 1.1: Ableitung des Skalars $(p_J \cdot J)$**
Wir nutzen die Produktregel für Skalare und die Kettenregel mit $\frac{\partial J}{\partial \boldsymbol{C}} = \frac{1}{2}J\boldsymbol{C}^{-1}$:

$$
\frac{\partial (p_J J)}{\partial \boldsymbol{C}} = \left( \frac{\partial p_J}{\partial J} J + p_J \right) \frac{\partial J}{\partial \boldsymbol{C}}
$$

Mit $W_{JJ} := \frac{\partial p_J}{\partial J} = K$ und $p_J = K(J-1)$ wird das zu:

$$
\frac{\partial (p_J J)}{\partial \boldsymbol{C}} = \left( K \cdot J + K(J-1) \right) \frac{1}{2} J \boldsymbol{C}^{-1} = K \left( J - \frac{1}{2} \right) J \boldsymbol{C}^{-1}
$$

**Schritt 1.2: Ableitung des Tensors $\boldsymbol{C}^{-1}$**
Die Ableitung eines inversen Tensors nach sich selbst ergibt immer den negativen Identitätstensor 4. Stufe $\mathbb{I}_{C^{-1}}$:

$$
\frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}} = -\mathbb{I}_{C^{-1}} =-\frac{1}{2} \left( C_{ik}^{-1} C_{jl}^{-1} + C_{il}^{-1} C_{jk}^{-1} \right)
$$

**Schritt 1.3: Volumetrischen Teil zusammensetzen**
Wir setzen beides in unsere Produktregel ein und multiplizieren mit 2 für den Steifigkeitstensor $C_{T, vol} = 2 \frac{\partial P_{vol}}{\partial C}$:

$$
\mathbb{C}_{T, vol} = 2 \left[ K \left( J - \frac{1}{2} \right) J C^{-1} \otimes C^{-1} - (p_J J) \mathbb{I}_{C^{-1}} \right]
$$

---

## 2. Der deviatorische Teil

Jetzt schauen wir uns $\boldsymbol{P}_{dev}$ an:

$$
\boldsymbol{P}_{dev} = 2 W_1 J^{-2/3} \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)
$$

Hier nutzen wir einen massiven Vorteil: Wir wissen aus der vorherigen Herleitung bereits, dass $J^{-2/3} \left( \boldsymbol{I} - \frac{1}{3} I_1 \boldsymbol{C}^{-1} \right)$ exakt die Ableitung $\frac{\partial \bar{I}_1}{\partial \boldsymbol{C}}$ ist! 

Wir können die Gleichung also extrem kompakt schreiben:

$$
\boldsymbol{P}_{dev} = 2 W_1 \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}}
$$

Auch hier wenden wir wieder die Produktregel $(f \cdot \boldsymbol{T})' = f' \otimes \boldsymbol{T} + f \cdot \boldsymbol{T}'$ an:

$$
\frac{\partial \boldsymbol{P}_{dev}}{\partial \boldsymbol{C}} = 2 \frac{\partial W_1}{\partial \boldsymbol{C}} \otimes \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} + 2 W_1 \frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C} \partial \boldsymbol{C}}
$$

**Schritt 2.1: Ableitung von $W_1$**
Mit der Kettenregel leiten wir die skalare Funktion $W_1$ ab:

$$
\frac{\partial W_1}{\partial \boldsymbol{C}} = \frac{\partial W_1}{\partial \bar{I}_1} \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} = W_{11} \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}}
$$

**Schritt 2.2: Die zweite Ableitung von $\bar{I}_1$**
Wir brauchen den Tensor 4. Stufe $\frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C} \partial \boldsymbol{C}}$. Dieser wird nach exakt denselben Regeln gebildet:

$$
\begin{aligned}
\frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C} \partial \boldsymbol{C}} &= J^{-2/3} \Bigg( -\frac{1}{3}\boldsymbol{C}^{-1} \otimes \boldsymbol{I} - \frac{1}{3}\boldsymbol{I} \otimes \boldsymbol{C}^{-1} \\
&\quad + \frac{1}{9}I_1\boldsymbol{C}^{-1} \otimes \boldsymbol{C}^{-1} + \frac{1}{3}I_1 \mathbb{I}_{C^{-1}} \Bigg)
\end{aligned}
$$

**Schritt 2.3: Deviatorischen Teil zusammensetzen**
Einsetzen liefert uns die Ableitung der deviatorischen Spannung. Um den Steifigkeitstensor zu erhalten multiplizieren wir es noch mit 2:

$$
\mathbb{C}_{T, dev} = 4 W_{11} \left( \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} \otimes \frac{\partial \bar{I}_1}{\partial \boldsymbol{C}} \right) + 4 W_1 \left( \frac{\partial^2 \bar{I}_1}{\partial \boldsymbol{C} \partial \boldsymbol{C}} \right)
$$

---

## 3. Das finale Ergebnis

Wenn du das in deinem Code implementierst, addierst du einfach die beiden berechneten Tensoren 4. Stufe:

$$
\mathbb{C}_T = \mathbb{C}_{T, dev} + \mathbb{C}_{T, vol}
$$


### Einschub: Herleitung von $\frac{\partial (p_J J)}{\partial \boldsymbol{C}}$

Bevor wir ableiten, müssen wir uns klar machen, wie die Variablen überhaupt zusammenhängen:
* $J$ (das Volumenverhältnis) ist eine Funktion des Tensors $\boldsymbol{C}$.
* $p_J$ (die Ableitung des volumetrischen Potentials) ist eine Funktion von $J$. 
* Weil $p_J$ von $J$ abhängt und $J$ wiederum von $\boldsymbol{C}$ abhängt, hängt $p_J$ indirekt auch von $\boldsymbol{C}$ ab.

Wir wollen nun das Produkt aus den beiden skalaren Werten $(p_J \cdot J)$ nach dem Tensor $\boldsymbol{C}$ ableiten. Das funktioniert in drei einfachen mathematischen Schritten.

#### Schritt 1: Die Produktregel anwenden
Die Produktregel besagt: $(u \cdot v)' = u' \cdot v + u \cdot v'$.
Wir setzen $u = p_J$ und $v = J$. Die formale Ableitung nach $\boldsymbol{C}$ lautet dann:

$$
\frac{\partial (p_J \cdot J)}{\partial \boldsymbol{C}} = \frac{\partial p_J}{\partial \boldsymbol{C}} \cdot J + p_J \cdot \frac{\partial J}{\partial \boldsymbol{C}}
$$

#### Schritt 2: Die Kettenregel für $p_J$ anwenden
Betrachten wir den vorderen Term $\frac{\partial p_J}{\partial \boldsymbol{C}}$. Da die Funktion $p_J$ nur $J$ als direkte Variable kennt, müssen wir den Umweg über die Kettenregel ("Äußere Ableitung mal innere Ableitung") gehen:

$$
\frac{\partial p_J}{\partial \boldsymbol{C}} = \underbrace{\frac{\partial p_J}{\partial J}}_{\text{äußere Abl.}} \cdot \underbrace{\frac{\partial J}{\partial \boldsymbol{C}}}_{\text{innere Abl.}}
$$

#### Schritt 3: Einsetzen und Ausklammern
Nun setzen wir unser Ergebnis für die Kettenregel aus Schritt 2 zurück in die Produktregel aus Schritt 1 ein:

$$
\frac{\partial (p_J \cdot J)}{\partial \boldsymbol{C}} = \left( \frac{\partial p_J}{\partial J} \cdot \frac{\partial J}{\partial \boldsymbol{C}} \right) \cdot J + p_J \cdot \frac{\partial J}{\partial \boldsymbol{C}}
$$

Wir können die Reihenfolge bei der Multiplikation von Skalaren tauschen, also ziehen wir das $J$ direkt an den Bruch heran:

$$
\frac{\partial (p_J \cdot J)}{\partial \boldsymbol{C}} = \left( \frac{\partial p_J}{\partial J} \cdot J \right) \frac{\partial J}{\partial \boldsymbol{C}} + p_J \cdot \frac{\partial J}{\partial \boldsymbol{C}}
$$

Jetzt sehen wir, dass beide Summanden auf der rechten Seite mit dem Tensor-Term $\frac{\partial J}{\partial \boldsymbol{C}}$ multipliziert werden. Wir können diesen gemeinsamen Term am Ende der Gleichung ausklammern und erhalten exakt die gesuchte Form:

$$
\frac{\partial (p_J \cdot J)}{\partial \boldsymbol{C}} = \left( \frac{\partial p_J}{\partial J} J + p_J \right) \frac{\partial J}{\partial \boldsymbol{C}}
$$
