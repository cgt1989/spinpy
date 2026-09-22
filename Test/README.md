# Test — réplica de los resultados publicados

Esta carpeta existe para poder responder a una pregunta concreta:

> *«¿De dónde sale este método y cómo sé que vuestra implementación lo hace
> bien?»*

La respuesta no es «confía en nuestras pruebas». Es: **aquí está la figura del
artículo original, y aquí al lado la misma figura regenerada con nuestro
código, con las comprobaciones numéricas que tendrían que fallar si la
implementación se hubiera desviado.**

No estamos reinventando el método. Lo estamos **implementando y verificando**.

---

## Qué hay aquí

```
Test/
  referencia/                    figuras publicadas, tal cual (ver FUENTES.md)
    Kumar2020_Fig1.png
    Kumar2020_Fig2.png
    Zheng2021_Fig1.png
    Guo2024_Fig7.png
    FUENTES.md                   cita, licencia y procedencia de cada imagen
  replicar_kumar2020.py          el método: las cuatro clases de anisotropía
  replicar_zheng2021.py          las cotas: Voigt, Hashin-Shtrikman, ortotropía
  replicar_guo2024.py            las curvaturas: el perfil (k1, k2)
  resultados/                    lo que producen (figura + informe JSON)
```

Son **tres réplicas de tres artículos**, y cada una comprueba algo que las
otras no:

| Réplica | Qué pone a prueba | Podría fallar si… |
|---|---|---|
| Kumar et al. (2020) | que el método generador es el suyo: unión de conos, conjunto de nivel, las cuatro clases | el muestreo de direcciones no fuera su ecuación (2) |
| Zheng et al. (2021) | **magnitudes absolutas**: cada superficie elástica dentro de Voigt y de Hashin-Shtrikman; la ortotropía del tensor; por qué ρ ≥ 0.3 | la homogeneización estuviera escalada por un factor, o el orden de Voigt fuera otro |
| Guo et al. (2024) | las **curvaturas** de la interfaz, contra dos respuestas cerradas | el estimador de curvatura o la normalización de longitudes estuvieran mal |

La primera comprueba **proporciones** (E₁/E₃); la segunda, **cotas**, que un
error de escala no puede pasar; la tercera, una métrica que las otras dos no
miran y que es la única que distingue una red conectada de un montón de islas.

Se ejecutan desde la aplicación —menú **Validación**, una pestaña por
artículo— o desde la consola:

```bash
python replicar_kumar2020.py            # completo, ~4 min
python replicar_zheng2021.py            # completo, ~90 min (homogeneiza a 40³)
python replicar_guo2024.py              # completo, ~4 min
python replicar_kumar2020.py --rapido   # mallas menores, ~1.5 min
```

La de Zheng es larga porque hace cuatro cosas caras: homogeneiza cuatro clases a
40³ (~25 min), barre trece densidades a 128³ (~35), promedia cuatro
realizaciones para la ortotropía (~8) y mide la dispersión con cinco semillas a
N = 2000 (~15). Su modo rápido baja las tres primeras pero **no** la cuarta, que
tiene configuración fija:

```bash
```

El modo `--rapido` sirve para comprobar que la tubería funciona. Para dar por
buena la implementación, ejecución completa: en la réplica de Guo el criterio
de error del estimador queda **justo en el límite** con mallas pequeñas, y
pasar por una centésima es pasar por suerte.

---

## La réplica de Kumar et al. (2020)

### Referencia replicada

> Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
> spinodoid metamaterials. *npj Computational Materials, 6*(1), 1–10.
> Artículo 73. https://doi.org/10.1038/s41524-020-0341-6

Licencia del artículo: **CC BY 4.0**
(https://creativecommons.org/licenses/by/4.0/), © The Author(s) 2020. Ese es
el motivo por el que sus figuras pueden reproducirse aquí, con atribución y
sin modificaciones. Ver `referencia/FUENTES.md`.

### Qué implementa nuestro código, y por qué es *el mismo* método

El artículo define el espinodoide en tres ecuaciones (pp. 3–4):

**(1) Campo gaussiano aleatorio**

```
phi(x) = sqrt(2/N) · SUM_i  cos( beta · n_i · x + gamma_i )
         n_i ~ U(S²),   gamma_i ~ U([0, 2*pi))
```

**(2) Anisotropía**, restringiendo las direcciones de onda a la **unión** de
tres conos alrededor de los ejes cartesianos:

```
n_i ~ U{ k in S² : |k·e1| > cos(th1)  o  |k·e2| > cos(th2)  o  |k·e3| > cos(th3) }
      th1, th2, th3  in  {0} U [th_min, 90°],   th_min = 15°
```

**(3) Conjunto de nivel**: sólido donde `phi(x) <= phi0`, con

```
phi0 = sqrt(2) · erf⁻¹(2*rho - 1)
```

Nuestra implementación **es esa**, no una parecida:

| Artículo | Nuestro código |
|---|---|
| ec. (2), unión de conos | `spinpy/grf.py::_waves_rechazo` — acepta un candidato isótropo si el ángulo respecto de **algún** eje es menor que el θ de ese eje |
| ec. (3), umbral | `spinpy/grf.py::level_set` — `sqrt(2)*erfinv(2*rho-1)` |
| θ_min = 15° | límite inferior de los presets de ángulo del ajuste |

El esquema alternativo `equitativo` (del *TPMS-Scaffolds-generator*) reparte
las ondas **en partes iguales** entre los conos y **no** es la ecuación (2).
La réplica usa `rechazo` porque es el único que corresponde al artículo — y la
comprobación C6 lo demuestra midiéndolo.

### Parámetros: lo citado y lo deducido

Una réplica que mete valores inventados no replica nada, así que se distingue:

| Parámetro | Valor | Procedencia |
|---|---|---|
| ρ | 0.5 | **verbatim**, leyenda Fig. 2 |
| β | 15π | **verbatim**, leyenda Fig. 5 («para comparación visual») |
| θ_min | 15° | **verbatim**, p. 4 |
| lamelar | (30, 0, 0)° | **verbatim**, leyenda Fig. 3 |
| columnar | (0, 30, 30)° | **verbatim**, leyenda Fig. 3 |
| cúbica | (30, 30, 30)° | *deducido* de la regla de clases (p. 4); el artículo no publica la terna |
| isótropa | (90, 90, 90)° | **verbatim** del texto (p. 4): «θ₃ = π/2, reduciéndose a la esfera unidad» |
| terna desigual | (15, 45, 90)° | **nuestra**, sólo para contrastar los dos esquemas de muestreo |

Las ternas de los seis paneles de la Fig. 2 están en la capa vectorial del PDF
y no son extraíbles, así que **no** se reproduce panel por panel: se reproducen
las cuatro clases que el texto nombra.

### Una diferencia visual que conviene declarar

Nuestras estructuras salen **más finas** que las de su Fig. 2: donde ellos
muestran unas seis láminas, aquí se ven muchas más. No es un error de
implementación sino de escala. El número de onda β fija el tamaño
característico, y el artículo **no publica el β de la Fig. 2**; el 15π que
usamos viene de la leyenda de la Fig. 5. Con un β menor las estructuras salen
más gruesas y el parecido visual es mayor, sin que cambie ninguna de las seis
comprobaciones: la anisotropía depende de los conos, no del tamaño de poro.

Se deja el valor citable en lugar del que «queda mejor». Una réplica que ajusta
un parámetro no publicado hasta que la figura se parece deja de ser una réplica.

### Qué se comprueba, y qué podría fallar

Una figura parecida no demuestra nada. El guion evalúa predicciones
**falsables** que salen del propio artículo (p. 4): *«las interfaces de las
topologías se alinean preferentemente perpendiculares a esos vectores nᵢ»*.
De ahí se sigue que un cono alrededor de un eje produce láminas apiladas a lo
largo de ese eje, es decir que **ese eje queda blando**:

| | Predicción | Criterio |
|---|---|---|
| C1 | Lamelar (cono sólo en e₁): e₁ es el eje **blando** | E₁/E₃ < 0.6 |
| C2 | Columnar (conos en e₂ y e₃): e₁ es el eje **rígido** | E₁/E₃ > 1.6 |
| C3 | Isótropa: E₁ = E₂ = E₃ | E_max/E_min < 1.25 |
| C4 | Cúbica: mucho **menos** anisótropa en los ejes que la lamelar y la columnar | < ½ de la mayor anisotropía |
| C5 | Lamelar y columnar ordenan al revés el eje e₁ | razón lamelar < razón columnar |
| C6 | Con conos **desiguales**, `equitativo` no es la ecuación (2) | difieren > 5 % en DA |

Los márgenes se declaran antes de medir y son holgados a propósito: el
generador es estocástico y la homogeneización va en malla pequeña. Lo que se
comprueba es el **signo y el orden de magnitud** de la anisotropía, que es lo
que el artículo afirma — no una cifra.

La superficie elástica E(d) se autocomprueba antes de dibujarse: con un tensor
isótropo tiene que salir constante. Sin eso, un error en los factores ½ y ¼ de
la notación de Voigt daría superficies convincentes y falsas.

### Dos comprobaciones que hubo que rehacer

Se dejan escritas porque son exactamente la clase de cosa que una carpeta de
validación debe poner a la vista, no barrer.

**C4 pasaba por una centésima.** La primera versión exigía E_max/E_min < 1.25
para la cúbica y obtenía **1.249**. Eso es pasar por suerte, no por buena
implementación. La causa es real y vale la pena entenderla: el artículo habla
del **ensemble** —con tres conos iguales la distribución de direcciones tiene
simetría cúbica y E₁ = E₂ = E₃ *en media*—, pero **una sola realización con
una semilla concreta no es simétrica**. Por eso C4 se formula ahora en
relativo, que es la afirmación robusta y sigue siendo falsable. La cifra cruda
se reporta igual: quien lea ve el 1.249, no una comprobación maquillada.

### Una nota sobre C6, que empezó fallando

La primera versión de C6 comparaba los dos esquemas sobre la terna **columnar**
(0, 30, 30) y **fallaba**: rechazo daba DA 1.804 y equitativo 1.843, un 2.2 %.
No era un fallo del código sino de la prueba: con dos conos **iguales** los dos
esquemas reparten las ondas mitad y mitad y coinciden por construcción. La
diferencia sólo aparece con conos de tamaños distintos, donde el reparto
equitativo deja de ser proporcional al ángulo sólido. Con (15, 45, 90) la
diferencia es del 27 %.

Queda escrito porque es exactamente la clase de error que una carpeta de
validación debe dejar a la vista, no barrer.

---

## La réplica de Zheng et al. (2021)

### Referencia replicada

> Zheng, L., Kumar, S., & Kochmann, D. M. (2021). Data-driven topology
> optimization of spinodoid metamaterials with seamlessly tunable anisotropy.
> *Computer Methods in Applied Mechanics and Engineering, 383*, 113894.
> https://doi.org/10.1016/j.cma.2021.113894

Licencia del artículo: **CC BY**, © 2021 The Author(s), Elsevier B.V.

### Qué se replica, y qué no

El grueso del artículo —optimización topológica multiescala con IPOPT y una red
neuronal sustituta de la homogeneización— **no** se replica: no hay problema
macroescala en esta aplicación ni motivo para tenerlo. Se replica su
**Sección 2 y su Figura 1(b–e)**, que es donde el artículo hace afirmaciones
comprobables sobre el espinodoide en sí.

### Las cuatro ternas, ahora verbatim

Los pies de su Fig. 1 publican las cuatro, **incluida la cúbica** que la Fig. 2
de Kumar et al. (2020) no da y que `replicar_kumar2020.py` tuvo que deducir:

| Clase | ρ | (θ₁, θ₂, θ₃) | Procedencia |
|---|---|---|---|
| lamelar | 0.5 | (0, 0, 15)° | **verbatim**, pie Fig. 1(b) |
| columnar | 0.5 | (15, 15, 0)° | **verbatim**, pie Fig. 1(c) |
| cúbica | 0.5 | (15, 15, 15)° | **verbatim**, pie Fig. 1(d) |
| isótropa | 0.5 | (90, 90, 90)° | **verbatim**, pie Fig. 1(e) |

Están **traspuestas** respecto de las de la réplica de Kumar —allí el cono
activo de la lamelar era e₁ y aquí es e₃—, así que comprobar el mapa «cono
activo → eje blando» con estas ternas es una comprobación independiente y no la
misma dos veces.

### Qué se comprueba

| | Predicción | Criterio |
|---|---|---|
| Z1 | Ninguna superficie elástica sale de la esfera de **Voigt** (la gris clara de su figura, de radio ρ·E_s = 0.5) | E(d) ≤ ρ·E_s en las cuatro clases |
| Z2 | La isótropa no supera la cota superior de **Hashin-Shtrikman** (la gris oscura del panel e, E = 0.333 E_s) | E(d) ≤ E_HS⁺ |
| Z3 | El tensor homogeneizado es **ortótropo** (p. 11): al promediar realizaciones, el acoplamiento normal–cortante baja hacia cero | el promedio baja respecto de las realizaciones sueltas y queda ≤ 5 % de la diagonal |
| Z4 | Con sus ternas, el eje del cono activo es el **blando** en la lamelar y el **rígido** en la columnar | lamelar E₃/E₁ < 0.6; columnar > 1.6 |
| Z5 | A ρ = 0.3 —su ρ_min— el material sigue formando caminos que atraviesan, con el θ_min = π/6 que fija la **misma frase** | ≥ 95 % de material portante **en todas las semillas** de las clases con θ ≥ 30° |
| Z6 | Por debajo de ρ_min el material desconectado crece en todas las clases | crece en las cuatro |
| Z7 | **Hallazgo**, no una predicción del artículo: con los 15° del pie de su Fig. 1 la lamelar **no** cumple el propósito de ρ_min; con los 30° del texto sí | al menos dos realizaciones de cinco rotas a 15°, y ninguna a 30° |

**Z1 y Z2 son cotas de verdad**, con fórmula cerrada, y ahí está su valor: la
réplica de Kumar comprueba cocientes E₁/E₃, que sobrevivirían intactos si toda
nuestra homogeneización estuviera multiplicada por un factor constante. Una
cota no sobrevive a eso.

**Z5 mide la fracción PORTANTE, no la «fuera de la mayor componente».** Las dos
cosas no son lo mismo y para una clase importa mucho. `1 − mayor_componente` es
la lectura literal de «dominios disjuntos» y sirve para la lamelar o la
isótropa; para la **columnar** no, porque esa clase son columnas paralelas que
atraviesan la probeta, legítimamente separadas unas de otras, y **todas cargan**.
Medida con una semilla desafortunada llegó a dar 13.6 % «desconectado» teniendo
el 98 % de material portante. Penalizarla por eso sería penalizarla por no ser
una sola pieza, que no es un defecto mecánico. El informe lleva las dos cifras.

Y se exige en la **peor semilla**, no en la media: lo que se afirma es que la
regla del artículo funciona, y una regla que funciona de media no funciona.

**Z3 tuvo que reformularse, y queda escrito por qué.** La primera versión
exigía un umbral pequeño a **una** realización y fallaba: a 40³ la clase
columnar da un 2.6 % de acoplamiento pero la lamelar da un **11.1 %**. No es un
error de implementación —lo que tiene simetría es la *distribución* de
direcciones de onda, no una realización suelta, y con las láminas un poco
torcidas el acoplamiento normal–cortante es real— sino un error de la prueba:
el artículo no afirma eso de una realización. Z3 comprueba ahora lo que el
artículo sí dice, que en el **ensemble** el acoplamiento se anula: se promedian
cuatro realizaciones y se mira si baja. Medido a 24³ en la clase lamelar:
4.8 %, 5.2 %, 1.6 % y 3.4 % sueltas → **1.8 %** el tensor promedio. Es la misma
corrección que hubo que hacerle a C4 en la réplica de Kumar, y por la misma
razón.

### Por qué ρ ≥ 0.3, medido — y la sorpresa

El artículo restringe el espacio de diseño a ρ ≥ 0.3 «para evitar dominios
sólidos disjuntos» (p. 5).

**Esta primera tabla es UNA realización por casilla**, a 128³ con N = 1000: sirve
para ver la *forma* de la curva y **no** para leer el valor de un punto, porque
en esa configuración la dispersión entre realizaciones es de decenas de puntos
porcentuales (ver más abajo). Las cifras con las que se deciden Z5 y Z7 están en
la tabla siguiente, medidas con cinco semillas. Fracción de hueso fuera de la
mayor componente:

| ρ | lamelar (15°) | columnar | cúbica | isótropa | lamelar (30°) |
|---|---|---|---|---|---|
| 0.15 | 93.9 % | 62.6 % | 78.4 % | 81.2 % | 90.4 % |
| 0.20 | 92.4 % | 9.4 % | 44.0 % | 13.8 % | 77.5 % |
| 0.25 | 86.2 % | 1.7 % | 3.6 % | 4.7 % | 33.9 % |
| **0.30** | **48.4 %** | **0.8 %** | **0.7 %** | **0.4 %** | **2.0 %** |
| 0.35 | 22.6 % | 0.4 % | 0.0 % | 0.2 % | 0.1 % |
| 0.40 | 0.5 % | 0.0 % | 0.1 % | 0.1 % | 0.0 % |
| 0.50 | 0.0 % | 0.0 % | 0.0 % | 0.0 % | 0.0 % |

Para columnar, cúbica e isótropa la regla del artículo funciona limpiamente: a
ρ = 0.3 queda menos del 1 % de hueso suelto y por debajo se desploma. **Para la
lamelar no**: a ρ = 0.3 casi la mitad del hueso está en componentes distintas, y
no baja del 5 % hasta ρ ≈ 0.40.

### Las cifras con las que se decide: cinco semillas, configuración fija

El bloque de dispersión **no sigue al modo rápido**: usa siempre 96³ y **N = 2000
ondas**, porque el fenómeno que miden Z5 y Z7 sólo existe en ese régimen — a
N = 400 la lamelar da 1.7 ± 2.5 % y no hay nada que ver, de modo que las dos
comprobaciones fallarían en modo rápido por mirar donde el efecto no está.

| clase | ρ | desconectado (5 semillas) | portante mínimo |
|---|---|---|---|
| lamelar (15°) | 0.30 | **48.7 ± 11.4 %** | **0.0 %** |
| lamelar (15°) | 0.35 | 8.9 ± 16.3 % | 0.0 % |
| columnar | 0.30 | 1.8 ± 1.7 % | 97.9 % |
| cúbica | 0.30 | 0.3 ± 0.3 % | 99.3 % |
| isótropa | 0.30 | 0.7 ± 0.7 % | 98.2 % |
| **lamelar (30°)** | 0.30 | **0.8 ± 1.3 %** | **96.9 %** |

Las cinco realizaciones de la lamelar a 15° salen rotas (30.6, 45.4, 51.4, 57.9
y 58.3 % desconectado) y en alguna la fracción portante cae a **cero**: no
atraviesa en absoluto. Con los 30° del texto, ninguna.

### Y ahí la discrepancia interna del artículo deja de ser inocua

El texto de la p. 5 dice `θ_min = π/6` (30°). El pie de su Fig. 3 dice
`θ_min = 15°`, los pies de la Fig. 1 usan 15° y el Benchmark I arranca de
(15, 0, 0). Los dos valores **no son intercambiables** para la otra mitad de esa
misma frase, la que fija ρ ≥ 0.3. Medido a ρ = 0.30:

| lamelar | desconectado |
|---|---|
| (0, 0, 15)° — el de los pies de figura | **48.4 %** |
| (0, 0, 30)° — el π/6 del texto | **2.0 %** |
| (0, 0, 45)° | 0.2 % |

La regla ρ_min = 0.3 cumple su propósito **si se aplica junto al θ_min = π/6 que
la acompaña en el texto**, y no lo cumple con los 15° de los pies de figura. La
razón es geométrica y no tiene misterio: un solo cono estrecho produce láminas
casi planas, y unas láminas planas y separadas son componentes distintas por
definición; con dos o tres conos, o con un cono ancho, las láminas se ondulan y
se tocan. Por eso Z5 comprueba la regla **donde el artículo la enuncia** y Z7
deja escrito el hallazgo.

### Cuidado al medir esto: la dispersión entre realizaciones es enorme

La primera versión de esta réplica medía la desconexión con **una** realización
por punto y presentaba series como si fueran tendencias. No lo eran. Medido
después con cinco semillas por punto sobre la lamelar de 15° a ρ = 0.30
(`Estudio_Percolacion`):

| N | media ± sd | valores |
|---|---|---|
| 200 | 19.3 ± 37.7 % | 0.0 · 86.5 · 0.7 · 8.8 · 0.8 |
| 400 | 1.7 ± 2.5 % | 2.3 · 0.2 · 0.0 · 0.0 · 5.9 |
| 1000 | 25.3 ± 34.5 % | 48.4 · 1.4 · 0.0 · 75.0 · 2.0 |
| **2000** | **46.0 ± 10.1 %** | 45.3 · 51.6 · 57.8 · 30.6 · 44.7 |

Hasta N = 1000 la dispersión es de **decenas de puntos porcentuales** y cualquier
cifra suelta de esa zona es una tirada, no una medida. A N = 2000 la dispersión
se desploma y la estructura se autopromedia: ahí sí hay señal. La afirmación
«con un cono estrecho la lamelar no cumple el propósito de ρ_min» se sostiene;
lo que no se sostenía era la evidencia de tres puntos con una semilla que se
publicó primero.

En **resolución**, en cambio, la cosa está tranquila: 128³, 160³ y 192³ dan 48.4,
48.3 y 48.4 % sobre la misma realización.

Queda escrito porque es el tipo de error que esta carpeta existe para no barrer:
una serie monótona de tres puntos, cada uno con una semilla distinta, parece una
tendencia y no lo es.

**Y hay una explicación para tanta dispersión.** `Estudio_Percolacion` mide el
umbral de percolación de esa clase en **ρ_c ≈ 0.31**: ρ = 0.30 cae justo encima,
que es exactamente donde cada realización se decide por un lado o por otro.

**Lo que esto significa para el ajuste de este proyecto.** El espinodoide
ajustado a nuestro VOI tiene ρ = 0.339 y un 13.8 % de hueso en islas: encaja
justo en la columna lamelar de la tabla (22.6 % a ρ = 0.35), y no en ninguna de
las otras tres. El ajuste no está aterrizando en un punto cualquiera del espacio
de diseño, está aterrizando en el rincón frágil. Conviene además revisar con
cuántas ondas se midió: `fit.ajustar_spinodoide` usa 700 por omisión, y entre
400 y 1000 esta métrica cambia de 2 % a 48 %.

### Las condiciones de contorno no son las mismas, y eso tiene signo

Ellos homogeneizan con condiciones **afines** y lo dicen expresamente (p. 5):
*«the computed response provides an upper bound to the actual effective
stiffness»*. Nosotros usamos **periódicas**, que dan el valor más bajo de los
dos. Por eso las comprobaciones están escritas como desigualdades y no como
comparaciones cifra a cifra: cifra a cifra no serían comparables.

---

## La réplica de Guo et al. (2024)

### Referencia replicada

> Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface curvatures
> by deep learning. *Advanced Intelligent Systems, 6*(6), 2300789.
> https://doi.org/10.1002/aisy.202300789

Licencia del artículo: **CC BY 4.0**, © 2024 The Authors, Wiley-VCH GmbH.

### Qué se replica, y qué no

No se replica su diseño inverso por redes neuronales —no tenemos su marco de
campo de fase ni sus 18.000 casos de entrenamiento— ni su muestra de hueso
trabecular, que es de Tozzi et al. y no la tenemos. Se replican las dos
columnas **«Target»** de su Figura 7, que son geometrías completamente
especificadas:

| Objetivo | Parámetros | Procedencia |
|---|---|---|
| espinodoide | β = 15π, Q = 1000, ρ = 0.3, θ = (60, 30, 10)° | **verbatim**, recuadro de la Fig. 7b |
| superficie nodal periódica | sin x·sin 1.8y + sin y·sin 1.8z + sin z·sin 1.8x = 0.5 | **verbatim**, ecuación (10) |

La PNS es la pieza clave: al ser implícita y analítica, sus curvaturas tienen
**fórmula cerrada**, y esa fórmula no es nuestra. Es la respuesta contra la que
se valida el estimador discreto que luego hay que aplicar al hueso real, donde
no hay campo analítico.

### La métrica nueva: `spinpy/curvatura.py`

Las curvaturas principales (κ₁ ≥ κ₂) de la interfaz hueso-vacío, por dos
caminos que se comparan entre sí:

- **exacto**, del gradiente y el hessiano del campo (el GRF de un espinodoide
  tiene derivadas analíticas: `grf.derivadas_grf`);
- **discreto**, ajustando la segunda forma fundamental por triángulo
  (Rusinkiewicz 2004) sobre la malla del campo continuo.

Convención de signo, declarada en el módulo y verificada en el bloque 09 de
`tests/`: la normal apunta **hacia el vacío**, de modo que una bola de hueso
tiene κ = +1/R, un poro κ = −1/R y una trabécula es silla de montar
(κ₁ > 0 > κ₂).

### Qué se comprueba

| | Predicción | Criterio |
|---|---|---|
| R1 | El perfil del espinodoide cae donde cae el suyo: silla, con el centro dentro de la caja leída de su panel (b) | κ₁ en [10, 60], κ₂ en [−40, 10] |
| R2 | Doblar β dobla **todas** las curvaturas | razón = 2.00 ± 5 % |
| R3 | A ρ = 0.5, sólido y vacío son intercambiables: la curvatura media del ensemble es cero | ⟨H⟩ ≤ 15 % de su desviación, en valor absoluto |
| R4 | En el espinodoide, el estimador discreto coincide con la fórmula cerrada | error mediano ≤ 10 % |
| R5 | En la PNS, ídem contra una fórmula que **no es nuestra** | error mediano ≤ 10 % |
| R6 | La PNS **no** es superficie mínima; el artículo la elige por eso | ⟨H⟩ > 10 % de la escala, en valor absoluto |
| R7 | Marching cubes sobre la máscara **binaria** sobreestima el área frente a la isosuperficie del campo continuo | razón entre 1.02 y 1.20 |

R2 y R3 son predicciones **exactas** que no dependen de leer ninguna figura:
salen de que el GRF con β' = 2β es el mismo campo con las longitudes a la
mitad, y de que el campo y su opuesto tienen la misma distribución. R7 es la
corrección C4 de este proyecto vista desde otro lado: el sesgo del binario se
midió en +8.5 % sobre una esfera y aquí reaparece en **+8.8 %** sobre un
espinodoide, que es una confirmación independiente y sobre otra geometría.

### La unidad de la curvatura: esto es deducido, no citado

Los ejes de sus perfiles van de −100 a +100. El artículo dice que el dominio es
[0, 100]³ y que las longitudes están normalizadas respecto de su tamaño, lo que
admite dos lecturas. La que cuadra con la física es que κ está en unidades de
1/(lado del dominio): la longitud de onda del espinodoide es 2π/β = 2/15 del
lado y el radio de una trabécula es un cuarto de eso → κ ~ 30; la celda de la
PNS es 10π y sus radios son de orden 1 → κ ~ 30; y el límite del eje, κ = 100,
es un radio de 1/100 del lado, o sea **un voxel** en una rejilla de 100³. Las
tres cosas cuadran a la vez. Si la lectura fuera la otra, nuestras cifras
saldrían multiplicadas por 100 y R1 fallaría a gritos, no en silencio.

### Una diferencia que conviene declarar

Nuestras manchas salen **más anchas** que las suyas. El centro cae donde tiene
que caer —que es lo que comprueba R1— pero la dispersión es mayor, y hay dos
causas: ellos miden sobre un campo de fase con interfaz **difusa**, más lisa
que nuestro conjunto de nivel muestreado en rejilla, y el estimador discreto
añade la suya. La segunda se puede separar y se separa: el informe lleva la
dispersión de H por los dos caminos sobre los mismos vértices.

---

## Todas las referencias que sostienen la aplicación

Citadas en APA 7 con DOI. Las que aparecen en la interfaz o en los informes
llevan además la cita en el propio código.

### Generación de espinodoides

- Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
  spinodoid metamaterials. *npj Computational Materials, 6*(1), 1–10.
  Artículo 73. https://doi.org/10.1038/s41524-020-0341-6
- Cahn, J. W. (1965). Phase separation by spinodal decomposition in isotropic
  systems. *The Journal of Chemical Physics, 42*(1), 93–99.
  https://doi.org/10.1063/1.1695731
- Moerman, K. M. (2018). GIBBON: The geometry and image-based bioengineering
  add-On. *Journal of Open Source Software, 3*(22), 506.
  https://doi.org/10.21105/joss.00506
- Zheng, L., Kumar, S., & Kochmann, D. M. (2021). Data-driven topology
  optimization of spinodoid metamaterials with seamlessly tunable anisotropy.
  *Computer Methods in Applied Mechanics and Engineering, 383*, 113894.
  https://doi.org/10.1016/j.cma.2021.113894

### Morfometría ósea

- Bouxsein, M. L., Boyd, S. K., Christiansen, B. A., Guldberg, R. E., Jepsen,
  K. J., & Müller, R. (2010). Guidelines for assessment of bone microstructure
  in rodents using micro-computed tomography. *Journal of Bone and Mineral
  Research, 25*(7), 1468–1486. https://doi.org/10.1002/jbmr.141
- Parfitt, A. M., Drezner, M. K., Glorieux, F. H., Kanis, J. A., Malluche, H.,
  Meunier, P. J., Ott, S. M., & Recker, R. R. (1987). Bone histomorphometry:
  Standardization of nomenclature, symbols, and units. *Journal of Bone and
  Mineral Research, 2*(6), 595–610. https://doi.org/10.1002/jbmr.5650020617
- Harrigan, T. P., & Mann, R. W. (1984). Characterization of microstructural
  anisotropy in orthotropic materials using a second rank tensor. *Journal of
  Materials Science, 19*(3), 761–767. https://doi.org/10.1007/BF00540446
- Odgaard, A., & Gundersen, H. J. G. (1993). Quantification of connectivity in
  cancellous bone, with special emphasis on 3-D reconstructions. *Bone, 14*(2),
  173–182. https://doi.org/10.1016/8756-3282(93)90245-6
- Hildebrand, T., & Rüegsegger, P. (1997). Quantification of bone
  microarchitecture with the structure model index. *Computer Methods in
  Biomechanics and Biomedical Engineering, 1*(1), 15–23.
  https://doi.org/10.1080/01495739708936692
- Salmon, P. L., Ohlsson, C., Shefelbine, S. J., & Doube, M. (2015). Structure
  model index does not measure rods and plates in trabecular bone. *Frontiers
  in Endocrinology, 6*, 162. https://doi.org/10.3389/fendo.2015.00162
- Doube, M., Kłosowski, M. M., Arganda-Carreras, I., Cordelières, F. P.,
  Dougherty, R. P., Jackson, J. S., Schmid, B., Hutchinson, J. R., &
  Shefelbine, S. J. (2010). BoneJ: Free and extensible bone image analysis in
  ImageJ. *Bone, 47*(6), 1076–1079. https://doi.org/10.1016/j.bone.2010.08.023

### Curvatura de la interfaz

- Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface curvatures
  by deep learning. *Advanced Intelligent Systems, 6*(6), 2300789.
  https://doi.org/10.1002/aisy.202300789
- Rusinkiewicz, S. (2004). Estimating curvatures and their derivatives on
  triangle meshes. En *2nd International Symposium on 3D Data Processing,
  Visualization and Transmission* (pp. 486–493). IEEE.
  https://doi.org/10.1109/TDPVT.2004.1335277
- Callens, S. J. P., Tourolle né Betts, D. C., Müller, R., & Zadpoor, A. A.
  (2021). The local and global geometry of trabecular bone. *Acta
  Biomaterialia, 130*, 343–361. https://doi.org/10.1016/j.actbio.2021.06.013
- do Carmo, M. P. (1976). *Differential geometry of curves and surfaces*.
  Prentice-Hall.

### Mecánica y homogeneización

- Andreassen, E., & Andreasen, C. S. (2014). How to determine composite
  material properties using numerical homogenization. *Computational Materials
  Science, 83*, 488–495. https://doi.org/10.1016/j.commatsci.2013.09.006
- Hill, R. (1952). The elastic behaviour of a crystalline aggregate.
  *Proceedings of the Physical Society. Section A, 65*(5), 349–354.
  https://doi.org/10.1088/0370-1298/65/5/307
- Hashin, Z., & Shtrikman, S. (1963). A variational approach to the theory of
  the elastic behaviour of multiphase materials. *Journal of the Mechanics and
  Physics of Solids, 11*(2), 127–140.
  https://doi.org/10.1016/0022-5096(63)90060-7
- Gibson, L. J., & Ashby, M. F. (1997). *Cellular solids: Structure and
  properties* (2.ª ed.). Cambridge University Press.
  https://doi.org/10.1017/CBO9781139878326
- Pistoia, W., van Rietbergen, B., Lochmüller, E.-M., Lill, C. A., Eckstein,
  F., & Rüegsegger, P. (2002). Estimation of distal radius failure load with
  micro-finite element analysis models based on three-dimensional peripheral
  quantitative computed tomography images. *Bone, 30*(6), 842–848.
  https://doi.org/10.1016/S8756-3282(02)00736-6

### Protocolo micro-CT/FEA de referencia

- Tapia, D., González, A., Vidal, F., & Salinas, P. (2026). A specimen-based
  comparative micro-CT–FEA analysis of vertebral trabecular bone
  microarchitecture and mechanical response in two South American cervids.
  *Biology, 15*, 722.

### Herramientas

- Sullivan, C. B., & Kaszynski, A. A. (2019). PyVista: 3D plotting and mesh
  analysis through a streamlined interface for the Visualization Toolkit (VTK).
  *Journal of Open Source Software, 4*(37), 1450.
  https://doi.org/10.21105/joss.01450
- Si, H. (2015). TetGen, a Delaunay-based quality tetrahedral mesh generator.
  *ACM Transactions on Mathematical Software, 41*(2), Artículo 11.
  https://doi.org/10.1145/2629697
- Attene, M. (2010). A lightweight approach to repairing digitized polygon
  meshes. *The Visual Computer, 26*(11), 1393–1406.
  https://doi.org/10.1007/s00371-010-0416-3
- Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0:
  Fundamental algorithms for scientific computing in Python. *Nature Methods,
  17*(3), 261–272. https://doi.org/10.1038/s41592-019-0686-2
- Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array
  programming with NumPy. *Nature, 585*(7825), 357–362.
  https://doi.org/10.1038/s41586-020-2649-2
- van der Walt, S., Schönberger, J. L., Nunez-Iglesias, J., Boulogne, F.,
  Warner, J. D., Yager, N., Gouillart, E., & Yu, T. (2014). scikit-image: Image
  processing in Python. *PeerJ, 2*, e453. https://doi.org/10.7717/peerj.453

---

## Qué NO es esta carpeta

No sustituye a la suite de pruebas (`tests/`, bloques 01–09 contra
soluciones cerradas) ni a la validación contra la app de MATLAB
(`validar_ajuste.py`) ni a la comparación con BoneJ. Aquellas comprueban que
el código hace bien **sus** cuentas; ésta comprueba que **el método es el
publicado**. Son preguntas distintas y hacen falta las dos.
