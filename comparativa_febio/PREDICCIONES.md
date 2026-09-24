# Predicciones antes de medir — spinpy frente a FEBio

**Este documento se escribe ANTES de correr FEBio sobre ningún VOI, y no se toca
después.** Mismo protocolo que `comparativa_bonej/PREDICCIONES.md`: una
discrepancia fuera de lo declarado aquí es un hallazgo, no algo que se explica
a posteriori.

La única corrida previa fue la prueba de humo del escritor sobre un bloque
macizo de 4×4×6 a 1 MPa, necesaria para comprobar que FEBio lee el archivo.
Dio E_app = 20 000 MPa con una no linealidad geométrica de 3·10⁻⁵ (tensión de
Cauchy −0,99997 MPa: el área lateral crece con el Poisson). Ese número fija la
carga que se usa abajo; no es un resultado de la comparación.

---

## Qué se compara, y qué NO

Es la comparación de la **mecánica**, que BoneJ no toca. Se enfrenta
`resistencia.ensayo_compresion` (el ensayo de compresión en z de la app) con
FEBio 4.5 resolviendo el **mismo problema discreto**:

| | spinpy | FEBio |
|---|---|---|
| malla | un hexaedro por vóxel portante | la misma, escrita por `escribe.escribir_febio` |
| elemento | hex8 trilineal, Gauss 2×2×2 (`hex8_ke`) | hex8, GAUSS8 (su valor por omisión) |
| material | elástico lineal, E_s = 20 GPa, ν = 0,3 | `isotropic elastic` (St. Venant–Kirchhoff) |
| carga | fuerzas nodales, ¼ por nodo y cara del techo | presión no seguidora integrada por FEBio |
| apoyo | deslizante: uz = 0 en la base, 3 GDL en dos esquinas | los mismos nodos, los mismos GDL |
| solver | AMG con modos rígidos + CG, residuo 1e-8 | Newton completo con Pardiso (directo) |
| tensión | en el centro del elemento | media de los 8 puntos de Gauss |

Lo que **no** se compara: la homogeneización periódica (`elastic.homogeneizar`).
FEBio no trae condiciones periódicas listas para una malla de vóxeles, y
escribirlas a mano sería comprobar nuestra propia implementación de las
restricciones. Queda fuera, y se declara.

**Por qué el problema es el mismo y no uno parecido.** Tres cosas podrían
diferir y se han comprobado sobre el papel antes de medir:

- La media de las tensiones en los 8 puntos de Gauss de un hexaedro
  rectangular es **exactamente** el valor en el centro: la deformación de un
  elemento trilineal es multilineal en las otras coordenadas, y sus términos
  impares se anulan en la cuadratura simétrica.
- La presión constante sobre una cara bilineal da, integrada de forma
  consistente, un cuarto de la fuerza a cada nodo: el mismo reparto por área
  tributaria que spinpy aplica a mano.
- St. Venant–Kirchhoff se reduce a la elasticidad lineal con error del orden
  de la deformación del tejido. Por eso se carga con **σ = 1 kPa** en FEBio y
  se reescala ×1000 (el problema es lineal): la no linealidad baja a ~3·10⁻⁸.
  spinpy corre con su 1 MPa por omisión.

Así que **no se espera ninguna diferencia de modelo**. Lo que queda es precisión
numérica: el residuo 1e-8 del solver iterativo de spinpy (todos los casos
reales superan `UMBRAL_DIRECTO` = 6000 GDL, así que ninguno va por LU), la no
linealidad residual de FEBio y las 12 cifras con que FEBio escribe su logfile.

---

## Casos

| caso | por qué |
|---|---|
| bloque macizo 8×8×12 | solución exacta: E_app = E_s, σ_zz uniforme |
| VOI H4 proximal, medio y distal, a n = 32 y n = 48 | los mismos de la comparación con BoneJ; BV/TV 0,28 / 0,54 / 0,77; 48 es la resolución que el estudio de convergencia pide para citar |
| espinodoide ajustado a H4 proximal, n = 48 | el objeto del problema abierto (β = 15π, N = 800, θ = (15, 15, 60), ρ = 0,30848, semilla 20260720, escalado a Lc = 4,994 mm) |
| H4 proximal n = 48 con apoyo empotrado | la otra condición de contorno que ofrece la app |

Remuestreo con `elastic.remuestrear_bw` (vecino más próximo, el de la app).

---

## Predicciones

### Bloque macizo — exacto en los dos

E_app = 20 000 MPa en spinpy y en FEBio, **|Δ| < 1e-6 relativo** frente al
valor exacto. σ_zz = −σ_app en todos los elementos con la misma tolerancia. Si
esto falla, el escritor o la lectura están mal y no se sigue.

### Equilibrio — la carga que aplica FEBio es la de spinpy

Σ Rz en la base de FEBio = F_total de spinpy (σ_app · A_bruta), **|Δ| < 1e-9
relativo**. Comprueba que la presión integrada por FEBio sobre el área ósea del
techo suma exactamente la fuerza de la sección bruta.

### E_app — coincidencia a la precisión del solver

**|Δ| / E_app < 1e-6** en los siete casos. E_app es la media de uz en el techo:
una magnitud global, poco sensible al error local del iterativo.

**Qué lo falsaría:** un |Δ| > 1e-3 no es precisión numérica, es otro problema
—un GDL mal fijado, una esquina de apoyo distinta, una cara del techo de más o
de menos—. Sería un error del escritor o de spinpy, y habría que encontrarlo.

### Campo de desplazamientos

**max |Δu| / max |u| < 1e-5**, nodo a nodo. Los nodos se escriben en el orden
de `ensayo_compresion`; si no coincidiera, este número saldría de orden 1, así
que también verifica la correspondencia.

### Tensiones por elemento

**max |Δσ| / max(von Mises) < 1e-4**, componente a componente sobre todos los
elementos. Es la tolerancia más holgada porque el error local del iterativo se
concentra en los puntales finos, precisamente donde está el pico. Si E_app
coincide y esto no, el hallazgo es que el residuo 1e-8 de spinpy no basta para
el campo local, no un error de modelo.

### Estadísticos citables

`vm_p99_superficie`, factor de Pistoia (y con él σ_fallo y F_fallo),
calculados por las MISMAS funciones de spinpy sobre los dos campos:
**|Δ| relativo < 1e-5**. Se usan las mismas funciones a propósito: lo que se
contrasta es la solución de EF, no el cálculo de un percentil.

`vm_max` se reporta pero **no se predice**: es un único elemento, y es
justamente el que el estudio de convergencia declaró no citable.

---

## Lo que esta comparación NO puede demostrar

Coincidir con FEBio prueba que spinpy resuelve bien el problema que plantea. No
prueba que ese problema sea el correcto: el mismo mallado de vóxeles, el mismo
material homogéneo y los mismos apoyos en los dos lados comparten cualquier
error de modelado (el escalonado de la superficie, la falta de
convergencia del máximo, el 1e-6 del vacío que aquí ni siquiera aparece porque
solo se malla el hueso). Esa validación es la del estudio de convergencia y la
de la solución cerrada de Goodier, no esta.
