# Ladrillos mejorados como alternativa a TET10: estudio de viabilidad

Fecha: 2026-09-30. spinpy V1.0.2, solver propio (AMG con modos rígidos,
`resistencia`) con rigidez por elemento. Código: `hex_mejorado.py`,
`giroide.py`, `estudio.py`. Registros: `resultados/`. No se usó FEBio (no está
en el contenedor); los tiempos se midieron con procesos concurrentes y solo
sirven para comparar variantes entre sí.

## 1. Planteamiento

TET10 resuelve dos defectos del ladrillo de vóxel (superficie escalonada y
rigidez excesiva a flexión del elemento trilineal), pero multiplica los grados
de libertad: en el VOI proximal de H4 a 48³ son 1,0 M GDL y unos 7 GB, y a 32³
la malla suave sale 2,3 veces más blanda que la de ladrillos porque el
suavizado estrecha los puntales (INFORME de esta carpeta, §5). TET4 no es
alternativa: es aún más rígido a flexión que el hex8.

Los dos defectos se pueden atacar por separado sin cambiar la malla de
ladrillos ni su número de GDL:

| variante | geometría | elemento | GDL |
|---|---|---|---|
| hex8 | vóxeles | trilineal, Gauss 2×2×2 | los del vóxel |
| hex8I | vóxeles | modos incompatibles (Wilson-Taylor, 9 modos condensados) | los mismos |
| taubin10 | nodos de superficie suavizados (Taubin, 10 it., tope 0,5 h) | trilineal | los mismos |
| proy | nodos de superficie proyectados sobre la isosuperficie continua (Newton, tope 0,5 h) | trilineal | los mismos |
| proyI | proyección | modos incompatibles | los mismos |

La proyección necesita un campo continuo cuya isosuperficie sea la frontera: en
un spinodoide es el propio GRF (`grf.derivadas_grf` da el gradiente exacto); en
un VOI de micro-CT sería la imagen en grises interpolada al umbral de
segmentación. Esto último **no se ha probado aquí**.

## 2. Verificaciones previas

| prueba | resultado |
|---|---|
| `ke_lote` sobre un vóxel = `elastic.hex8_ke` | 1,0·10⁻¹² |
| Bloque macizo, E_app = E_s (hex8 y proyección) | exacto a 10⁻¹⁵ |
| hex8I: modos de energía nula | 6 (los rígidos) |
| hex8I: prueba de la parcela en un elemento distorsionado (energía de un campo lineal) | 3·10⁻¹⁶ |
| Cavidad, solver propio hex8 frente a FEBio hex8 (E_app y pico, n = 16) | 18 645,66 frente a 18 645,7 MPa; 1,8426 frente a 1,8426 |

Voladizo L/t = 10 con carga en la punta, flecha frente a Timoshenko
(`estudio.py viga`):

| elementos en el espesor | hex8 | hex8I |
|---|---|---|
| 1 | 0,643 | 0,985 |
| 2 | 0,869 | 0,987 |
| 4 | 0,957 | 0,990 |

Con dos elementos en el espesor, la situación de un puntal de hueso a 32³, el
ladrillo estándar es un 13 % más rígido a flexión; el de modos incompatibles,
un 1,3 %.

## 3. Giroide esquelético: rigidez con puntales finos

Giroide de tres celdas, BV/TV continuo 0,28, desfasado respecto al borde. Su
topología no cambia con la resolución (un spinodoide con estos parámetros sí
cambiaba de conectividad entre 48³ y 64³ y se descartó como referencia).
Referencia: hex8I a 128³ (Tb.Th/h ≈ 10), E_app = 688,6 MPa; hex8 a 128³ da
+0,7 %. Error de E_app frente a la referencia:

| n | Tb.Th/h | BV/TV vóxel | hex8 | hex8I | taubin10 | proy | proyI | GDL |
|---|---|---|---|---|---|---|---|---|
| 16 | 1,3 | 0,277 | +13,9 % | −0,3 % | +1,6 % | −1,2 % | −12,9 % | 9 639 |
| 20 | 1,6 | 0,286 | +2,5 % | −7,1 % | +13,2 % | +6,8 % | −2,2 % | 16 486 |
| 24 | 1,9 | 0,273 | +19,2 % | +10,8 % | +18,6 % | +9,6 % | **+2,1 %** | 24 273 |
| 32 | 2,6 | 0,283 | +5,6 % | +0,1 % | +14,3 % | +5,9 % | **+1,4 %** | 51 863 |
| 40 | 3,2 | 0,279 | +2,6 % | −1,3 % | +8,7 % | +5,3 % | **+2,2 %** | 90 636 |
| 48 | 3,9 | 0,277 | −6,4 % | −9,3 % | +0,2 % | +2,0 % | **−0,3 %** | 146 026 |
| 64 | 5,2 | 0,281 | +0,7 % | −1,1 % | +5,1 % | +3,8 % | **+2,4 %** | 314 609 |
| 96 | 7,7 | 0,283 | +3,1 % | +2,2 % | | | | 956 460 |
| 128 | 10,3 | 0,280 | +0,7 % | ref. | | | | 2 127 649 |

Lectura:

1. **El ladrillo de vóxel oscila entre −6 % y +19 %** hasta 48³, y la
   oscilación no sigue al BV/TV de los vóxeles (a 24³ tiene el menor BV/TV y es
   el más rígido): la dominan los puentes y estrangulamientos que la
   voxelización crea o destruye en los nudos. Los modos incompatibles quitan el
   sesgo de flexión (unos 4 a 5 puntos) pero no esa oscilación.
2. **Suavizar (Taubin) no ayuda**: rellena rincones cóncavos y rigidiza
   (+9 a +19 % entre 20³ y 40³), y a 16³ pierde un 8 % de volumen.
3. **Proyectar sobre la superficie verdadera elimina la oscilación**; lo que
   queda es un sesgo rígido de +4 a +7 % que coincide con el del elemento
   trilineal a flexión.
4. **Las dos correcciones juntas (proyI) quedan dentro de ±2,4 % desde 20³**
   (Tb.Th/h ≈ 1,6), con los mismos GDL que el ladrillo. A 16³ (Tb.Th/h ≈ 1,3)
   falla: el tope de 0,5 h y el control de jacobiano impiden alcanzar la
   superficie y la malla pierde un 12 % de volumen.

La pérdida de volumen de proyI frente al BV/TV continuo es −3,7 % a 20³,
−4,5 % a 24³, −1,1 % a 32³ y < 1 % desde 40³. Aun así la rigidez queda dentro
de la banda: parte del material que falta es el de los puentes espurios.

Coste: mismos GDL, iteraciones de AMG entre −4 % y +18 % respecto al hex8, y un mallado de 0,1 s a
24³ y 1,2 s a 64³. Sólo los elementos de la capa superficial necesitan matriz
propia; en vóxeles regulares hex8I sigue siendo una única matriz 24×24.

## 4. Cavidad esférica: pico de tensión

El mismo problema del §4 del INFORME de TET10 (referencia FEBio TET10 sobre
la esfera analítica: pico 2,128 σ₀, E_app 18 644,6 MPa). Pico de von Mises
en la pared, en el centro del elemento, error frente a la referencia:

| n | hex8 | hex8I | proy | proyI | TET10 (FEBio, §4) |
|---|---|---|---|---|---|
| 16 | −13,4 % | −11,9 % | −14,4 % | −12,5 % | +1,7 % |
| 24 | −8,7 % | −7,1 % | −10,4 % | −9,1 % | +5,2 % |
| 32 | −11,8 % | −11,0 % | −7,1 % | −6,3 % | −1,1 % |
| 40 | +1,5 % | +2,5 % | −1,6 % | −0,1 % | +12,1 % |
| 48 | −4,6 % | −4,3 % | −3,4 % | −2,7 % | |

E_app de proyI queda a menos de 0,05 % de la referencia desde 24³ (hex8:
−0,45 % a 24³, por la esfera vóxelizada). El pico de proyI se acerca a la
referencia de forma casi monótona, pero desde abajo y lentamente: en el
centro de un elemento lineal se promedia el gradiente de la pared. Ninguna
variante, TET10 incluido, hace converger el pico a ±3 % en todas las
resoluciones; la conclusión del INFORME de TET10 sobre el máximo sigue en pie.

## 5. Conclusión y límites

Para rigidez y cargas de fallo basadas en percentiles, un ladrillo con la
superficie proyectada sobre la isosuperficie continua y con modos
incompatibles da, en este caso de prueba, la precisión que se buscaba en TET10
con el coste del ladrillo. Queda por demostrar:

1. **Hueso real.** La proyección sobre la imagen en grises del micro-CT no se
   ha probado; el ruido de la imagen puede exigir un filtrado previo, y ese
   filtrado es un parámetro de modelado como lo era Taubin.
2. **FEBio.** La geometría proyectada es un hex8 estándar con nodos movidos y
   FEBio la acepta tal cual. Los modos incompatibles no los encontré en la
   documentación de FEBio (sólo reglas de integración para hex8); sí existen
   en Abaqus (C3D8I) y ANSYS (SOLID185 con deformación mejorada). Sin ellos,
   en FEBio se obtiene la columna «proy», no «proyI».
3. **Un caso.** Un giroide y una cavidad; hace falta al menos el VOI de H4 y
   un spinodoide ajustado antes de integrarlo en `spinpy`.
