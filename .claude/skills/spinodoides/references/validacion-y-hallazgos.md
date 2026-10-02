# Validación, hallazgos y errores ya corregidos

Las tolerancias se declaran antes de medir y no se ajustan después. Un fallo es un hallazgo y se lista primero.

## Réplica de Kumar 2020 (`Test/replicar_kumar2020.py`, ~4 min; `--rapido` ~1.5 min)

Predicciones falsables derivadas del artículo (rechazo, rho = 0.5, beta = 15pi):

- C1 lamelar: e1 blando, E1/E3 < 0.6. C2 columnar: e1 rígido, E1/E3 > 1.6. C3 isótropa: Emax/Emin < 1.25.
- C4 cúbica mucho menos anisótropa que lamelar y columnar (< 1/2 de la mayor anisotropía). C5 lamelar y columnar ordenan al revés el eje e1. C6 con conos desiguales, `equitativo` no es la ec. (2) (> 5 % en DA).

Correcciones documentadas: C4 pasaba por 1.249 frente a 1.25 (suerte); se reformuló en relativo porque una realización no tiene la simetría del ensemble. C6 falló con la columnar (2.2 %) porque con conos iguales los esquemas coinciden; se movió a conos desiguales (27 %).

## Réplica de Zheng 2021 (`Test/replicar_zheng2021.py`, ~90 min completa)

Se replica la Sección 2 y la Fig. 1(b-e), no la optimización topológica ni la red neuronal.

- Z1 superficies elásticas dentro de la esfera de Voigt (rho·Es). Z2 isótropa bajo Hashin-Shtrikman superior (E = 0.333 Es en la figura). Son cotas con fórmula cerrada: sobreviven a lo que un cociente E1/E3 no detecta (un factor de escala global).
- Z3 ortotropía: se verifica en el ensemble (promedio de cuatro realizaciones). Medido a 24³ lamelar: 4.8, 5.2, 1.6, 3.4 % sueltas y 1.8 % promediado. Una realización lamelar a 40³ dio 11.1 %.
- Z4 con sus ternas, el eje del cono activo es blando en la lamelar y rígido en la columnar.
- Z5 y Z7 percolación (abajo).

Zheng homogeniza con condiciones afines (cota superior de la rigidez); spinpy con periódicas (valor más bajo). Las comprobaciones son desigualdades.

## Percolación y regla rho >= 0.3

Zheng restringe rho >= 0.3 "para evitar dominios sólidos disjuntos" (p. 5). Medido (5 semillas, 96³, N = 2000, rho = 0.30):

| clase | desconectado | portante mínimo |
|---|---|---|
| lamelar 15° | 48.7 ± 11.4 % | 0.0 % |
| columnar | 1.8 ± 1.7 % | 97.9 % |
| cúbica | 0.3 ± 0.3 % | 99.3 % |
| isótropa | 0.7 ± 0.7 % | 98.2 % |
| lamelar 30° | 0.8 ± 1.3 % | 96.9 % |

Hallazgo: la regla se cumple con theta_min = pi/6 (texto) y falla con 15° (pies de figura) en la lamelar: a rho = 0.30, 48.4 % con 15°, 2.0 % con 30°, 0.2 % con 45°. Umbral de percolación de la lamelar ~ rho_c 0.31. Un cono estrecho produce láminas casi planas y separadas.

Usar la **fracción portante** y no "1 − mayor componente": la columnar son columnas paralelas legítimamente separadas que cargan (llegó a 13.6 % desconectado con 98 % portante). Exigir en la peor semilla, no en la media.

Dispersión: a rho = 0.30 lamelar 15° la desconexión fue 19.3 ± 37.7 % (N = 200), 1.7 ± 2.5 % (N = 400), 25.3 ± 34.5 % (N = 1000), 46.0 ± 10.1 % (N = 2000). Una serie de puntos con una semilla distinta cada uno parece tendencia y no lo es. La resolución es estable (128³, 160³, 192³ dan 48.4, 48.3, 48.4 %). El efecto solo existe en el régimen N ~ 2000, por eso el bloque de dispersión no sigue al modo rápido.

## Réplica de Guo 2024 (`Test/replicar_guo2024.py`, ~4 min)

Curvaturas principales (k1 >= k2) de la interfaz, ponderadas por **área**, no por elemento (ec. 3 de Guo). Convención: normal hacia el vacío; bola de hueso k = +1/R, poro −1/R, trabécula silla (k1 > 0 > k2), superficie mínima H = 0. Dos caminos: exacto (gradiente y hessiano del GRF, `grf.derivadas_grf`) y discreto (segunda forma fundamental por triángulo, Rusinkiewicz 2004).

- R1 perfil espinodal en silla dentro de la caja leída del panel b (k1 en [10, 60], k2 en [−40, 10]). R2 doblar beta dobla todas las curvaturas (2.00 ± 5 %). R3 a rho = 0.5 la curvatura media del ensemble es nula. R4 y R5 estimador discreto frente a fórmula cerrada (espinodoide y PNS sin x sin 1.8y + sin y sin 1.8z + sin z sin 1.8x = 0.5), error mediano <= 10 %. R6 la PNS no es superficie mínima. R7 marching cubes sobre máscara binaria sobreestima el área (1.02 a 1.20; +8.8 % en espinodoide, +8.5 % en esfera).
- La unidad de curvatura es 1/(lado del dominio): **deducida**, no citada. El criterio de R4 queda justo en el límite con mallas pequeñas; usar la ejecución completa.
- No se replica el diseño inverso por redes neuronales ni la muestra de hueso de Tozzi et al.
- Las curvaturas distinguen una red conectada de islas con igual BV/TV, Tb.Th y DA.

## Verificación contra soluciones cerradas (`docs/validacion_literatura/INFORME.md`)

61 de 66 comprobaciones dentro de tolerancia. No pasaron: espesor local en losas (+33 % a +7 %, sesgo de medio vóxel) y el exponente de Gibson-Ashby n = 3.88 fuera de [1, 3] (medido con tres densidades, rho = 0.262, 0.35 y 0.509; el informe lo lista como hallazgo y no establece su causa). Invariantes comprobados: Euler y Conn.D (bola, toro, retículo), SMI de losa, cilindro y esfera, Hooke en cubo macizo, Pistoia (sigma_fallo = eps_crit·Es), Voigt y Hashin-Shtrikman en spinodoides, estanqueidad de malla.

## Cómo reportar una réplica

Distinguir parámetros citados, deducidos y propios. Declarar diferencias visuales y su causa. Dar la cifra cruda aunque pase por poco. Si una comprobación falló y se reformuló, decirlo y explicar por qué la prueba y no el código era el error.
