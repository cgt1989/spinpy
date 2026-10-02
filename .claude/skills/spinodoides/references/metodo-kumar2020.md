# Método de Kumar et al. (2020) y su implementación

Fuente primaria: Kumar, Tan, Zheng y Kochmann (2020), npj Computational Materials 6:73. Ver `bibliografia.md`.

## Tres ecuaciones

1. Campo gaussiano: `phi(x) = sqrt(2/N) SUM cos(beta n_i . x + gamma_i)`, `n_i ~ U(S2)`, `gamma_i ~ U[0, 2pi)`.
2. Anisotropía: `n_i` uniforme sobre `{k : |k.e1| > cos(th1) o |k.e2| > cos(th2) o |k.e3| > cos(th3)}`, con `th_j` en `{0} U [th_min, 90°]`.
3. Conjunto de nivel: sólido donde `phi <= phi0`, `phi0 = sqrt(2) erfinv(2 rho - 1)`.

Las interfaces se alinean preferentemente perpendiculares a las `n_i`. De ahí: cono alrededor de un eje, láminas apiladas a lo largo de ese eje, eje blando.

## Correspondencia con el código

| Artículo | spinpy |
|---|---|
| ec. (2), unión de conos | `grf.py::_waves_rechazo` (candidato isótropo aceptado si el ángulo a algún eje es menor que su theta) |
| ec. (3), umbral | `grf.py::level_set` (erfinv exacta de scipy; la aproximación de Winitzki, error ~2e-3, existe solo para reproducir el repositorio del profesor) |
| theta_min = 15° | límite inferior de los presets de ángulo del ajuste |

Convención de ejes: la dimensión 1 de la máscara es X, como `readVTKVOI` en AppFinal_V2.m.

## Dos convenciones de muestreo (declarar siempre)

- `rechazo` (GIBBON `spinodoid.m:178-206`): reparto proporcional al ángulo sólido. Con (15, 45, 0) el cono ancho recibe ~90 % de las ondas.
- `equitativo` (TPMS-Scaffolds-generator, `TPMS.py:223`): `num_waves // n_conos` por cono, 50/50 aunque los conos difieran.

Coinciden solo con conos iguales, un cono o caso isótropo. Con la terna columnar (0, 30, 30) la diferencia en DA fue 2.2 % (coinciden por construcción); con (15, 45, 90) fue 27 %.

## Ternas de ángulos y su procedencia

| Clase | Terna (grados) | Fuente |
|---|---|---|
| lamelar | (30, 0, 0) | Kumar 2020, leyenda Fig. 3 (citado) |
| columnar | (0, 30, 30) | Kumar 2020, leyenda Fig. 3 (citado) |
| cúbica | (30, 30, 30) | deducida de la regla de clases (Kumar no la publica) |
| isótropa | (90, 90, 90) | Kumar 2020, texto p. 4 (citado) |
| lamelar | (0, 0, 15) | Zheng 2021, pie Fig. 1b (citado) |
| columnar | (15, 15, 0) | Zheng 2021, pie Fig. 1c (citado) |
| cúbica | (15, 15, 15) | Zheng 2021, pie Fig. 1d (citado) |
| isótropa | (90, 90, 90) | Zheng 2021, pie Fig. 1e (citado) |
| espinodal de Guo | (60, 30, 10), beta = 15pi, Q = 1000, rho = 0.3 | Guo 2024, Fig. 7b (citado) |

Las ternas de Zheng están traspuestas respecto de las de Kumar (cono activo e3 en vez de e1), lo que hace independiente la comprobación "cono activo, eje blando".

Parámetros citados: rho = 0.5 (leyenda Fig. 2), beta = 15pi (leyenda Fig. 5, "comparación visual"). El beta de la Fig. 2 no está publicado, por eso las réplicas salen más finas que la figura; se mantiene el valor citable.

## Otros generadores del proyecto

- **dual-lattice** (`dual_lattice.py`, Vafaeefar et al.): red dual de una triangulación de Delaunay perturbada, conectividad nodal 4, puntales de sección redonda (GIBBON usa sección triangular). El radio sale de la densidad (cuantil rho de la distancia al esqueleto), como el umbral erfinv en el spinodoide. Vafaeefar comparan gyroide, spinodoide y dual-lattice como modelos de hueso trabecular y reportan para E3 (MPa): hueso 1188, gyroide 887, spinodoide 714, dual-lattice 1313, con BV/TV, Tb.Th y DA igualados. Verificar la cifra en el artículo antes de citarla.
- **TPMS** (giroide, Schwarz): periódicas, a diferencia del spinodoide.
