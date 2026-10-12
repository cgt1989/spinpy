# Tiempos de ejecución de los motores de elementos finitos de spinpy

Fecha: 2026-10-12. spinpy V2.1.1. Banco de pruebas: `comparativa_motores/correr.py` (opción `--salida`). Datos: `comparativa_motores/tiempos_sesion/*.jsonl`. Tablas completas: `tablas_tiempos.md`. Figuras: `t_escala_hex.png`, `t_escala_tet.png`, `t_etapas_hex64.png`, `t_etapas_tet32.png`. Generador: `comparativa_motores/informe_tiempos.py`.

## Resumen

**Objetivo.** Medir cuánto tarda cada motor de elementos finitos integrado en spinpy (la app, NGSolve, FEniCSx, scikit-fem y SfePy) en los mismos problemas discretos de la comparativa publicada en `comparativa_motores/INFORME.md`, y contrastar la ordenación con esa campaña.

**Métodos.** Se repitió la campaña completa (95 corridas) con las mismas mallas y resolvedores: bloque macizo, cavidad esférica, espinodoide con hexaedros (24³ a 80³, de 30 672 a 661 965 grados de libertad, GDL), espinodoide con TET10 (314 880 y 790 599 GDL) y problemas no lineales con St. Venant-Kirchhoff y neo-Hookeano. Cada corrida se ejecutó en un proceso hijo, con 4 hilos y de forma secuencial, y se registró el tiempo de reloj de pared por etapas y la memoria residente de pico.

**Resultados.** Con malla de hexaedros, los resolvedores directos de NGSolve y FEniCSx fueron los más rápidos (9,2 y 9,9 s a 365 655 GDL, frente a 47,5 s de la app y 130 s de scikit-fem con su resolvedor iterativo). Con TET10 a 790 599 GDL, FEniCSx con MUMPS resolvió en 18,5 s y NGSolve con Cholesky propio en 32,0 s, mientras que scikit-fem y SfePy necesitaron 673 y 806 s. En el espinodoide no lineal a 24³, FEniCSx tardó 3,2 s (St. Venant-Kirchhoff), NGSolve 6,3 s y scikit-fem 110 s; SfePy no convergió. Respecto de la campaña publicada, los tiempos de esta sesión fueron menores (mediana del cociente 0,73; percentiles 10 y 90 de 0,62 y 0,80, n = 46 pares), y la ordenación entre motores se mantuvo salvo en un punto: con hexaedros, el directo de NGSolve y el de FEniCSx quedaron prácticamente empatados.

**Conclusión.** La ordenación de motores por tiempo es estable entre campañas y equipos del mismo tipo; las cifras absolutas no lo son. Para malla de hexaedros lineal, el resolvedor propio de la app es competitivo con los iterativos de las bibliotecas, pero queda por detrás de cualquier factorización directa multihilo mientras esta quepa en memoria.

## 1. Diseño

### 1.1 Equipo y versiones

Equipo: Intel Xeon a 2,10 GHz, 4 CPU virtuales, 15 GB de memoria, Linux 6.18, Python 3.13 (FEniCSx con el Python 3.x de su entorno de conda-forge). Versiones: spinpy 2.1.1, NGSolve 6.2.2608, scikit-fem 12.0.2, SfePy 2026.3, FEniCSx 0.11.0, numpy 2.5.3, scipy 1.18.1, pyamg 5.3.0. La configuración de hilos fue `OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 4`.

### 1.2 Problemas

Son los casos de `comparativa_motores/casos.py`, sin modificaciones, descritos en la sección 2.2 del informe de motores: bloque macizo (hexaedros y TET10), cavidad esférica (hexaedros a 32³ y TET10), espinodoide con hexaedros a 24³, 32³, 48³, 64³ y 80³, espinodoide con TET10 a 32³ y 48³, y espinodoide no lineal a 24³ con plato rígido hasta el 2 % de deformación. Todos los motores reciben el mismo problema discreto (`spinpy.motores.problema_de_malla`).

### 1.3 Medida

El tiempo informado es el `tiempo_total_s` de cada motor: malla propia del motor, compilación JIT, montaje, resolución y extracción de resultados. No incluye la importación de la biblioteca ni la construcción de la malla de spinpy, que es común a todos. La memoria es el pico residente (`VmHWM`) del proceso hijo. Un vigilante detenía cada corrida por encima de 13 GB o de una hora. Cada combinación de caso, motor y resolvedor se ejecutó **una sola vez**.

Resolvedores comparados: la app usa su CG con multigrid algebraico (pyamg) y modos rígidos. Los demás motores se ejecutaron con un resolvedor directo (NGSolve: Cholesky propio; FEniCSx: MUMPS; scikit-fem y SfePy: SuperLU) y con uno iterativo (CG con multigrid algebraico o BDDC). Los directos de scikit-fem y SfePy solo se ejecutaron hasta 1,5·10⁵ GDL, y el directo de los demás hasta 4·10⁵ GDL con hexaedros y 1,2·10⁶ con TET10, que son los límites del banco de pruebas.

## 2. Resultados

### 2.1 Espinodoide con hexaedros

Tiempo total del motor en segundos. `n/d`: no ejecutado por el límite de GDL del resolvedor directo.

| Malla | GDL | App | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 24³ | 30 672 | 2,40 | 0,69 | 3,78 | 1,13 | 2,49 | 4,97 | 5,91 | 4,04 | 3,10 |
| 32³ | 62 595 | 5,33 | 1,01 | 6,61 | 1,72 | 5,20 | 11,3 | 12,5 | 10,0 | 6,68 |
| 48³ | 173 484 | 24,9 | 3,80 | 32,0 | 4,35 | 21,4 | n/d | 61,9 | n/d | 32,5 |
| 64³ | 365 655 | 47,5 | 9,21 | 50,9 | 9,91 | 43,3 | n/d | 130,0 | n/d | 55,1 |
| 80³ | 661 965 | 102,3 | n/d | 107,3 | n/d | 78,9 | n/d | 291,9 | n/d | 115,6 |

![Espinodoide con hexaedros: tiempo frente a tamaño](t_escala_hex.png)

Con los resolvedores iterativos, NGSolve, FEniCSx y SfePy tardaron entre 0,8 y 1,6 veces lo que la app, y scikit-fem entre 2,3 y 2,9 veces, y el exponente de escala frente a los GDL estuvo entre 1,1 y 1,3 en todos. Con resolvedor directo, a 365 655 GDL NGSolve y FEniCSx tardaron la quinta parte que la app. El exponente del directo de NGSolve fue 1,09 (4 puntos) y el de FEniCSx 0,88 (4 puntos); con tan pocos puntos y con un coste fijo que domina en los tamaños pequeños, estos exponentes no deben leerse como una ley asintótica.

Memoria de pico a 365 655 GDL: app 1,4 GB; NGSolve 1,5 GB (directo); FEniCSx 1,7 GB (directo); SfePy 1,5 GB (iterativo); scikit-fem 3,0 GB (iterativo).

### 2.2 Espinodoide con TET10

| Malla | GDL | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem ite. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|
| 32³ | 314 880 | 8,73 | 11,6 | 5,71 | fallo | 194,4 | 184,9 |
| 48³ | 790 599 | 32,0 | 43,1 | 18,5 | 225,0 | 672,7 | 806,1 |

![Espinodoide con TET10: tiempo frente a tamaño](t_escala_tet.png)

La app no resuelve TET10. Con 790 599 GDL, el directo de FEniCSx fue 1,7 veces más rápido que el de NGSolve y 36 veces más rápido que el iterativo de scikit-fem (673 s), con una memoria de pico de 3,6 GB frente a 3,2 GB de NGSolve y 6,8 GB de scikit-fem. El iterativo de FEniCSx (CG con GAMG) no convergió en la malla de 314 880 GDL (PETSc KSP, motivo -8) y necesitó 225 s en la de 790 599 GDL, lo cual es coherente con la divergencia descrita en el informe de motores.

![Etapas del tiempo en el espinodoide TET10 a 32³](t_etapas_tet32.png)

La figura de etapas muestra que, en scikit-fem, el montaje (unos 40 s) pesa tanto como el 20 % del total, y que en los motores iterativos la resolución domina el tiempo (más del 75 %).

### 2.3 Problemas pequeños y cavidad esférica

| Caso | GDL | App | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bloque (hex) | 3 159 | 0,20 | 0,32 | 0,33 | 0,72 | 0,74 | 0,89 | 0,90 | 0,89 | 0,91 |
| bloque (TET10) | 13 446 | n/d | 0,42 | 0,51 | 0,85 | 1,36 | 1,84 | 1,93 | 2,91 | 1,71 |
| cavidad (hex 32³) | 105 558 | 7,59 | 26,7 | 5,97 | 5,17 | 6,32 | 316,8 | 28,2 | 591,2 | 7,71 |
| cavidad (TET10) | 206 142 | n/d | 27,4 | 31,7 | 7,80 | 26,0 | n/d | 46,3 | n/d | 35,8 |

En los problemas de unos pocos miles de GDL, la app fue el motor más rápido (0,20 s), porque el resto paga un coste fijo de arranque (montaje de formas, compilación JIT en FEniCSx). En la cavidad con hexaedros, el directo de NGSolve (26,7 s) fue más lento que su propio iterativo (5,97 s), y los directos SuperLU de scikit-fem y SfePy (317 y 591 s) fueron entre 11 y 77 veces más lentos que los iterativos de esos motores, un efecto del relleno de SuperLU en 3D sin ordenación de disección anidada.

### 2.4 Problemas no lineales

Newton con resolvedor directo, plato rígido. Tiempo total en segundos.

| Caso | GDL | NGSolve | FEniCSx | scikit-fem | SfePy |
|---|---:|---:|---:|---:|---:|
| bloque, St. Venant-Kirchhoff | 525 | 0,31 | 2,81 | 2,42 | 1,77 |
| bloque, neo-Hookeano | 525 | 0,33 | 2,51 | 4,72 | n/d |
| espinodoide 24³, St. Venant-Kirchhoff | 30 672 | 6,26 | 3,23 | 110,4 | fallo |
| espinodoide 24³, neo-Hookeano | 30 672 | 7,57 | 3,78 | 175,0 | n/d |

En el bloque de 525 GDL, el tiempo de FEniCSx está dominado por la compilación JIT de las formas (1,8 a 2,1 s de 2,5 a 2,8 s totales); esa compilación solo ocurre la primera vez que se usa cada forma, y después queda en caché en disco. Con 30 672 GDL, FEniCSx fue el más rápido (3,2 y 3,8 s) y scikit-fem el más lento (110 y 175 s), 34 y 46 veces más que FEniCSx. El Newton de SfePy no convergió en el espinodoide a la carga del 2 % (30 iteraciones), y en el bloque necesitó entre 22 y 27 iteraciones por paso de carga frente a 3 de los demás motores, lo que confirma que su tangente no es consistente con el residuo.

### 2.5 Reproducibilidad frente a la campaña publicada

En los 46 pares de corridas del espinodoide (hexaedros y TET10) que existen en ambas campañas, el tiempo de esta sesión fue menor que el publicado, con una mediana del cociente de 0,73 (percentiles 10 y 90: 0,62 y 0,80). El cociente más bajo (0,44 a 0,49) correspondió al directo de NGSolve, y es esperable: la versión publicada usaba PARDISO (MKL) como directo por omisión, mientras que NGSolve instalado desde PyPI ya no incluye MKL y esta sesión usó su Cholesky propio. No se midió la causa del resto de la diferencia (otro procesador, otras versiones de las bibliotecas). Dado que se trata de una única corrida por combinación, no se estimó la variabilidad entre repeticiones.

## 3. Discusión

La comparación de tiempos depende más del resolvedor que del motor. Con el mismo tipo de resolvedor, las diferencias entre NGSolve, FEniCSx y SfePy son modestas (entre 0,8 y 1,6 veces el tiempo de la app en los iterativos con hexaedros); las diferencias grandes aparecen entre factorizaciones directas multihilo (NGSolve, FEniCSx) y los resolvedores secuenciales de SciPy (SuperLU en scikit-fem y SfePy), y con TET10 también entre el iterativo de SfePy y scikit-fem (un orden de magnitud más lentos) y los directos de NGSolve y FEniCSx. La ventaja práctica de los directos tiene como contrapartida la memoria, que crece más deprisa que el tiempo, de modo que con la malla de hexaedros a 80³ el banco de pruebas no ejecutó el directo.

Para la aplicación, estos datos respaldan lo ya recomendado en el informe de motores: NGSolve para TET10 y no lineal (instalable con pip en Windows, 32 s a 790 599 GDL), y el resolvedor propio para hexaedros lineales. FEniCSx fue más rápido que NGSolve en TET10 y en el no lineal de tamaño medio, pero requiere conda y un compilador de C en tiempo de ejecución, lo que limita su distribución.

## 4. Limitaciones

- Una sola corrida por combinación y un solo equipo: no se estimó la variabilidad entre repeticiones ni entre máquinas. Las diferencias menores de aproximadamente 20 % entre motores no deben interpretarse como ordenación firme.
- El tiempo de NGSolve con el directo corresponde a su Cholesky propio, no a PARDISO. La comparación con la campaña publicada mezcla, por tanto, dos resolvedores directos distintos para este motor.
- Durante la corrida de SfePy con SuperLU en la cavidad de hexaedros (591 s) se ejecutó en paralelo, durante unos 5 minutos, una prueba corta de la réplica del protocolo de Tapia (la app y NGSolve sobre un VOI de 40³). Esa medición puede estar algo inflada y no cambia la conclusión (SuperLU quedó entre 11 y 77 veces por detrás de los iterativos).
- No se ejecutaron los directos de scikit-fem y SfePy por encima de 1,5·10⁵ GDL, ni los de NGSolve y FEniCSx con hexaedros a 80³ (límite de 4·10⁵ GDL del banco de pruebas), por lo que no hay comparación en esas celdas.
- El tiempo no incluye la importación de las bibliotecas ni el mallado de spinpy. En FEniCSx, la compilación JIT de la primera ejecución se informa aparte y no se repitió.
- Esta es una medición computacional y no un estudio observacional ni clínico; no existe una guía de reporte (STROBE, CONSORT, PRISMA, etc.) que le sea aplicable, por lo que no se incluye lista de chequeo.

## 5. Reproducir

```
pip install scipy scikit-image pyamg scikit-fem ngsolve sfepy tetgen pymeshfix pyvista pandas joblib matplotlib
# FEniCSx (conda-forge): fenics-dolfinx, numpy, scipy, petsc4py, mpich, y los paquetes anteriores
export FENICSX_PYTHON=/ruta/al/entorno/fx/bin/python
cd comparativa_motores
python correr.py exactos cavidad nl espinodoide_hex espinodoide_tet --salida tiempos_sesion
python informe_tiempos.py
```
