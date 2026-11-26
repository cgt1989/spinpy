# Arquitectura de spinpy

Mapa para quien abre el código por primera vez. Todo el código está en español
(identificadores, comentarios, docstrings) y cada módulo empieza con una
cabecera que explica **por qué** existe y qué decisiones lo moldearon, no solo
qué hace. Léela antes de tocar el módulo: muchas líneas que parecen
redundantes son correcciones ya pagadas.

## Capas

Los módulos de cálculo de `spinpy/` no importan Qt. La interfaz (`visor.py` y
los `dialogo_*.py`) importa el paquete; el paquete nunca importa la interfaz.
Se puede usar el núcleo desde un script o un cuaderno sin arrastrar la GUI.
Hay dos excepciones declaradas, que son código de interfaz y viven en el
paquete por historia: `pantalla.py` (diálogos que se ajustan a la pantalla) e
`idioma.py` (traducción de la interfaz, que importa Qt dentro de una función).

```
 visor.py, dialogo_*.py            interfaz (PyQt5, GPL)         ← depende de todo
 ─────────────────────────────────────────────────────────────
 informe.py, figuras*.py, estudios.py, lista_chequeo.py, md2pdf.py
                                   informe para publicación
 ─────────────────────────────────────────────────────────────
 fit.py, fit_dual.py, metodos.py, lote.py, simulacion.py, fem.py
                                   flujos: ajustar, lote, simular, ensayar
 ─────────────────────────────────────────────────────────────
 morphometry.py, curvatura.py, espesor.py, elipsoide.py, elastic.py,
 resistencia.py, cargas.py, error.py, estadistica.py, incertidumbre.py
                                   medir y resolver (numérico puro)
 ─────────────────────────────────────────────────────────────
 grf.py, dual_lattice.py           generadores de microestructura
 voi.py, io.py                     entrada: pilas TIFF, .vtk, .mat, recorte
 motores/                          motores FEM intercambiables
```

## Qué hay en la raíz del repositorio

| Carpeta o archivo | Qué es |
|---|---|
| `spinpy/` | el paquete (núcleo, informe, motores FEM) |
| `visor.py`, `dialogo_*.py` | la interfaz gráfica; se lanza con `python visor.py` |
| `validacion/` | guiones que contrastan el código con MATLAB y con soluciones analíticas (`python validacion/validar_error.py`); escriben en `resultados/` |
| `tests/` | la suite de verificación |
| `Test/` | réplicas de las figuras de los artículos en que se basa el método |
| `comparativa_*/` | campañas de comparación con FEBio y entre motores (informes y guiones) |
| `ejemplos/`, `datos/` | cuaderno e inicio rápido, y los diez VOIs porcinos abiertos |
| `instalador/` | construcción de los ejecutables de Windows, macOS y Linux |
| `docs/` | manual, guías y medios del README |
| `matlab/`, `resultados/` | lado MATLAB de la validación y sus salidas de referencia |
| `idioma_revisar.py` | comprobador del diccionario de traducción de la interfaz |

## Dónde está cada cosa

| Si quieres tocar… | Mira en | Prueba |
|---|---|---|
| Leer VOIs, umbralizar, recortar (cubo o cilindro) | `voi.py`, `io.py` | 08, 31, 33, 34, 35 |
| Generador espinodal (campo aleatorio gaussiano) | `grf.py` | 05, 11, 16-18, 30 |
| Generador dual-lattice | `dual_lattice.py` | 12, 13 |
| BV/TV, Tb.Th, DA, Conn.D, SMI | `morphometry.py` | 02, 03, 10, 18, 35 |
| Espesor local, curvatura, Ellipsoid Factor | `espesor.py`, `curvatura.py`, `elipsoide.py` | 04, 09, 11 |
| Homogeneización elástica | `elastic.py` | 05, 17 |
| Ensayo de compresión, von Mises, Pistoia | `resistencia.py` | 05, 06, 20, 26 |
| Cargas y apoyos a elegir, barrido de ángulo | `cargas.py` (solver y regiones), `escena_cargas.py` (3D y clic), `dialogo_cargas.py` (ventana), `informe_cargas.py` (sección del informe) | 37 |
| **Motores FEM** (NGSolve, scikit-fem…) | `motores/` | 29, 36 |
| Mallado TET10, exportar a Abaqus/ANSYS/FEBio | `fem.py`, `solido.py`, `escribe.py` | 07, 25, 27 |
| Ajuste de parámetros a un VOI | `fit.py`, `fit_dual.py`, `error.py` | 13, 17 |
| Informe, citabilidad, lista de chequeo | `informe*.py`, `lista_chequeo.py` | 19, 22, 23 |
| Textos de la interfaz en inglés | `idioma_textos.py` | `idioma_revisar.py` |

## Convenciones que no se negocian

- **Ejes**: las máscaras salen con la dimensión 1 = X. `imread` devuelve
  (fila, columna) = (y, x), así que cada rebanada se traspone al apilar.
  Mezclar las dos convenciones es el error B1/B2 del historial del proyecto.
- **Una sola medición**: `morphometry.morfometria` es el único sitio donde se
  calculan BV/TV, BS/BV, Tb.Th, Tb.Sp, Tb.N y DA, porque el sesgo del método
  solo se cancela si el VOI y el candidato recorren el mismo camino
  (corrección C4). Una métrica nueva pasa por ahí. Para un VOI que no es un
  cubo, la medida pasa por `morfometria_dominio`, que la restringe a su
  dominio.
- **Unidades internas**: mm, mm², mm³, MPa, N.
- **El signo de los ejes PCA es arbitrario.** Nada debe depender de él; los
  nombres proximal y distal los pone quien conoce la pieza.
- **Los números llevan su procedencia** (`procedencia.py`): versión, familia,
  semilla y número de onda, para que un resultado se pueda regenerar.

## Cómo se prueba

`tests/` tiene un bloque por tema (`test_NN_tema.py`). No son pruebas de
regresión ("el número no ha cambiado"): son verificaciones contra respuestas
conocidas de antemano (invariantes topológicos, soluciones cerradas,
definiciones publicadas), con la **tolerancia declarada antes de medir** y su
justificación en la cabecera del módulo. Un fallo es un hallazgo: se explica y,
si se acepta, se marca `xfail(strict=True)`. Ver `tests/conftest.py`.

## La interfaz

`visor.py` es una sola ventana grande por razones históricas (nació como
traducción de una app de MATLAB). Sus partes, en orden de aparición:
utilidades de render, diálogos (`Dialogo*`), la ventana principal `Visor` y su
panel de controles. Toda la lógica numérica vive en el paquete; la interfaz
solo orquesta, lanza el trabajo en un `Trabajador` (hilo) y pinta. Si vas a
añadir una función de cálculo, escríbela en el paquete con su prueba y deja en
`visor.py` solo el botón.
