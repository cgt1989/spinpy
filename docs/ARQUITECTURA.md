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
 fit.py, fit_dual.py, fit_levelset.py, metodos.py, lote.py, simulacion.py,
 fem.py
                                   flujos: ajustar, lote, simular, ensayar
 ─────────────────────────────────────────────────────────────
 morphometry.py, curvatura.py, espesor.py, elipsoide.py, elastic.py,
 resistencia.py, cargas.py, error.py, estadistica.py, incertidumbre.py
                                   medir y resolver (numérico puro)
 ─────────────────────────────────────────────────────────────
 grf.py, dual_lattice.py, level_set.py
                                   generadores de microestructura
 familias.py                       las tres familias: códigos y claves
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
| Generador level-set foam (morphoGen, opcional) y su ajuste | `level_set.py`, `fit_levelset.py` | 38 |
| Qué familias hay y cómo se llaman sus claves de resultados | `familias.py` | 38 |
| BV/TV, Tb.Th, DA, Conn.D, SMI | `morphometry.py` | 02, 03, 10, 18, 35 |
| Espesor local, curvatura, Ellipsoid Factor | `espesor.py`, `curvatura.py`, `elipsoide.py` | 04, 09, 11 |
| Homogeneización elástica | `elastic.py` | 05, 17 |
| Ensayo de compresión, von Mises, Pistoia | `resistencia.py` | 05, 06, 20, 26 |
| Cargas y apoyos a elegir, barrido de ángulo | `cargas.py` (solver y regiones), `escena_cargas.py` (3D y clic), `dialogo_cargas.py` (ventana), `informe_cargas.py` (sección del informe) | 37 |
| **Motores FEM** (NGSolve, scikit-fem…) | `motores/` | 29, 36 |
| Mallado TET10, exportar a Abaqus/ANSYS/FEBio | `fem.py`, `solido.py`, `escribe.py` | 07, 25, 27 |
| Ajuste de parámetros a un VOI | `fit.py`, `fit_dual.py`, `fit_levelset.py`, `error.py` | 13, 17, 38 |
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

## Añadir una familia de microestructura

La tercera familia (level-set foam) se añadió sin tocar la lógica de las dos
primeras; esa es la receta. Lo que recorre las familias usa `spinpy.familias`
(`FAMILIAS`, `sufijo`, `clave`, `clave_morfometria`), no tuplas escritas a mano.

1. **Generador** `generar_X(resolution, ..., seed) -> (BW, campo, info)`, con
   `info["rho_objetivo"]` y `info["rho_obtenida"]`. `BW` es booleana, `dim 1 = X`.
2. **Ajuste** `ajustar_X(VOI, spacing, modo, m_voi, seed, ...)` que devuelve el
   mismo diccionario que `fit_dual.ajustar_dual_lattice` (claves en su
   docstring de retorno; la máscara y las métricas del candidato van en
   `mascara` y `metricas_spin`).
3. **Registro de métodos** (`metodos.REGISTROS`), `gobierna_de`, y las réplicas
   (`generar_replicas`).
4. **Traducción parámetros → generador**, en un único sitio:
   `procedencia.kwargs_generador`, `procedencia.bloque` y
   `procedencia.desde_parametros`.
5. **Lote y línea de órdenes**: `lote.PARAMETROS_FAMILIA`, `generar_desde_parametros`,
   `dispersion_semillas`, `correr_lote` y `cli --familia`.
6. **Familia en `spinpy/familias.py`** (`INFO`).
7. **Informe y figuras**: `informe.ESTRUCTURAS`, `_regenerar`, el párrafo de
   métodos, la sección de modelos de `informe_modelos.py`, y en `figuras.py` el
   color (`COLOR`, validado con pares adyacentes), `ORDEN` y `NOMBRES`.
8. **Interfaz** (`visor.py`): el grupo de deslizadores (prefijo propio en
   `Visor.sl`, ver `PREFIJO_SL`), `_params`, `_generar`, `_ajuste_listo` y la
   rama de `ajustar`; el selector sale de `MODOS_FAMILIA`.
9. **Traducciones**: `python idioma_revisar.py` y `--esqueleto`.

Trampa conocida: un nombre de familia nuevo que llega a un `else:` cae **sin
error** en el spinodoide (genera con `generar_mascara` y pide `wave_number`).
Búscalos con `grep -n '"dual-lattice"'` antes de dar la familia por integrada.

Si la librería de la familia tiene otra licencia o exige otro Python, hazla
opcional (extra en `pyproject.toml`, importación perezosa y una función
`disponible()`), y que el selector la deje desactivada en vez de fallar.
