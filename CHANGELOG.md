# Registro de cambios

Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/).
Las versiones siguen [SemVer](https://semver.org/lang/es/): un cambio en el
tercer número no altera ni los resultados ni la interfaz de programación.

## [1.0.1] - 2026-09-28

Versión de orden interno. No cambia ningún cálculo, ninguna salida numérica,
ninguna opción de la interfaz ni la API pública de `spinpy`. Lo único que
cambia en los resultados es el número de versión que viaja en el bloque de
procedencia y en los informes.

### Reorganización

- Los diálogos secundarios de la aplicación pasan de la raíz a `interfaz/`:
  `dialogo_febio.py`, `dialogo_metodos.py` y `dialogo_validacion.py`.
  `visor.py` sigue en la raíz porque es el punto de entrada de PyInstaller,
  de `instalador/construir.ps1` y de la orden `python visor.py`.
- El contraste con la implementación MATLAB queda reunido en
  `validacion_matlab/`: los guiones Python (`comparar.py`, `validar_*.py`,
  `probar_voi_h4.py`), los guiones MATLAB que antes estaban en `matlab/` y
  las figuras del contraste que antes estaban en `figs/`. Las salidas de
  referencia siguen en `resultados/`, que también usan la suite de pruebas y
  el anexo de validación.
- `idioma_revisar.py` pasa a `herramientas/`. Se ejecuta con
  `python herramientas/idioma_revisar.py`.
- `Test/`, `comparativa_febio/` y `comparativa_febio_tet/` conservan su nombre:
  están citadas en los informes PDF ya publicados, en textos visibles de la
  interfaz y en las claves del diccionario de traducción.

### Corregido

- Los guiones MATLAB calculaban mal la raíz del proyecto desde que se
  movieron a `matlab/`: buscaban `Validacion_Anexo/` y `H4/` dentro del
  repositorio y escribían en `matlab/resultados/`, donde los guiones Python
  no leen. Ahora escriben en `resultados/` de la raíz, como indica su propia
  cabecera.

### Limpieza de código

- Importaciones sin uso eliminadas en `visor.py`, `spinpy/resistencia.py`,
  `instalador/febio_minimo.py` y dos pruebas.
- Variables locales sin uso eliminadas en `visor.py` y
  `spinpy/resistencia.py`.
- Entrada duplicada (con la misma traducción) eliminada en
  `spinpy/idioma_textos.py`.
- `morfometria_malla`, ya importada en `spinpy/__init__.py`, figura ahora en
  `__all__`.
- `.gitignore` excluye `Test/resultados/`, la salida regenerable de las
  réplicas.

### Documentación

- Los README (español e inglés) incluyen la estructura del repositorio y las
  rutas nuevas.
- Versión 1.0.1 en `spinpy/__init__.py`, `pyproject.toml`, `CITATION.cff`,
  los README y el instalador de Windows.

## [1.0.0] - 2026-09-27

Primera versión estable.
