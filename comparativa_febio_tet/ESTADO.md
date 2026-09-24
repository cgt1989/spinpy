# Estado de la integración con FEBio (2026-09-24, segunda sesión)

Plan: `../../PLAN_Analizar_con_FEBio.md`. Resultados medidos: `INFORME.md`.

## Hecho

| pieza | estado |
|---|---|
| `spinpy/febio.py` (localizar, correr cancelable, leer, mallar, `mallar_aislado`, ensayo, `analizar`, `homogeneizar`, `tamano_previsto`, CLI con `--conv-malla`) | hecho y probado |
| `escribe.escribir_febio_ensayo` (hex8/tet10, `adaptativo`) | bloque 27 |
| estadísticos ponderados por volumen (`resistencia`) | bloque 27, bit a bit |
| `informe.comprobar` → `_items_febio` (motivos `febio_*`, `detalle`) | bloque 27 |
| `informe_modelos`: subsección «Resolución en FEBio» + `maas2012` (Crossref) | bloques 19/22/23 pasan |
| `tiempos.py`: modelo FEBio PROVISIONAL; memoria y tamaño de malla con test | bloque 27 |
| GUI: `dialogo_febio.py` (FEM automático según la maqueta, panel corto, parámetros, resultados, mapas, figura de 5 columnas), `visor.py` (`act_fem_auto`, `btn_febio`, `_febio_una`, `_febio_listo`, `fem_auto`, `_fem_auto_fin`, detener mata FEBio) | prueba de humo sin pantalla OK; `idioma_revisar.py` 0/0/0 |
| Empaquetado: `instalador/febio_minimo.py` (19 archivos, 103 MB, probado aislado), `spinpy.spec` (carpeta `febio/` + EULA + `LICENCIA_FEBIO.txt`), autocomprobación con ensayo FEBio hex8/TET10, `actualizar.ps1` (dialogo_febio, comparativa_febio_tet), LEEME del instalador y sección del README de GitHub | escrito; **`actualizar.ps1` NO lanzado** |
| Suite completa | 205 correctas, 2 fallos conocidos (bloques 04 y 05) |

Licencia: los binarios de FEBio Studio no se pueden redistribuir (EULA 4.0).
El usuario decidió incluirlos en Instalador y RAR para uso interno del
laboratorio, con la licencia de redistribución EN TRÁMITE; GitHub no lleva
binarios.

## Pendiente

1. Lanzar `"Para compartir/actualizar.ps1"` (avisar antes: compite por CPU)
   y verificar que la autocomprobación DEL EJECUTABLE pase la prueba de FEBio
   (proceso hijo de `mallar_aislado` dentro del ejecutable congelado).
2. Calibrar `tiempos.py` con el equipo tranquilo: `calibrar_tiempos.py`.
3. TET10 del VOI real a 40-48³ (~5-7 GB) con el equipo libre; ver §5 del
   INFORME para lo medido a 32³.
4. Probar a mano las dos ventanas con pantalla (la prueba de humo fue sin
   pantalla).
5. Cuando la otra sesión termine: `comparativa_febio/*.py` pueden importar de
   `spinpy.febio` (hoy se copió, no se movió).
