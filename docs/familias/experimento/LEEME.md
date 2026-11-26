# Scripts del experimento de level-set foam

Repiten, con el código de spinpy, la comparación de familias descrita en
[`../LEVEL_SET_FOAM.md`](../LEVEL_SET_FOAM.md). Se ejecutan **desde esta
carpeta** (usan rutas relativas) y necesitan morphoGen (`pip install
"spinpy[level-set]"`, Python ≥ 3.12) y los VOI de `datos/porcino_koria2020/`.

| Script | Qué hace |
|---|---|
| `objetivo.py` | Parches EN MEMORIA sobre spinpy: registra Conn.D, mide los términos de forma del VOI al vóxel del candidato y da el mismo peso de DA₂ a todas las familias. No modifica el repositorio. |
| `exp.py VOI FAMILIA` | Ronda 1: objetivo base. FAMILIA = `spinodoide`, `dual`, `open_foam`, `levelset`, `levelset_tpl`. |
| `exp3.py VOI FAMILIA` | Objetivo ampliado con la densidad fijada al BV/TV medido. |
| `exp4.py VOI` | `level_set_foam` con sus parámetros de forma (K_c, c₁, c₂) ajustados. |
| `stage2.py VOI FAMILIA` | Propiedades que no entran en el objetivo: Conn.D, SMI, espesor, curvatura, poros y rigidez elástica periódica. |
| `resumen3.py`, `resumen4.py` | Tablas de resultados. |

Variables de entorno: `EXP_OUT` (carpeta de salida; por defecto `./out` o
`./out3`/`./out4` según el script), `NUMBA_NUM_THREADS=1` y
`OMP_NUM_THREADS=1` si se lanzan varios procesos a la vez (uno por núcleo).

Ejemplo (un VOI, una familia):

```bash
cd docs/familias/experimento
EXP_OUT=./out3 python exp3.py C1 levelset
EXP_OUT=./out3 python stage2.py C1 levelset
```

Los resultados (`out*/res_*.json`) no se versionan: pesan poco pero se
regeneran. La búsqueda completa de una familia tarda de 4 a 10 minutos por VOI en
un núcleo.
