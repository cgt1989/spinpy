## Tablas de tiempos

Tiempo total del motor en segundos (malla propia del motor, montaje, resolucion y extraccion; sin importar la biblioteca). `dir.` = factorizacion directa, `ite.` = CG con multigrid algebraico; la app usa su propio CG + AMG.


### exactos

| Caso | GDL | App pro. | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bloque_hex | 3 159 | 0.20 | 0.32 | 0.33 | 0.72 | 0.74 | 0.89 | 0.90 | 0.89 | 0.91 |
| bloque_tet | 13 446 | n/d | 0.42 | 0.51 | 0.85 | 1.36 | 1.84 | 1.93 | 2.91 | 1.71 |


### cavidad

| Caso | GDL | App pro. | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cavidad_hex_n32 | 105 558 | 7.59 | 26.7 | 5.97 | 5.17 | 6.32 | 316.8 | 28.2 | 591.2 | 7.71 |
| cavidad_tet_esfera64 | 206 142 | n/d | 27.4 | 31.7 | 7.80 | 26.0 | n/d | 46.3 | n/d | 35.8 |


### espinodoide_hex

| Caso | GDL | App pro. | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| espinodoide_hex_n24 | 30 672 | 2.40 | 0.69 | 3.78 | 1.13 | 2.49 | 4.97 | 5.91 | 4.04 | 3.10 |
| espinodoide_hex_n32 | 62 595 | 5.33 | 1.01 | 6.61 | 1.72 | 5.20 | 11.3 | 12.5 | 10.0 | 6.68 |
| espinodoide_hex_n48 | 173 484 | 24.9 | 3.80 | 32.0 | 4.35 | 21.4 | n/d | 61.9 | n/d | 32.5 |
| espinodoide_hex_n64 | 365 655 | 47.5 | 9.21 | 50.9 | 9.91 | 43.3 | n/d | 130.0 | n/d | 55.1 |
| espinodoide_hex_n80 | 661 965 | 102.3 | n/d | 107.3 | n/d | 78.9 | n/d | 291.9 | n/d | 115.6 |


### espinodoide_tet

| Caso | GDL | App pro. | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| espinodoide_tet_n32 | 314 880 | n/d | 8.73 | 11.6 | 5.71 | n/d | n/d | 194.4 | n/d | 184.9 |
| espinodoide_tet_n48 | 790 599 | n/d | 32.0 | 43.1 | 18.5 | 225.0 | n/d | 672.7 | n/d | 806.1 |


### nl (Newton, resolvedor directo)

| Caso | GDL | NGSolve | FEniCSx | scikit-fem | SfePy |
|---|---:|---:|---:|---:|---:|
| bloque_hex_nl_svk | 525 | 0.31 | 2.81 | 2.42 | 1.77 |
| bloque_hex_nl_neohookeano | 525 | 0.33 | 2.51 | 4.72 | n/d |
| espinodoide_hex_n24_nl_svk | 30 672 | 6.26 | 3.23 | 110.4 | FALLO |
| espinodoide_hex_n24_nl_neohookeano | 30 672 | 7.57 | 3.78 | 175.0 | n/d |


## Configuracion mas rapida de cada motor


### espinodoide_hex

| Caso | GDL | App (spinpy) | NGSolve | FEniCSx | scikit-fem | SfePy |
|---|---:|---:|---:|---:|---:|---:|
| espinodoide_hex_n24 | 30 672 | 2.40 (pro.; 0.2 GB) | 0.69 (dir.; 0.2 GB) | 1.13 (dir.; 0.2 GB) | 4.97 (dir.; 0.3 GB) | 3.10 (ite.; 0.2 GB) |
| espinodoide_hex_n32 | 62 595 | 5.33 (pro.; 0.3 GB) | 1.01 (dir.; 0.3 GB) | 1.72 (dir.; 0.3 GB) | 11.3 (dir.; 0.7 GB) | 6.68 (ite.; 0.3 GB) |
| espinodoide_hex_n48 | 173 484 | 24.9 (pro.; 0.7 GB) | 3.80 (dir.; 0.7 GB) | 4.35 (dir.; 0.8 GB) | 61.9 (ite.; 1.2 GB) | 32.5 (ite.; 0.7 GB) |
| espinodoide_hex_n64 | 365 655 | 47.5 (pro.; 1.4 GB) | 9.21 (dir.; 1.5 GB) | 9.91 (dir.; 1.7 GB) | 130.0 (ite.; 3.0 GB) | 55.1 (ite.; 1.5 GB) |
| espinodoide_hex_n80 | 661 965 | 102.3 (pro.; 2.5 GB) | 107.3 (ite.; 4.0 GB) | 78.9 (ite.; 1.2 GB) | 291.9 (ite.; 5.6 GB) | 115.6 (ite.; 2.7 GB) |


### espinodoide_tet

| Caso | GDL | App (spinpy) | NGSolve | FEniCSx | scikit-fem | SfePy |
|---|---:|---:|---:|---:|---:|---:|
| espinodoide_tet_n32 | 314 880 | n/d | 8.73 (dir.; 1.0 GB) | 5.71 (dir.; 1.4 GB) | 194.4 (ite.; 2.9 GB) | 184.9 (ite.; 1.4 GB) |
| espinodoide_tet_n48 | 790 599 | n/d | 32.0 (dir.; 3.2 GB) | 18.5 (dir.; 3.6 GB) | 672.7 (ite.; 6.8 GB) | 806.1 (ite.; 3.2 GB) |


## Exponente de escala

Pendiente de la regresion log-log del tiempo total frente a los GDL (hexaedros; todos los tamanos del espinodoide con resultado). Un exponente de 1 es escala lineal.


| Motor | Resolvedor | Puntos | Exponente |
|---|---:|---:|---:|
| App (spinpy) | propio | 5 | 1.23 |
| NGSolve | directo | 4 | 1.09 |
| NGSolve | iterativo | 5 | 1.11 |
| FEniCSx | directo | 4 | 0.88 |
| FEniCSx | iterativo | 5 | 1.15 |
| scikit-fem | iterativo | 5 | 1.29 |
| SfePy | iterativo | 5 | 1.19 |


## Contraste con la campana publicada

Cociente entre el tiempo medido ahora y el publicado (`resultados/`) para la misma combinacion de caso, motor y resolvedor (hexaedros del espinodoide). Mide la reproducibilidad entre campanas, no el error del motor.

| Caso | App pro. | NGSolve dir. | NGSolve ite. | FEniCSx dir. | FEniCSx ite. | scikit-fem dir. | scikit-fem ite. | SfePy dir. | SfePy ite. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| espinodoide_hex_n24 | 0.62 | 0.62 | 0.75 | 1.16 | 0.67 | 0.72 | 0.71 | 0.76 | 0.62 |
| espinodoide_hex_n32 | 0.66 | 0.48 | 0.61 | 0.88 | 0.71 | 0.60 | 0.65 | 0.62 | 0.59 |
| espinodoide_hex_n48 | 0.87 | 0.49 | 0.75 | 0.78 | 0.72 | n/d | 0.76 | n/d | 0.74 |
| espinodoide_hex_n64 | 0.79 | 0.46 | 0.72 | 0.73 | 0.69 | n/d | 0.78 | n/d | 0.79 |
| espinodoide_hex_n80 | 0.80 | n/d | 0.76 | n/d | 0.67 | n/d | 0.75 | n/d | 0.77 |
| espinodoide_tet_n32 | n/d | 0.44 | 0.63 | 0.72 | n/d | n/d | 0.75 | n/d | 0.78 |
| espinodoide_tet_n48 | n/d | 0.46 | 0.73 | 0.81 | 0.69 | n/d | 0.74 | n/d | 0.80 |


## Corridas fallidas o omitidas

| Caso | Motor | Resolvedor | Motivo |
|---|---:|---:|---:|
| espinodoide_tet_n32 | fenicsx | iterativo | spinpy.motores._comun.ErrorMotor: PETSc KSP no convergio (-8) |
| espinodoide_hex_n24_nl_svk | sfepy | directo | spinpy.motores._comun.ErrorMotor: Newton de SfePy no convergio en la carga 0.02 (30 iteraciones) |
