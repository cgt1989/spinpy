# Tablas de la comparativa de motores (generadas por figuras.py)

#### Casos con solución exacta (bloque macizo)

| Caso | GDL | Motor | Resolvedor | Tiempo (s) | JIT (s) | Memoria (MB) | Residuo | Iter. | Δu/ref | ΔvM/ref | E_app (MPa) | p99 sup. (MPa) | σ fallo (MPa) | dF |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bloque_hex | 3159 | App (spinpy) | LU directo | 0.27 | n/d | 107 | n/d | n/d | 3.4e-13 | 1.8e-13 | 20000 | n/d | 140 | 1.1e-14 |
| bloque_hex | 3159 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 1.08 | n/d | 124 | 7.0e-11 | 12 | 2.9e-11 | 3.7e-11 | 20000 | n/d | 140 | 8.6e-13 |
| bloque_hex | 3159 | scikit-fem | SuperLU (scipy) | 0.97 | n/d | 131 | 2.8e-14 | n/d | 2.2e-14 | 1.3e-14 | 20000 | n/d | 140 | 1.9e-15 |
| bloque_hex | 3159 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 0.85 | n/d | 169 | n/d | n/d | 2.9e-11 | 3.7e-11 | 20000 | n/d | 140 | 8.6e-13 |
| bloque_hex | 3159 | SfePy | SuperLU (ls.scipy_direct) | 1.10 | n/d | 170 | n/d | n/d | 1.3e-13 | 4.8e-14 | 20000 | n/d | 140 | 3.1e-15 |
| bloque_hex | 3159 | NGSolve | CG + pyamg (matriz de NGSolve) | 0.50 | n/d | 202 | 3.8e-11 | 12 | 1.1e-11 | 2.6e-11 | 20000 | n/d | 140 | 1.2e-13 |
| bloque_hex | 3159 | NGSolve | PARDISO (MKL) | 0.41 | n/d | 226 | 7.1e-15 | n/d | 0 | 0 | 20000 | n/d | 140 | 6.9e-16 |
| bloque_hex | 3159 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 0.58 | 1.55 | 173 | 9.1e-11 | 16 | 2.7e-11 | 5.4e-11 | 20000 | n/d | 140 | 1.3e-12 |
| bloque_hex | 3159 | FEniCSx | Cholesky MUMPS (PETSc) | 0.14 | 0.00 | 184 | 1.8e-14 | 1 | 5.6e-13 | 2.8e-13 | 20000 | n/d | 140 | 9.5e-15 |
| bloque_tet | 13446 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 2.46 | n/d | 208 | 9.5e-11 | 45 | 7.7e-11 | 1.0e-07 | 20000 | n/d | 140 | 7.6e-14 |
| bloque_tet | 13446 | scikit-fem | SuperLU (scipy) | 2.66 | n/d | 279 | 1.3e-13 | n/d | 2.5e-13 | 2.3e-11 | 20000 | n/d | 140 | 4.3e-15 |
| bloque_tet | 13446 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 2.13 | n/d | 208 | n/d | n/d | 5.8e-11 | 8.8e-08 | 20000 | n/d | 140 | 3.7e-13 |
| bloque_tet | 13446 | SfePy | SuperLU (ls.scipy_direct) | 4.60 | n/d | 276 | n/d | n/d | 1.1e-13 | 2.9e-11 | 20000 | n/d | 140 | 9.0e-15 |
| bloque_tet | 13446 | NGSolve | CG + BDDC (NGSolve) | 0.82 | n/d | 226 | 3.6e-14 | 2 | 2.5e-13 | 1.6e-11 | 20000 | n/d | 140 | 2.8e-15 |
| bloque_tet | 13446 | NGSolve | PARDISO (MKL) | 0.82 | n/d | 305 | 1.8e-14 | n/d | 0 | 0 | 20000 | n/d | 140 | 3.5e-15 |
| bloque_tet | 13446 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 1.81 | 1.83 | 188 | 9.9e-11 | 53 | 1.6e-10 | 1.0e-07 | 20000 | n/d | 140 | 3.3e-13 |
| bloque_tet | 13446 | FEniCSx | Cholesky MUMPS (PETSc) | 0.32 | 0.00 | 216 | 7.9e-14 | 1 | 1.2e-13 | 1.3e-10 | 20000 | n/d | 140 | 1.9e-15 |

#### Cavidad esférica

| Caso | GDL | Motor | Resolvedor | Tiempo (s) | JIT (s) | Memoria (MB) | Residuo | Iter. | Δu/ref | ΔvM/ref | E_app (MPa) | p99 sup. (MPa) | σ fallo (MPa) | dF |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cavidad_hex_n32 | 105558 | App (spinpy) | AMG (con modos rigidos) + CG | 10.24 | n/d | 552 | 2.8e-09 | n/d | 7.7e-10 | 9.0e-10 | 18690 | 1.8775 | 106.04 | 3.8e-12 |
| cavidad_hex_n32 | 105558 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 43.74 | n/d | 1303 | 1.7e-11 | 18 | 3.7e-12 | 6.5e-12 | 18690 | 1.8775 | 106.04 | 1.0e-13 |
| cavidad_hex_n32 | 105558 | scikit-fem | SuperLU (scipy) | 516.19 | n/d | 7311 | 2.5e-13 | n/d | 1.5e-13 | 8.2e-14 | 18690 | 1.8775 | 106.04 | 1.7e-14 |
| cavidad_hex_n32 | 105558 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 11.39 | n/d | 658 | n/d | n/d | 3.8e-12 | 6.5e-12 | 18690 | 1.8775 | 106.04 | 6.8e-14 |
| cavidad_hex_n32 | 105558 | SfePy | SuperLU (ls.scipy_direct) | 922.71 | n/d | 6998 | n/d | n/d | 5.1e-13 | 5.1e-13 | 18690 | 1.8775 | 106.04 | 2.5e-14 |
| cavidad_hex_n32 | 105558 | NGSolve | CG + pyamg (matriz de NGSolve) | 11.80 | n/d | 915 | 9.2e-11 | 17 | 1.9e-11 | 5.8e-11 | 18690 | 1.8775 | 106.04 | 1.7e-13 |
| cavidad_hex_n32 | 105558 | NGSolve | PARDISO (MKL) | 25.01 | n/d | 3570 | 3.4e-14 | n/d | 0 | 0 | 18690 | 1.8775 | 106.04 | 1.1e-15 |
| cavidad_hex_n32 | 105558 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 8.14 | 0.01 | 346 | 5.7e-11 | 17 | 1.6e-11 | 3.5e-11 | 18690 | 1.8775 | 106.04 | 3.1e-13 |
| cavidad_hex_n32 | 105558 | FEniCSx | Cholesky MUMPS (PETSc) | 7.30 | 0.01 | 1525 | 1.2e-13 | 1 | 1.4e-11 | 1.7e-11 | 18690 | 1.8775 | 106.04 | 1.7e-13 |
| cavidad_tet_esfera64 | 206142 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 73.80 | n/d | 2438 | 7.6e-11 | 66 | 2.0e-10 | 2.0e-08 | 18647 | 2.0965 | 105.69 | 2.3e-13 |
| cavidad_tet_esfera64 | 206142 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 49.14 | n/d | 1071 | n/d | n/d | 2.4e-10 | 1.6e-08 | 18647 | 2.0965 | 105.69 | 4.6e-13 |
| cavidad_tet_esfera64 | 206142 | NGSolve | CG + BDDC (NGSolve) | 45.53 | n/d | 3898 | 2.0e-13 | 2 | 5.3e-12 | 5.5e-12 | 18647 | 2.0965 | 105.69 | 2.0e-14 |
| cavidad_tet_esfera64 | 206142 | NGSolve | PARDISO (MKL) | 42.89 | n/d | 5385 | 1.7e-14 | n/d | 0 | 0 | 18647 | 2.0965 | 105.69 | 5.6e-15 |
| cavidad_tet_esfera64 | 206142 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 33.33 | 0.01 | 548 | 8.2e-11 | 59 | 3.7e-10 | 3.5e-08 | 18647 | 2.0965 | 105.69 | 1.8e-12 |
| cavidad_tet_esfera64 | 206142 | FEniCSx | Cholesky MUMPS (PETSc) | 10.29 | 0.01 | 2452 | 9.6e-14 | 1 | 2.5e-13 | 1.9e-11 | 18647 | 2.0965 | 105.69 | 2.9e-15 |

#### Espinodoide, malla de ladrillos (hex8)

| Caso | GDL | Motor | Resolvedor | Tiempo (s) | JIT (s) | Memoria (MB) | Residuo | Iter. | Δu/ref | ΔvM/ref | E_app (MPa) | p99 sup. (MPa) | σ fallo (MPa) | dF |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| espinodoide_hex_n24 | 30672 | App (spinpy) | AMG (con modos rigidos) + CG | 3.88 | n/d | 183 | 9.6e-09 | n/d | 3.9e-10 | 7.2e-10 | 326.534 | 36.551 | 4.7185 | 4.4e-11 |
| espinodoide_hex_n24 | 30672 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 8.26 | n/d | 241 | 8.6e-11 | 109 | 8.7e-12 | 6.2e-12 | 326.534 | 36.551 | 4.7185 | 1.4e-13 |
| espinodoide_hex_n24 | 30672 | scikit-fem | SuperLU (scipy) | 6.87 | n/d | 349 | 6.2e-13 | n/d | 2.9e-13 | 2.3e-13 | 326.534 | 36.551 | 4.7185 | 1.8e-14 |
| espinodoide_hex_n24 | 30672 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 4.78 | n/d | 249 | n/d | n/d | 8.7e-12 | 6.1e-12 | 326.534 | 36.551 | 4.7185 | 4.7e-14 |
| espinodoide_hex_n24 | 30672 | SfePy | SuperLU (ls.scipy_direct) | 5.11 | n/d | 348 | n/d | n/d | 4.8e-13 | 4.2e-13 | 326.534 | 36.551 | 4.7185 | 8.3e-14 |
| espinodoide_hex_n24 | 30672 | NGSolve | CG + pyamg (matriz de NGSolve) | 5.05 | n/d | 325 | 9.1e-11 | 108 | 1.0e-11 | 7.2e-12 | 326.534 | 36.551 | 4.7185 | 3.4e-13 |
| espinodoide_hex_n24 | 30672 | NGSolve | PARDISO (MKL) | 1.09 | n/d | 364 | 3.2e-13 | n/d | 0 | 0 | 326.534 | 36.551 | 4.7185 | 2.6e-15 |
| espinodoide_hex_n24 | 30672 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 3.39 | 0.01 | 201 | 7.7e-11 | 81 | 8.6e-12 | 4.9e-12 | 326.534 | 36.551 | 4.7185 | 1.2e-12 |
| espinodoide_hex_n24 | 30672 | FEniCSx | Cholesky MUMPS (PETSc) | 0.64 | 0.01 | 237 | 5.0e-13 | 1 | 3.5e-12 | 5.3e-12 | 326.534 | 36.551 | 4.7185 | 1.4e-12 |
| espinodoide_hex_n32 | 62595 | App (spinpy) | AMG (con modos rigidos) + CG | 8.01 | n/d | 312 | 8.9e-09 | n/d | 8.4e-10 | 6.3e-10 | 261.146 | 49.69 | 3.7951 | 1.6e-11 |
| espinodoide_hex_n32 | 62595 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 19.34 | n/d | 446 | 9.2e-11 | 105 | 3.7e-12 | 7.2e-12 | 261.146 | 49.69 | 3.7951 | 8.5e-14 |
| espinodoide_hex_n32 | 62595 | scikit-fem | SuperLU (scipy) | 18.93 | n/d | 716 | 1.4e-12 | n/d | 2.1e-13 | 1.7e-13 | 261.146 | 49.69 | 3.7951 | 4.1e-14 |
| espinodoide_hex_n32 | 62595 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 11.06 | n/d | 364 | n/d | n/d | 3.5e-12 | 7.2e-12 | 261.146 | 49.69 | 3.7951 | 1.5e-13 |
| espinodoide_hex_n32 | 62595 | SfePy | SuperLU (ls.scipy_direct) | 16.03 | n/d | 628 | n/d | n/d | 4.0e-13 | 3.6e-13 | 261.146 | 49.69 | 3.7951 | 5.1e-14 |
| espinodoide_hex_n32 | 62595 | NGSolve | CG + pyamg (matriz de NGSolve) | 10.85 | n/d | 488 | 8.2e-11 | 105 | 4.0e-12 | 5.1e-12 | 261.146 | 49.69 | 3.7951 | 1.8e-13 |
| espinodoide_hex_n32 | 62595 | NGSolve | PARDISO (MKL) | 2.10 | n/d | 581 | 6.2e-13 | n/d | 0 | 0 | 261.146 | 49.69 | 3.7951 | 2.5e-14 |
| espinodoide_hex_n32 | 62595 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 7.03 | 0.01 | 249 | 8.2e-11 | 80 | 1.2e-11 | 1.4e-11 | 261.146 | 49.69 | 3.7951 | 1.1e-12 |
| espinodoide_hex_n32 | 62595 | FEniCSx | Cholesky MUMPS (PETSc) | 1.60 | 0.01 | 337 | 1.0e-12 | 1 | 1.1e-11 | 1.4e-11 | 261.146 | 49.69 | 3.7951 | 1.0e-12 |
| espinodoide_hex_n48 | 173484 | App (spinpy) | AMG (con modos rigidos) + CG | 28.62 | n/d | 684 | 7.6e-09 | n/d | 5.8e-10 | 7.7e-10 | 298.888 | 53.529 | 3.8215 | 4.7e-11 |
| espinodoide_hex_n48 | 173484 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 81.36 | n/d | 1266 | 7.1e-11 | 154 | 5.2e-12 | 8.3e-12 | 298.888 | 53.529 | 3.8215 | 6.8e-14 |
| espinodoide_hex_n48 | 173484 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 43.79 | n/d | 769 | n/d | n/d | 5.3e-12 | 8.2e-12 | 298.888 | 53.529 | 3.8215 | 5.9e-13 |
| espinodoide_hex_n48 | 173484 | NGSolve | CG + pyamg (matriz de NGSolve) | 42.91 | n/d | 1103 | 8.2e-11 | 153 | 5.3e-12 | 8.7e-12 | 298.888 | 53.529 | 3.8215 | 3.4e-13 |
| espinodoide_hex_n48 | 173484 | NGSolve | PARDISO (MKL) | 7.59 | n/d | 1546 | 1.1e-12 | n/d | 0 | 0 | 298.888 | 53.529 | 3.8215 | 1.4e-14 |
| espinodoide_hex_n48 | 173484 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 29.27 | 0.01 | 415 | 7.0e-11 | 114 | 2.2e-11 | 4.0e-11 | 298.888 | 53.529 | 3.8215 | 3.7e-12 |
| espinodoide_hex_n48 | 173484 | FEniCSx | Cholesky MUMPS (PETSc) | 5.21 | 0.01 | 787 | 2.3e-12 | 1 | 1.9e-11 | 3.4e-11 | 298.888 | 53.529 | 3.8215 | 3.7e-12 |
| espinodoide_hex_n64 | 365655 | App (spinpy) | AMG (con modos rigidos) + CG | 60.10 | n/d | 1413 | 8.0e-09 | n/d | 9.2e-10 | 1.3e-09 | 280.474 | 60.867 | 3.5535 | 1.5e-11 |
| espinodoide_hex_n64 | 365655 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 165.88 | n/d | 2982 | 9.5e-11 | 84 | 7.0e-12 | 1.3e-11 | 280.474 | 60.867 | 3.5535 | 2.2e-13 |
| espinodoide_hex_n64 | 365655 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 69.64 | n/d | 1565 | n/d | n/d | 6.3e-12 | 1.2e-11 | 280.474 | 60.867 | 3.5535 | 2.3e-13 |
| espinodoide_hex_n64 | 365655 | NGSolve | CG + pyamg (matriz de NGSolve) | 70.41 | n/d | 2248 | 8.6e-11 | 84 | 5.7e-12 | 1.2e-11 | 280.474 | 60.867 | 3.5535 | 5.1e-14 |
| espinodoide_hex_n64 | 365655 | NGSolve | PARDISO (MKL) | 19.62 | n/d | 3719 | 1.8e-12 | n/d | 0 | 0 | 280.474 | 60.867 | 3.5535 | 5.1e-14 |
| espinodoide_hex_n64 | 365655 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 62.53 | 0.03 | 726 | 8.2e-11 | 104 | 2.0e-11 | 3.6e-11 | 280.474 | 60.867 | 3.5535 | 7.7e-13 |
| espinodoide_hex_n64 | 365655 | FEniCSx | Cholesky MUMPS (PETSc) | 13.15 | 0.03 | 1722 | 4.1e-12 | 1 | 1.9e-11 | 3.7e-11 | 280.474 | 60.867 | 3.5535 | 5.8e-13 |
| espinodoide_hex_n80 | 661965 | App (spinpy) | AMG (con modos rigidos) + CG | 127.98 | n/d | 2561 | 7.6e-09 | n/d | 5.0e-10 | 6.7e-10 | 264.85 | 66.723 | 3.4328 | 3.0e-12 |
| espinodoide_hex_n80 | 661965 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 388.47 | n/d | 5661 | 8.8e-11 | 91 | 2.3e-11 | 3.3e-11 | 264.85 | 66.723 | 3.4328 | 1.8e-13 |
| espinodoide_hex_n80 | 661965 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 149.31 | n/d | 2793 | n/d | n/d | 2.6e-11 | 3.0e-11 | 264.85 | 66.723 | 3.4328 | 2.3e-13 |
| espinodoide_hex_n80 | 661965 | NGSolve | CG + pyamg (matriz de NGSolve) | 141.23 | n/d | 4085 | 8.1e-11 | 91 | 2.7e-11 | 3.3e-11 | 264.85 | 66.723 | 3.4328 | 5.3e-13 |
| espinodoide_hex_n80 | 661965 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 117.99 | 0.05 | 1214 | 7.8e-11 | 100 | 0 | 0 | 264.85 | 66.723 | 3.4328 | 8.6e-12 |
| espinodoide_hex_n80 | 661965 | NGSolve | PARDISO (MKL) | 51.71 | n/d | 7920 | 2.6e-12 | n/d | 0 | 0 | 264.85 | 66.723 | 3.4328 | 2.4e-14 |
| espinodoide_hex_n80 | 661965 | FEniCSx | Cholesky MUMPS (PETSc) | 27.68 | 0.05 | 3532 | 6.4e-12 | 1 | 2.7e-11 | 3.6e-11 | 264.85 | 66.723 | 3.4328 | 8.6e-12 |
| espinodoide_hex_n48 | 173484 | NGSolve | sparsecholesky (NGSolve) | 5.26 | n/d | 670 | 3.9e-12 | n/d | 0 | 0 | 298.888 | 53.529 | 3.8215 | 1.1e-13 |
| espinodoide_hex_n64 | 365655 | NGSolve | sparsecholesky (NGSolve) | 14.43 | n/d | 1535 | 9.1e-12 | n/d | 0 | 0 | 280.474 | 60.867 | 3.5535 | 7.2e-13 |
| espinodoide_hex_n80 | 661965 | NGSolve | sparsecholesky (NGSolve) | 32.67 | n/d | 3352 | 1.6e-11 | n/d | 0 | 0 | 264.85 | 66.723 | 3.4328 | 1.2e-12 |

#### Espinodoide, malla suave (TET10)

| Caso | GDL | Motor | Resolvedor | Tiempo (s) | JIT (s) | Memoria (MB) | Residuo | Iter. | Δu/ref | ΔvM/ref | E_app (MPa) | p99 sup. (MPa) | σ fallo (MPa) | dF |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| espinodoide_tet_n32 | 314880 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 258.37 | n/d | 3017 | 9.3e-11 | 398 | 8.9e-03 | 7.2e-02 | 157.039 | 85.891 | 2.7254 | 4.5e-13 |
| espinodoide_tet_n32 | 314880 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 236.16 | n/d | 1403 | n/d | n/d | 1.8e-11 | 3.0e-10 | 156.683 | 85.767 | 2.7196 | 3.5e-13 |
| espinodoide_tet_n32 | 314880 | NGSolve | CG + BDDC (NGSolve) | 18.31 | n/d | 1299 | 6.4e-12 | 2 | 5.6e-12 | 6.8e-12 | 156.683 | 85.767 | 2.7196 | 3.3e-13 |
| espinodoide_tet_n32 | 314880 | NGSolve | PARDISO (MKL) | 19.46 | n/d | 2935 | 4.4e-12 | n/d | 0 | 0 | 156.683 | 85.767 | 2.7196 | 4.4e-14 |
| espinodoide_tet_n32 | 314880 | FEniCSx | iterativo | **falló**: spinpy.motores._comun.ErrorMotor: PETSc KSP no convergio (-8 | | | | | | | | | | |
| espinodoide_tet_n32 | 314880 | FEniCSx | Cholesky MUMPS (PETSc) | 7.48 | 0.01 | 1413 | 1.7e-11 | 1 | 4.1e-12 | 1.4e-11 | 156.683 | 85.767 | 2.7196 | 4.7e-13 |
| espinodoide_tet_n48 | 790599 | scikit-fem | CG + AMG agregacion suavizada (pyamg, modos rigidos) | 905.07 | n/d | 7045 | 1.0e-10 | 541 | 1.6e-03 | 5.2e-02 | 219.046 | 79.024 | 3.0737 | 4.2e-13 |
| espinodoide_tet_n48 | 790599 | SfePy | CG (ls.scipy_iterative) + pyamg con modos rigidos | 1009.66 | n/d | 3259 | n/d | n/d | 8.2e-12 | 1.4e-11 | 219.024 | 79.029 | 3.0729 | 5.5e-13 |
| espinodoide_tet_n48 | 790599 | NGSolve | CG + BDDC (NGSolve) | 58.71 | n/d | 3904 | 9.2e-12 | 2 | 1.5e-12 | 1.9e-12 | 219.024 | 79.029 | 3.0729 | 1.0e-12 |
| espinodoide_tet_n48 | 790599 | NGSolve | PARDISO (MKL) | 69.21 | n/d | 8310 | 5.4e-12 | n/d | 0 | 0 | 219.024 | 79.029 | 3.0729 | 2.5e-13 |
| espinodoide_tet_n48 | 790599 | FEniCSx | CG + GAMG (PETSc, modos rigidos) | 325.04 | 0.02 | 1484 | 1.1e-10 | 244 | 5.7e-12 | 4.6e-11 | 219.024 | 79.029 | 3.0729 | 1.4e-13 |
| espinodoide_tet_n48 | 790599 | FEniCSx | Cholesky MUMPS (PETSc) | 22.34 | 0.02 | 3770 | 2.1e-11 | 1 | 6.9e-12 | 6.4e-12 | 219.024 | 79.029 | 3.0729 | 3.9e-13 |
| espinodoide_tet_n32 | 314880 | NGSolve | sparsecholesky (NGSolve) | 13.51 | n/d | 1036 | 6.4e-12 | n/d | 0 | 0 | 156.683 | 85.767 | 2.7196 | 5.2e-14 |
| espinodoide_tet_n48 | 790599 | NGSolve | sparsecholesky (NGSolve) | 44.81 | n/d | 3243 | 9.2e-12 | n/d | 0 | 0 | 219.024 | 79.029 | 3.0729 | 3.8e-13 |

#### No lineal (control por desplazamiento del plato)

| Caso | GDL | Motor | Tiempo (s) | Memoria (MB) | Iteraciones de Newton por paso | F por paso (N) | Error frente a la solución cerrada | Δu/ref |
|---|---:|---|---:|---:|---|---|---:|---:|
| bloque_hex_nl_svk | 525 | scikit-fem | 3.67 | 93 | [3, 3, 3, 3] | 37.05, 68.4, 94.35, 115.2 | 4.1e-12 | 2.4e-15 |
| bloque_hex_nl_svk | 525 | SfePy | 2.30 | 150 | [22, 24, 24, 27] | 37.05, 68.4, 94.35, 115.2 | 4.4e-16 | 5.0e-10 |
| bloque_hex_nl_svk | 525 | NGSolve | 0.67 | 203 | [3, 3, 3, 3] | 37.05, 68.4, 94.35, 115.2 | 4.1e-12 | 0 |
| bloque_hex_nl_svk | 525 | FEniCSx | 0.16 | 176 | [3, 3, 3, 3] | 37.05, 68.4, 94.35, 115.2 | 4.1e-12 | 2.7e-14 |
| bloque_hex_nl_neohookeano | 525 | scikit-fem | 7.91 | 93 | [3, 3, 3, 3] | 41.6387, 86.9777, 136.776, 192.004 | 8.1e-13 | 4.0e-16 |
| bloque_hex_nl_neohookeano | 525 | NGSolve | 0.69 | 203 | [3, 3, 3, 3] | 41.6387, 86.9777, 136.776, 192.004 | 8.1e-13 | 0 |
| bloque_hex_nl_neohookeano | 525 | FEniCSx | 0.17 | 176 | [3, 3, 3, 3] | 41.6387, 86.9777, 136.776, 192.004 | 8.1e-13 | 1.2e-15 |
| espinodoide_hex_n24_nl_svk | 30672 | scikit-fem | 188.85 | 399 | [3, 3, 4, 4] | 33.2523, 66.0776, 130.454, 254.156 | n/d | 3.8e-15 |
| espinodoide_hex_n24_nl_svk | 30672 | SfePy | **falló**: spinpy.motores._comun.ErrorMotor: Newton de SfePy no convergio en la c | | | | | |
| espinodoide_hex_n24_nl_svk | 30672 | NGSolve | 18.44 | 580 | [3, 3, 4, 4] | 33.2523, 66.0776, 130.454, 254.156 | n/d | 0 |
| espinodoide_hex_n24_nl_svk | 30672 | FEniCSx | 4.23 | 298 | [3, 3, 3, 4] | 33.2523, 66.0776, 130.454, 254.156 | n/d | 5.4e-15 |
| espinodoide_hex_n24_nl_neohookeano | 30672 | scikit-fem | 284.19 | 400 | [3, 3, 4, 4] | 33.352, 66.4736, 132.016, 260.223 | n/d | 2.5e-15 |
| espinodoide_hex_n24_nl_neohookeano | 30672 | NGSolve | 16.65 | 593 | [3, 3, 4, 4] | 33.352, 66.4736, 132.016, 260.223 | n/d | 0 |
| espinodoide_hex_n24_nl_neohookeano | 30672 | FEniCSx | 4.09 | 297 | [3, 3, 3, 4] | 33.352, 66.4736, 132.016, 260.223 | n/d | 5.6e-15 |
| bloque_hex_nl_svk | 525 | NGSolve | 0.43 | 169 | [3, 3, 3, 3] | 37.05, 68.4, 94.35, 115.2 | 4.1e-12 | 0 |
| bloque_hex_nl_neohookeano | 525 | NGSolve | 0.50 | 170 | [3, 3, 3, 3] | 41.6387, 86.9777, 136.776, 192.004 | 8.1e-13 | 0 |
| espinodoide_hex_n24_nl_svk | 30672 | NGSolve | 9.12 | 246 | [3, 3, 4, 4] | 33.2523, 66.0776, 130.454, 254.156 | n/d | 0 |
| espinodoide_hex_n24_nl_neohookeano | 30672 | NGSolve | 10.14 | 259 | [3, 3, 4, 4] | 33.352, 66.4736, 132.016, 260.223 | n/d | 0 |

#### Convergencia con solución exacta (refinamiento h)

| Elemento | Motor | n | GDL | error L2 rel. | error H1 rel. | orden local L2 | orden local H1 | Δu/ref | tiempo (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| hex8 (Q1) | App (spinpy) | 4 | 375 | 2.5455e-02 | 1.6470e-01 | n/d | n/d | 4.7e-16 | 0.29 |
| hex8 (Q1) | App (spinpy) | 8 | 2187 | 6.4374e-03 | 8.2649e-02 | 1.983 | 0.995 | 2.5e-15 | 0.03 |
| hex8 (Q1) | App (spinpy) | 16 | 14739 | 1.6134e-03 | 4.1362e-02 | 1.996 | 0.999 | 9.0e-10 | 0.66 |
| hex8 (Q1) | App (spinpy) | 32 | 107811 | 4.0356e-04 | 2.0686e-02 | 1.999 | 1.000 | 2.1e-08 | 4.63 |
| hex8 (Q1) | NGSolve | 4 | 375 | 2.5455e-02 | 1.6470e-01 | n/d | n/d | 0 | 0.25 |
| hex8 (Q1) | NGSolve | 8 | 2187 | 6.4374e-03 | 8.2649e-02 | 1.983 | 0.995 | 0 | 0.10 |
| hex8 (Q1) | NGSolve | 16 | 14739 | 1.6134e-03 | 4.1362e-02 | 1.996 | 0.999 | 0 | 1.14 |
| hex8 (Q1) | NGSolve | 32 | 107811 | 4.0356e-04 | 2.0686e-02 | 1.999 | 1.000 | 0 | 34.74 |
| hex8 (Q1) | FEniCSx | 4 | 375 | 2.5455e-02 | 1.6470e-01 | n/d | n/d | 1.1e-15 | 0.51 |
| hex8 (Q1) | FEniCSx | 8 | 2187 | 6.4374e-03 | 8.2649e-02 | 1.983 | 0.995 | 2.3e-15 | 0.55 |
| hex8 (Q1) | FEniCSx | 16 | 14739 | 1.6134e-03 | 4.1362e-02 | 1.996 | 0.999 | 5.6e-15 | 1.02 |
| hex8 (Q1) | FEniCSx | 32 | 107811 | 4.0356e-04 | 2.0686e-02 | 1.999 | 1.000 | 1.5e-13 | 9.92 |
| hex8 (Q1) | scikit-fem | 4 | 375 | 2.5455e-02 | 1.6470e-01 | n/d | n/d | 4.7e-16 | 0.06 |
| hex8 (Q1) | scikit-fem | 8 | 2187 | 6.4374e-03 | 8.2649e-02 | 1.983 | 0.995 | 1.9e-15 | 0.17 |
| hex8 (Q1) | scikit-fem | 16 | 14739 | 1.6134e-03 | 4.1362e-02 | 1.996 | 0.999 | 6.2e-15 | 2.73 |
| hex8 (Q1) | scikit-fem | 32 | 107811 | 4.0356e-04 | 2.0686e-02 | 1.999 | 1.000 | 1.2e-12 | 25.11 |
| hex8 (Q1) | SfePy | 4 | 375 | 2.5455e-02 | 1.6470e-01 | n/d | n/d | 4.0e-16 | 0.51 |
| hex8 (Q1) | SfePy | 8 | 2187 | 6.4374e-03 | 8.2649e-02 | 1.983 | 0.995 | 1.8e-15 | 0.09 |
| hex8 (Q1) | SfePy | 16 | 14739 | 1.6134e-03 | 4.1362e-02 | 1.996 | 0.999 | 6.3e-15 | 2.62 |
| hex8 (Q1) | SfePy | 32 | 107811 | 4.0356e-04 | 2.0686e-02 | 1.999 | 1.000 | 1.2e-12 | 7.70 |
| TET4 (P1) | NGSolve | 4 | 375 | 5.2970e-02 | 3.0234e-01 | n/d | n/d | 0 | 0.06 |
| TET4 (P1) | NGSolve | 8 | 2187 | 1.3673e-02 | 1.5277e-01 | 1.954 | 0.985 | 0 | 0.46 |
| TET4 (P1) | NGSolve | 16 | 14739 | 3.4528e-03 | 7.6585e-02 | 1.986 | 0.996 | 0 | 3.98 |
| TET4 (P1) | NGSolve | 32 | 107811 | 8.6554e-04 | 3.8317e-02 | 1.996 | 0.999 | 0 | 42.48 |
| TET4 (P1) | FEniCSx | 4 | 375 | 5.2970e-02 | 3.0234e-01 | n/d | n/d | 4.0e-16 | 0.51 |
| TET4 (P1) | FEniCSx | 8 | 2187 | 1.3673e-02 | 1.5277e-01 | 1.954 | 0.985 | 3.0e-15 | 0.54 |
| TET4 (P1) | FEniCSx | 16 | 14739 | 3.4528e-03 | 7.6585e-02 | 1.986 | 0.996 | 7.1e-15 | 0.91 |
| TET4 (P1) | FEniCSx | 32 | 107811 | 8.6554e-04 | 3.8317e-02 | 1.996 | 0.999 | 1.8e-13 | 5.86 |
| TET4 (P1) | scikit-fem | 4 | 375 | 5.2970e-02 | 3.0234e-01 | n/d | n/d | 5.5e-16 | 0.02 |
| TET4 (P1) | scikit-fem | 8 | 2187 | 1.3673e-02 | 1.5277e-01 | 1.954 | 0.985 | 3.6e-15 | 0.18 |
| TET4 (P1) | scikit-fem | 16 | 14739 | 3.4528e-03 | 7.6585e-02 | 1.986 | 0.996 | 7.0e-15 | 2.25 |
| TET4 (P1) | scikit-fem | 32 | 107811 | 8.6554e-04 | 3.8317e-02 | 1.996 | 0.999 | 7.0e-13 | 25.25 |
| TET4 (P1) | SfePy | 4 | 375 | 5.2970e-02 | 3.0234e-01 | n/d | n/d | 4.4e-16 | 0.02 |
| TET4 (P1) | SfePy | 8 | 2187 | 1.3673e-02 | 1.5277e-01 | 1.954 | 0.985 | 3.4e-15 | 0.08 |
| TET4 (P1) | SfePy | 16 | 14739 | 3.4528e-03 | 7.6585e-02 | 1.986 | 0.996 | 7.1e-15 | 2.25 |
| TET4 (P1) | SfePy | 32 | 107811 | 8.6554e-04 | 3.8317e-02 | 1.996 | 0.999 | 1.1e-12 | 7.53 |
| TET10 (P2) | NGSolve | 4 | 2187 | 2.9529e-03 | 2.8582e-02 | n/d | n/d | 0 | 0.23 |
| TET10 (P2) | NGSolve | 8 | 14739 | 3.6950e-04 | 7.2446e-03 | 2.998 | 1.980 | 0 | 1.29 |
| TET10 (P2) | NGSolve | 16 | 107811 | 4.6041e-05 | 1.8169e-03 | 3.005 | 1.995 | 0 | 23.26 |
| TET10 (P2) | NGSolve | 32 | 823875 | 5.7477e-06 | 4.5456e-04 | 3.002 | 1.999 | 0 | 129.75 |
| TET10 (P2) | FEniCSx | 4 | 2187 | 2.9529e-03 | 2.8582e-02 | n/d | n/d | 3.0e-15 | 0.63 |
| TET10 (P2) | FEniCSx | 8 | 14739 | 3.6950e-04 | 7.2446e-03 | 2.998 | 1.980 | 5.1e-15 | 0.99 |
| TET10 (P2) | FEniCSx | 16 | 107811 | 4.6041e-05 | 1.8169e-03 | 3.005 | 1.995 | 3.9e-14 | 5.20 |
| TET10 (P2) | FEniCSx | 32 | 823875 | 5.7477e-06 | 4.5456e-04 | 3.002 | 1.999 | 2.2e-12 | 77.29 |
| TET10 (P2) | scikit-fem | 4 | 2187 | 2.9529e-03 | 2.8582e-02 | n/d | n/d | 2.2e-15 | 0.30 |
| TET10 (P2) | scikit-fem | 8 | 14739 | 3.6950e-04 | 7.2446e-03 | 2.998 | 1.980 | 8.5e-15 | 4.77 |
| TET10 (P2) | scikit-fem | 16 | 107811 | 4.6041e-05 | 1.8169e-03 | 3.005 | 1.995 | 1.0e-12 | 30.97 |
| TET10 (P2) | SfePy | 4 | 2187 | 2.9529e-03 | 2.8582e-02 | n/d | n/d | 2.1e-15 | 0.06 |
| TET10 (P2) | SfePy | 8 | 14739 | 3.6950e-04 | 7.2446e-03 | 2.998 | 1.980 | 6.6e-15 | 3.90 |
| TET10 (P2) | SfePy | 16 | 107811 | 4.6041e-05 | 1.8169e-03 | 3.005 | 1.995 | 1.2e-12 | 8.32 |

#### Convergencia con solución exacta (refinamiento p, malla fija)

| Motor | n | grado p | GDL | error L2 rel. | error H1 rel. | tiempo (s) |
|---|---:|---:|---:|---:|---:|---:|
| NGSolve | 4 | 1 | 375 | 4.4368e-02 | 3.0431e-01 | 0.09 |
| NGSolve | 4 | 2 | 2187 | 2.8738e-03 | 2.8493e-02 | 0.17 |
| NGSolve | 4 | 3 | 6591 | 1.3019e-04 | 1.8394e-03 | 0.49 |
| NGSolve | 4 | 4 | 14739 | 6.0648e-06 | 9.8100e-05 | 0.90 |
| NGSolve | 4 | 5 | 27783 | 2.1128e-07 | 4.1294e-06 | 2.24 |
| NGSolve | 4 | 6 | 46875 | 5.8160e-09 | 1.3846e-07 | 7.23 |
| NGSolve | 4 | 7 | 73167 | 1.5621e-10 | 4.0932e-09 | 18.58 |
| NGSolve | 4 | 8 | 107811 | 3.8196e-12 | 1.0998e-10 | 49.55 |
| FEniCSx | 4 | 1 | 375 | 5.2970e-02 | 3.0234e-01 | 0.03 |
| FEniCSx | 4 | 2 | 2187 | 2.9529e-03 | 2.8582e-02 | 0.05 |
| FEniCSx | 4 | 3 | 6591 | 1.3226e-04 | 1.8118e-03 | 0.21 |
| FEniCSx | 4 | 4 | 14739 | 5.9788e-06 | 9.7022e-05 | 0.82 |
| FEniCSx | 4 | 5 | 27783 | 2.1024e-07 | 4.1098e-06 | 2.53 |
| FEniCSx | 4 | 6 | 46875 | 5.7683e-09 | 1.3806e-07 | 8.13 |
| FEniCSx | 4 | 7 | 73167 | 1.5575e-10 | 4.0877e-09 | 23.43 |
| FEniCSx | 4 | 8 | 107811 | 4.8119e-12 | 1.3012e-10 | 112.06 |
