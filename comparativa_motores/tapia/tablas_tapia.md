## Resultados por motor


### C2 · VOI

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 3219.601 | 91.233 | 203.07 | 0.02722 | 9.3 | 0.55 | 5.9e-09 | 5.8e-09 |
| NGSolve | sparsecholesky (NGSolve) | 3219.601 | 91.233 | 203.07 | 0.02722 | 6.9 | 1.30 | 9.4e-14 | 1.2e-13 |
| FEniCSx | Cholesky MUMPS (PETSc) | 3219.601 | 91.233 | 203.07 | 0.02722 | 5.4 | 0.96 | 3.0e-12 | 3.1e-12 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 3219.601 | 91.233 | 203.07 | 0.02722 | 25.8 | 0.95 | 4.2e-11 | 2.1e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 3219.601 | 91.233 | 203.07 | 0.02722 | 11.3 | 0.62 | 4.1e-11 | 2.1e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 3219.601 | 91.233 | 203.07 | 0.02722 | 7.0 | 1.34 | — | — |

### C2 · Spinodoide

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 1361.423 | 162.621 | 874.71 | 0.12137 | 12.1 | 0.55 | 7.2e-09 | 2.1e-09 |
| NGSolve | sparsecholesky (NGSolve) | 1361.423 | 162.621 | 874.71 | 0.12137 | 3.1 | 0.71 | 1.5e-13 | 1.0e-13 |
| FEniCSx | Cholesky MUMPS (PETSc) | 1361.423 | 162.621 | 874.71 | 0.12137 | 3.6 | 0.67 | 2.2e-11 | 1.8e-11 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 1361.423 | 162.621 | 874.71 | 0.12137 | 29.3 | 0.94 | 6.0e-11 | 6.0e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 1361.423 | 162.621 | 874.71 | 0.12137 | 15.8 | 0.62 | 6.1e-11 | 6.1e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 1361.423 | 162.621 | 874.71 | 0.12137 | 3.0 | 0.72 | — | — |

### C2 · Dual-lattice

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 2570.985 | 109.065 | 509.93 | 0.14528 | 11.8 | 0.55 | 9.1e-09 | 3.7e-09 |
| NGSolve | sparsecholesky (NGSolve) | 2570.985 | 109.065 | 509.93 | 0.14528 | 3.7 | 0.89 | 9.0e-14 | 5.0e-14 |
| FEniCSx | Cholesky MUMPS (PETSc) | 2570.985 | 109.065 | 509.93 | 0.14528 | 3.7 | 0.72 | 6.0e-12 | 9.7e-12 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 2570.985 | 109.065 | 509.93 | 0.14528 | 29.2 | 0.95 | 7.4e-11 | 8.3e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 2570.985 | 109.065 | 509.93 | 0.14528 | 13.0 | 0.62 | 7.6e-11 | 8.5e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 2570.985 | 109.065 | 509.93 | 0.14528 | 3.6 | 0.94 | — | — |

### C4 · VOI

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 3017.587 | 94.585 | 275.38 | 0.04079 | 8.6 | 0.53 | 5.1e-09 | 2.1e-09 |
| NGSolve | sparsecholesky (NGSolve) | 3017.587 | 94.585 | 275.38 | 0.04079 | 4.8 | 1.10 | 5.8e-14 | 1.7e-14 |
| FEniCSx | Cholesky MUMPS (PETSc) | 3017.587 | 94.585 | 275.38 | 0.04079 | 3.6 | 0.78 | 2.8e-12 | 1.7e-12 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 3017.587 | 94.585 | 275.38 | 0.04079 | 23.8 | 0.82 | 5.9e-11 | 3.4e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 3017.587 | 94.585 | 275.38 | 0.04079 | 10.7 | 0.58 | 6.2e-11 | 3.6e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 3017.587 | 94.585 | 275.38 | 0.04079 | 4.6 | 1.15 | — | — |

### C4 · Spinodoide

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 857.553 | 231.901 | 482.88 | 0.12603 | 12.7 | 0.51 | 3.7e-09 | 3.9e-09 |
| NGSolve | sparsecholesky (NGSolve) | 857.553 | 231.901 | 482.88 | 0.12603 | 2.1 | 0.44 | 1.2e-13 | 1.4e-13 |
| FEniCSx | Cholesky MUMPS (PETSc) | 857.553 | 231.901 | 482.88 | 0.12603 | 3.1 | 0.52 | 1.4e-11 | 6.4e-12 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 857.553 | 231.901 | 482.88 | 0.12603 | 28.1 | 0.79 | 3.7e-11 | 4.5e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 857.553 | 231.901 | 482.88 | 0.12603 | 14.1 | 0.57 | 3.7e-11 | 4.5e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 857.553 | 231.901 | 482.88 | 0.12603 | 2.0 | 0.44 | — | — |

### C4 · Dual-lattice

| Motor | Resolvedor | E_app (MPa) | p99 vM (MPa) | máx vM (MPa) | |u| máx (mm) | Tiempo motor (s) | Memoria (GB) | máx Δσ | máx Δ|u| |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| App (spinpy) | AMG (con modos rigidos) + CG | 2392.614 | 117.291 | 270.78 | 0.03232 | 9.0 | 0.52 | 5.2e-09 | 2.1e-09 |
| NGSolve | sparsecholesky (NGSolve) | 2392.614 | 117.291 | 270.78 | 0.03232 | 3.1 | 0.75 | 4.3e-14 | 3.7e-14 |
| FEniCSx | Cholesky MUMPS (PETSc) | 2392.614 | 117.291 | 270.78 | 0.03232 | 3.4 | 0.66 | 3.3e-12 | 2.0e-12 |
| scikit-fem | CG + AMG agregacion suavizada (pya | 2392.614 | 117.291 | 270.78 | 0.03232 | 22.1 | 0.82 | 1.4e-10 | 4.6e-11 |
| SfePy | CG (ls.scipy_iterative) + pyamg co | 2392.614 | 117.291 | 270.78 | 0.03232 | 10.3 | 0.58 | 1.3e-10 | 4.5e-11 |
| NGSolve (directo, ref.) | sparsecholesky (NGSolve) | 2392.614 | 117.291 | 270.78 | 0.03232 | 3.1 | 0.74 | — | — |
