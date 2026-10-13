<p align="right">
  🌐 <a href="README.md">Español</a> · <b>English</b>
</p>

<p align="center">
  <img src="docs/media/banner_en.png" alt="spinpy — spinodal microstructures fitted to micro-CT trabecular bone" width="100%">
</p>

<p align="center">
  <img alt="version V2.2.0" src="https://img.shields.io/badge/version-V2.2.0-e8a33d">
  <img alt="Python 3.10+" src="https://img.shields.io/badge/python-3.10%2B-3776ab?logo=python&logoColor=white">
  <img alt="Windows" src="https://img.shields.io/badge/executable-Windows%2064--bit-0078d6?logo=windows&logoColor=white">
  <img alt="Code MIT" src="https://img.shields.io/badge/code-MIT-2ea44f">
  <img alt="Executable GPL-3.0" src="https://img.shields.io/badge/executable-GPL--3.0-8a8a8a">
  <img alt="Verification: 31 blocks" src="https://img.shields.io/badge/verification-31%20blocks-5b6b7f">
  <img alt="Interface ES/EN" src="https://img.shields.io/badge/interface-ES%20%7C%20EN-5b6b7f">
</p>

<p align="center">
  <b>From micro-CT to mechanical model, saying how much each number is worth.</b><br>
  <a href="#what-it-does">What it does</a> ·
  <a href="#the-application-in-motion">Demo</a> ·
  <a href="#installation">Installation</a> ·
  <a href="#usage">Usage</a> ·
  <a href="#validation-against-published-results">Validation</a> ·
  <a href="#citing">Citing</a>
</p>

---

**spinpy** fits **spinodal** (and *dual-lattice*) microstructures to
**trabecular bone** volumes of interest obtained by micro-computed
tomography, and then measures and tests them: morphometry, elastic
homogenization, finite-element compression testing and a publication-ready
report. It is the Python port of the numerical half of `AppFinal_V2.m`,
verified against analytical solutions and against the original MATLAB
implementation.

Developed on equine sesamoid bones scanned at 51.489 µm, with support also
for porcine, mouse-vertebra and femur volumes.

<table>
<tr>
<td>

```
micro-CT TIFF stack
 → segmentation, PCA framing
 → cubic VOI
 → morphometry
     BV/TV, Tb.Th, DA, Conn.D,
     SMI, Ellipsoid Factor
 → fit: spinodoid
     or dual-lattice
 → elastic homogenization
 → FE compression test
     von Mises, Pistoia
 → report and export
     PDF, STL, VTU, Abaqus, ANSYS, FEBio
```

</td>
<td align="center">
<img src="docs/media/giro.gif" alt="Fitted spinodoid rotating" width="300">
</td>
</tr>
</table>

---

## The application in motion

The animations are screen captures of the interface itself pressing its own
buttons, and they are regenerated with
`python docs/media/hacer_medios.py --idioma en`: they are not mock-ups.

### 1 · Exploring the design space

Density, wave number and the three cone angles define the structure. Every
change is regenerated on the spot, from the isotropic class to the columnar,
the lamellar or one fitted to a real bone.

<p align="center"><img src="docs/media/01_explorar_en.gif" alt="Exploring spinodoid classes in the viewer" width="100%"></p>

### 2 · Fitting to the micro-CT VOI

The VOI is loaded (`.vtk`, `.mat`, a folder of TIFF slices or a multipage
TIFF), the microstructure is fitted by a staged search, and the table
compares bone and candidate metric by metric. Both views rotate in sync.

<p align="center"><img src="docs/media/02_ajuste_en.gif" alt="Fitting a spinodoid to an equine VOI" width="100%"></p>

### 3 · Compression testing

Finite-element test along X, Y, Z or all three axes, with effective-strain
and von Mises stress fields, **on the same colour scale for the bone and the
candidate**, and failure load by the Pistoia criterion.

<p align="center"><img src="docs/media/03_mecanica_en.gif" alt="Von Mises field of the compression test" width="100%"></p>

---

## What it does

| Stage | What it includes |
|---|---|
| **Reads** | `.vtk` and `.mat` VOIs, and micro-CT TIFF stacks (folder of slices or multipage) with the scale read from the SkyScan *log*. Thresholds, frames by PCA and crops the cube. |
| **Generates** | Spinodoids from a Gaussian random field (Kumar et al. 2020), with both wave-sampling conventions found in the literature, and *dual-lattice* (Vafaeefar et al. 2022). |
| **Measures** | BV/TV, BS/BV, BS/PV, Tb.Th, Tb.Sp, Tb.N, MIL tensor and DA, Conn.D, SMI, local thickness, pore size, **Ellipsoid Factor** and principal curvatures of the interface. |
| **Fits** | The parameters that best reproduce a real VOI: staged search, replicas with fresh seeds, self-consistent noise floor and optional mechanical tie-break. |
| **Homogenizes** | The elastic tensor by periodic unit cell on the voxel grid, with algebraic multigrid and a residual check. |
| **Tests** | Compression along X, Y or Z with von Mises, effective-strain and displacement fields; Pistoia failure; mesh-convergence study. On bricks or on a smooth mesh of quadratic tetrahedra, linear or nonlinear, with built-in FE engines (the app, NGSolve and, if installed, FEniCSx, scikit-fem and SfePy), one or several at once. E and p99 **corrected for boundary artefacts** (rigid platen and measurement in the VOI core), validated against the same bone embedded in bone. |
| **Simulates** | *In silico* bone loss (thinning, thin-trabeculae loss, disuse, recovery) and progressive failure on the VOI's digital twin. |
| **Reports** | A publication report in Spanish and English (Markdown and PDF): methods written from the values actually used, 600 dpi figures, a citability table, a reporting checklist (Bouxsein et al. 2010 for micro-CT, Erdemir et al. 2012 for finite elements) and a SHA-256 fingerprint to reproduce each mask. |
| **Exports** | Hexahedral or TET10 solid to Abaqus, ANSYS APDL, VTU and STL, and the full compression test to FEBio 4 (`.feb`). |
| **Batch processing** | Decomposes variance between and within specimens, with ICC, effective N and equivalence by TOST. |

All of it with a graphical interface (`visor.py`) or from Python.

---

## What sets it apart

Most bone-morphometry tools give a number. This one also says **how much
that number is worth**:

- **Verification against manufactured solutions**, not only against another
  program: spheres, cylinders and plates with a known analytical answer,
  Backus averaging for the elastic tensor, and the MIL tensor against
  geometries of exact anisotropy.
- **Biases are measured and reported next to the number.** Local thickness
  underestimates by ~30 % at 3-4 voxels; SMI carries a −15 % discretization
  bias and is confounded by concavity; MIL has a noise floor of DA ≈ 1.07.
  The interface shows this where the value is read.
- **What is not citable is said so.** The maximum of a von Mises field does
  not converge with the mesh, so the report cites the p99 of the surface
  layer and marks every value as *citable*, *with reservations* or *not
  citable*.
- **The generator is stochastic, and that is measured**, not assumed: the
  same parameters with K seeds give the noise floor below which no difference
  means anything.
- **Replicas do not inflate N.** The statistics module computes the
  effective N and puts it before the number of rows, because treating K
  replicas as K specimens is pseudoreplication.
- **Every saved number carries its provenance**: version, family, scheme,
  seed and the wave number with both its readings, so that an exported result
  can be regenerated bit for bit.

---

## Installation

### Without Python: Windows executable

For those who only want to **use** the application there is a 64-bit
executable that needs neither Python nor any dependency. It installs into the
user's folder and **does not ask for an administrator password**, which is
what a university computer requires.

| Executable | |
|---|---|
| Installed | ~650 MB |
| Installer / zip | ~220-300 MB |
| Saves to | `Documents\spinpy` (not deleted on uninstall) |
| Licence | **GPL-3.0** (see below) |

If something does not work on a given computer, the start menu has a
shortcut **"spinpy - autocomprobación"** (self-test): it exercises every
computation path one by one and leaves a report in
`Documents/spinpy/autocomprobacion.txt`. That is what to send for diagnosis.
It can also be launched from the console:

```
spinpy.exe --autocomprobacion
```

To rebuild it, `instalador/construir.ps1` runs the whole process, and
`instalador/LEEME_DESARROLLO.md` explains the packaging pitfalls (in
Spanish).

### macOS test build (Apple Silicon)

Since V2.1.3 each release also attaches `spinpy-V*-macos-arm64-PRUEBA.zip`
with `spinpy.app`, for Macs with an Apple Silicon chip (M1 or later). **It is a
test build:** it builds and passes the `.app` self-test on a GitHub Actions
Apple Silicon runner (including the three built-in FE engines), but the
graphical interface has not been tried by hand on a Mac with a screen. It is
neither signed nor notarised: the first time, right-click `spinpy.app` and
choose "Open" (or `xattr -dr com.apple.quarantine spinpy.app`). It uses Qt
5.15.17 instead of the 5.15.2 used on Windows, which does not exist for Apple
Silicon. There is no build for Intel Macs. Instructions are in
`instalador/LEEME_MAC.txt` (in Spanish); to rebuild it, use
`instalador/construir_mac.sh`.

### With Python: as a library

```
pip install -e .
```

The core only needs `numpy`, `scipy`, `scikit-image` and `pyvista`. The
extras are installed as needed:

```
pip install -e ".[elastico]"   # multigrid: homogenization and FE test
pip install -e ".[malla]"      # TET10 meshing
pip install -e ".[lote]"       # batch and statistics
pip install -e ".[gui]"        # graphical interface
pip install -e ".[informe]"    # report PDF and figures
pip install -e ".[todo]"       # everything
```

`pyamg` is listed as an extra, but it is not optional in practice: without
it the solver falls back to conjugate gradient with Jacobi, which does not
converge with the 10⁻⁶ contrast between bone and void.

---

## Usage

Graphical interface:

```
python visor.py
```

The panels follow the order of the work: reference VOI, microstructure,
visualization, morphometry, fitting, mechanical analysis, *in silico*
simulations, export and batch. The **"Publication report (auto)…"** button
chains every stage and estimates beforehand how long it will take.

### Language

The interface is in **Spanish and English** and is switched from the
*Language* menu without rebuilding the window: the loaded VOI, the metrics
and the fields are kept. Metric names (BV/TV, Tb.Th, DA, SMI…) are
deliberately not translated: they are international notation (Bouxsein et
al. 2010).

After touching any text, run the dictionary checker:

```
python idioma_revisar.py
```

### From Python

The API is in Spanish, like the rest of the code:

```python
from spinpy import leer_voi, morfometria, ajustar_spinodoide

VOI, spacing = leer_voi("VOI_medio_cubico.vtk")
m = morfometria(VOI, spacing, extra=True)
print(m["BVTV"], m["TbTh"], m["DA"], m["ConnD"])

r = ajustar_spinodoide(VOI, spacing, modo="completo")
print(r["parametros"], r["error"])
```

### Command line

```
spinpy-lote folder_with_VOIs --replicas 20 --salida results
spinpy-informe resultados_sesion.json        # rebuilds the report from a JSON
spinpy-validar-exp predicciones.csv          # comparison with physical tests
```

### Quick start with open data

The repository ships ten porcine subchondral bone VOIs in
[`datos/porcino_koria2020/`](datos/porcino_koria2020/) (C1 to C5 and V1 to V5;
3.0 mm cube, 188³ voxels of 16 µm; Koria, Mengoni and Brockett 2020,
[doi:10.5518/787](https://doi.org/10.5518/787), CC BY 4.0), so the installation
can be checked without your own data. Each file has its fingerprint in
`SHA256SUMS.txt`.

```
python ejemplos/inicio_rapido.py                    # VOI_V1, central 64³ cube
python ejemplos/inicio_rapido.py --voi C3 --lado 96
python ejemplos/inicio_rapido.py --voi V1 --lado 187   # whole VOI
```

The script verifies the file fingerprint, measures the morphometry and runs a
compression test along z (sliding support, E_s = 20 GPa, 1 MPa), reporting
the apparent modulus and the Pistoia failure load. With the 40³ cube
(`--rapido`) it takes a few seconds. The central crop is meant for learning
the workflow and checking the installation; use the whole VOI to quote values.

---

## Validation against published results

The method `spinpy` implements **is not ours**: it is that of Kumar et al.
(2020). The `Test/` folder exists to prove it, placing the published figure
next to the same figure regenerated with this code, and evaluating
falsifiable predictions from each article:

| Script | What it checks |
|---|---|
| `replicar_kumar2020.py` | The method: cone union, level set and the four classes. |
| `replicar_zheng2021.py` | Bounds: every elastic surface within Voigt and Hashin-Shtrikman, and the ρ ≥ 0.3 curve. |
| `replicar_guo2024.py` | Curvature: the (κ₁, κ₂) profile against a nodal surface with exact curvatures. |

They are opened from the application in the **Validation** menu, or from
the console:

```
python Test/replicar_kumar2020.py            # full, ~7 min
python Test/replicar_kumar2020.py --rapido   # ~1.5 min
```

> Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
> spinodoid metamaterials. *npj Computational Materials, 6*, 73.
> https://doi.org/10.1038/s41524-020-0341-6 — open access under CC BY 4.0,
> which is what allows its figure to be reproduced here with attribution.
> Provenance of every reference figure in `Test/referencia/FUENTES.md`.

---

## Verification

```
python -m pytest tests/ -q
```

The suite runs on every pull request and every change to `main`
(`.github/workflows/pruebas.yml`: Linux with Python 3.11 and 3.13, Windows and
macOS with 3.13), and building the installer requires it to pass. Each run
stores as an artifact the tolerance log and the exact environment
(`pip freeze`) it was measured with.

31 blocks with tolerances **declared before measuring**: topology, SMI,
thickness, mechanics, Pistoia, TIFF stacks, curvature, Ellipsoid Factor,
*dual-lattice*, provenance, objective function, MIL sampling, report, von
Mises surface layer… A failure here is a finding, not a bug in the suite.
The two findings already documented (the slab of block 04 comes out one voxel
thicker, and the Gibson-Ashby exponent of block 05 comes out ≈ 3.9) are
marked as strict `xfail`: the suite stays green, the log still records them
as failed and, should either start passing, the suite fails to force a review
of its documentation.

Against the original MATLAB implementation (`validar_*.py`, reference
outputs in `resultados/`):

| What is checked | Result |
|---|---|
| Morphometry on real VOIs | 5.9 × 10⁻¹⁴ % in everything that does not go through marching cubes |
| Elastic tensor against Backus averaging | 2.7 × 10⁻¹⁵ relative difference |
| Error function, over 1173 pairs | 3.55 × 10⁻¹⁶ |
| Fit search grids | exact |
| von Mises under uniaxial stress | 1.000000000000 |

---

## Comparison with BoneJ

Checked against BoneJ 7.2.2 on the same volumes:

- **BV/TV agrees to machine precision**: both tools see the same mask.
- **BS differs by up to 46 %, and it is a difference of convention**: BoneJ
  closes the surface against the volume boundary and spinpy does not,
  because the cube faces are the VOI's artificial cut. With the convention
  matched, the difference drops to ±6 %.
- **DA uses different scales** in the two tools, although both agree on the
  ranking of the specimens.
- **SMI cannot be compared**: BoneJ withdrew it in version 2, for the same
  criticism of confounding with concavity that spinpy documents.

And something that changes how to read all of the above: **segmentation
uncertainty is an order of magnitude larger than the differences between
implementations.** Moving the threshold by ±15 % moves Tb.Th by 62 %, against
the ±6 % separating the two tools.

---

## Comparison with FEBio

> Validation done with V1.0.x. Since V2.0.1 the application no longer uses
> FEBio (see "Built-in FE engines" below); these campaigns remain the external
> validation of the app's solver, and the engine comparison reproduced their
> brick-mesh cavity figures to 10⁻⁸.

Every mechanical analysis of the application has been solved again,
independently, with **FEBio 4.5** (Maas et al. 2012), the reference
finite-element program in biomechanics. Two campaigns:

| | structures | what was checked | report (Spanish) |
|---|---|---|---|
| **Equine H4** | 3 VOIs (BV/TV 0.28–0.77) at 32³ and 48³ + a spinodoid | the compression test | [`comparativa_febio/INFORME.pdf`](comparativa_febio/INFORME.pdf) |
| **Porcine V1** | the VOI and **both candidates fitted to it** in one session (spinodoid and *dual-lattice*) | **every** mechanical analysis of the app | [`comparativa_febio/porcino/INFORME.pdf`](comparativa_febio/porcino/INFORME.pdf) |

### The same problem, two programs

For the comparison to measure anything, both programs must solve **the same
discrete problem**, not a similar one. The app's exporter
(`escribe.escribir_febio`) writes the same mesh —one hexahedron per
load-bearing voxel, node for node—, the same support, the same load over the
gross section and the same material; FEBio solves it with Newton and a direct
factorization (Pardiso), and spinpy with algebraic multigrid and conjugate
gradient. Derived quantities (von Mises p99 on the surface layer, Pistoia
factor) are computed with the **same spinpy functions** on both fields.

FEBio is geometrically nonlinear; spinpy is linear. The linear problem is
compared by extrapolating FEBio to zero load, `u = 2·u(1 kPa) − u(2 kPa)`,
which cancels the first-order nonlinear term. And the difference is put to
use: FEBio is also run **at the test load and at the failure load**, which
measures how far the app's linear hypothesis holds.

### Porcine VOI and its two candidates: every analysis

Subchondral bone of the porcine talus (specimen V1 of Koria, Mengoni and
Brockett 2020, [doi:10.5518/787](https://doi.org/10.5518/787), CC BY 4.0;
188³ voxels of 16 µm, BV/TV 0.40). The three structures are **those of the
automatic-report session**: the VOI is checked by SHA-256, and each candidate
is regenerated from `reproduccion.json` with its fingerprint required bit for
bit.

| app analysis | configuration | worst spinpy/FEBio difference |
|---|---|---|
| Compression test + Pistoia | 40³, sliding, 20 GPa, 1 MPa | 4·10⁻⁸ |
| Compared analysis (Tapia et al.) | 40³, fixed, 18 GPa, 100 N | 4·10⁻⁸ |
| Periodic elastic tensor | 32³, 6 load cases | 2·10⁻¹⁰ |
| Mesh convergence | 22, 28, 34 and 40³ | 2·10⁻⁹ |
| Progressive failure | 32³, 10 steps, full loop with FEBio as the solver | **0 differing voxels** in 33 steps |

<p align="center"><img src="comparativa_febio/porcino/figs/en/fig_coincidencia.png" alt="Relative difference between spinpy and FEBio for each quantity and analysis" width="80%"></p>

All differences lie between 10⁻¹¹ and 10⁻⁷, within the tolerances declared
before the first comparison with FEBio. Three details make the result
stronger than it looks:

- **The periodic tensor too.** FEBio has no periodic conditions for a voxel
  mesh; they were written with its linear constraints
  (`u(x') = u(x) + E·(x' − x)`, eliminated from the system, not penalized),
  and the stiffness comes out of FEBio by **another route** —volume average
  of stress— than in spinpy —cell energy—.
- **Progressive failure was repeated in full.** It was not compared step by
  step on spinpy's damage: at each step FEBio decides which tissue breaks from
  its own field and builds the next step on its own damage. In the three
  structures and the 33 steps, both loops break **exactly the same
  elements**, even though the failure rule compares each element with a
  percentile, and a small error in the wrong place would have separated the
  two histories for good.
- **The session reproduces.** The regenerated structures give the numbers
  the automatic report saved, with differences ≤ 3·10⁻¹¹.

#### The colour maps

The maps are drawn with the same recipe as figure 8 of the app's report and
with **a single colour scale** for app and FEBio. Columns: the app; linear
FEBio; their difference (log scale); **nonlinear** FEBio at the test load;
and its difference from the app.

<p align="center"><img src="comparativa_febio/porcino/figs/en/mapa_vm_comparado.png" alt="Von Mises from the app and from FEBio on the VOI, the spinodoid and the dual-lattice, Tapia protocol, 100 N" width="100%"></p>

The first two columns are the same image: the difference does not exceed
3·10⁻⁷ of the p99 in any element (the ceiling is set by storing in
`float32`; in double precision it is 2·10⁻⁸). The fifth is the only one that
changes: there FEBio lets the part actually deform, and what shows in red are
trabeculae on the spinodoid's loaded face that bend more than a linear
calculation predicts.

<p align="center"><img src="comparativa_febio/porcino/figs/en/mapa_desp_comparado.png" alt="Total displacement from the app and from FEBio on the three structures" width="100%"></p>

<p align="center">
  <img src="comparativa_febio/porcino/figs/en/fig_paridad_comparado.png" alt="Per-element von Mises, spinpy versus FEBio" width="100%">
  <img src="comparativa_febio/porcino/figs/en/fig_colas_superficie.png" alt="Distribution of von Mises on the surface layer, spinpy and FEBio" width="100%">
</p>

#### Tensor, convergence and progressive failure

<p align="center"><img src="comparativa_febio/porcino/figs/en/fig_tensor.png" alt="Engineering constants of the periodic tensor, spinpy and FEBio" width="85%"></p>

<p align="center">
  <img src="comparativa_febio/porcino/figs/en/fig_convergencia.png" alt="Mesh convergence of E_app in spinpy and FEBio" width="55%">
</p>
<p align="center">
  <img src="comparativa_febio/porcino/figs/en/fig_fallo.png" alt="Progressive failure with spinpy and with FEBio as the solver" width="100%">
</p>

#### How far linearity holds

| structure | app E_app (MPa) | nonlinear deviation at 1 MPa | Pistoia failure load (MPa) | nonlinear deviation at that load | E_app rigid platen / force |
|---|---|---|---|---|---|
| VOI | 4021 | −0.07 % | 21.8 | −1.6 % | 1.06 |
| Spinodoid | 1725 | −0.70 % | 12.5 | **−8.0 %** | **1.69** |
| Dual-lattice | 2699 | −0.22 % | 15.4 | −3.4 % | 1.23 |

<p align="center"><img src="comparativa_febio/porcino/figs/en/fig_rigidez.png" alt="Apparent modulus in spinpy and FEBio, linear, nonlinear and with a rigid platen" width="85%"></p>

- **On the VOI, the app's linear test holds up to its failure load**
  (−1.6 %; the Pistoia criterion is met at 2.07 % of the tissue instead of
  2 %).
- **On the spinodoid, not entirely.** Near its failure load it loses 8 % of
  its stiffness and its surface p99 rises by 8.7 %; with the 100 N of the
  Tapia protocol (11 MPa over a 3 mm side), the linear p99 falls 7.7 % short.
- **The cause is the loading condition, as in H4.** With force imposed on
  the bone of the top face, the trabeculae cut by the loaded face act as
  cantilevers. A rigid platen reduces the spinodoid's deviation to −2.6 % and
  raises its linear modulus by 69 %.

#### What the comparison says about the candidates

With both programs in agreement, the differences between structures belong
to the structures. With BV/TV matched to within 1 %, **the *dual-lattice*
has 0.67 times the VOI's stiffness and the spinodoid 0.43**; the spinodoid
doubles the bone's surface p99 and loses 78 % of its stiffness when the first
2 % of the tissue breaks (the VOI, 24 %). And a warning the app did not give:
**the candidates are not mesh-converged the way the VOI is** —between 34³ and
40³ the VOI moves −0.8 %, the *dual-lattice* −1.9 % and the spinodoid
−5.5 %—, so their stiffness at 40³ must be cited with reservations.

### Equine H4 VOIs

The compression test on three bone VOIs (BV/TV 0.28 to 0.77, at 32³ and
48³) and a fitted spinodoid:

- **The computation agrees to seven or eight digits**, within tolerances
  declared before measuring.
- **On the most porous VOI (BV/TV 0.28) the nonlinear response at 1 MPa
  departs by 16 % at 32³ and by 48 % at 48³**; on those with BV/TV 0.55 and
  0.77, by 0.06 %. With a rigid platen the deviation drops to 0.8 %, and the
  linear modulus itself is 2.1 times higher. The three VOIs come from a single
  specimen: it is one measured case, not a general rate. On very porous VOIs,
  the apparent modulus and the failure load must be cited stating the loading
  condition.

### What this comparison does not prove

Agreeing with FEBio proves that spinpy **solves the problems it poses
correctly**, not that those problems represent real bone. The voxel mesh,
the homogeneous isotropic tissue, the boundary conditions and the Pistoia
criterion are assumptions both programs share; only a physical test can
check them.

### Reproducing

```
python comparativa_febio/comparar_febio.py          # H4, ~90 min
python comparativa_febio/porcino/validar_porcino.py # porcine, every analysis, ~2 h
python comparativa_febio/porcino/figuras_porcino.py # figures, only reads results
python -m pytest tests/test_25_febio.py             # small version, ~10 s
```

Requires FEBio 4.5 (FEBio Studio 2). Results are left in
`comparativa_febio/**/resultados/*.jsonl` with their provenance block; the
`.feb` files and FEBio outputs are regenerated and not versioned.
`tests/test_25_febio.py` is skipped if FEBio is not installed.

> Maas, S. A., Ellis, B. J., Ateshian, G. A., & Weiss, J. A. (2012). FEBio:
> finite elements for biomechanics. *Journal of Biomechanical Engineering,
> 134*(1), 011005. https://doi.org/10.1115/1.4005694

---

## Built-in FE engines: no external programs

Since V2.0.1 every test is solved **inside spinpy**. FEBio is no longer used
or searched for; `.feb` export is kept as an exchange format. In
**Mechanical analysis**, "FE engine" chooses who solves, and "Mesh" and
"Analysis" what is solved:

| engine | meshes | analyses | in the Windows installer |
|---|---|---|---|
| **App (spinpy)** | bricks (hex8) | linear | yes (the usual, validated solver) |
| **NGSolve** (recommended) | hex8 and TET10 | linear and nonlinear (St. Venant-Kirchhoff, neo-Hookean) | yes |
| scikit-fem | hex8 and TET10 | linear and nonlinear | yes |
| FEniCSx | hex8 and TET10 | linear and nonlinear | no (conda-forge only; compiles C at run time) |
| SfePy | hex8 and TET10 | linear; nonlinear SVK only | no (no Windows wheels) |

"Compare engines…" (panel) and **"FEM report (Auto)…"** (top bar) solve the
same mesh with **several engines at once** and give a comparison table
between engines and against the app, which also goes into the publication
report. The automatic report has its own "FE engines (comparison table)" row.
Each solve runs in a child process: "Stop" kills it and running out of memory
does not close the window. The same computation runs from the command line:

```bash
python -m spinpy.fem VOI.vtk --protocolo app tapia2026 --malla hex8 tet10 --motores app ngsolve
```

### Report: engine comparison

📄 **[Full report (PDF, Spanish)](comparativa_motores/INFORME.pdf)** ·
[Markdown version](comparativa_motores/INFORME.md) ·
[tables](comparativa_motores/tablas.md) ·
[data](comparativa_motores/resultados/)

| result | value |
|---|---|
| solid block, E_app = E_s (hex8 and TET10, all five engines) | ≤ 10⁻¹¹ |
| nonlinear compression up to 20 %, against the closed-form solution | ≤ 4·10⁻¹² |
| cavity on bricks, engines against FEBio 4.5 (E_app, p99) | 7.7·10⁻⁹ · 1.7·10⁻⁹ |
| spinodoid, displacements between engines (direct solvers) | 10⁻¹³ to 10⁻¹¹ |
| TET10 at 790,599 DOF: NGSolve (Cholesky) · FEniCSx (MUMPS) · scikit-fem · SfePy | 45 s · 22 s · 905 s · 1010 s |
| NGSolve: its own Cholesky against MKL PARDISO | 1.4-1.6× faster, 2.3-2.8× less memory |

<p align="center"><img src="comparativa_motores/figs/fig3_tiempo_tet10.png" alt="Solution time of the spinodoid with a smooth TET10 mesh against degrees of freedom, per engine and solver" width="90%"></p>

Three measured reliability findings drove the configuration: SfePy's St.
Venant-Kirchhoff tangent is not the derivative of its residual (37 % error,
linear Newton convergence); FEniCSx CG + GAMG diverges on the TET10 mesh of a
real spinodoid; and scikit-fem CG + pyamg stops at a 10⁻¹⁰ residual with a
displacement error of up to 0.9 % on those meshes. Hence NGSolve solves the
smooth mesh with BDDC or its Cholesky, not with a generic algebraic
multigrid. The work also found and fixed two bugs in the app's TET10 meshing
(flat tetrahedra on the cube faces that aborted the mesh, and mid-side nodes
off their edge midpoint after snapping to the planes).

### Background: smooth mesh (TET10) against bricks, with FEBio

📄 [Report (PDF)](comparativa_febio_tet/INFORME.pdf) ·
[Markdown version](comparativa_febio_tet/INFORME.md) ·
[data](comparativa_febio_tet/resultados/cavidad.jsonl)

With FEBio 4.5 (V1.0.x) it was established that the smooth mesh does not make
the von Mises peak converge on the spherical cavity (stiffness is stable) and
that on real bone at 32³ two-voxel struts come out too soft: the useful
comparison is at 48³ or finer. Those conclusions concern the mesh, not the
program, and still hold with the built-in engines.

## Compute times and the Tapia protocol in every engine

📄 **[Full report (PDF, in Spanish)](comparativa_motores/INFORME_TIEMPOS_Y_TAPIA.pdf)** ·
[timings (Markdown, in Spanish)](comparativa_motores/tiempos_sesion/INFORME_TIEMPOS.md) ·
[Tapia protocol (Markdown, in Spanish)](comparativa_motores/tapia/INFORME_TAPIA_SOLVERS.md) ·
[timing data](comparativa_motores/tiempos_sesion/) ·
[Tapia protocol data](comparativa_motores/tapia/resultados/tapia.jsonl)

### Run time of each engine

The engine-comparison campaign was repeated (95 runs, 4 threads, one child
process per run, a single repetition per combination) on a 2.10 GHz Xeon with
4 CPUs and 15 GB. Total engine time, in seconds:

| problem | DOF | App | NGSolve | FEniCSx | scikit-fem | SfePy |
|---|---:|---:|---:|---:|---:|---:|
| spinodoid hex8 64³, direct | 365,655 | 47.5 (CG + AMG) | 9.2 | 9.9 | n/a | n/a |
| spinodoid hex8 64³, iterative | 365,655 | 47.5 | 50.9 | 43.3 | 130.0 | 55.1 |
| spinodoid TET10 48³, direct | 790,599 | n/a | 32.0 | 18.5 | n/a | n/a |
| spinodoid TET10 48³, iterative | 790,599 | n/a | 43.1 | 225.0 | 672.7 | 806.1 |
| spinodoid 24³ nonlinear (SVK) | 30,672 | n/a | 6.3 | 3.2 | 110.4 | failed |

n/a: not run, or not solved by that engine (the app does not solve TET10; the
direct solvers of scikit-fem and SfePy were not run above 1.5·10⁵ DOF). The
ranking of the engines matches the published campaign; this session's times
were lower (median 0.73 times the published ones, n = 46 pairs). NGSolve used
its own Cholesky, not PARDISO, because the PyPI build no longer ships MKL.

<p align="center"><img src="comparativa_motores/tiempos_sesion/t_escala_hex.png" alt="Time against degrees of freedom for the brick-mesh spinodoid, by engine and solver" width="90%"></p>

### Tapia et al. (2026) protocol solved with every engine

The protocol (brick mesh at 40³, fixed base, E_s = 18 GPa, ν = 0.3, 100 N axial
load, linear elasticity) was applied to two real porcine micro-CT VOIs (C2 and
C4) and to their fitted spinodoid and dual-lattice, six structures in all, with
the five engines. What is replicated is the test protocol and the type of
output, not the original study (cervid vertebrae, 0.05 mm tetrahedra in
ANSYS): absolute values are not comparable with the published ones, but they
are comparable between engines.

| structure | E_app (MPa) | p99 von Mises (MPa) | max. von Mises (MPa) | max. ‖u‖ (mm) |
|---|---:|---:|---:|---:|
| C2 · VOI | 3219.6 | 91.2 | 203.1 | 0.0272 |
| C2 · spinodoid | 1361.4 | 162.6 | 874.7 | 0.1214 |
| C2 · dual-lattice | 2571.0 | 109.1 | 509.9 | 0.1453 |
| C4 · VOI | 3017.6 | 94.6 | 275.4 | 0.0408 |
| C4 · spinodoid | 857.6 | 231.9 | 482.9 | 0.1260 |
| C4 · dual-lattice | 2392.6 | 117.3 | 270.8 | 0.0323 |

The five engines gave these same figures. Element by element, the largest von
Mises difference against NGSolve with a direct solver (as a fraction of the
p99) was 3.7·10⁻⁹ to 9.1·10⁻⁹ with the app, 3.7·10⁻¹¹ to 1.4·10⁻¹⁰ with
scikit-fem and SfePy, 2.8·10⁻¹² to 2.2·10⁻¹¹ with FEniCSx and 4.3·10⁻¹⁴ to
1.5·10⁻¹³ with NGSolve (which uses the same algorithm as the reference, so it
only measures round-off). The ordering matches that of the solver tolerances.
Total deformation gives differences of similar size. The maps from the five
engines are indistinguishable by eye.

<p align="center"><img src="comparativa_motores/tapia/figs/mapa_vm_solvers.png" alt="Von Mises stress from each engine in the six structures, Tapia protocol" width="100%"></p>

<p align="center"><img src="comparativa_motores/tapia/figs/mapa_u_solvers.png" alt="Total deformation from each engine in the six structures, Tapia protocol" width="100%"></p>

<p align="center"><img src="comparativa_motores/tapia/figs/resumen_diferencias.png" alt="Largest per-element difference of each engine against direct NGSolve" width="100%"></p>

With a common mesh and converged solvers, the choice of engine does not change
the fields of the protocol; the differences between structures exceed those
between engines by several orders of magnitude. Limits: two VOIs out of four,
one fit and one seed per VOI and family, brick meshes at 40³ only, no
nonlinearity or mesh-convergence study, and solver tolerances not equalised.
The VOIs are not distributed with the repository.

---

## Documentation

The mechanical validation against FEBio is in
[`comparativa_febio/porcino/INFORME.pdf`](comparativa_febio/porcino/INFORME.pdf),
and the comparison of the built-in FE engines in
[`comparativa_motores/INFORME.pdf`](comparativa_motores/INFORME.pdf).

`docs/MANUAL_spinpy.pdf` documents every module and every function, and is
generated from the code itself (`python docs/generar_manual.py`). The
docstrings of this project are not a summary: they carry the design
decisions, the measured biases and the pitfalls that took effort to find.
Code, docstrings and reports are in Spanish.

---

## Experimental validation: status and procedure

**Current status.** The absolute mechanical values from spinpy (apparent
modulus, failure load, stresses) are simulation results verified numerically
against analytical solutions, FEBio and five finite element engines. **They
have not been compared with mechanical tests of the bone the VOIs come from.**
Until that comparison exists, they should be read as relative comparisons
between structures under one loading condition, not as predictions of the
strength or stiffness of a real specimen. The morphometry, in contrast, is
checked against manufactured solutions and against BoneJ (previous section).

**Available procedure.** `spinpy.validacion_exp` (and the `spinpy-validar-exp`
command) compares predicted with measured values for one or more quantities
and computes: root mean square error, mean absolute relative error and bias;
Lin's concordance correlation coefficient with its confidence interval;
regression with confidence intervals for slope and intercept; Bland-Altman on a logarithmic scale (the
measured/predicted ratio, because trabecular mechanics errors are
multiplicative); Spearman correlation; and a calibration of the tissue modulus
E_s with leave-one-out cross-validation, which separates a scale error (a
wrong E_s, correctable) from a shape error (the geometry does not explain the
differences between specimens).

```
python ejemplos/predecir_para_ensayo.py --salida predicciones.csv   # empty 'medido' column
# fill 'medido' with the physical test of the same specimen
spinpy-validar-exp predicciones.csv --ccc-min 0.9 --mape-max 0.25 \
    --pendiente-delta 0.15 --Es-usado 20000 --json contraste.json
```

Rules of use. Acceptance criteria are declared **before** looking at the
measured data; without criteria the verdict is "no evaluado" (not evaluated),
and a criterion that cannot be computed gives "no evaluable", never a silent
"does not meet". The CCC criterion compares the lower limit of its interval,
and the slope criterion requires **equivalence** (the whole interval within
1 ± δ), because a wide interval from few data does not demonstrate agreement.
The same criteria apply to every quantity in the CSV.
The test conditions (tissue modulus, support, tested region) must be
reproduced in the prediction. With several VOIs per animal the pairs are not
independent and the intervals come out too narrow: use one pair per animal or
a mixed model. With fewer than 8 pairs the module warns that the intervals
cannot support a claim of equivalence. The template is at
[`datos/porcino_koria2020/ensayos_experimentales_PLANTILLA.csv`](datos/porcino_koria2020/ensayos_experimentales_PLANTILLA.csv).

---

## Limitations worth knowing

- **No experimental validation of the absolute mechanical values** (see the
  previous section).
- **Mesh-based measurement on Apple Silicon macOS.** In the automated
  verification, `tests/test_07_malla.py::test_losa_tapas_excluidas` fails on
  macOS (BS +2.34 % against the 1 % tolerance) and passes on Linux and Windows.
  On a symmetric slab, two of the four lateral faces lose cap area only on ARM,
  which points to the mesh repair (PyMeshFix) or the VTK clipping. Until this is
  clarified, do not use the mesh mode on macOS for values you intend to quote;
  the voxel mode is not affected.
- **Segmentation is a global threshold.** It is the largest source of
  uncertainty in the whole process, which is why the application declares
  which one was applied instead of hiding it.
- **The spinodal family has a limit of its own**: a single characteristic
  length gives a uniform trabecular thickness, and bone is not uniform. A
  spinodoid matching BV/TV, Tb.Th and DA can be considerably softer than the
  bone it was fitted to; spinpy measures it and says so.
- The Pareto, Bayesian and MOBO optimizers remain in MATLAB only.
- Below ρ ≈ 0.25 (isotropic class) periodic homogenization does not reach the
  declared tolerance: it is a property of the regime near the rigidity
  threshold, and the application reports it.

---

## Citing

If you use this software, cite it with the metadata in
[`CITATION.cff`](CITATION.cff).

**Authors:** Carlos González-Torres
([ORCID](https://orcid.org/0000-0002-7765-6388)), David Ortiz-Puerta
([ORCID](https://orcid.org/0000-0001-6285-3066)) and Mauricio A.
Sarabia-Vallejos ([ORCID](https://orcid.org/0000-0001-5128-796X)) —
Universidad de Valparaíso.

Methods implemented: Kumar et al. 2020 (spinodal generator), Vafaeefar et
al. 2022 (*dual-lattice*), Andreassen & Andreasen 2014 (homogenization),
Harrigan & Mann 1984 (MIL tensor), Pistoia et al. 2002 (failure criterion),
Odgaard & Gundersen 1993 (Conn.D), Hildebrand & Rüegsegger 1997 (local
thickness and SMI), Doube 2015 (Ellipsoid Factor), Parfitt (plate model).

## Licence

**The source code is MIT.** See [`LICENSE`](LICENSE). Anyone importing the
`spinpy` package from a script or a notebook takes on no copyleft
obligation: the core does not depend on Qt.

**The Windows executable is GPL-3.0**, because it includes PyQt5, which is
GPL v3. It is a combined work and its redistribution is subject to that
licence; use, on the other hand, is unrestricted, and what you produce with it
is yours. The details, with the list of third-party components, are in
`instalador/LICENCIA_BINARIO.txt`.

Migrating the interface to PySide6 (LGPL) would not by itself be enough for a
binary without copyleft: `pymeshfix` is GPL v3 and `tetgen` derives from AGPL
code.
