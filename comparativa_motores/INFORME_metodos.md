## 1. Objetivo

Hasta la versión 1.0.2, spinpy resolvía sus ensayos de compresión con dos programas: su propio resolvedor sobre hexaedros de un vóxel (`resistencia.ensayo_compresion`) y FEBio 4.5, un ejecutable externo que se lanzaba como proceso aparte para la malla suave de tetraedros cuadráticos (TET10) y para los análisis no lineales. Esa segunda vía hacía a la aplicación dependiente de un programa cuya licencia no permite redistribuirlo. Este informe evalúa cuatro bibliotecas de elementos finitos que pueden correr dentro del proceso de spinpy, sin programas externos, y las compara con el resolvedor actual de la aplicación en cinco ejes: exactitud frente a soluciones cerradas, concordancia entre motores sobre la misma malla, tiempo de cálculo, memoria de pico y viabilidad de distribución (licencia, disponibilidad en Windows, dependencias en tiempo de ejecución).

Los cinco motores comparados son:

| Motor | Versión | Licencia | Lenguaje del núcleo | Resolvedores usados |
|---|---|---|---|---|
| App (spinpy) | 1.0.2 | MIT | Python (numpy, scipy, pyamg) | LU de SuperLU por debajo de 6000 GDL; CG + multigrid algebraico (pyamg) con modos rígidos por encima |
| NGSolve | 6.2.2607 | LGPL-2.1 | C++ con interfaz Python | PARDISO (MKL); CG + BDDC (P2); CG + pyamg sobre su matriz (hex8) |
| FEniCSx (DOLFINx) | 0.11.0 | LGPL-3.0 | C++ con UFL y compilación JIT | Cholesky de MUMPS; CG + GAMG (PETSc) con modos rígidos |
| scikit-fem | 12.0.2 | BSD-3 | Python (numpy, scipy) | SuperLU; CG + pyamg con modos rígidos |
| SfePy | 2026.3 | BSD-3 | Python con extensiones C/Cython | SuperLU (`ls.scipy_direct`); CG (`ls.scipy_iterative`) + pyamg con modos rígidos |

## 2. Diseño de la comparación

### 2.1 Un mismo problema discreto para todos

La comparación aísla el motor. Cada caso es un problema discreto completo que construye spinpy (malla, condiciones de contorno, material y carga) y que se entrega idéntico a los cinco motores (`spinpy.motores.problema_de_malla`). Cada motor devuelve el desplazamiento nodal en el orden de nodos de spinpy; el emparejamiento se hace por coordenadas y se aborta si algún nodo no tiene pareja a 10⁻⁹ del tamaño de la malla. El postproceso (tensión por elemento, módulo aparente, capa superficial, criterio de Pistoia) es único y se aplica al desplazamiento de cada motor, de modo que cualquier diferencia entre motores procede de su montaje, su cuadratura, sus condiciones de contorno o su resolvedor, nunca de la malla ni de la definición de las magnitudes.

Los motores que construyen su propio espacio P2 sobre las esquinas del tetraedro (todos salvo la app, que no resuelve TET10) necesitan aristas rectas. En las mallas TET10 de spinpy, el ajuste de los nodos a los planos del cubo (`fem.preparar_tet10`) mueve esquinas situadas a menos de 10⁻³ h de un plano sin mover el nodo intermedio de sus aristas; se midió un desplazamiento de hasta 8,4·10⁻⁴ h respecto del punto medio. Antes de resolver, cada nodo intermedio se recoloca en el punto medio de su arista (ecuación 1), y el desplazamiento máximo aplicado se registra con cada problema.

```math
%eq 1
\mathbf{x}_{ij}=\frac{1}{2}\left(\mathbf{x}_{i}+\mathbf{x}_{j}\right),\quad (i,j)\in\{(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)\}
```

### 2.2 Casos

1. **Bloque macizo** de 8 × 8 × 12 vóxeles (h = 0,05 mm), con hexaedros y con TET10 (tetgen sobre el mismo cubo). El ensayo lineal tiene solución exacta del problema discreto: deformación uniforme y E_app = E_s.
2. **Bloque en compresión uniaxial no lineal** (4 × 4 × 6 vóxeles) con plato rígido sin fricción hasta el 20 % de acortamiento, con St. Venant-Kirchhoff y con el neo-Hookeano compresible de FEBio. La deformación es homogénea y la fuerza tiene expresión cerrada (sección 3.6).
3. **Cavidad esférica**: cubo de 1 mm con una cavidad centrada de radio 0,2 mm. Con hexaedros a 32³ (vóxeles) y con TET10 sobre la esfera analítica (la malla de referencia de `comparativa_febio_tet/cavidad.py`, la misma función y los mismos argumentos), para contrastar con el resultado de FEBio guardado en el repositorio.
4. **Espinodoide** trabecular (ρ = 0,30, número de onda 12π, 700 ondas, conos de 30°, 30° y 90°, semilla 1, lado de 5 mm) con hexaedros a 24³, 32³, 48³, 64³ y 80³. Con la misma semilla el campo aleatorio es el mismo a cualquier resolución; al subir n cambia la discretización, no la estructura.
5. **Espinodoide con malla suave TET10** a 32³ y 48³ (`fem.mallar`, opciones por omisión: 20 iteraciones de Taubin, decimado 0,5, corrección de volumen).
6. **Espinodoide no lineal** a 24³ con hexaedros, plato rígido al 0,25 %, 0,5 %, 1 % y 2 %, con los dos materiales.

Material: E_s = 20 GPa, ν = 0,30. Ensayo de la app: tensión aparente de 1 MPa sobre la sección bruta, apoyo deslizante.

### 2.3 Medida

Cada resolución corre en un proceso hijo (`comparativa_motores/correr.py`), con 4 hilos permitidos (`OMP_NUM_THREADS = MKL_NUM_THREADS = 4`) en un equipo de 4 núcleos y 15 GB, sin otra carga durante la campaña. El tiempo es de reloj de pared por etapas (malla del motor, compilación JIT, montaje, resolución, extracción) y excluye la importación del motor. En FEniCSx, la compilación JIT de las formas se informa aparte y solo ocurre la primera vez que se usa cada forma (después queda en caché en disco). La memoria es el pico de memoria residente del proceso hijo leído de `VmHWM` (no de `ru_maxrss`, que en Linux se hereda a través de `exec` y daría el del proceso lanzador). Un vigilante detiene el hijo por encima de 13 GB o de una hora.

## 3. Formulación matemática

### 3.1 Elasticidad lineal

Con desplazamiento **u**, deformación infinitesimal y ley de Hooke isótropa (Lamé λ, μ):

```math
%eq 2
\varepsilon(\mathbf{u})=\frac{1}{2}\left(\nabla\mathbf{u}+\nabla\mathbf{u}^{T}\right),\qquad \sigma=\lambda\,\mathrm{tr}(\varepsilon)\,\mathbf{I}+2\mu\,\varepsilon,\qquad \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\quad \mu=\frac{E}{2(1+\nu)}
```

En notación de Voigt (orden xx, yy, zz, yz, xz, xy, con distorsiones de ingeniería γ = 2ε), σ = **D** ε con

```math
%eq 3
\mathbf{D}=\left[\begin{array}{cccccc}\lambda+2\mu&\lambda&\lambda&0&0&0\\ \lambda&\lambda+2\mu&\lambda&0&0&0\\ \lambda&\lambda&\lambda+2\mu&0&0&0\\ 0&0&0&\mu&0&0\\ 0&0&0&0&\mu&0\\ 0&0&0&0&0&\mu\end{array}\right]
```

La forma débil del ensayo de compresión es: hallar **u** que cumpla las condiciones esenciales tal que, para todo **v** admisible,

```math
%eq 4
\int_{\Omega}\sigma(\mathbf{u}):\varepsilon(\mathbf{v})\,d\Omega=\int_{\Gamma_{t}}\mathbf{t}\cdot\mathbf{v}\,d\Gamma,\qquad \mathbf{t}=-p\,\mathbf{e}_{z},\qquad p=\frac{\sigma_{0}A_{b}}{A_{t}}
```

donde σ₀ es la tensión aparente sobre la sección bruta A_b y A_t es el área ósea de las caras del techo, de modo que la fuerza total es σ₀ A_b, la convención de la app.

**Condiciones de contorno.** Apoyo deslizante: u_z = 0 en todos los nodos de la base (intermedios incluidos con TET10); u_x = u_y = 0 en el nodo de la base que minimiza x + y y u_y = 0 en el que maximiza x − y, lo mínimo para eliminar los movimientos de sólido rígido. Con TET10 las dos anclas se eligen entre las esquinas. Apoyo empotrado: los tres componentes en toda la base. Plato rígido sin fricción: además, u_z = −ε_p H en todos los nodos del techo.

### 3.2 Elementos

**Hexaedro trilineal (hex8).** Sobre el cubo de referencia ξ, η, ζ ∈ [−1, 1], con los vértices en el orden C3D8 (ξ_a, η_a, ζ_a ∈ {−1, 1}):

```math
%eq 5
N_{a}(\xi,\eta,\zeta)=\frac{1}{8}(1+\xi_{a}\xi)(1+\eta_{a}\eta)(1+\zeta_{a}\zeta),\qquad a=1,\dots,8
```

Integración de Gauss 2 × 2 × 2 (puntos ±1/√3, pesos 1). En un vóxel el jacobiano es constante, J = diag(h_x, h_y, h_z)/2.

**Tetraedro cuadrático (TET10).** Con las coordenadas baricéntricas L_i (i = 0, …, 3) y el orden C3D10 (aristas 01, 12, 02, 03, 13, 23):

```math
%eq 6
N_{i}=L_{i}(2L_{i}-1),\qquad N_{ij}=4L_{i}L_{j}
```

Con aristas rectas el mapeo es afín, la deformación es lineal en el elemento y la matriz de rigidez (integrando de grado 2) se integra exactamente con cualquier regla de grado ≥ 2. Cada motor usa la suya: scikit-fem, FEniCSx y SfePy, grado 2; NGSolve, grado 4. Todas son exactas para este integrando, así que el sistema lineal es el mismo. Con hexaedros de un vóxel el mapeo también es afín y el integrando es de grado 2 en cada dirección; todos los motores usan Gauss con dos puntos por dirección, exacto hasta grado 3.

**Matriz de rigidez y carga.** Con **B** la matriz de derivadas de forma en Voigt,

```math
%eq 7
\mathbf{K}=\sum_{e}\int_{\Omega_{e}}\mathbf{B}^{T}\mathbf{D}\,\mathbf{B}\,d\Omega,\qquad \mathbf{f}_{a}=\int_{\Gamma_{t}}N_{a}\,\mathbf{t}\,d\Gamma
```

En una cara cuadrada del techo cada uno de sus cuatro nodos recibe p A/4 (el reparto por área tributaria de la app); en una cara tri6 las esquinas reciben cero y cada nodo intermedio p A/3.

### 3.3 Sistema reducido y resolvedores

Con los grados de libertad particionados en libres (f) y prescritos (d), el sistema que se resuelve es

```math
%eq 8
\mathbf{K}_{ff}\,\mathbf{u}_{f}=\mathbf{f}_{f}-\mathbf{K}_{fd}\,\mathbf{u}_{d}
```

Los resolvedores directos factorizan K_ff (LU de SuperLU con ordenación COLAMD; Cholesky de MUMPS; LDLᵀ de PARDISO con disección anidada). Los iterativos usan gradiente conjugado precondicionado por multigrid: agregación suavizada (pyamg, GAMG de PETSc) con los seis modos de sólido rígido como espacio casi nulo, o BDDC (NGSolve, con TET10). Los modos rígidos son las tres traslaciones y los tres giros infinitesimales:

```math
%eq 9
\mathbf{r}_{1..3}=\mathbf{e}_{x},\,\mathbf{e}_{y},\,\mathbf{e}_{z};\qquad \mathbf{r}_{4}=(-y,\,x,\,0),\quad \mathbf{r}_{5}=(0,\,-z,\,y),\quad \mathbf{r}_{6}=(z,\,0,\,-x)
```

El criterio de parada es el residuo relativo no precondicionado (ecuación 10), 10⁻¹⁰ en los motores externos y 10⁻⁸ en la app (su valor de siempre). En PETSc se fijó explícitamente la norma no precondicionada, que no es la que usa por omisión.

```math
%eq 10
r=\frac{\|\mathbf{f}_{f}-\mathbf{K}_{ff}\mathbf{u}_{f}\|_{2}}{\|\mathbf{f}_{f}\|_{2}}
```

### 3.4 Postproceso común

La deformación de cada elemento se evalúa en su centro (hex8, la convención de la app, igual a la media de volumen en un vóxel) o en su centroide (TET10, igual a la media de volumen porque la deformación es lineal); la tensión es σ_e = **D** ε_e. De ella salen la tensión de von Mises y la deformación efectiva de Pistoia:

```math
%eq 11
\sigma_{vM}=\sqrt{\frac{1}{2}\left[(\sigma_{xx}-\sigma_{yy})^{2}+(\sigma_{yy}-\sigma_{zz})^{2}+(\sigma_{zz}-\sigma_{xx})^{2}\right]+3\left(\tau_{yz}^{2}+\tau_{xz}^{2}+\tau_{xy}^{2}\right)}
```

```math
%eq 12
U_{e}=\frac{1}{2}\,\sigma_{e}\cdot\varepsilon_{e},\qquad \varepsilon_{eff,e}=\sqrt{\frac{2\,U_{e}}{E_{s}}}
```

El módulo aparente usa el desplazamiento vertical medio del techo: la media simple de sus nodos (la definición de la app con vóxeles) o la media ponderada por área de sus caras cargadas (la que se usa con TET10):

```math
%eq 13
E_{app}=\frac{\sigma_{0}}{|\bar{u}_{z}|/H},\qquad \bar{u}_{z}^{nodal}=\frac{1}{n_{t}}\sum_{k\in\mathrm{techo}}u_{z,k},\qquad \bar{u}_{z}^{\acute{a}rea}=\frac{\sum_{c}A_{c}\,\bar{u}_{z,c}}{\sum_{c}A_{c}}
```

El equilibrio se comprueba, con independencia del motor, por la integral de volumen de la tensión (exacta en el problema discreto lineal, con el desplazamiento virtual v = z e_z):

```math
%eq 14
F=-\frac{1}{H}\sum_{e}\sigma_{zz,e}\,V_{e},\qquad dF=\frac{|F-\sigma_{0}A_{b}|}{\sigma_{0}A_{b}}
```

El criterio de Pistoia et al. (2002) busca el factor k que lleva el percentil 98 (en volumen) de ε_eff a la deformación crítica, y el p99 de von Mises se toma en la capa superficial (elementos con una cara hacia el poro, sin contar las caras del cubo), ambos ponderados por volumen con TET10:

```math
%eq 15
k=\frac{\varepsilon_{crit}}{P_{98}(\varepsilon_{eff})},\qquad \sigma_{fallo}=k\,\sigma_{0},\qquad \varepsilon_{crit}=0.007
```

### 3.5 Análisis no lineal (lagrangiano total)

Con el gradiente de deformación **F**, el tensor de Cauchy-Green derecho **C** y el de Green-Lagrange **E**:

```math
%eq 16
\mathbf{F}=\mathbf{I}+\nabla_{0}\mathbf{u},\qquad \mathbf{C}=\mathbf{F}^{T}\mathbf{F},\qquad \mathbf{E}=\frac{1}{2}(\mathbf{C}-\mathbf{I}),\qquad J=\det\mathbf{F}
```

Los dos materiales hiperelásticos comparados son St. Venant-Kirchhoff (el `isotropic elastic` de FEBio) y el neo-Hookeano compresible de FEBio:

```math
%eq 17
W_{SVK}=\frac{\lambda}{2}(\mathrm{tr}\,\mathbf{E})^{2}+\mu\,\mathbf{E}:\mathbf{E},\qquad W_{NH}=\frac{\mu}{2}(\mathrm{tr}\,\mathbf{C}-3)-\mu\ln J+\frac{\lambda}{2}(\ln J)^{2}
```

Sus segundos tensores de Piola-Kirchhoff **S** = ∂W/∂**E** son

```math
%eq 18
\mathbf{S}_{SVK}=\lambda\,\mathrm{tr}(\mathbf{E})\,\mathbf{I}+2\mu\,\mathbf{E},\qquad \mathbf{S}_{NH}=\mu(\mathbf{I}-\mathbf{C}^{-1})+\lambda\ln J\,\mathbf{C}^{-1}
```

El residuo (fuerzas internas menos externas) con el primer tensor de Piola-Kirchhoff **P** = **F S** es

```math
%eq 19
\mathbf{R}(\mathbf{u})\cdot\mathbf{v}=\int_{\Omega_{0}}\mathbf{F}\mathbf{S}:\nabla_{0}\mathbf{v}\,d\Omega_{0}-\int_{\Gamma_{t}}\mathbf{t}\cdot\mathbf{v}\,d\Gamma
```

y su linealización consistente en la dirección δ**u**, con δ**E** = sym(**F**ᵀ∇₀δ**u**), es

```math
%eq 20
D\mathbf{R}[\delta\mathbf{u}]\cdot\mathbf{v}=\int_{\Omega_{0}}\left(\nabla_{0}\delta\mathbf{u}\,\mathbf{S}+\mathbf{F}\,\delta\mathbf{S}\right):\nabla_{0}\mathbf{v}\,d\Omega_{0}
```

```math
%eq 21
\delta\mathbf{S}_{SVK}=\lambda\,\mathrm{tr}(\delta\mathbf{E})\,\mathbf{I}+2\mu\,\delta\mathbf{E},\qquad \delta\mathbf{S}_{NH}=\lambda(\mathbf{C}^{-1}:\delta\mathbf{E})\,\mathbf{C}^{-1}+2(\mu-\lambda\ln J)\,\mathbf{C}^{-1}\delta\mathbf{E}\,\mathbf{C}^{-1}
```

El primer término de la ecuación 20 es la rigidez geométrica. scikit-fem recibe las ecuaciones 18 a 21 escritas a mano; NGSolve y FEniCSx derivan la energía W automáticamente (`Variation` y `derivative` de UFL); SfePy usa su término `dw_tl_he_svk`. Cada paso de carga se resuelve por Newton-Raphson,

```math
%eq 22
\mathbf{K}_{T}(\mathbf{u}^{k})\,\Delta\mathbf{u}=-\mathbf{R}(\mathbf{u}^{k}),\qquad \mathbf{u}^{k+1}=\mathbf{u}^{k}+\Delta\mathbf{u},\qquad \frac{\|\mathbf{R}_{f}\|}{\|\mathbf{f}_{int,techo}\|}\leq10^{-10}
```

**Predictor consistente del plato.** Con control por desplazamiento, imponer el incremento Δu_d solo en los nodos del techo deja una primera iteración muy distorsionada. En el bloque al 20 % con St. Venant-Kirchhoff, cuya energía no es convexa, el Newton que se escribió así para NGSolve convergió en 22 iteraciones a otra rama de equilibrio (fuerza de −133,4 N en lugar de 115,2 N). La primera iteración de cada paso se hace, por eso, con el incremento del plato dentro del sistema linealizado, el mismo esquema que aplica SNES de PETSc con las condiciones de contorno por elevación:

```math
%eq 23
\mathbf{K}_{T,ff}\,\Delta\mathbf{u}_{f}=-\mathbf{R}_{f}-\mathbf{K}_{T,fd}\,\Delta\mathbf{u}_{d}
```

**Fuerza de reacción.** Es la suma de las fuerzas internas en z sobre los nodos del techo, positiva en compresión:

```math
%eq 24
F_{reac}=-\sum_{k\in\mathrm{techo}}f_{int,z,k}
```

**Carga muerta frente a presión seguidora.** Con control por fuerza, los motores internos aplican una tracción muerta (fija en la configuración de referencia); la presión de FEBio sigue a la cara deformada. La diferencia es del orden de la rotación de las caras del techo y se declara en el informe de la aplicación.

### 3.6 Soluciones cerradas del bloque uniaxial

Con plato sin fricción, base deslizante y laterales libres, la deformación es homogénea, **F** = diag(a, a, λ_z), y S_xx = S_yy = 0 fija el estiramiento lateral a. Para St. Venant-Kirchhoff la solución es explícita:

```math
%eq 25
E_{xx}=-\nu\,E_{zz},\qquad E_{zz}=\frac{\lambda_{z}^{2}-1}{2},\qquad P_{zz}=\lambda_{z}\,E\,\frac{\lambda_{z}^{2}-1}{2},\qquad F=-P_{zz}A_{0}
```

Para el neo-Hookeano de FEBio, a es la raíz de la ecuación 26 (resuelta por el método de Brent a 10⁻¹⁵) y la tensión nominal es la ecuación 27:

```math
%eq 26
\mu(a^{2}-1)+\lambda\ln(a^{2}\lambda_{z})=0
```

```math
%eq 27
P_{zz}=\lambda_{z}\left[\mu\left(1-\lambda_{z}^{-2}\right)+\lambda\,\ln J\,\lambda_{z}^{-2}\right],\qquad J=a^{2}\lambda_{z}
```

### 3.7 Métricas de comparación

Entre dos soluciones del mismo problema discreto (un motor frente a la referencia del caso, que es una solución por resolvedor directo; la de NGSolve si existe, si no la de FEniCSx, SfePy o scikit-fem):

```math
%eq 28
\Delta u=\frac{\max_{k}|\mathbf{u}_{k}-\mathbf{u}_{k}^{ref}|}{\max_{k}|\mathbf{u}_{k}^{ref}|},\qquad \Delta\sigma_{vM}=\frac{\max_{e}|\sigma_{vM,e}-\sigma_{vM,e}^{ref}|}{\max_{e}|\sigma_{vM,e}^{ref}|},\qquad \Delta E=\frac{E_{app}}{E_{app}^{ref}}-1
```

Para comprobar la tangente de un motor se comparó su matriz tangente con la derivada numérica centrada de su propio residuo, en un estado aleatorio **x** y una dirección aleatoria **d**:

```math
%eq 29
e_{T}=\frac{\left\|\mathbf{K}_{T}(\mathbf{x})\,\mathbf{d}-\frac{\mathbf{R}(\mathbf{x}+h\mathbf{d})-\mathbf{R}(\mathbf{x}-h\mathbf{d})}{2h}\right\|}{\left\|\frac{\mathbf{R}(\mathbf{x}+h\mathbf{d})-\mathbf{R}(\mathbf{x}-h\mathbf{d})}{2h}\right\|},\qquad h=10^{-7}H
```

Una tangente consistente da e_T del orden del error de truncamiento de la diferencia centrada (≲ 10⁻⁸); un valor de orden uno indica una tangente que no es la derivada del residuo, y con ella Newton pierde la convergencia cuadrática.
