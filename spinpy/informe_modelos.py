"""
informe_modelos.py — Seccion «Modelos matematicos y metodos numericos» del
informe para publicacion.

POR QUE UNA SECCION APARTE
--------------------------
El parrafo de metodos (`informe.parrafos_metodos`) es el que se pega en un
articulo: dice QUE se hizo con los valores usados. Un revisor pide ademas COMO:
que ecuacion da la tension de von Mises, sobre que conjunto de elementos se
toma el percentil, que condiciones de contorno tiene el ensayo, que convencion
de percentil se usa. Esta seccion responde a eso con las ecuaciones que
EJECUTA el codigo, no con las de un libro.

REGLAS PARA QUE NO DIVERJA DEL CODIGO
-------------------------------------
  * Toda constante numerica que aparece en el texto se LEE del modulo que la
    usa (umbrales, pesos, tolerancias, numero de direcciones del MIL...). Si
    alguien cambia un umbral, el informe cambia con el.
  * Las ecuaciones estan comprobadas contra el codigo en el bloque 22 de la
    suite (`tests/test_22_modelos_informe.py`): von Mises, deformacion
    efectiva, E_{app}, Pistoia, percentiles, capa superficial, identidades de
    Parfitt, DA, umbral del GRF, constantes elasticas y E(n).
  * Cada subseccion aparece solo si el documento contiene ese calculo: no se
    describe un ensayo que no se hizo.

Las ecuaciones van en bloques ```math (LaTeX del subconjunto de mathtext de
matplotlib) con una linea `%eq N`; `md2pdf` las compone numeradas y en GitHub
se ven como LaTeX.
"""

from __future__ import annotations

import inspect

import numpy as np

from . import avisos, error as err_mod, incertidumbre
from .elastic import UMBRAL_DIRECTO, VOID_SCALE, homogeneizar
from .espesor import espesor_local  # noqa: F401  (referencia del texto)
from .morphometry import indice_smi, tensor_mil
from .resistencia import (EPS_CRITICA, FRAC_CRITICA, SIGMA0_DEF,
                          ensayo_compresion)


# Numeros de onda (beta/pi) de las etapas A del ajuste del spinodoide. Son las
# expresiones literales de `fit.ajustar_spinodoide`, que no las expone como
# constantes; el bloque 22 comprueba que el ajuste recorre exactamente estas.
ONDAS_RAPIDO = np.unique(np.clip([10.0, 15.0, 20.0], 8, 25))
ONDAS_COMPLETO = np.unique(np.round(np.linspace(8, 25, 6)))


def _defecto(fn, nombre):
    return inspect.signature(fn).parameters[nombre].default


class _Eq:
    """Numeracion de ecuaciones y bloques ```math."""

    def __init__(self):
        self.n = 0

    def __call__(self, L, tex):
        self.n += 1
        L += ["```math", f"%eq {self.n}", tex, "```", ""]
        return self.n


def seccion_modelos(doc, idioma, r, numero):
    """Lineas de Markdown de la seccion, citando con el redactor `r`.

    `numero` es el numero de seccion. Devuelve [] si el documento no tiene
    ningun calculo que describir.
    """
    from . import informe as I          # import tardio: informe importa esto

    es = idioma == "es"

    def T(a, b):
        return a if es else b

    def n(x, espec=".4g"):
        return r.n(x, espec)

    res = doc.get("resultados") or {}
    m_voi = doc.get("morfometria_voi") or {}
    fams = [(f, suf) for f, suf in (("spinodoide", ""), ("dual-lattice",
                                                         "_dual"))
            if isinstance(res.get("ajuste" + suf), dict)]
    ela = [res.get("elastico" + s) for _f, s in (("", ""), ("", "_dual"))]
    ela = [x for x in ela if isinstance(x, dict)]
    fe = [res.get("resistencia" + s) for s in ("", "_dual")]
    fe = [x for x in fe if isinstance(x, dict)]
    comp = [res.get("analisis_comparado" + s) for s in ("", "_dual")]
    comp = [x for x in comp if isinstance(x, dict)]
    if not (m_voi or fams or ela or fe or comp):
        return []

    eq = _Eq()
    L = []
    k = [0]

    def sub(es_t, en_t):
        k[0] += 1
        L.append(f"### {numero}.{k[0]} " + T(es_t, en_t))
        L.append("")

    def par(es_t, en_t):
        L.append(T(es_t, en_t))
        L.append("")

    from . import __version__
    L.append(f"## {numero}. " + T("Modelos matemáticos y métodos numéricos",
                                  "Mathematical models and numerical methods"))
    L.append("")
    par(f"Esta sección describe cada cálculo tal como lo ejecuta spinpy "
        f"{__version__}: las ecuaciones son las del código, no las de una "
        "referencia genérica, y las constantes numéricas del texto se leen "
        "del propio código al escribir el informe. Solo se describen los "
        "cálculos que contiene esta sesión.",
        f"This section describes every computation exactly as spinpy "
        f"{__version__} executes it: the equations are those of the code, "
        "not of a generic reference, and the numerical constants in the text "
        "are read from the code itself when the report is written. Only the "
        "computations contained in this session are described.")

    # ------------------------------------------------------------------
    sub("Convenciones y unidades", "Conventions and units")
    par("Las longitudes están en mm, las tensiones y módulos en Pa (las "
        "tablas los dan en MPa) y las fuerzas en N. La imagen es una máscara "
        "binaria B(i, j, k) ∈ {0, 1} (1 = hueso), con el primer índice a lo "
        "largo de X, sobre un cubo de n vóxeles de arista h y lado L = n·h. "
        "Deformaciones y tensiones se escriben en notación de Voigt con "
        "deformaciones angulares de ingeniería:",
        "Lengths are in mm, stresses and moduli in Pa (tables give MPa) and "
        "forces in N. The image is a binary mask B(i, j, k) ∈ {0, 1} (1 = "
        "bone), with the first index along X, on a cube of n voxels of edge h "
        "and side L = n·h. Strains and stresses are written in Voigt notation "
        "with engineering shear strains:")
    eq(L, r"\mathbf{\varepsilon}=[\varepsilon_{xx},\,"
                    r"\varepsilon_{yy},\,\varepsilon_{zz},\,\gamma_{yz},\,"
                    r"\gamma_{xz},\,\gamma_{xy}]^{T},\quad \gamma_{ij}="
                    r"2\,\varepsilon_{ij},\quad \mathbf{\sigma}="
                    r"[\sigma_{xx},\,\sigma_{yy},\,\sigma_{zz},\,\tau_{yz},\,"
                    r"\tau_{xz},\,\tau_{xy}]^{T}")
    par("El tejido es elástico lineal e isótropo, con σ = D ε y",
        "The tissue is linear elastic and isotropic, with σ = D ε and")
    eq(L, r"\sigma_{ii}=\lambda\,(\varepsilon_{xx}+\varepsilon_{yy}+"
          r"\varepsilon_{zz})+2\mu\,\varepsilon_{ii},\quad \tau_{ij}=\mu\,"
          r"\gamma_{ij},\quad\lambda=\frac{E_s\,\nu_s}{(1+\nu_s)(1-2\nu_s)},"
          r"\quad\mu=\frac{E_s}{2(1+\nu_s)}")
    par("Cuando la carga se da como fuerza, la tensión aparente se obtiene "
        "sobre la sección bruta del cubo con el factor de unidades de un área "
        "en mm²:",
        "When the load is given as a force, the apparent stress is taken over "
        "the gross cross-section of the cube with the unit factor of an area "
        "in mm²:")
    eq(L, r"\sigma_0\,[\mathrm{Pa}]=\frac{10^{6}\,F\,[\mathrm{N}]}{A\,"
          r"[\mathrm{mm}^{2}]},\qquad A=L_x\,L_y")
    par("Se usan dos definiciones de percentil, y cada resultado declara "
        "cuál. Sobre la muestra ordenada x(1) ≤ … ≤ x(n), con p ∈ [0, 1], "
        f"la de tipo 7 de Hyndman y Fan {r.c('hyndman1996')} (interpolación "
        "lineal, la de numpy por omisión) y la de tipo 5 (posiciones "
        "(i − ½)/n, la de `prctile` de MATLAB) son",
        "Two percentile definitions are used, and every result states which. "
        "On the sorted sample x(1) ≤ … ≤ x(n), with p ∈ [0, 1], Hyndman and "
        f"Fan's type 7 {r.c('hyndman1996')} (linear interpolation, numpy's "
        "default) and type 5 (positions (i − ½)/n, MATLAB's `prctile`) are")
    e_q7 = eq(L, r"Q_7(p):\ \ h=(n-1)\,p+1,\qquad\qquad Q_5(p):\ \ h="
                 r"n\,p+\frac{1}{2}")
    eq(L, r"Q(p)=x_{(\lfloor h\rfloor)}+(h-\lfloor h\rfloor)\,\left("
          r"x_{(\lfloor h\rfloor+1)}-x_{(\lfloor h\rfloor)}\right)")

    # ------------------------------------------------------------------
    if m_voi:
        sub("Morfometría", "Morphometry")
        modo = (res.get("morfometria") or {}).get("modo") or m_voi.get("modo")
        par("BV/TV es la fracción de vóxeles de hueso. La superficie ósea BS "
            "es el área de la isosuperficie de nivel ½ de la máscara, "
            f"extraída con marching cubes {r.c('lorensen1987', 'lewiner2003')}"
            " sobre los centros de vóxel y sumada triángulo a triángulo. La "
            "isosuperficie se deja ABIERTA en las seis caras del cubo: las "
            "tapas de corte no son interfaz ósea y no cuentan.",
            "BV/TV is the fraction of bone voxels. The bone surface BS is the "
            "area of the level-½ isosurface of the mask, extracted with "
            f"marching cubes {r.c('lorensen1987', 'lewiner2003')} on the "
            "voxel centres and summed triangle by triangle. The isosurface is "
            "left OPEN at the six faces of the cube: the cut caps are not a "
            "bone interface and are not counted.")
        eq(L, r"\frac{BV}{TV}=\frac{1}{n_x n_y n_z}\sum_{ijk}B_{ijk},\qquad "
              r"BS=\sum_{t}\frac{1}{2}\left\|(\mathbf{b}_t-\mathbf{a}_t)"
              r"\times(\mathbf{c}_t-\mathbf{a}_t)\right\|")
        par("Espesor, separación y número trabeculares siguen el modelo de "
            f"placas de Parfitt {r.c('parfitt1987')}, como CTAn y Scanco "
            f"{r.c('bouxsein2010')}:",
            "Trabecular thickness, separation and number follow Parfitt's "
            f"plate model {r.c('parfitt1987')}, as CTAn and Scanco do "
            f"{r.c('bouxsein2010')}:")
        eq(L, r"Tb.Th=\frac{2\,BV}{BS},\qquad Tb.Sp=Tb.Th\left(\frac{TV}"
              r"{BV}-1\right),\qquad Tb.N=\frac{BV/TV}{Tb.Th}")
        n_dirs = _defecto(tensor_mil, "n_dirs")
        n_lin = _defecto(tensor_mil, "n_lines")
        cuerda = _defecto(tensor_mil, "min_chord_frac")
        muestreo = (doc.get("muestreo_mil")
                    or (res.get("morfometria") or {}).get("muestreo_mil")
                    or "voxel")
        par(f"La anisotropía sale del tensor de longitud media de intercepción "
            f"(MIL) {r.c('harrigan1984')}. Se toman {n_dirs} direcciones "
            "cuasiuniformes sobre el hemisferio (espiral de Fibonacci: "
            f"cos θₖ = (k + ½)/{n_dirs}, φₖ = π(1 + √5)(k + ½)). Por cada "
            f"dirección d se lanza una rejilla de {n_lin}×{n_lin} rectas "
            "paralelas por el centro del cubo, muestreadas cada h/2; se "
            f"descartan las cuerdas de longitud útil menor que {n(cuerda)}·L. "
            "Con la longitud total de hueso atravesada, ℓ(d), y el número de "
            "transiciones hueso↔poro, N(d), el MIL y el tensor de fábrica M "
            "(ajustado por mínimos cuadrados) son",
            "Anisotropy comes from the mean intercept length (MIL) tensor "
            f"{r.c('harrigan1984')}. {n_dirs} quasi-uniform directions over "
            "the hemisphere are used (Fibonacci spiral: cos θₖ = (k + ½)/"
            f"{n_dirs}, φₖ = π(1 + √5)(k + ½)). For each direction d a grid "
            f"of {n_lin}×{n_lin} parallel lines through the centre of the cube "
            "is cast, sampled every h/2; chords with a useful length below "
            f"{n(cuerda)}·L are discarded. With the total bone length crossed, "
            "ℓ(d), and the number of bone↔pore transitions, N(d), the MIL and "
            "the fabric tensor M (least-squares fit) are")
        eq(L, r"\mathrm{MIL}(\mathbf{d})=\frac{\ell(\mathbf{d})}"
                      r"{\frac{1}{2}N(\mathbf{d})},\qquad \frac{1}"
                      r"{\mathrm{MIL}(\mathbf{d})^{2}}=\mathbf{d}^{T}\mathbf{M}"
                      r"\,\mathbf{d},\qquad DA=\sqrt{\frac{\lambda_3}"
                      r"{\lambda_1}},\qquad DA_2=\sqrt{\frac{\lambda_2}"
                      r"{\lambda_1}}")
        par("con λ₁ ≤ λ₂ ≤ λ₃ los autovalores de M. El autovector de λ₁ (MIL "
            "máximo) es la dirección principal. Si el ajuste da algún "
            "autovalor no positivo (estructuras laminares) se acota a la "
            "mitad del menor 1/MIL² medido y el DA se declara cota inferior. "
            f"Muestreo de las rectas en esta sesión: «{muestreo}»"
            + (", con la máscara suavizada por una gaussiana de 0,7 vóxeles e "
               "interpolación trilineal con umbral ½ (corrección M1)"
               if muestreo == "suavizado" else
               ", leyendo el vóxel que contiene cada muestra") + ". "
            f"Diferencias de DA menores que ~0,07 están por debajo del suelo "
            f"de ruido medido del estimador, y con DA < "
            f"{n(avisos.DA_MIN_EJE, '.3g')} la dirección principal no se "
            "interpreta.",
            "with λ₁ ≤ λ₂ ≤ λ₃ the eigenvalues of M. The eigenvector of λ₁ "
            "(largest MIL) is the principal direction. If the fit yields a "
            "non-positive eigenvalue (laminar structures) it is clamped to "
            "half the smallest measured 1/MIL² and DA is declared a lower "
            f"bound. Line sampling in this session: \"{muestreo}\""
            + (", with the mask smoothed by a 0.7-voxel Gaussian and "
               "trilinear interpolation with threshold ½ (correction M1)"
               if muestreo == "suavizado" else
               ", reading the voxel that contains each sample") + ". "
            "DA differences below ~0.07 are under the measured noise floor "
            f"of the estimator, and with DA < {n(avisos.DA_MIN_EJE, '.3g')} "
            "the principal direction is not interpreted.")
        if np.isfinite(float(m_voi.get("ConnD", np.nan) or np.nan)):
            dr = _defecto(indice_smi, "dr_rel")
            par("La densidad de conectividad usa la característica de Euler χ "
                f"de la fase sólida con vecindad 26 {r.c('odgaard1993')}, y el "
                f"índice de modelo estructural {r.c('hildebrand1997smi')} "
                "desplaza cada vértice de la isosuperficie a lo largo de su "
                f"normal una distancia δ = {n(dr, '.2g')}·h:",
                "Connectivity density uses the Euler characteristic χ of the "
                f"solid phase with 26-neighbourhood {r.c('odgaard1993')}, and "
                f"the structure model index {r.c('hildebrand1997smi')} "
                "displaces every isosurface vertex along its normal by a "
                f"distance δ = {n(dr, '.2g')}·h:")
            eq(L, r"Conn.D=\frac{1-\chi}{TV},\qquad SMI=\frac{6\,BV}{BS^{2}}"
                  r"\,\frac{BS(\delta)-BS(0)}{\delta}")
            par(f"El SMI está confundido por la concavidad {r.c('salmon2015')}"
                ": con BV/TV alto puede salir negativo y perder la lectura "
                "placa/barra.",
                "The SMI is confounded by concavity "
                f"{r.c('salmon2015')}: at high BV/TV it can be negative and "
                "lose its plate/rod reading.")
        if True:
            par("El espesor local (y el tamaño de poro, aplicado a la fase "
                "complementaria) es el diámetro de la mayor esfera inscrita "
                f"que contiene cada punto {r.c('hildebrand1997thickness')}. "
                "Se calcula con transformadas de distancia euclídea: con r(x) "
                "la distancia al fondo, se recorren los radios desde el máximo "
                "hacia abajo en pasos de h y cada vóxel toma el primer radio "
                "cuya bola lo cubre. Po.Dm es la media de ese espesor sobre "
                "los vóxeles de poro. Es también la magnitud de la figura de "
                "distribuciones.",
                "Local thickness (and pore size, applied to the complementary "
                "phase) is the diameter of the largest inscribed sphere that "
                f"contains each point {r.c('hildebrand1997thickness')}. It is "
                "computed with Euclidean distance transforms: with r(x) the "
                "distance to the background, radii are swept from the maximum "
                "downwards in steps of h and every voxel takes the first "
                "radius whose ball covers it. Po.Dm is the mean of that "
                "thickness over pore voxels. It is also the quantity of the "
                "distributions figure.")
            eq(L, r"\tau(\mathbf{x})=2\,\max\left\lbrace r:\ \exists\,"
                  r"\mathbf{c},\ \mathbf{x}\in S(\mathbf{c},r)\subseteq\Omega"
                  r"\right\rbrace,\qquad Po.Dm=\frac{1}{|\Omega^{c}|}\sum_{"
                  r"\mathbf{x}\in\Omega^{c}}\tau_{\Omega^{c}}(\mathbf{x})")
        if np.isfinite(float(m_voi.get("EF", np.nan) or np.nan)):
            par("El Ellipsoid Factor busca en cada punto el mayor elipsoide "
                "contenido en el hueso que lo contiene (semiejes a ≤ b ≤ c) "
                f"{r.c('doube2015')}; la fracción de placas es la del tejido "
                "con EF < −0,2 y la de barras con EF > 0,2.",
                "The Ellipsoid Factor finds, at each point, the largest "
                "ellipsoid inside the bone that contains it (semi-axes a ≤ b "
                f"≤ c) {r.c('doube2015')}; the plate fraction is the tissue "
                "with EF < −0.2 and the rod fraction that with EF > 0.2.")
            eq(L, r"EF=\frac{a}{b}-\frac{b}{c}\in[-1,\,1]")
        if modo == "malla":
            par("Esta sesión midió en modo «malla»: BV y BS salen de la "
                "isosuperficie cerrada, suavizada con Taubin, reparada y "
                "recortada medio vóxel por dentro del cubo, en lugar del "
                "conteo de vóxeles y la superficie abierta. DA, DA2 y Conn.D "
                "siguen saliendo de los vóxeles.",
                "This session measured in \"mesh\" mode: BV and BS come from "
                "the closed isosurface, Taubin-smoothed, repaired and cropped "
                "half a voxel inside the cube, instead of the voxel count and "
                "the open surface. DA, DA2 and Conn.D still come from voxels.")
        par("El VOI y cada candidato recorren exactamente el mismo código, a "
            "su propia resolución, con el candidato escalado al lado físico "
            "del VOI: el sesgo de marching cubes sobre datos binarios "
            "(~8,5 % de área medido sobre una esfera) afecta por igual a "
            "ambos.",
            "The VOI and every candidate go through exactly the same code, "
            "each at its own resolution, with the candidate scaled to the "
            "physical side of the VOI: the bias of marching cubes on binary "
            "data (~8.5 % of area measured on a sphere) affects both "
            "equally.")

    # ------------------------------------------------------------------
    nombres = [f for f, _s in fams]
    if fams:
        sub("Microestructuras sintéticas", "Synthetic microstructures")
    if "spinodoide" in nombres:
        par(f"Un espinodoide {r.c('kumar2020', 'soyarslan2018')} es un "
            "conjunto de nivel de un campo aleatorio gaussiano. El campo se "
            "evalúa en n puntos equiespaciados de [0, 1] por eje:",
            f"A spinodoid {r.c('kumar2020', 'soyarslan2018')} is a level set "
            "of a Gaussian random field. The field is evaluated at n equally "
            "spaced points of [0, 1] per axis:")
        eq(L, r"\phi(\mathbf{x})=\sqrt{\frac{2}{N}}\sum_{i=1}^{N}\cos\left("
              r"\beta\,\mathbf{n}_i\cdot\mathbf{x}+\gamma_i\right),\qquad "
              r"\gamma_i\sim\mathcal{U}[0,2\pi)")
        par("Las direcciones de onda nᵢ se restringen a una unión de conos "
            "alrededor de los ejes rotados R e_{j}, de semiángulos θ_{j}. Con el "
            "esquema «rechazo» se muestrea un candidato uniforme en la esfera "
            "y se acepta si arccos|nᵢ·R e_{j}| < θ_{j} para algún j (el reparto "
            "entre conos queda proporcional a su ángulo sólido); con "
            "«equitativo» se toman N/k ondas por cono activo, uniformes dentro "
            "de cada cono. La fase sólida es la región bajo el umbral que da "
            "la densidad relativa ρ para un campo N(0, 1):",
            "Wave directions nᵢ are restricted to a union of cones about the "
            "rotated axes R e_{j}, of half-angles θ_{j}. With the \"rechazo\" "
            "(rejection) scheme a uniform candidate on the sphere is drawn and "
            "accepted if arccos|nᵢ·R e_{j}| < θ_{j} for some j (cones receive "
            "waves in proportion to their solid angle); with \"equitativo\" "
            "N/k waves are drawn per active cone, uniformly inside each cone. "
            "The solid phase is the region below the threshold that gives the "
            "relative density ρ for an N(0, 1) field:")
        eq(L, r"\Omega=\left\lbrace\mathbf{x}:\ \phi(\mathbf{x})\leq"
                       r"\phi_0\right\rbrace,\qquad \phi_0=\sqrt{2}\,"
                       r"\mathrm{erf}^{-1}(2\rho-1)")
        par("El número de onda se guarda con sus dos lecturas, β/π y β en "
            "rad. La máscara se escala al lado del VOI con vóxel h = L/n.",
            "The wave number is stored with both readings, β/π and β in rad. "
            "The mask is scaled to the VOI side with voxel h = L/n.")
    if "dual-lattice" in nombres:
        par(f"El dual-lattice {r.c('vafaeefar2022')} se construye sobre una "
            "rejilla entera perturbada, g = z + ι·u con u ∼ U(−½, ½)³ y ι la "
            f"irregularidad, triangulada por Delaunay {r.c('virtanen2020')} "
            "ANTES de deformarla. Los vértices se llevan al cubo con el "
            "espaciado d = 1/c (c celdas por lado) y el estiramiento "
            "normalizado e (e_{x} e_{y} e_{z} = 1); el esqueleto son los "
            "semipuntales que unen el centroide de cada tetraedro con los de "
            "sus cuatro caras, y la fase sólida son los puntos a menos de r "
            "del esqueleto, con r el cuantil que da la densidad exacta:",
            f"The dual-lattice {r.c('vafaeefar2022')} is built on a perturbed "
            "integer grid, g = z + ι·u with u ∼ U(−½, ½)³ and ι the "
            f"irregularity, Delaunay-triangulated {r.c('virtanen2020')} "
            "BEFORE deforming it. Vertices are mapped to the cube with spacing "
            "d = 1/c (c cells per side) and normalised stretch e (e_{x} e_{y} e_{z} "
            "= 1); the skeleton is the set of half-struts joining each "
            "tetrahedron centroid to the centroids of its four faces, and the "
            "solid phase is the set of points closer than r to the skeleton, "
            "with r the quantile that gives the exact density:")
        eq(L, r"\mathbf{X}=\mathbf{R}\,\mathrm{diag}(d\,\mathbf{e})\,\mathbf{g}"
              r"+\frac{L}{2},\qquad \Omega=\left\lbrace\mathbf{x}:\ D(\mathbf{x}"
              r")\leq r\right\rbrace,\qquad r=D_{(k)},\ \ k=\mathrm{round}("
              r"\rho\,n^{3})")

    # ------------------------------------------------------------------
    if fams:
        sub("Ajuste al VOI", "Fit to the VOI")
        pesos = dict(zip(err_mod.CAMPOS, err_mod.PESOS_BASE))
        par("Cada candidato se mide con la morfometría anterior y se compara "
            "con el VOI con el error ponderado de diferencias relativas al "
            "cuadrado, normalizado por el peso REALMENTE usado (un término "
            "sin valor no cuenta como error nulo):",
            "Every candidate is measured with the morphometry above and "
            "compared with the VOI through the weighted error of squared "
            "relative differences, normalised by the weight ACTUALLY used (a "
            "term without a value does not count as a zero error):")
        e_err = eq(L, r"e=\frac{\sum_{k\in\mathcal{K}}w_k\left(\frac{m_k^{c}-"
                      r"m_k^{V}}{|m_k^{V}|}\right)^{2}}{\sum_{k\in\mathcal{K}}"
                      r"w_k},\qquad \mathcal{K}=\left\lbrace k:\ m_k^{c},\,"
                      r"m_k^{V}\ \mathrm{" + T("finitos", "finite")
                      + r"},\ m_k^{V}\neq0,\ w_k>0\right\rbrace")
        etq = {"BVTV": "BV/TV", "DA": "DA", "BSBV": "BS/BV", "TbTh": "Tb.Th",
               "TbSp": "Tb.Sp", "TbN": "Tb.N", "PoTot": "Po.tot",
               "Ez_rel": "E_{z}/E_{s}", "Ez_Ex": "E_{z}/E_{x}"}
        filas = [f"| {etq[c]} | {n(w, '.3g')} |" for c, w in pesos.items()
                 if c not in ("Ez_rel", "Ez_Ex")]
        L.append("| " + T("Término", "Term") + " | "
                 + T("Peso w", "Weight w") + " |")
        L.append("|---|---|")
        L += filas
        if "dual-lattice" in nombres:
            from .fit_dual import PESO_DA2
            L.append(f"| DA₂ ({T('solo dual-lattice', 'dual-lattice only')})"
                     f" | {n(PESO_DA2, '.3g')} |")
        L.append("")
        objetivos = [(res["ajuste" + s].get("objetivo") or {})
                     for _f, s in fams]
        if any(o.get("pesos") for o in objetivos):
            par("Atención: esta sesión usó pesos distintos de los de la "
                "tabla; constan en `resultados_sesion.json` "
                "(ajuste → objetivo → pesos).",
                "Note: this session used weights different from those in the "
                "table; they are recorded in `resultados_sesion.json` (ajuste "
                "→ objetivo → pesos).")
        if any((o.get("peso_mecanico") or 0) > 0 for o in objetivos):
            par("Se añadieron los términos mecánicos E_{z}/E_{s} y E_{z}/E_{x}, con "
                f"pesos {n(pesos['Ez_rel'])}·w_{{mec}} y {n(pesos['Ez_Ex'])}·w_{{mec}},"
                " evaluados solo entre los finalistas.",
                "The mechanical terms E_{z}/E_{s} and E_{z}/E_{x} were added, with "
                f"weights {n(pesos['Ez_rel'])}·w_{{mec}} and "
                f"{n(pesos['Ez_Ex'])}·w_{{mec}}, evaluated among finalists only.")
        from .fit import RHO_MIN, alinear_marco_fabrica
        # Las mismas expresiones que `fit.ajustar_spinodoide` (no hay constante
        # que leer); el bloque 22 comprueba que siguen coincidiendo.
        ondas_r = ", ".join(n(x, ".3g") for x in ONDAS_RAPIDO)
        ondas_c = ", ".join(n(x, ".3g") for x in ONDAS_COMPLETO)
        rmin = n(RHO_MIN, ".2g")
        par("La búsqueda es escalonada: (A) rejilla de densidad × número de "
            "onda con los ángulos fijos, (B) conjunto de ángulos de cono "
            "predefinidos con los mejores parámetros de A y (C) refinado local "
            "de densidad y número de ondas. El método «rápido» usa 3 "
            "densidades equiespaciadas entre 0,85 y 1,15 veces el BV/TV del "
            f"VOI (nunca por debajo de {rmin}), β/π ∈ {{{ondas_r}}} y 4 juegos "
            "de ángulos, sin refinado; el «completo», 5 densidades entre 0,7 y "
            f"1,3 veces el BV/TV del VOI, β/π ∈ {{{ondas_c}}}, 14 juegos de "
            "ángulos y refinado; el «equitativo» es el completo con el otro "
            "reparto de ondas. El dual-lattice usa densidad × celdas alrededor "
            "de la escala estimada con el Tb.N del VOI y conjuntos de "
            "estiramiento axiales y triaxiales.",
            "The search is staged: (A) grid of density × wave number with "
            "fixed angles, (B) set of predefined cone angles with the best "
            "parameters from A and (C) local refinement of density and number "
            "of waves. The \"fast\" method uses 3 densities equally spaced "
            f"between 0.85 and 1.15 times the VOI BV/TV (never below {rmin}), "
            f"β/π ∈ {{{ondas_r}}} and 4 angle sets, without refinement; the "
            "\"full\" method, 5 densities between 0.7 and 1.3 times the VOI "
            f"BV/TV, β/π ∈ {{{ondas_c}}}, 14 angle sets and refinement; "
            "\"equitable\" is the full method with the other wave allocation. "
            "The dual-lattice uses density × cells around the scale estimated "
            "from the VOI Tb.N and axial and triaxial stretch sets.")
        if any(isinstance(res.get("ajuste_comparado" + s), dict)
               for _f, s in fams):
            par("De los métodos que comparten la definición de error de la "
                "ecuación ({e}) se devuelve el de menor error. El desempate "
                "mecánico queda fuera porque su error suma términos que los "
                "demás no tienen.".format(e=e_err),
                "Among the methods that share the error definition of "
                "equation ({e}), the one with the lowest error is returned. "
                "The mechanical tie-break is excluded because its error adds "
                "terms the others lack.".format(e=e_err))
        plano = n(_defecto(alinear_marco_fabrica, "da2_plano"), ".3g")
        itmax = _defecto(alinear_marco_fabrica, "iter_max")
        toldeg = n(_defecto(alinear_marco_fabrica, "tol_deg"), ".3g")
        par("La orientación se impone después: se alinea el marco de "
            "autovectores del MIL del candidato con el del VOI (la normal del "
            f"plano si el candidato es plano, DA₂ < {plano}), se regenera y se "
            f"vuelve a medir, hasta {itmax} iteraciones o {toldeg}° de "
            "tolerancia, conservando la mejor orientación medida.",
            "Orientation is imposed afterwards: the MIL eigenvector frame of "
            "the candidate is aligned with that of the VOI (the plane normal "
            f"if the candidate is planar, DA₂ < {plano}), regenerated and "
            f"measured again, for up to {itmax} iterations or a {toldeg}° "
            "tolerance, keeping the best measured orientation.")
        Ks = [(res["ajuste" + s].get("incertidumbre") or {}).get("K")
              for _f, s in fams]
        if any(Ks):
            par("Como el generador es estocástico, el error de la semilla de "
                "búsqueda es el mínimo de muchas tiradas. La incertidumbre se "
                f"mide con K réplicas (mínimo {incertidumbre.K_MIN}) de "
                "semillas nuevas s + 1, …, s + K, nunca la de búsqueda; el "
                "suelo autoconsistente es el error medio entre pares "
                "ordenados de réplicas del mismo candidato, es decir, lo que "
                "difieren dos realizaciones idénticas:",
                "Because the generator is stochastic, the error of the search "
                "seed is the minimum of many draws. Uncertainty is measured "
                f"with K replicas (at least {incertidumbre.K_MIN}) of fresh "
                "seeds s + 1, …, s + K, never the search seed; the "
                "self-consistent floor is the mean error between ordered "
                "pairs of replicas of the same candidate, i.e. how much two "
                "identical realisations differ:")
            eq(L, r"\bar{e}=\frac{1}{K}\sum_{j=1}^{K}e_{j},\qquad s_e=\sqrt{"
                  r"\frac{1}{K-1}\sum_{j=1}^{K}(e_{j}-\bar{e})^{2}},\qquad "
                  r"e_{\mathrm{suelo}}=\frac{1}{K(K-1)}\sum_{i\neq j}e(m_i,"
                  r"\,m_j)")

    # ------------------------------------------------------------------
    if ela:
        rec = ela[0]
        tol = _defecto(homogeneizar, "tol")
        sub("Homogeneización periódica", "Periodic homogenisation")
        par(f"El tensor de rigidez efectivo se calcula por homogeneización "
            f"numérica periódica {r.c('andreassen2014')} sobre la rejilla de "
            "vóxeles, remuestreada al vecino más próximo a "
            f"{_i(rec.get('resolucion'))}³ (el vóxel nuevo i' toma el original "
            "i = round(i'·(n − 1)/(n' − 1)); arista h' = h·n/n'). Cada vóxel "
            "es un hexaedro "
            "trilineal de 8 nodos; su matriz de rigidez se integra con "
            "cuadratura de Gauss 2×2×2:",
            "The effective stiffness tensor is computed by periodic numerical "
            f"homogenisation {r.c('andreassen2014')} on the voxel grid, "
            "nearest-neighbour resampled to "
            f"{_i(rec.get('resolucion'))}³ (the new voxel i' takes the "
            "original i = round(i'·(n − 1)/(n' − 1)); edge h' = h·n/n'). Every "
            "voxel is an 8-node trilinear "
            "hexahedron; its stiffness matrix is integrated with 2×2×2 Gauss "
            "quadrature:")
        eq(L, r"N_a=\frac{1}{8}(1+\xi\xi_a)(1+\eta\eta_a)(1+\zeta\zeta_a),"
              r"\qquad \mathbf{k}_e=\int_{V_e}\mathbf{B}^{T}\mathbf{D}\,"
              r"\mathbf{B}\,dV")
        par(f"Hueso con E_{{s}} = {n((_f(rec, 'E_s_Pa') or np.nan) / 1e9, '.4g')} "
            f"GPa y ν_{{s}} = {n(rec.get('nu_s'), '.3g')}; los poros NO se quitan: "
            f"se rellenan con rigidez s_{{e}} = {n(VOID_SCALE, '.0e')}·E_{{s}} para "
            "que la celda no sea singular. Con condiciones de contorno "
            "periódicas y los seis estados de deformación macroscópica "
            "unitaria ε⁰(k) (desplazamientos nodales u⁰(k) = ε⁰(k)·x), el "
            "campo corrector χ(k) y el tensor homogeneizado son",
            f"Bone with E_{{s}} = {n((_f(rec, 'E_s_Pa') or np.nan) / 1e9, '.4g')} "
            f"GPa and ν_{{s}} = {n(rec.get('nu_s'), '.3g')}; pores are NOT "
            f"removed: they are filled with stiffness s_{{e}} = "
            f"{n(VOID_SCALE, '.0e')}·E_{{s}} so that the cell is not singular. "
            "With periodic boundary conditions and the six unit macroscopic "
            "strain states ε⁰(k) (nodal displacements u⁰(k) = ε⁰(k)·x), the "
            "corrector field χ(k) and the homogenised tensor are")
        eq(L, r"\mathbf{K}\,\mathbf{\chi}^{(k)}=\sum_e s_e\,"
                        r"\mathbf{k}_e\,\mathbf{u}^{0(k)}_e,\qquad C_{ij}="
                        r"\frac{1}{|V|}\sum_e s_e\left(\mathbf{u}^{0(i)}_e-"
                        r"\mathbf{\chi}^{(i)}_e\right)^{T}\mathbf{k}_e\left("
                        r"\mathbf{u}^{0(j)}_e-\mathbf{\chi}^{(j)}_e\right)")
        par("K es singular solo por las traslaciones rígidas y se fija el "
            f"nodo 0. Por debajo de {UMBRAL_DIRECTO} grados de libertad se "
            "resuelve con LU directa; por encima, con gradiente conjugado "
            "precondicionado con multimalla algebraica de agregación "
            f"suavizada {r.c('bell2022')} (tolerancia {n(tol, '.0e')}, 500 "
            "iteraciones como máximo). Un solver iterativo no avisa cuando no "
            "converge, así que en ese caso se comprueba el residuo relativo "
            "de los seis casos y el tensor se rechaza si",
            "K is singular only through rigid translations and node 0 is "
            f"fixed. Below {UMBRAL_DIRECTO} degrees of freedom a direct LU "
            "factorisation is used; above, conjugate gradients preconditioned "
            f"with smoothed-aggregation algebraic multigrid {r.c('bell2022')} "
            f"(tolerance {n(tol, '.0e')}, at most 500 iterations). An "
            "iterative solver does not warn when it fails to converge, so in "
            "that case the relative residual of the six cases is checked and "
            "the tensor is rejected if")
        eq(L, r"\max_{k}\ \frac{\Vert\mathbf{f}^{(k)}-\mathbf{K}\,\mathbf{\chi}"
              r"^{(k)}\Vert}{\Vert\mathbf{f}^{(k)}\Vert}>\max\left(10^{-6},\ "
              r"100\,\mathrm{tol}\right)")
        par("Las constantes de ingeniería salen de la flexibilidad S = C⁻¹, y "
            "el módulo de Young en una dirección cualquiera n (figura de "
            "anisotropía) de la misma S:",
            "Engineering constants come from the compliance S = C⁻¹, and "
            "Young's modulus along any direction n (anisotropy figure) from "
            "the same S:")
        e_En = eq(L, r"E_x=\frac{1}{S_{11}},\ \ G_{yz}=\frac{1}{S_{44}},\ \ "
                     r"\nu_{xy}=-\frac{S_{12}}{S_{11}},\qquad \frac{1}{E("
                     r"\mathbf{n})}=\mathbf{a}^{T}\mathbf{S}\,\mathbf{a},\ \ "
                     r"\mathbf{a}=[n_1^{2},n_2^{2},n_3^{2},n_2n_3,n_1n_3,"
                     r"n_1n_2]^{T}")

    # ------------------------------------------------------------------
    if fe or comp:
        rec = (fe or comp)[0]
        tol = _defecto(ensayo_compresion, "tol")
        sub("Ensayo de compresión por elementos finitos",
            "Finite-element compression test")
        par("Al contrario que la homogeneización, el ensayo malla SOLO el "
            f"hueso {r.c('vanrietbergen1995')}: rellenar los poros con un "
            "material blando arruina el condicionamiento en estructuras poco "
            "densas. La máscara se remuestrea al vecino más próximo a la "
            "resolución del ensayo y se conservan únicamente las componentes "
            "conexas por caras (vecindad 6) que tocan a la vez la capa inferior "
            "y la superior: un fragmento que no une los dos platos no "
            "transmite carga. La fracción portante f_{port} es la parte del "
            "hueso que sobrevive a ese filtro.",
            "Unlike homogenisation, the test meshes ONLY the bone "
            f"{r.c('vanrietbergen1995')}: filling pores with a soft material "
            "ruins the conditioning of low-density structures. The mask is "
            "nearest-neighbour resampled to the test resolution and only the "
            "face-connected components (6-neighbourhood) that touch both the "
            "bottom and the top layers are kept: a fragment that does not "
            "join the two platens carries no load. The load-bearing fraction "
            "f_{port} is the part of the bone that survives that filter.")
        par("La carga es una compresión en −z repartida por área tributaria: "
            "cada uno de los N_{top} elementos con la cara superior en el techo "
            "da un cuarto de su parte a cada uno de sus cuatro nodos "
            "superiores, que es el equivalente nodal consistente de una "
            "presión uniforme. Apoyos: «empotrado» fija los tres "
            "desplazamientos de todos los nodos de la base; «deslizante» fija "
            "u_{z} = 0 en la base y, para eliminar solo los movimientos de "
            "sólido rígido restantes, u_{x} = u_{y} = 0 en la esquina (x mín, "
            "y mín) y u_{y} = 0 en la esquina (x máx, y mín).",
            "The load is a compression along −z distributed by tributary "
            "area: each of the N_{top} elements whose upper face lies on the top "
            "gives a quarter of its share to each of its four upper nodes, "
            "which is the consistent nodal equivalent of a uniform pressure. "
            "Supports: \"empotrado\" (fixed) clamps the three displacements of "
            "every base node; \"deslizante\" (sliding) sets u_{z} = 0 on the "
            "base and, to remove only the remaining rigid-body motions, u_{x} = "
            "u_{y} = 0 at the (min x, min y) corner and u_{y} = 0 at the (max x, "
            "min y) corner.")
        eq(L, r"f_{z,a}=-\frac{F}{4\,N_{\mathrm{top}}}\quad\mathrm{"
              + T(r"por\ elemento\ del\ techo\ y\ nodo\ superior",
                  r"per\ top\ element\ and\ upper\ node")
              + r"},\qquad F=\sigma_0\,A")
        par(f"Resolución: LU directa por debajo de {UMBRAL_DIRECTO} grados de "
            "libertad; por encima, multimalla algebraica con los SEIS modos de "
            "sólido rígido (tres traslaciones y tres rotaciones infinitesimales)"
            " como espacio casi nulo, acelerada con gradiente conjugado "
            f"(tolerancia {n(tol, '.0e')}); con el solver iterativo se rechaza "
            "la solución si el residuo relativo supera max(10⁻⁶, 100·tol). "
            "Deformación y tensión se evalúan en "
            "el centro de cada elemento (ξ = η = ζ = 0), y de ahí la densidad "
            "de energía, la deformación efectiva de Pistoia y la tensión de "
            "von Mises:",
            f"Solution: direct LU below {UMBRAL_DIRECTO} degrees of freedom; "
            "above, algebraic multigrid with the SIX rigid-body modes (three "
            "translations and three infinitesimal rotations) as near-null "
            "space, accelerated with conjugate gradients (tolerance "
            f"{n(tol, '.0e')}); with the iterative solver the solution is "
            "rejected if the relative residual exceeds max(10⁻⁶, 100·tol). "
            "Strain and stress are evaluated at the "
            "centre of each element (ξ = η = ζ = 0), and from them the energy "
            "density, Pistoia's effective strain and the von Mises stress:")
        e_sig = eq(L, r"\mathbf{\varepsilon}_e=\mathbf{B}(\mathbf{0})\,"
                      r"\mathbf{u}_e,\qquad \mathbf{\sigma}_e=\mathbf{D}\,"
                      r"\mathbf{\varepsilon}_e,\qquad U_e=\frac{1}{2}\,"
                      r"\mathbf{\sigma}_e^{T}\mathbf{\varepsilon}_e,"
                      r"\qquad \varepsilon_{\mathrm{eff}}=\sqrt{\frac{2\,U_e}"
                      r"{E_s}}")
        e_vm = eq(L, r"\sigma_{\mathrm{vM}}=\sqrt{\frac{1}{2}\left[(\sigma_{xx}"
                     r"-\sigma_{yy})^{2}+(\sigma_{yy}-\sigma_{zz})^{2}+(\sigma"
                     r"_{zz}-\sigma_{xx})^{2}\right]+3\left(\tau_{yz}^{2}+\tau_"
                     r"{xz}^{2}+\tau_{xy}^{2}\right)}=\sqrt{3\,J_2}")
        par("La ecuación ({v}) usa las tensiones tangenciales τ, no las "
            "deformaciones angulares γ, y es el invariante J₂ del desviador: "
            "no depende de la parte hidrostática. Por eso no es un reescalado "
            "de ε_{{eff}}, que sale de la energía total; coinciden salvo un factor "
            "solo en estado uniaxial. La rigidez aparente sale del "
            "desplazamiento medio del techo, con H la altura de la probeta:"
            .format(v=e_vm),
            "Equation ({v}) uses the shear stresses τ, not the engineering "
            "shear strains γ, and is the J₂ invariant of the deviator: it "
            "does not depend on the hydrostatic part. That is why it is not a "
            "rescaling of ε_{{eff}}, which comes from the total energy; they "
            "coincide up to a factor only in a uniaxial state. The apparent "
            "stiffness comes from the mean displacement of the top, with H "
            "the height of the specimen:".format(v=e_vm))
        eq(L, r"\varepsilon_{\mathrm{app}}=\frac{\left|\langle u_{z}"
                       r"\rangle_{\mathrm{top}}\right|}{H},\qquad E_{\mathrm"
                       r"{app}}=\frac{\sigma_0}{\varepsilon_{\mathrm{app}}},"
                       r"\qquad |\mathbf{u}|=\sqrt{u_x^{2}+u_{y}^{2}+u_{z}^{2}}")
        par("Los ensayos en X e Y aplican una permutación CÍCLICA de los ejes "
            "(determinante +1, la terna sigue siendo dextrógira), resuelven el "
            "mismo problema en z y devuelven los campos al marco original. "
            f"Sin carga impuesta, la tensión de referencia es σ₀ = "
            f"{n(SIGMA0_DEF / 1e6, '.3g')} MPa; al ser el problema lineal, "
            "E_{app} no depende de ella.",
            "Tests along X and Y apply a CYCLIC permutation of the axes "
            "(determinant +1, the frame stays right-handed), solve the same "
            "problem along z and return the fields to the original frame. "
            f"Without an imposed load, the reference stress is σ₀ = "
            f"{n(SIGMA0_DEF / 1e6, '.3g')} MPa; since the problem is linear, "
            "E_{app} does not depend on it.")

        sub("Estadísticos de la tensión de von Mises",
            "Von Mises stress statistics")
        par("El máximo de la ecuación ({v}) sobre una malla de vóxeles NO "
            "converge con la resolución: en una cavidad esférica periódica "
            "con solución cerrada el pico binario supera el valor exacto de "
            "Goodier en +5,6 % en la malla más fina y oscila ~10 puntos al "
            "refinar. Se informa en cambio un percentil alto sobre la CAPA "
            "SUPERFICIAL S, los elementos de hueso con al menos un vecino por "
            "cara en el vacío. Fuera del dominio se supone material, de modo "
            "que las seis caras del cubo (incluidas la cargada y la apoyada) "
            "no cuentan como superficie libre:".format(v=e_vm),
            "The maximum of equation ({v}) on a voxel mesh does NOT converge "
            "with resolution: on a periodic spherical cavity with a "
            "closed-form answer the binary peak exceeds Goodier's exact value "
            "by +5.6 % on the finest mesh and oscillates by ~10 points under "
            "refinement. A high percentile over the SURFACE LAYER S is "
            "reported instead: the bone elements with at least one face "
            "neighbour in the void. Outside the domain counts as material, so "
            "the six faces of the cube (including the loaded and supported "
            "ones) do not count as free surface:".format(v=e_vm))
        eq(L, r"S=\left\lbrace e\in\Omega_h:\ \exists\,e'\in\mathcal"
              r"{N}_6(e),\ B_{e'}=0\right\rbrace,\qquad B\equiv1\ \mathrm{"
              + T(r"fuera\ del\ cubo", r"outside\ the\ cube") + r"}")
        eq(L, r"\sigma_{\mathrm{vM}}^{p99,S}=Q_5\left(0.99;\ \left\lbrace"
              r"\sigma_{\mathrm{vM},e}\right\rbrace_{e\in S}\right),\qquad "
              r"\sigma_{\mathrm{vM}}^{p99}=Q_7\left(0.99;\ \left\lbrace\sigma"
              r"_{\mathrm{vM},e}\right\rbrace_{e\in\Omega_h}\right)")
        par("Sobre esa cavidad el percentil 99 de la capa queda a −0,3 % del "
            f"valor exacto. Se declaran tres sabores del pico: el máximo (no "
            "citable), el percentil 99 de todo el tejido con la definición "
            "Q₇ de la ecuación ({q}) (con reservas: en estructuras densas lo "
            "domina el material a granel a la tensión nominal) y el de la "
            "capa con Q₅ (citable si la capa tiene al menos "
            f"{I.N_SUPERFICIE_MIN} elementos, un umbral propio de este "
            "proyecto).".format(q=e_q7),
            "On that cavity the surface-layer 99th percentile lies within "
            "−0.3 % of the exact value. Three flavours of the peak are "
            "declared: the maximum (not citable), the 99th percentile over "
            "all tissue with the Q₇ definition of equation ({q}) (with "
            "caveats: in dense structures it is dominated by bulk material at "
            "the nominal stress) and the surface-layer one with Q₅ (citable "
            f"if the layer has at least {I.N_SUPERFICIE_MIN} elements, a "
            "threshold of this project).".format(q=e_q7))

        if fe:
            sub("Carga de fallo: criterio de Pistoia",
                "Failure load: Pistoia criterion")
            par(f"El hueso se declara roto cuando el {n(100 * FRAC_CRITICA, '.3g')}"
                f" % del tejido supera una deformación efectiva de "
                f"{n(100 * EPS_CRITICA, '.3g')} % {r.c('pistoia2002')}. Como "
                "el problema es lineal, todos los campos escalan con la carga "
                "y el factor de escala sale de un cociente:",
                f"Bone is declared failed when {n(100 * FRAC_CRITICA, '.3g')} "
                f"% of the tissue exceeds an effective strain of "
                f"{n(100 * EPS_CRITICA, '.3g')} % {r.c('pistoia2002')}. Since "
                "the problem is linear, every field scales with the load and "
                "the scale factor is a quotient:")
            eq(L, r"\varepsilon^{*}=Q_7\left(1-f;\ \left\lbrace\varepsilon_{"
                  r"\mathrm{eff},e}\right\rbrace\right),\qquad k=\frac{"
                  r"\varepsilon_{\mathrm{crit}}}{\varepsilon^{*}},\qquad "
                  r"\sigma_{\mathrm{fallo}}=k\,\sigma_0,\qquad F_{\mathrm{fallo}"
                  r"}=k\,F,\qquad \sigma_{\mathrm{vM,fallo}}=k\,\sigma_{\mathrm"
                  r"{vM}}")
            par(f"Con f = {n(FRAC_CRITICA, '.3g')} y ε_{{crit}} = "
                f"{n(EPS_CRITICA, '.3g')}. Los dos parámetros se calibraron en "
                "radio distal humano y no son constantes físicas: la carga de "
                "fallo es una estimación comparativa entre estructuras "
                "ensayadas igual, no una predicción absoluta para este hueso.",
                f"With f = {n(FRAC_CRITICA, '.3g')} and ε_{{crit}} = "
                f"{n(EPS_CRITICA, '.3g')}. Both parameters were calibrated on "
                "the human distal radius and are not physical constants: the "
                "failure load is a comparative estimate between structures "
                "tested alike, not an absolute prediction for this bone.")

        if comp:
            c0 = comp[0]
            sub("Protocolo del análisis comparado", "Compared-analysis protocol")
            par(f"Aplica a VOI y candidato el mismo ensayo de la ecuación "
                f"({e_sig}) con el protocolo de {c0.get('protocolo') or '—'}: "
                f"tejido de E_{{s}} = {n((_f(c0, 'E_s_Pa') or np.nan) / 1e9, '.4g')}"
                f" GPa y ν_{{s}} = {n(c0.get('nu_s'), '.3g')}, fuerza total F = "
                f"{n(c0.get('carga_N'), '.4g')} N en −z, apoyo "
                f"«{I.nombre_apoyo(c0.get('apoyo'), 'es')}», a {_i(c0.get('resolucion'))}³. "
                "Se informan la tensión de von Mises (media, máximo y "
                "percentiles de la sección anterior, en MPa) y la deformación "
                "total |u|, cuyo máximo nodal está siempre en la cara cargada. "
                "Para saber si las dos estructuras se cargan por el mismo eje "
                "se mide el MIL sobre la MISMA malla ensayada y se da el ángulo "
                "entre direcciones principales:",
                "It applies to VOI and candidate the same test of equation "
                f"({e_sig}) with the protocol of {c0.get('protocolo') or '—'}:"
                f" tissue E_{{s}} = {n((_f(c0, 'E_s_Pa') or np.nan) / 1e9, '.4g')} "
                f"GPa and ν_{{s}} = {n(c0.get('nu_s'), '.3g')}, total force F = "
                f"{n(c0.get('carga_N'), '.4g')} N along −z, support "
                f"\"{I.nombre_apoyo(c0.get('apoyo'), 'en')}\", at "
                f"{_i(c0.get('resolucion'))}³. The von Mises stress (mean, "
                "maximum and the percentiles of the previous section, in MPa) "
                "and the total deformation |u|, whose nodal maximum is always "
                "on the loaded face, are reported. To know whether both "
                "structures are loaded along the same axis, the MIL is "
                "measured on the SAME tested mesh and the angle between "
                "principal directions is given:")
            eq(L, r"\alpha=\arccos\left|\mathbf{d}_V\cdot\mathbf{d}_C\right|,"
                  r"\qquad \alpha>" + f"{avisos.ANGULO_MECANICO_MAX_DEG:.0f}"
                  + r"^{\circ}\ \Rightarrow\ \mathrm{"
                  + T(r"comparaci\acute{o}n\ mec\acute{a}nica\ no\ "
                      r"interpretable",
                      r"mechanical\ comparison\ not\ interpretable") + r"}")

    # ------------------------------------------------------------------
    from . import figuras as F
    en_es = (f" y E(n) de la ecuación ({e_En})" if ela else "")
    en_en = (f" and E(n) of equation ({e_En})" if ela else "")
    sub("Magnitudes de las figuras", "Quantities shown in the figures")
    par("Figura de anisotropía: MIL(n) = (nᵀ M n)^(−½)"
        + en_es + ", en los planos XY, XZ e YZ. Distribuciones: histogramas de "
        "densidad del espesor local de hueso y de poro, con clases del ancho "
        "del vóxel más grueso (el espesor local toma valores discretos, "
        "múltiplos de h). Cortes: planos centrales de la máscara binaria, "
        "sin suavizar. Renders 3D: isosuperficie cerrada contra las caras del "
        f"cubo, opcionalmente suavizada con Taubin {r.c('taubin1995')} "
        f"({F.TAUBIN_ITER} iteraciones, banda de paso "
        f"{n(F.TAUBIN_BANDA, '.2g')}), que no encoge la pieza; el "
        "suavizado es SOLO de la figura y ninguna medida lo usa. Proyección "
        "paralela con una cámara que depende solo del cubo, de modo que todas "
        "las estructuras salen a la misma escala. Mapa de von Mises: cada "
        "vértice toma el valor del vóxel sólido más próximo (sin filtro de "
        "máximo, que inflaría los picos) y la escala de color común va de "
        "los percentiles 1 a 99 de todo el tejido de todas las "
        "estructuras.",
        "Anisotropy figure: MIL(n) = (nᵀ M n)^(−½)"
        + en_en + ", on the XY, XZ and YZ planes. Distributions: density "
        "histograms of the local thickness of bone and pore, with bins as "
        "wide as the coarsest voxel (local thickness takes discrete values, "
        "multiples of h). Sections: central planes of the binary mask, "
        "unsmoothed. 3D renders: isosurface closed against the cube faces, "
        f"optionally smoothed with Taubin {r.c('taubin1995')} "
        f"({F.TAUBIN_ITER} iterations, pass band "
        f"{n(F.TAUBIN_BANDA, '.2g')}), which does not shrink the part; "
        "smoothing is for "
        "the FIGURE only and no measurement uses it. Parallel projection with "
        "a camera that depends only on the cube, so every structure is drawn "
        "at the same scale. Von Mises map: each vertex takes the value of the "
        "nearest solid voxel (no maximum filter, which would inflate peaks) "
        "and the common colour scale spans the 1st to 99th percentiles of all "
        "the tissue of all structures.")

    # ------------------------------------------------------------------
    sub("Criterios de citabilidad", "Citability criteria")
    par("Cada resultado de la sección de citabilidad se gradúa con las reglas "
        "siguientes. «Medido» indica que el umbral sale de un estudio del "
        "proyecto; «propio», que es una convención declarada.",
        "Every result in the citability section is graded with the following "
        "rules. \"Measured\" means the threshold comes from a project study; "
        "\"own\", that it is a declared convention.")
    G = I.GRAVEDAD
    E = {k: I.ESTADOS[v][0 if es else 1] for k, v in G.items()}
    med, pro = T("medido", "measured"), T("propio", "own")
    reglas = [
        (T("Residuo relativo del solver", "Solver relative residual"),
         f"> {n(I.RESIDUO_MAX, '.0e')}", E["residuo"],
         T("corrección F10", "correction F10")),
        (T("Resolución del cálculo mecánico", "Mechanical resolution"),
         f"n < {I.N_MECANICA_MIN}", E["resolucion_mecanica"],
         med + " (Estudio_Convergencia)"),
        (T("Hueso desconectado", "Disconnected bone"),
         f"> {n(100 * I.DESCONEXION_MAX, '.3g')} %", E["desconexion_alta"],
         pro),
        (T("Hueso desconectado", "Disconnected bone"),
         f"> {n(100 * I.DESCONEXION_RESERVAS, '.3g')} %", E["desconexion"],
         pro),
        (T("Máximo de von Mises", "Von Mises maximum"), T("siempre", "always"),
         E["vm_maximo"], med + " (Estudio_Convergencia)"),
        (T("Percentil 99 de todo el tejido", "Whole-tissue 99th percentile"),
         T("siempre", "always"), E["vm_p99_global"],
         med + " (Estudio_Convergencia)"),
        (T("Percentil 99 de la capa superficial",
           "Surface-layer 99th percentile"),
         f"n_{{S}} < {I.N_SUPERFICIE_MIN}", E["vm_superficie_pocos"], pro),
        (T("Carga de fallo (Pistoia)", "Failure load (Pistoia)"),
         T("siempre", "always"), E["pistoia_calibracion"],
         T("calibración ajena", "external calibration")),
        (T("Dirección principal", "Principal direction"),
         f"DA < {n(avisos.DA_MIN_EJE, '.3g')}", E["da_sin_eje"], med),
        (T("VOI no trabecular", "Non-trabecular VOI"),
         f"BV/TV ≥ {n(avisos.BVTV_NO_TRABECULAR, '.3g')}",
         E["voi_no_trabecular"], med + " (Estudio_Familias)"),
        (T("Candidato desalineado", "Misaligned candidate"),
         f"α > {avisos.ANGULO_MECANICO_MAX_DEG:.0f}°", E["desalineado"], med),
        (T("Réplicas del error de ajuste", "Fit error replicas"),
         f"K < {incertidumbre.K_MIN}", E["sin_replicas"], pro),
        (T("Muestreo MIL suavizado (M1)", "Smoothed MIL sampling (M1)"),
         f"Tb.Th/h < {n(I.TBTH_H_M1_MIN, '.2g')}", E["m1_resolucion"],
         med + " (Estudio_MIL)"),
        (T("Ángulos de cono en la región degenerada",
           "Cone angles in the degenerate region"),
         T("unión de conos = esfera", "cone union = sphere"),
         E["thetas_degenerados"], T("analítico", "analytic")),
        (T("Número de onda sin sus dos lecturas",
           "Wave number without both readings"),
         T("clave «wave» o lecturas incoherentes",
           "\"wave\" key or inconsistent readings"),
         E["onda_inconsistente"], med + " (Estudio_Percolacion)"),
    ]
    L.append("| " + " | ".join((T("Comprobación", "Check"),
                                T("Condición", "Condition"),
                                T("Estado", "Status"),
                                T("Origen", "Origin"))) + " |")
    L.append("|---|---|---|---|")
    for fila in reglas:
        L.append("| " + " | ".join(fila) + " |")
    L.append("")

    L.append("> " + T(
        "**En palabras sencillas.** Aquí está la «receta» exacta de cada "
        "número del informe: qué ecuación lo produce, sobre qué parte de la "
        "estructura se calcula y con qué tolerancias. Si un revisor pregunta "
        "cómo se obtuvo la tensión de von Mises o por qué no se da su valor "
        "máximo, la respuesta está en esta sección, y las ecuaciones están "
        "comprobadas contra el propio programa.",
        "**In plain words.** This is the exact \"recipe\" of every number in "
        "the report: which equation produces it, on which part of the "
        "structure it is computed and with which tolerances. If a reviewer "
        "asks how the von Mises stress was obtained or why its maximum is not "
        "reported, the answer is in this section, and the equations are "
        "checked against the program itself."))
    L.append("")
    return L


def _i(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return "—"


def _f(d, clave):
    try:
        x = float((d or {}).get(clave))
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None
