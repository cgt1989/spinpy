# -*- coding: utf-8 -*-
"""
idioma_textos.py — Diccionario espanol -> ingles de la interfaz.

Solo datos. La logica esta en `idioma.py`.

REGLAS AL ANADIR O CAMBIAR UN TEXTO
  * La clave es la cadena espanola EXACTA, con sus tildes, sus comillas
    angulares y sus puntos suspensivos de un solo caracter. Una diferencia de
    un caracter deja el texto sin traducir, en silencio.
  * Los marcadores {} tienen que aparecer en los DOS idiomas y con el mismo
    nombre. `idioma_revisar.py` lo comprueba.
  * Los nombres de metrica no se traducen: BV/TV, Tb.Th, Tb.Sp, Tb.N, BS/BV,
    DA, SMI, Conn.D son notacion internacional (Parfitt / Bouxsein et al.
    2010). Traducirlos haria la tabla ilegible para un lector de hueso.
  * Las unidades tampoco: mm, mm2, MPa, GPa.
  * Las citas se dejan como estan. 'Salmon et al. 2015' no tiene version
    espanola.

Comprobar el estado:  python idioma_revisar.py
"""

EN = {

    # -- procedencia, objetivo e incertidumbre -----------------------------
    "thetas en la region degenerada: la estructura es la isotropa (90, 90, "
    "90) y mover estos angulos no cambia nada":
        "thetas in the degenerate region: the structure is the isotropic "
        "(90, 90, 90) one and moving these angles changes nothing",
    "thetas junto a la region degenerada: {p}% de probabilidad de que esta "
    "realizacion sea exactamente la isotropa":
        "thetas next to the degenerate region: {p}% probability that this "
        "realisation is exactly the isotropic one",
    "orientacion impuesta sin llegar a la tolerancia: eje a {ang}°":
        "orientation imposed without reaching the tolerance: axis at {ang}°",
    "error con {K} semillas nuevas: {m} ± {s} (suelo {f})":
        "error over {K} fresh seeds: {m} ± {s} (floor {f})",
    "Sesion en un formato anterior, sin procedencia: no consta con que "
    "version de spinpy se guardo. El numero de onda se lee como multiplo de "
    "pi, que es lo que siempre guardo el deslizador.":
        "Session in an older format, without provenance: the spinpy version "
        "that saved it is unknown. The wave number is read as a multiple of "
        "pi, which is what the slider always stored.",
    "El numero de onda de la sesion no cuadra con su procedencia: {e}":
        "The session's wave number does not match its provenance: {e}",

    # -- menus y ventana ---------------------------------------------------
    "&Sesion": "&Session",
    "Guardar sesion…": "Save session…",
    "Cargar sesion…": "Load session…",
    "Exportar resultados (CSV + JSON)…": "Export results (CSV + JSON)…",
    "Cargar VOI…": "Load VOI…",
    "Cargar pila TIFF (micro-CT)…": "Load TIFF stack (micro-CT)…",
    "&Idioma": "&Language",
    "Español": "Español",
    "English": "English",
    "Listo. Carga un VOI o pulsa Generar.":
        "Ready. Load a VOI or press Generate.",

    # -- cabeceras de seccion ----------------------------------------------
    "VOI de referencia": "Reference VOI",
    "Spinodoide": "Spinodoid",
    "Visualizacion": "Visualisation",
    "Morfometria": "Morphometry",
    "Ajuste al VOI": "Fit to the VOI",
    "Analisis mecanico": "Mechanical analysis",
    "Exportacion": "Export",
    "Lote y varianza": "Batch and variance",

    # -- pies de seccion: librerias ----------------------------------------
    "Librerias de Python que hacen el calculo de esta seccion.":
        "Python libraries behind this section's computation.",
    "Librerias: numpy · scipy.io para los VOI en .mat · scikit-image y Pillow "
    "para las pilas TIFF":
        "Libraries: numpy · scipy.io for .mat VOIs · scikit-image and Pillow "
        "for TIFF stacks",
    "Carpeta de rebanadas o TIFF multipagina, con su tamano de voxel real":
        "Slice folder or multipage TIFF, with its real voxel size",
    "Librerias: numpy · scipy.special (el erf⁻¹ del umbral) · pyvista/VTK "
    "para la isosuperficie de la vista":
        "Libraries: numpy · scipy.special (the erf⁻¹ of the threshold) · "
        "pyvista/VTK for the preview isosurface",
    "Librerias: pyvista y pyvistaqt (los dos paneles) · matplotlib "
    "(histograma) · scipy.ndimage (espesor local por transformada de "
    "distancia)":
        "Libraries: pyvista and pyvistaqt (both panels) · matplotlib "
        "(histogram) · scipy.ndimage (local thickness by distance transform)",
    "Librerias: numpy · scipy.ndimage (componentes conexas y distancias) · "
    "scikit-image (marching cubes y numero de Euler) · pyvista y pymeshfix "
    "en el modo malla":
        "Libraries: numpy · scipy.ndimage (connected components and "
        "distances) · scikit-image (marching cubes and Euler number) · "
        "pyvista and pymeshfix in mesh mode",
    "Librerias: numpy · scipy.ndimage · scipy.sparse solo si se activa el "
    "desempate mecanico":
        "Libraries: numpy · scipy.ndimage · scipy.sparse only if the "
        "mechanical tie-break is enabled",
    "Librerias: numpy · scipy.sparse y scipy.sparse.linalg (cg, splu) · "
    "matplotlib (distribuciones y convergencia)":
        "Libraries: numpy · scipy.sparse and scipy.sparse.linalg (cg, splu) "
        "· matplotlib (distributions and convergence)",
    "Librerias: pyvista (.vtu y .stl) · tetgen (TET10) · pymeshfix "
    "(estanqueidad de la superficie)":
        "Libraries: pyvista (.vtu and .stl) · tetgen (TET10) · pymeshfix "
        "(surface watertightness)",
    "Librerias: numpy · pandas · joblib (paralelo) · scipy.stats (TOST)":
        "Libraries: numpy · pandas · joblib (parallel) · scipy.stats (TOST)",

    # -- VOI ---------------------------------------------------------------
    "Cargar VOI (.vtk / .mat)…": "Load VOI (.vtk / .mat)…",

    # -- carga de pila TIFF ------------------------------------------------
    "Cargar pila de imagenes (micro-CT)": "Load image stack (micro-CT)",
    "Cargar pila TIFF": "Load TIFF stack",
    "Origen": "Source",
    "Patron de archivo": "File pattern",
    "Tamano de voxel": "Voxel size",
    "Umbral": "Threshold",
    "automatico": "automatic",
    "sin comprobar": "not checked yet",
    "Carpeta…": "Folder…",
    "Archivo…": "File…",
    "Carpeta con las rebanadas": "Folder with the slices",
    "TIFF multipagina": "Multipage TIFF",
    "TIFF (*.tif *.tiff);;Todos (*)": "TIFF (*.tif *.tiff);;All (*)",
    "carpeta de rebanadas o TIFF multipagina":
        "slice folder or multipage TIFF",
    "Falta el origen": "No source selected",
    "Elige la carpeta de rebanadas o el TIFF multipagina.":
        "Choose the slice folder or the multipage TIFF.",
    "Leyendo rebanadas…": "Reading slices…",
    "Leyendo rebanada {i} de {n}…": "Reading slice {i} of {n}…",
    "Volumen muy grande": "Very large volume",
    "<b>La imagen debe estar en formato TIFF (.tif / .tiff)</b>, que es lo que "
    "exportan los micro-CT.<br><br>Se admiten dos formas:<br>&nbsp;&nbsp;• una "
    "<b>carpeta</b> con una rebanada por archivo (se ordenan por el numero del "
    "nombre), o<br>&nbsp;&nbsp;• un unico <b>TIFF multipagina</b>.<br><br>Todas "
    "las rebanadas deben tener el <b>mismo tamano</b>: reescalar una cambiaria "
    "la morfometria sin avisar, asi que la carga se detiene si alguna difiere. "
    "Las imagenes pueden estar ya binarizadas o en escala de grises.":
        "<b>Images must be TIFF (.tif / .tiff)</b>, which is what micro-CT "
        "scanners export.<br><br>Two layouts are accepted:<br>&nbsp;&nbsp;• a "
        "<b>folder</b> with one slice per file (sorted by the number in the "
        "name), or<br>&nbsp;&nbsp;• a single <b>multipage TIFF</b>.<br><br>All "
        "slices must have the <b>same size</b>: rescaling one would change the "
        "morphometry without warning, so loading stops if any differs. Images "
        "may already be binarised or in greyscale.",
    "Solo para carpetas. Util cuando conviven varias series en el mismo sitio: "
    "en H4/Segmentadas hay BW_*.tif y SEG_*.tif mezcladas.":
        "Folders only. Useful when several series share a location: "
        "H4/Segmentadas holds BW_*.tif and SEG_*.tif mixed together.",
    "Si la pila ya es binaria toma > minimo. Si viene en escala de grises "
    "aplica Otsu y lo deja anotado.":
        "If the stack is already binary it takes > minimum. If it comes in "
        "greyscale it applies Otsu and records it.",
    "No se encontro ningun <tt>*_rec.log</tt> del escaner. <b>Escribe el "
    "tamano de voxel a mano</b> — no hay valor por defecto que sea seguro.":
        "No scanner <tt>*_rec.log</tt> was found. <b>Type the voxel size "
        "manually</b> — there is no safe default.",
    "Leido de <tt>{log}</tt> — <tt>{linea}</tt>":
        "Read from <tt>{log}</tt> — <tt>{linea}</tt>",
    "La pila tiene {nv:,} voxeles ({f}).\n\nSe puede visualizar y medir la "
    "morfometria, pero el analisis mecanico y el ajuste NO son viables a este "
    "tamano: el coste va con el cubo del lado.\n\nPara eso hay que recortar un "
    "VOI cubico (spinpy.voi.vois_por_tercios / extraer_cubo).":
        "The stack has {nv:,} voxels ({f}).\n\nIt can be displayed and its "
        "morphometry measured, but mechanical analysis and fitting are NOT "
        "viable at this size: cost grows with the cube of the side.\n\nFor "
        "that, crop a cubic VOI (spinpy.voi.vois_por_tercios / extraer_cubo).",

    # -- recorte del VOI cubico --------------------------------------------
    "Recortar VOI cubico": "Crop cubic VOI",
    "Recortar VOI cubico…": "Crop cubic VOI…",
    "Orientando el hueso por PCA…": "Orienting the bone by PCA…",
    "Orienta el hueso por PCA y recorta un cubo: tres tramos sugeridos o "
    "posicion libre":
        "Orients the bone by PCA and crops a cube: three suggested segments "
        "or a free position",
    "El hueso se orienta por PCA y se corta en tres tramos a lo largo de su "
    "eje mayor. <b>El sentido del eje es arbitrario</b>: los tramos van "
    "ordenados, pero cual es proximal lo decides tu. Pista: en estos "
    "sesamoideos la densidad crece de proximal a distal.":
        "The bone is oriented by PCA and cut into three segments along its "
        "major axis. <b>The axis direction is arbitrary</b>: the segments are "
        "ordered, but which one is proximal is your call. Hint: in these "
        "sesamoids density grows from proximal to distal.",
    "Lado del cubo (mm)": "Cube side (mm)",
    "Recalcular tramos": "Recompute segments",
    "El tramo 1 es": "Segment 1 is",
    "proximal": "proximal",
    "centro": "middle",
    "distal": "distal",
    "libre": "free",
    "Elegir libremente": "Choose freely",
    "usar esta posicion": "use this position",
    "posicion en el eje (%)": "position along the axis (%)",
    "desplazar eje 2 (mm)": "offset on axis 2 (mm)",
    "desplazar eje 3 (mm)": "offset on axis 3 (mm)",
    "BV/TV {b:.4f} · lado {lv} vox": "BV/TV {b:.4f} · side {lv} vox",
    "  ·  <b style='color:#c0392b;'>{p:.1f} % del cubo cae fuera del "
    "volumen</b>":
        "  ·  <b style='color:#c0392b;'>{p:.1f} % of the cube falls "
        "outside the volume</b>",
    "  ·  <b style='color:#c0392b;'>{p:.1f} % fuera</b>":
        "  ·  <b style='color:#c0392b;'>{p:.1f} % outside</b>",
    "{p:.1f} % del cubo cae fuera del volumen: el BV/TV esta subestimado.":
        "{p:.1f} % of the cube falls outside the volume: BV/TV is "
        "underestimated.",
    "Nada seleccionado": "Nothing selected",
    "Elige un tramo o el modo libre.": "Choose a segment or the free mode.",
    "Sin vista previa": "No preview yet",
    "Mueve algun mando para calcular el cubo libre.":
        "Move a control to compute the free cube.",
    "VOI recortado: {et}, {lv}³ vox, BV/TV {b:.4f}. El volumen completo "
    "ya no esta cargado.":
        "VOI cropped: {et}, {lv}³ vox, BV/TV {b:.4f}. The full volume is "
        "no longer loaded.",
    "3 clases: aire / medula / hueso": "3 classes: air / marrow / bone",
    "2 clases: Otsu clasico": "2 classes: classic Otsu",
    "En una reconstruccion cruda el campo es aire en su mayor parte y hay "
    "<b>tres</b> poblaciones: con 2 clases la medula acaba contada como "
    "hueso. Medido en H4 frente a la segmentacion manual (tercios 0.28 / "
    "0.55 / 0.77): con 2 clases salen 0.61 / 0.79 / 0.94, con 3 clases 0.35 / "
    "0.59 / 0.81. <b>Ningun automatico reproduce una segmentacion manual</b> "
    "— comprueba el valor.":
        "In a raw reconstruction the field is mostly air and there are "
        "<b>three</b> populations: with 2 classes the marrow ends up counted "
        "as bone. Measured on H4 against the manual segmentation (thirds "
        "0.28 / 0.55 / 0.77): 2 classes give 0.61 / 0.79 / 0.94, 3 classes "
        "0.35 / 0.59 / 0.81. <b>No automatic threshold reproduces a manual "
        "segmentation</b> — check the value.",
    "Umbral por {m} = {u:.1f}. La segmentacion es la mayor fuente de "
    "incertidumbre: comprueba el valor.":
        "Threshold by {m} = {u:.1f}. Segmentation is the largest source of "
        "uncertainty: check the value.",
    "Archivos apartados": "Files set aside",
    "No se pudo cargar la pila": "The stack could not be loaded",
    "{n} archivo(s) apartados por no ser rebanadas de la serie: {lista}":
        "{n} file(s) set aside for not being slices of the series: {lista}",
    "Se cargaron {ok} rebanadas de {tot}x{ancho} px.\n\nEstos {n} archivos se "
    "apartaron por tener otra forma — son las proyecciones, previsualizaciones "
    "y hojas que deja el escaner junto a la reconstruccion:\n\n{lista}":
        "{ok} slices of {tot}x{ancho} px were loaded.\n\nThese {n} files were "
        "set aside because their shape differs — they are the projections, "
        "previews and sheets the scanner leaves next to the "
        "reconstruction:\n\n{lista}",
    "Pila cargada: {n} rebanadas · {mm:.6f} mm/vox ({orig}) · umbral {u} [{met}]":
        "Stack loaded: {n} slices · {mm:.6f} mm/vox ({orig}) · "
        "threshold {u} [{met}]",
    "ninguno": "none",

    # -- parametros del spinodoide -----------------------------------------
    "Parametros del campo": "Field parameters",
    "Muestreo de ondas:": "Wave sampling:",
    "rechazo (GIBBON)": "rejection (GIBBON)",
    "equitativo (TPMS-Scaffolds)": "uniform (TPMS-Scaffolds)",
    "Con conos de angulos DESIGUALES los dos esquemas producen\n"
    "anisotropias distintas para los mismos thetas. No es un ajuste\n"
    "cosmetico: hay que declarar cual se uso.":
        "With UNEQUAL cone angles the two schemes give different\n"
        "anisotropies for the same thetas. This is not a cosmetic\n"
        "setting: the scheme used must be reported.",
    "Semilla:": "Seed:",
    "Otra": "New",
    "El generador es estocastico: la misma parametrizacion da\n"
    "realizaciones distintas. Cambiar la semilla muestra esa\n"
    "dispersion, que es el suelo de ruido de cualquier ajuste.":
        "The generator is stochastic: the same parameters give\n"
        "different realisations. Changing the seed shows that\n"
        "spread, which is the noise floor of any fit.",
    "Resolucion": "Resolution",
    "Actualizar vista automaticamente": "Update the view automatically",
    "Generar vista": "Generate view",

    # -- visualizacion -----------------------------------------------------
    "Camara y color": "Camera and colour",
    "Enlazar las camaras de los dos paneles": "Link both panel cameras",
    "Con encuadres distintos, comparar grosores a ojo entre los dos\n"
    "paneles no significa nada: una estructura puede parecer mas\n"
    "gruesa solo por estar mas cerca de la camara.":
        "With different framings, comparing thickness by eye between\n"
        "the two panels means nothing: a structure can look thicker\n"
        "merely by being closer to the camera.",
    "Mostrar ejes y cotas en mm": "Show axes and dimensions in mm",
    "Mostrar fabrica MIL (eje y elipsoide)": "Show MIL fabric",
    "Sombras de contacto (SSAO)": "Contact shadows (SSAO)",
    "Fondo blanco (para figuras)": "White background (for figures)",
    "En pantalla el degradado da a la silueta algo contra lo que\n"
    "recortarse y el volumen se lee mejor. Para una figura de tesis o\n"
    "una captura de un informe conviene el blanco puro.":
        "On screen the gradient gives the silhouette something to stand\n"
        "out against and the volume reads better. For a thesis figure or\n"
        "a report screenshot, pure white is the right choice.",
    "Colorear por:": "Colour by:",
    "Material (liso)": "Material (plain)",
    "Espesor local": "Local thickness",
    "Deformacion efectiva": "Effective strain",
    "Tension de von Mises": "von Mises stress",
    "Deformacion total (mm)": "Total deformation (mm)",
    "Histograma de espesor…": "Thickness histogram…",
    "Recorte (solo visual)": "Clipping (visual only)",
    "Eje:": "Axis:",
    "Fraccion recortada": "Clipped fraction",

    # -- morfometria -------------------------------------------------------
    "Superficie:": "Surface:",
    "voxeles (marching cubes)": "voxels (marching cubes)",
    "malla suavizada + PyMeshFix": "smoothed mesh + PyMeshFix",
    "Medir (morfometria completa)": "Measure (full morphometry)",
    "Dispersion…": "Spread…",
    "K:": "K:",
    "Anadir Conn.D y SMI a la medida": "Add Conn.D and SMI",
    "Metrica": "Metric",
    "VOI": "VOI",
    "Dif. relativa": "Relative diff.",
    "MIL acotado: estructura laminar, el DA es una COTA INFERIOR.":
        "MIL clamped: laminar structure, the DA is a LOWER BOUND.",
    "(sin calcular)": "(not computed)",

    # -- ajuste ------------------------------------------------------------
    "Ajustar al VOI": "Fit to the VOI",
    "Ajustar con TODOS los metodos…": "Fit with ALL methods…",
    "Desempate mecanico": "Mechanical tie-break",
    "Peso:": "Weight:",
    "Malla:": "Mesh:",

    # -- mecanica ----------------------------------------------------------
    "Tensor elastico": "Elastic tensor",
    "Resolucion del tensor:": "Tensor resolution:",
    "Ensayo de compresion": "Compression test",
    "Apoyo:": "Support:",
    "deslizante": "sliding",
    "empotrado": "fixed",
    "Resolucion:": "Resolution:",
    "Direccion:": "Direction:",
    "Z (axial)": "Z (axial)",
    "Los tres ejes": "All three axes",
    "Resolver y estimar el fallo": "Solve and estimate failure",
    "Distribuciones…": "Distributions…",
    "Convergencia…": "Convergence…",
    "Analisis comparado VOI / spinodoide…": "Compare VOI / spinodoid…",
    "Ajustar antes al VOI (elige el de menor error)":
        "Fit first (lowest error wins)",

    # -- exportacion -------------------------------------------------------
    "Hexaedrica (voxeles)": "Hexahedral (voxels)",
    "Tetraedrica TET10": "Tetrahedral TET10",
    "Formatos:": "Formats:",
    ".vtu": ".vtu",
    ".inp  (Abaqus)": ".inp  (Abaqus)",
    ".apdl  (ANSYS)": ".apdl  (ANSYS)",
    ".stl": ".stl",
    "Exportar solido…": "Export solid…",

    # -- lote --------------------------------------------------------------
    "Replicas por VOI:": "Replicates per VOI:",
    "Modo:": "Mode:",
    "rapido (A+B)": "fast (A+B)",
    "completo (A+B+C)": "full (A+B+C)",
    "Elegir VOIs y correr lote…": "Choose VOIs and run batch…",

    # -- rotulos de los paneles 3D ------------------------------------------
    "VOI real": "Real VOI",
    "Spinodoide (vista)": "Spinodoid (view)",

    # -- dialogos de figura -------------------------------------------------
    "Guardar figura…": "Save figure…",
    "Cerrar": "Close",
    "Guardar la figura": "Save the figure",
    "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)":
        "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)",

    # -- ayudas largas: son la documentacion de la interfaz ---------------
    'Abaqus / CalculiX, con material y los conjuntos de nodos BASE y\nTECHO. Es el archivo mas pesado con diferencia.':
        'Abaqus / CalculiX, with material and the BASE and TOP node\nsets. By far the heaviest file.',
    'Ajusta cada VOI, genera las replicas y descompone la varianza\nen entre-especimenes y dentro, con ICC, N efectivo y TOST.\n\nTARDA MUCHO: es un ajuste completo por VOI mas K generaciones.':
        'Fits every VOI, generates the replicates and decomposes the\nvariance into between- and within-specimen, with ICC, effective\nN and TOST.\n\nSLOW: one full fit per VOI plus K generations.',
    'Anade una etapa D: reordena los mejores candidatos incluyendo\nEz/Es y Ez/Ex, obtenidos por homogeneizacion.\n\nNO es una busqueda con el termino mecanico dentro: homogeneizar\nlas ~53 evaluaciones seria de tres a cuatro ordenes de magnitud\nmas caro. Solo se paga entre candidatos que ya son buenos\nmorfometricamente, asi que si el optimo mecanico esta en una\nregion que la busqueda descarto pronto, esto no lo encuentra.':
        'Adds a stage D: re-ranks the best candidates including Ez/Es\nand Ez/Ex, obtained by homogenisation.\n\nThis is NOT a search with the mechanical term inside:\nhomogenising all ~53 evaluations would be three to four orders\nof magnitude more expensive. It is paid only among candidates\nthat are already morphometrically good, so if the mechanical\noptimum lies in a region the search discarded early, this will\nnot find it.',
    'Aplica a las DOS estructuras el mismo ensayo del articulo de\nTapia, Gonzalez, Vidal y Salinas (Biology 2026;15:722):\ntejido lineal elastico isotropo de 18 GPa y nu = 0.30, carga\naxial de compresion de 100 N y apoyo empotrado en la cara\nopuesta.\n\nDevuelve las dos variables que ese trabajo reporta:\n  · distribucion espacial de la tension equivalente de von Mises\n  · distribucion espacial de la deformacion total (mm)\n\nLos dos paneles quedan coloreados con la MISMA escala, que es lo\nque permite compararlos a ojo. La malla aqui es hexaedrica de un\nvoxel, no TET10 de ANSYS: los patrones son comparables, las\ncifras absolutas no.':
        'Applies to BOTH structures the same test as the paper by\nTapia, Gonzalez, Vidal & Salinas (Biology 2026;15:722):\nlinear elastic isotropic tissue of 18 GPa and nu = 0.30, a\n100 N axial compressive load and a fixed support on the\nopposite face.\n\nIt returns the two variables that work reports:\n  · spatial distribution of von Mises equivalent stress\n  · spatial distribution of total deformation (mm)\n\nBoth panels are coloured on the SAME scale, which is what makes\nthem comparable by eye. The mesh here is one hexahedron per\nvoxel, not the TET10 of ANSYS: the patterns are comparable, the\nabsolute figures are not.',
    'Busqueda escalonada de fitSpinodoidToVOI: densidad y numero de\nonda, despues angulos conicos, despues refinado. Son 53\nevaluaciones y tarda varios minutos a la resolucion del VOI.':
        'Staged search from fitSpinodoidToVOI: density and wave number,\nthen cone angles, then refinement. 53 evaluations, several\nminutes at the VOI resolution.',
    'Con los tres ejes sale E_max/E_min, una anisotropia MECANICA\nindependiente del DA del tensor MIL, que es puramente geometrico,\ny contrastable con Ex/Ey/Ez de la homogeneizacion periodica.\n\nCuesta el triple de tiempo: son tres sistemas completos.':
        'All three axes give E_max/E_min, a MECHANICAL anisotropy that is\nindependent of the DA of the MIL tensor -which is purely\ngeometric- and can be checked against Ex/Ey/Ez from periodic\nhomogenisation.\n\nIt costs three times as long: three complete systems.',
    "Conn.D (Odgaard & Gundersen) mide las conexiones REDUNDANTES de\nla red: separa 'hueso mas fino' de 'hueso roto', que es algo que\nBV/TV y Tb.Th no distinguen.\n\nSMI (Hildebrand & Ruegsegger) mide si la estructura tiende a\nplaca o a barra. Leerlo con la salvedad de Salmon et al. 2015:\na BV/TV alto la concavidad lo hace negativo y pierde sentido.\n\nEstan fuera del bucle de ajuste a proposito: cuestan una\nmarching cubes de mas cada una y no entran en el error.":
        "Conn.D (Odgaard & Gundersen) measures the REDUNDANT connections\nof the network: it tells 'thinner bone' from 'broken bone',\nwhich BV/TV and Tb.Th cannot.\n\nSMI (Hildebrand & Ruegsegger) measures whether the structure\ntends to plate or to rod. Read it with the caveat of Salmon et\nal. 2015: at high BV/TV concavity turns it negative and it\nstops meaning anything.\n\nThey are deliberately outside the fitting loop: each costs an\nextra marching cubes and neither enters the error.",
    'Corre los tres metodos que comparten la misma definicion de\nerror —rapido, completo y equitativo—, se queda con el de menor\nerror y ensaya ESE. El desempate mecanico se excluye a proposito:\nsu error incluye dos terminos mas y no es comparable con el de\nlos otros.\n\nSin esto se ensaya el spinodoide que haya en pantalla, que puede\nno tener nada que ver con el VOI: a densidad igual, los angulos\nde cono cambian E_app en un factor de 15.\n\nTARDA: del orden de veinte minutos a la resolucion del VOI, mas\nel ensayo. Desmarcalo si ya has ajustado.':
        'Runs the three methods that share the same error definition\n-fast, full and uniform-, keeps the one with the lowest error\nand tests THAT one. The mechanical tie-break is deliberately\nexcluded: its error carries two extra terms and is not\ncomparable with the others.\n\nWithout this, whatever spinodoid is on screen gets tested, and\nit may have nothing to do with the VOI: at equal density, the\ncone angles change E_app by a factor of 15.\n\nSLOW: of the order of twenty minutes at the VOI resolution,\nplus the test. Untick it if you have already fitted.',
    'Criterio de Pistoia: la carga a la que el 2% del tejido supera el 0.7% de deformacion efectiva. Los dos valores son convenciones calibradas en radio distal humano, no constantes fisicas.\nEl ensayo deja ademas los campos de deformacion efectiva y de von Mises listos para colorear.':
        'Pistoia criterion: the load at which 2% of the tissue exceeds\n0.7% effective strain. Both values are conventions calibrated on\nthe human distal radius, not physical constants.\nThe test also leaves the effective strain and von Mises fields\nready for colouring.',
    'Dibuja el elipsoide de fabrica del tensor MIL y, en magenta, su\neje mayor: la direccion en la que la estructura es mas continua\ny, por tanto, mas rigida.\n\nEs el mismo dato que da el DA de la tabla, pero orientado. Dos\nestructuras pueden tener el MISMO DA y estar orientadas a 90\ngrados una de otra: eso la tabla no lo distingue y esto si.\n\nNecesita haber medido. El elipsoide del spinodoide corresponde a\nla ultima medida, no a los deslizadores: si cambias parametros,\ndesaparece hasta que vuelvas a medir.':
        'Draws the fabric ellipsoid of the MIL tensor and, in magenta,\nits major axis: the direction in which the structure is most\ncontinuous and therefore stiffest.\n\nIt is the same information as the DA in the table, but oriented.\nTwo structures can have the SAME DA and be 90 degrees apart:\nthe table cannot tell them apart, this can.\n\nIt needs a measurement. The spinodoid ellipsoid belongs to the\nLAST measurement, not to the sliders: if you change parameters\nit disappears until you measure again.',
    'Distribucion completa del espesor local, VOI y spinodoide\nsuperpuestos. Dos estructuras con el mismo Tb.Th pueden tener\ndistribuciones muy distintas, y eso el escalar no lo distingue.':
        'Full distribution of local thickness, VOI and spinodoid\noverlaid. Two structures with the same Tb.Th can have very\ndifferent distributions, and the scalar hides that.',
    "El coste crece con el CUBO del lado. MEDIDO aqui, generar y\nmedir tarda 4.7 s a 64 vox, 31 s a 128, 105 s a 192 y 263 s a\n256, con 1.5 GB de pico. Lo caro es generar, no medir.\n\nPor encima de ~96 vox conviene DESMARCAR 'Actualizar vista\nautomaticamente': si no, cada deslizador que muevas lanza un\ncalculo de minutos.\n\nLa resolucion de MEDIDA por defecto es 64 porque es la que usa\nla app de MATLAB para comparar. Cambiarla cambia la\ndiscretizacion del candidato, asi que hay que declararla.":
        "Cost grows with the CUBE of the side. MEASURED here, generating\nand measuring takes 4.7 s at 64 vox, 31 s at 128, 105 s at 192\nand 263 s at 256, peaking at 1.5 GB. Generating is the expensive\npart, not measuring.\n\nAbove ~96 vox it is worth UNTICKING 'Update the view\nautomatically': otherwise every slider you move launches a\ncomputation of minutes.\n\nThe default MEASUREMENT resolution is 64 because that is what\nthe MATLAB app uses for comparison. Changing it changes the\ndiscretisation of the candidate, so it must be reported.",
    'El espesor local muestra DONDE la estructura es gruesa o fina.\nLos dos paneles comparten escala de color, sin lo cual la\ncomparacion no significaria nada.\n\nAVISO: a 3-4 voxeles de grosor el metodo subestima un 30%, asi\nque los valores absolutos NO son citables como Tb.Th. Para eso\nesta el Tb.Th de Parfitt de la tabla.\n\nLa deformacion efectiva y von Mises salen del mismo ensayo de\ncompresion, pero NO son el mismo mapa reescalado: la primera\nviene de la energia total y la segunda solo de la parte\ndesviadora. Se separan en los nudos, donde el estado es triaxial.':
        'Local thickness shows WHERE the structure is thick or thin.\nBoth panels share the colour scale, without which the\ncomparison would mean nothing.\n\nWARNING: at 3-4 voxels thick the method underestimates by 30%,\nso the absolute values are NOT quotable as Tb.Th. The Parfitt\nTb.Th in the table is what you quote.\n\nEffective strain and von Mises come from the same compression\ntest, but they are NOT the same map rescaled: the first comes\nfrom the total energy and the second only from the deviatoric\npart. They separate at the nodes, where the state is triaxial.',
    'El numero de filas NO es el N del estudio. Las replicas son tecnicas: reducen la incertidumbre de cada espécimen pero no anaden grados de libertad para comparar entre animales.':
        'The number of rows is NOT the N of the study. Replicates are\ntechnical: they reduce the uncertainty of each specimen but add\nno degrees of freedom for comparing between animals.',
    'Ez/Es depende MUCHO de esta malla: sobre el mismo VOI salio\n0.108 a 12 voxeles y 0.064 a 16. No es citable como rigidez.\nSirve porque VOI y candidatos se miden con la MISMA malla y el\nsesgo se cancela al comparar (mismo argumento que C4).':
        'Ez/Es depends STRONGLY on this mesh: on the same VOI it came out\n0.108 at 12 voxels and 0.064 at 16. It is not quotable as a\nstiffness. It is useful because VOI and candidates are measured\non the SAME mesh and the bias cancels in the comparison (the\nsame argument as C4).',
    'Hexaedrica: un elemento por voxel. Exacta, instantanea y no puede\nfallar al mallar, pero el elemento es lineal y la superficie queda\nescalonada. Es el metodo estandar de micro-elementos finitos.\n\nTET10: superficie suavizada y tetraedralizada. Elemento cuadratico\ny superficie lisa, pero mas lenta y pierde algo de volumen.':
        'Hexahedral: one element per voxel. Exact, instantaneous and it\ncannot fail to mesh, but the element is linear and the surface\nis stair-stepped. It is the standard micro-finite-element\nmethod.\n\nTET10: smoothed and tetrahedralised surface. Quadratic element\nand smooth surface, but slower and it loses some volume.',
    'Histograma de la deformacion efectiva y de von Mises.\nPistoia no mira el maximo sino el percentil 98, que es un punto\nde la cola: para saber si es robusto hay que ver la cola entera.':
        'Histogram of effective strain and of von Mises.\nPistoia does not look at the maximum but at the 98th percentile,\na point on the tail: to know whether it is robust you have to\nsee the whole tail.',
    'Homogeneizacion periodica sobre la rejilla de voxeles.\nEl coste crece con el CUBO del lado, asi que la estructura se\nremuestrea antes al limite indicado.':
        'Periodic homogenisation on the voxel grid.\nCost grows with the CUBE of the side, so the structure is\nresampled first down to the limit given here.',
    'Idioma':
        'Language',
    'Idioma cambiado. Los resultados ya calculados se conservan.':
        'Language changed. Results already computed are kept.',
    'La vista se recalcula al mover un deslizador; la medida solo al pulsar Medir. Misma semilla, asi que la vista previa es una version gruesa de lo que se medira.\nEl coste crece con el CUBO del lado: medido aqui, 64³ tarda 5 s, 128³ 31 s y 256³ 263 s con 1.5 GB. Por encima de ~96 vox conviene quitar la actualizacion automatica.':
        'The view is recomputed when you move a slider; the measurement only when you press Measure. Same seed, so the preview is a coarse version of what will be measured.\nCost grows with the CUBE of the side: measured here, 64³ takes 5 s, 128³ 31 s and 256³ 263 s with 1.5 GB. Above ~96 vox it is worth turning off the automatic update.',
    'Lanza los cuatro metodos disponibles a la vez, con una barra de\nprogreso por metodo, y al terminar los compara lado a lado con\nsus metricas frente al VOI para que elijas con cual seguir.\n\nDesde ahi se pueden generar ademas N replicas del elegido,\nidenticas o con variacion en las metricas que marques.':
        'Runs the four available methods at once, with one progress bar\nper method, and compares them side by side at the end with\ntheir metrics against the VOI so you can choose which to keep.\n\nFrom there you can also generate N replicates of the chosen one,\nidentical or varying in the metrics you tick.',
    'Malla volumetrica para ParaView o pyvista. Lleva la calidad y el\nvolumen de cada elemento, para localizar los degenerados antes de\nmandar la malla a un solver.':
        'Volumetric mesh for ParaView or pyvista. It carries the quality\nand volume of each element, so degenerate ones can be found\nbefore the mesh goes to a solver.',
    'Oclusion ambiental: oscurece los huecos donde una trabecula pasa\npor detras de otra. Sin ella la marana se lee plana y no se sabe\nque esta delante.\n\nEs un efecto de RENDER. No toca ninguna medida. Desmarcala si el\nequipo va lento o si el driver de video la dibuja mal.':
        'Ambient occlusion: darkens the gaps where one trabecula passes\nbehind another. Without it the tangle reads flat and you cannot\ntell what is in front.\n\nIt is a RENDERING effect. It touches no measurement. Untick it\nif the machine is slow or the video driver draws it badly.',
    'Repite la MISMA parametrizacion con K semillas distintas y da\nmedia +- sd de cada metrica. Es el suelo de ruido por debajo del\ncual ninguna diferencia significa nada.\n\nCon un VOI cargado anade la z: cuantas sd separan al VOI de la\nmedia de las realizaciones. |z| < 2 es indistinguible del ruido.':
        'Repeats the SAME parameters with K different seeds and gives\nmean +- sd of each metric. It is the noise floor below which no\ndifference means anything.\n\nWith a VOI loaded it adds the z: how many sd separate the VOI\nfrom the mean of the realisations. |z| < 2 is indistinguishable\nfrom noise.',
    'Resuelve a varias resoluciones y traza E_app frente al tamano\nde elemento. Es lo que convierte un E_app en citable o no.\n\nTarda: son varios ensayos completos, el mas fino el mas caro.':
        'Solves at several resolutions and plots E_app against element\nsize. It is what makes an E_app quotable or not.\n\nSlow: several complete tests, the finest the most expensive.',
    'Se exporta a la resolucion de MEDIDA, no a la de vista, regenerando con la misma semilla. El .stl solo existe por la via TET10: es la superficie suavizada y reparada, cerrada y por tanto imprimible. La hexaedrica no la produce porque seria la piel escalonada de los voxeles.':
        'Export uses the MEASUREMENT resolution, not the view one, regenerating with the same seed. The .stl exists only through the TET10 route: it is the smoothed and repaired surface, closed and therefore printable. The hexahedral mesh does not produce it because it would be the stair-stepped skin of the voxels.',
    'Sesion':
        'Session',
    'Solo disponible con la malla TET10. La hexaedrica solo puede dar la piel escalonada de los voxeles, que no interesa imprimir.':
        'Only available with the TET10 mesh. The hexahedral one can only give the stair-stepped skin of the voxels, which is not worth printing.',
    'X':
        'X',
    'Y':
        'Y',
    'Z':
        'Z',
    'deslizante: uz=0 en la base, sin friccion. Estado uniaxial.\nempotrado: base totalmente fija, como un ensayo con friccion\nentre probeta y platos. Coacciona lateralmente y rigidiza.\n\nLa eleccion cambia E_app y la carga de fallo: hay que declararla.':
        'sliding: uz=0 at the base, frictionless. Uniaxial state.\nfixed: base fully clamped, like a test with friction between\nspecimen and platens. It constrains laterally and stiffens.\n\nThe choice changes E_app and the failure load: it must be\nreported.',
    'voxeles: marching cubes crudo sobre la mascara. Es el modo de la\napp de MATLAB y el que usa el ajuste. Sobreestima el area un\n~8.5 % (medido sobre una esfera), igual en el VOI y en el\ncandidato, asi que el sesgo se cancela al comparar.\n\nmalla: superficie cerrada, suavizada (Taubin) y reparada con\nPyMeshFix, recortada medio voxel por dentro (caja de centros de\nvoxel, la misma de la superficie abierta) y sin las seis tapas.\nBaja ese sesgo a <1 %. El BV/TV sale ~3-5 % menor que el de\nvoxeles: es la diferencia entre la superficie a nivel 0.5 y el\nconteo de cubos enteros, no una perdida. DA y fraccion portante\nsiguen saliendo de los voxeles.\n\nVOI y candidato se miden SIEMPRE con el mismo modo. Cambiarlo\nvuelve a medir el VOI. Cada tabla debe declarar cual se uso.':
        'voxels: raw marching cubes on the mask. It is the mode of the\nMATLAB app and the one the fit uses. It overestimates the area\nby ~8.5 % (measured on a sphere), equally in the VOI and in the\ncandidate, so the bias cancels in the comparison.\n\nmesh: closed surface, smoothed (Taubin) and repaired with\nPyMeshFix, clipped half a voxel inwards (the voxel-centre box,\nthe same one as the open surface) and without the six caps.\nThat brings the bias below 1 %. BV/TV comes out ~3-5 % lower\nthan the voxel one: that is the difference between the surface\nat level 0.5 and counting whole cubes, not a loss. DA and the\nload-bearing fraction still come from the voxels.\n\nVOI and candidate are ALWAYS measured in the same mode.\nChanging it re-measures the VOI. Every table must state which\none was used.',

    # -- deslizadores: solo el texto base; el valor y la unidad los pone
    # `Deslizador._pinta`. Theta y Rotacion llevan el eje pegado porque en la
    # interfaz son nueve controles distintos, no uno con un selector.
    "Densidad relativa": "Relative density",
    "Numero de onda": "Wave number",
    "Numero de ondas": "Number of waves",
    "Theta X": "Theta X",
    "Theta Y": "Theta Y",
    "Theta Z": "Theta Z",
    "Rotacion X": "Rotation X",
    "Rotacion Y": "Rotation Y",
    "Rotacion Z": "Rotation Z",
    "Vista (interactiva)": "View (interactive)",
    "Medida": "Measurement",

    # -- mensajes de marcha: barra de estado, dialogos y ficheros -------
    # Los marcadores {} tienen que sobrevivir a la traduccion; el orden si
    # puede cambiar, que es justamente para lo que sirven los nombres.
    ', tapas excluidas {tap} mm2, perdida vs superficie cruda {per} %':
        ', caps excluded {tap} mm2, loss vs raw surface {per} %',
    '<p>{nv} VOI(s) × {nr} replicas, modo <b>{modo}</b>.</p><p>Estimacion muy gruesa: <b>{a}–{b} minutos</b>, segun el tamano de los VOIs. La ventana queda inutilizable mientras tanto.</p><p>¿Seguir?</p>':
        '<p>{nv} VOI(s) × {nr} replicates, mode <b>{modo}</b>.</p><p>Very rough estimate: <b>{a}–{b} minutes</b>, depending on the size of the VOIs. The window is unusable meanwhile.</p><p>Continue?</p>',
    'A {res}³ y densidad {dens} % salen del orden de <b>{n} elementos</b>.<br><br>El .inp de Abaqus ronda los 200 bytes por elemento, asi que serian varios GB, y la via TET10 puede no terminar.<br><br>¿Exportar de todas formas?':
        'At {res}³ and density {dens} % this comes to about <b>{n} elements</b>.<br><br>The Abaqus .inp runs about 200 bytes per element, so it would be several GB, and the TET10 route may never finish.<br><br>Export anyway?',
    'Abrir VOI':
        'Open VOI',
    'Ajustando — etapa {etapa}: {i} de {n}':
        'Fitting — stage {etapa}: {i} of {n}',
    'Ajuste comparado cancelado.':
        'Comparative fit cancelled.',
    'Ajuste terminado: error={err}, {n} evaluaciones en {t} s':
        'Fit finished: error={err}, {n} evaluations in {t} s',
    'Analisis comparado listo en {t} s — von Mises y deformacion total disponibles en «Colorear por»':
        'Comparative analysis ready in {t} s — von Mises and total deformation available under «Colour by»',
    'Analisis mecanico comparado':
        'Comparative mechanical analysis',
    'Camaras enlazadas: las dos vistas comparten encuadre.':
        'Cameras linked: both views share the framing.',
    'Camaras independientes.':
        'Cameras independent.',
    'Carga de fallo aparente — {hay}   [{t} s; los campos que se pintan son los del eje {eje}: cada eje es un caso de carga distinto]':
        'Apparent failure load — {hay}   [{t} s; the fields drawn are those of axis {eje}: each axis is a different load case]',
    'Carga un VOI o genera un spinodoide.':
        'Load a VOI or generate a spinodoid.',
    'Cargar una sesion':
        'Load a session',
    'Coloreado por deformacion efectiva [escala comun {a}–{b}; el criterio de fallo usa 7.0e-03]':
        'Coloured by effective strain [common scale {a}–{b}; the failure criterion uses 7.0e-03]',
    'Coloreado por deformacion total [escala comun {a}–{b} mm a la carga de referencia; el maximo esta en la cara cargada porque acumula todo el desplazamiento]':
        'Coloured by total deformation [common scale {a}–{b} mm at the reference load; the maximum is on the loaded face because it accumulates the whole displacement]',
    'Coloreado por material.':
        'Coloured by material.',
    'Coloreado por von Mises [escala comun {a}–{b} MPa a la carga de referencia; escala con la carga porque el problema es lineal]':
        'Coloured by von Mises [common scale {a}–{b} MPa at the reference load; it scales with the load because the problem is linear]',
    'Con {n} voxeles no salen tres mallas distintas por encima de 12. Sube la resolucion del ensayo.':
        'With {n} voxels there are no three distinct meshes above 12. Raise the test resolution.',
    'Confirmar el lote':
        'Confirm the batch',
    'Convergencia en {t} s — {n} de {tot} mallas resueltas.':
        'Convergence in {t} s — {n} of {tot} meshes solved.',
    'Desempate mecanico':
        'Mechanical tie-break',
    'El ajuste necesita un VOI de referencia. Cargalo primero.':
        'The fit needs a reference VOI. Load one first.',
    'El analisis comparado necesita las DOS: carga un VOI y genera el spinodoide.':
        'The comparative analysis needs BOTH: load a VOI and generate the spinodoid.',
    'El analisis no se resolvio':
        'The analysis did not converge',
    "El archivo no lleva la marca 'spinpy.visor/sesion'.":
        "The file does not carry the 'spinpy.visor/sesion' marker.",
    'Error en el calculo':
        'Error during the computation',
    'Espesor local — {partes}   [escala comun {a}–{b} mm; sesgo ~30% a 3 voxeles: usar como medida RELATIVA]':
        'Local thickness — {partes}   [common scale {a}–{b} mm; ~30% bias at 3 voxels: use as a RELATIVE measure]',
    'Exportación completada':
        'Export completed',
    'Exportado {tipo}: {n} elementos, {pct}% del volumen, en {t} s':
        'Exported {tipo}: {n} elements, {pct}% of the volume, in {t} s',
    'Exportando — {etapa}':
        'Exporting — {etapa}',
    'Exportar el spinodoide como solido':
        'Export the spinodoid as a solid',
    'Exportar los resultados':
        'Export the results',
    'Falta el VOI':
        'No VOI loaded',
    'Falta el ensayo':
        'No test run yet',
    'Faltan estructuras':
        'Structures missing',
    'Figura guardada en {f}':
        'Figure saved to {f}',
    'Genera primero un spinodoide.':
        'Generate a spinodoid first.',
    'Guardar la sesion':
        'Save the session',
    'Homogeneizacion terminada en {t} s':
        'Homogenisation finished in {t} s',
    'JSON (*.json)':
        'JSON (*.json)',
    'JSON (*.json);;Todos (*)':
        'JSON (*.json);;All (*)',
    'La malla va a ser muy grande':
        'The mesh is going to be very large',
    'Lote terminado en {t} min — {n} filas, pero N efectivo = {nef}. CSV en resultados/lote_*_{sello}.csv':
        'Batch finished in {t} min — {n} rows, but effective N = {nef}. CSV in resultados/lote_*_{sello}.csv',
    'Marca al menos un formato que exportar.':
        'Tick at least one format to export.',
    'Medido en {t} s   [superficie: {modo}{extra}]   [{esq}, semilla {sem}]':
        'Measured in {t} s   [surface: {modo}{extra}]   [{esq}, seed {sem}]',
    'Mide, homogeneiza o ensaya algo primero.':
        'Measure, homogenise or test something first.',
    'Midiendo el VOI antes de ajustar…':
        'Measuring the VOI before fitting…',
    'Misma semilla: la estructura es la misma realizacion.':
        'Same seed: the structure is the same realisation.',
    'Nada que exportar':
        'Nothing to export',
    'Ningun formato marcado':
        'No format ticked',
    'No es una sesion':
        'Not a session file',
    'No se pudo leer':
        'Could not be read',
    'No se pudo leer el VOI':
        'The VOI could not be read',
    'Nombre base (sin extension) (*)':
        'Base name (no extension) (*)',
    'Pulsa antes «Resolver y estimar el fallo»: las distribuciones salen de los campos de ese calculo.':
        'Press «Solve and estimate failure» first: the distributions come from the fields of that computation.',
    'Pulsa antes «Resolver y estimar el fallo»: los campos de deformacion efectiva y de von Mises salen de ese calculo.':
        'Press «Solve and estimate failure» first: the effective strain and von Mises fields come from that computation.',
    'Replicas generadas':
        'Replicates generated',
    'Resistencia estimada':
        'Estimated strength',
    'Resolucion insuficiente':
        'Resolution too low',
    'Resultados en {csv} y {json} ({n} magnitudes; el JSON lleva ademas el tensor C completo y la semilla)':
        'Results in {csv} and {json} ({n} quantities; the JSON also carries the full C tensor and the seed)',
    'Sesion cargada a medias':
        'Session only partly loaded',
    'Sesion de {fecha} cargada.':
        'Session from {fecha} loaded.',
    'Sesion guardada en {f} — regenera la MISMA realizacion, porque la semilla viaja dentro.':
        'Session saved to {f} — it regenerates the SAME realisation, because the seed travels inside.',
    'Sin estructura':
        'No structure',
    'Tensor elastico homogeneizado':
        'Homogenised elastic tensor',
    'VOI (*.vtk *.mat);;Todos (*)':
        'VOI (*.vtk *.mat);;All (*)',
    'VOI (*.vtk *.mat);;VTK legacy (*.vtk);;MATLAB (*.mat);;Todos (*)':
        'VOI (*.vtk *.mat);;VTK legacy (*.vtk);;MATLAB (*.mat);;All (*)',
    'VOI cargado. Pulsa Medir para la morfometria (el candidato se escalara a {lado} mm de lado).':
        'VOI loaded. Press Measure for the morphometry (the candidate will be scaled to {lado} mm a side).',
    'VOIs del lote (se pueden elegir varios)':
        'Batch VOIs (several allowed)',
    'Vista a {n}^3 en {t} s — densidad pedida {ped}%, obtenida {obt}%   [{esq}]':
        'View at {n}^3 in {t} s — density requested {ped}%, obtained {obt}%   [{esq}]',
    '{n} replicas en {archivo} — recuerda que no son {n} especimenes':
        '{n} replicates in {archivo} — remember they are not {n} specimens',

    # -- informes HTML de las ventanas de resultados --------------------
    # Las etiquetas de simbolo (E_app, sigma, epsilon, BV/TV) se dejan
    # como estan: son notacion, no idioma.
    ' · base {b} / techo {t} nodos':
        ' · base {b} / top {t} nodes',
    '&sigma;<sub>fallo</sub> es tension APARENTE sobre la seccion bruta; &sigma;<sub>vM</sub> es la del TEJIDO, y compararla con el limite elastico del hueso mineralizado (~150–200 MPa) dice si el criterio de Pistoia esta prediciendo algo fisicamente coherente. Apoyo: {apoyo}.<br>El 2% del tejido y el 0.7% de deformacion son convenciones calibradas en radio distal humano, no constantes fisicas.<br><b>E<sub>app</sub> no converge con la malla</b> en estructuras poco densas: a BV/TV 0.28 se midio 700, 797 y 365 MPa a 32³, 48³ y 64³. Haz un estudio de convergencia antes de citar un valor.<br>Resuelto en {t} s.':
        '&sigma;<sub>fallo</sub> is APPARENT stress on the gross section; &sigma;<sub>vM</sub> is that of the TISSUE, and comparing it with the yield strength of mineralised bone (~150–200 MPa) tells whether the Pistoia criterion is predicting something physically coherent. Support: {apoyo}.<br>The 2% of tissue and the 0.7% strain are conventions calibrated on the human distal radius, not physical constants.<br><b>E<sub>app</sub> does not converge with the mesh</b> in low-density structures: at BV/TV 0.28 it measured 700, 797 and 365 MPa at 32³, 48³ and 64³. Run a convergence study before quoting a value.<br>Solved in {t} s.',
    '&sigma;<sub>vM</sub> max al fallo':
        '&sigma;<sub>vM</sub> max at failure',
    '<b>Analisis mecanico comparado</b><br>':
        '<b>Comparative mechanical analysis</b><br>',
    '<b>Aviso:</b> el mallado perdió más del 5% del volumen de la superficie cruda. Con trabéculas de uno o dos vóxeles el suavizado y la reparación se comen material: sube la resolución de medida, o usa la malla hexaédrica, que es exacta.':
        '<b>Warning:</b> meshing lost more than 5% of the volume of the raw surface. With trabeculae one or two voxels thick, smoothing and repair eat material: raise the measurement resolution, or use the hexahedral mesh, which is exact.',
    '<b>Aviso:</b> el spinodoide actual no procede de un ajuste al VOI. Aqui eso pesa mas que en la morfometria: a densidad IGUAL, mover solo los angulos de cono cambia E<sub>app</sub> de 93 a 1390 MPa —un factor de 15— porque el cono ancho marca el eje BLANDO. Ajusta al VOI antes de comparar, o estaras midiendo la orientacion que dejaron los deslizadores.':
        '<b>Warning:</b> the current spinodoid does not come from a fit to the VOI. That matters more here than in the morphometry: at EQUAL density, moving only the cone angles changes E<sub>app</sub> from 93 to 1390 MPa —a factor of 15— because the wide cone marks the SOFT axis. Fit to the VOI before comparing, or you will be measuring the orientation the sliders happened to leave.',
    '<b>Aviso:</b> la deformacion aparente llega al {e} %. El analisis es LINEAL, asi que los campos escalan con la carga y siguen siendo correctos, pero el hueso real habria fallado mucho antes: los 100 N son un estimulo estandarizado de comparacion, no una condicion fisiologica. Para leer la carga de rotura esta el criterio de Pistoia.':
        '<b>Warning:</b> the apparent strain reaches {e} %. The analysis is LINEAR, so the fields scale with the load and remain correct, but real bone would have failed long before: the 100 N are a standardised comparison stimulus, not a physiological condition. For the failure load there is the Pistoia criterion.',
    '<b>Ensayo de compresion</b> — criterio de Pistoia<br>':
        '<b>Compression test</b> — Pistoia criterion<br>',
    '<b>Etapa D — desempate mecanico</b><br>':
        '<b>Stage D — mechanical tie-break</b><br>',
    '<b>La comparacion mecanica no es interpretable tal cual:</b> los ejes rigidos estan a {ang}° uno de otro. La carga va en Z, asi que cada estructura responde por una direccion distinta de su propia fabrica y la diferencia de E<sub>app</sub> mide sobre todo ese desalineamiento, no la microarquitectura.<br>El DA coincide porque es un ESCALAR: no distingue «rigido en Z» de «rigido en Y», y el error de ajuste solo mira ese escalar.':
        '<b>The mechanical comparison is not interpretable as it stands:</b> the stiff axes are {ang}° apart. The load runs along Z, so each structure responds through a different direction of its own fabric and the difference in E<sub>app</sub> measures mostly that misalignment, not the microarchitecture.<br>The DA agrees because it is a SCALAR: it does not tell «stiff in Z» from «stiff in Y», and the fitting error looks only at that scalar.',
    '<b>Malla exportada</b> — sistema de unidades mm · N · MPa<br>':
        '<b>Mesh exported</b> — unit system mm · N · MPa<br>',
    '<b>Malla:</b> aqui es hexaedrica de un voxel; en el articulo son tetraedros SOLID187 de 0.05 mm generados en ANSYS. Los patrones espaciales y el orden de magnitud son comparables; las cifras absolutas no lo son.':
        '<b>Mesh:</b> here it is one hexahedron per voxel; in the paper they are SOLID187 tetrahedra of 0.05 mm generated in ANSYS. The spatial patterns and the order of magnitude are comparable; the absolute figures are not.',
    '<b>Modulos efectivos</b> (material base E<sub>s</sub> = 20 GPa, &nu;<sub>s</sub> = 0.30)<br>':
        '<b>Effective moduli</b> (base material E<sub>s</sub> = 20 GPa, &nu;<sub>s</sub> = 0.30)<br>',
    '<b>Superficie NO estanca</b>: {nb} bordes abiertos, {nc} componente(s). El .stl no es imprimible tal cual; sube la resolución de medida o usa la malla hexaédrica.':
        '<b>Surface NOT watertight</b>: {nb} open edges, {nc} component(s). The .stl is not printable as it stands; raise the measurement resolution or use the hexahedral mesh.',
    '<b>Superficie estanca</b>: 0 bordes abiertos, {nc} componente(s). Apta para imprimir y para tetraedralizar.':
        '<b>Watertight surface</b>: 0 open edges, {nc} component(s). Fit to print and to tetrahedralise.',
    '<p><b>Ajuste previo</b> — se ensaya el de menor error:</p>':
        '<p><b>Previous fit</b> — the lowest-error one is tested:</p>',
    '<p><b>Orientacion del eje principal</b> (MIL, medido sobre la malla que se ensaya): VOI a {av}° del eje Z (DA {dav}), spinodoide a {as_}° (DA {das}). <b>Angulo entre ambos: {ang}°.</b></p>':
        '<p><b>Orientation of the principal axis</b> (MIL, measured on the mesh that is tested): VOI at {av}° from the Z axis (DA {dav}), spinodoid at {as_}° (DA {das}). <b>Angle between them: {ang}°.</b></p>',
    '<p>Frente a la superficie cerrada cruda: <b>{pct}%</b> del volumen (la superficie cruda es ya un {dif}% respecto al conteo de vóxeles, por definición).</p>':
        '<p>Against the raw closed surface: <b>{pct}%</b> of the volume (the raw surface is already {dif}% relative to the voxel count, by definition).</p>',
    '<p>Tension aparente sobre la seccion bruta: <b>{sig} MPa</b> ({carga} N). Resolucion {res}³, {n} elementos.</p>':
        '<p>Apparent stress on the gross section: <b>{sig} MPa</b> ({carga} N). Resolution {res}³, {n} elements.</p>',
    '<p>VOI: E<sub>z</sub>/E<sub>s</sub> = <b>{ezes}</b> · E<sub>z</sub>/E<sub>x</sub> = <b>{ezex}</b> (malla {res}³, peso {peso})</p>':
        '<p>VOI: E<sub>z</sub>/E<sub>s</sub> = <b>{ezes}</b> · E<sub>z</sub>/E<sub>x</sub> = <b>{ezex}</b> (mesh {res}³, weight {peso})</p>',
    '<p>{ne} elementos de tipo <b>{tipo}</b>, {nn} nodos ({ng} grados de libertad)<br>':
        '<p>{ne} elements of type <b>{tipo}</b>, {nn} nodes ({ng} degrees of freedom)<br>',
    '<tr><th></th><th>BV/TV</th><th>von Mises media<br>[MPa]</th><th>von Mises p99 superficie<br>[MPa]</th><th>von Mises p99 global<br>[MPa]</th><th>von Mises max<br>[MPa]</th><th>Deform. total max<br>[mm]</th><th>&epsilon; aparente<br>[%]</th><th>E<sub>app</sub><br>[MPa]</th></tr>':
        '<tr><th></th><th>BV/TV</th><th>von Mises mean<br>[MPa]</th><th>von Mises p99 surface<br>[MPa]</th><th>von Mises p99 global<br>[MPa]</th><th>von Mises max<br>[MPa]</th><th>Total deform. max<br>[mm]</th><th>apparent &epsilon;<br>[%]</th><th>E<sub>app</sub><br>[MPa]</th></tr>',
    '&sigma;<sub>vM</sub> p99 superficie al fallo':
        '&sigma;<sub>vM</sub> p99 surface at failure',
    '<tr><th></th><th>rejilla</th><th>E<sub>x</sub></th><th>E<sub>y</sub></th><th>E<sub>z</sub></th><th>E<sub>z</sub>/E<sub>s</sub></th><th>E<sub>z</sub>/E<sub>x</sub></th></tr>':
        '<tr><th></th><th>grid</th><th>E<sub>x</sub></th><th>E<sub>y</sub></th><th>E<sub>z</sub></th><th>E<sub>z</sub>/E<sub>s</sub></th><th>E<sub>z</sub>/E<sub>x</sub></th></tr>',
    '<tr><th>metodo</th><th>error</th><th>tiempo</th></tr>':
        '<tr><th>method</th><th>error</th><th>time</th></tr>',
    '<tr><th>orden morfo.</th><th>error morfo.</th><th>error + mecanica</th><th>E<sub>z</sub>/E<sub>s</sub></th><th>desvio vs VOI</th></tr>':
        '<tr><th>morpho. rank</th><th>morpho. error</th><th>error + mechanics</th><th>E<sub>z</sub>/E<sub>s</sub></th><th>deviation vs VOI</th></tr>',
    'Anisotropia mecanica':
        'Mechanical anisotropy',
    'Carpeta: {d}':
        'Folder: {d}',
    'Con una sola estructura resuelta <b>no hay comparacion</b>. Si el solver no convergio, la causa habitual es una resolucion demasiado baja para el grosor trabecular: a 20³ esta estructura no converge y a 24³ si. Sube la resolucion del ensayo y repite.':
        'With a single structure solved <b>there is no comparison</b>. If the solver did not converge, the usual cause is a resolution too low for the trabecular thickness: at 20³ this structure does not converge and at 24³ it does. Raise the test resolution and repeat.',
    'Confirmo':
        'Confirmed',
    'E<sub>z</sub>/E<sub>s</sub> depende mucho de la malla de homogeneizacion y <b>no es citable como rigidez</b>. Vale aqui porque VOI y candidatos se miden con la misma malla y el sesgo se cancela al comparar.<br>Solo se homogeneizaron los finalistas: si el optimo mecanico estaba en una region que la busqueda morfometrica descarto pronto, esta etapa no lo encuentra.':
        'E<sub>z</sub>/E<sub>s</sub> depends strongly on the homogenisation mesh and <b>is not quotable as a stiffness</b>. It is valid here because VOI and candidates are measured on the same mesh and the bias cancels in the comparison.<br>Only the finalists were homogenised: if the mechanical optimum lay in a region the morphometric search discarded early, this stage will not find it.',
    'E<sub>z</sub>/E<sub>s</sub> es la rigidez axial normalizada y E<sub>z</sub>/E<sub>x</sub> la anisotropia elastica. Son los dos terminos mecanicos que la app puede sumar al error de ajuste.<br>La estructura se remuestreo antes de homogeneizar: el coste crece con el cubo del lado.':
        'E<sub>z</sub>/E<sub>s</sub> is the normalised axial stiffness and E<sub>z</sub>/E<sub>x</sub> the elastic anisotropy. They are the two mechanical terms the app can add to the fitting error.<br>The structure was resampled before homogenising: cost grows with the cube of the side.',
    'El .inp lleva los conjuntos de nodos BASE y TECHO para aplicar el apoyo y la carga con un clic.':
        'The .inp carries the BASE and TOP node sets so the support and the load can be applied in one click.',
    'El desempate mecanico no compite aqui: su error incluye Ez/Es y Ez/Ex, dos terminos que los otros no tienen, asi que su cifra no es comparable con estas.':
        'The mechanical tie-break does not compete here: its error includes Ez/Es and Ez/Ex, two terms the others do not have, so its figure is not comparable with these.',
    'El ensayo no se resolvio.':
        'The test did not solve.',
    'Es una anisotropia MECANICA y no tiene por que coincidir con el DA del tensor MIL, que es puramente geometrico. Contrastarla con E<sub>x</sub>/E<sub>y</sub>/E<sub>z</sub> del boton «Tensor elastico»: alli el contorno es PERIODICO y aqui son platos reales sobre una probeta finita, asi que solo convergen si la estructura es grande frente a la trabecula. Que difieran es efecto de tamano o de borde, no un fallo.':
        'This is a MECHANICAL anisotropy and need not agree with the DA of the MIL tensor, which is purely geometric. Check it against E<sub>x</sub>/E<sub>y</sub>/E<sub>z</sub> from the «Elastic tensor» button: there the boundary is PERIODIC and here they are real platens on a finite specimen, so they only converge if the structure is large compared with the trabecula. A difference is a size or boundary effect, not a bug.',
    'Generado en {t} s.':
        'Generated in {t} s.',
    'La <b>von Mises</b> localiza donde se concentra la tension en el tejido; la <b>deformacion total</b> es el modulo del desplazamiento y mide la rigidez global, por lo que su maximo esta siempre en la cara cargada. No son el mismo mapa.':
        'The <b>von Mises</b> stress locates where stress concentrates in the tissue; the <b>total deformation</b> is the magnitude of the displacement and measures the global stiffness, so its maximum is always on the loaded face. They are not the same map.',
    'No se resolvio':
        'Did not solve',
    'Reordeno':
        'Re-ranked',
    'Resuelto en {t} s ({solver}, residuo {res}). Los dos campos quedan en «Colorear por» con la MISMA escala en los dos paneles.':
        'Solved in {t} s ({solver}, residual {res}). Both fields are left under «Colour by» with the SAME scale in both panels.',
    'Spinodoide ajustado ({m})':
        'Fitted spinodoid ({m})',
    'Volumen mallado {v} mm³ — <b>{pct}%</b> del volumen en vóxeles</p>':
        'Meshed volume {v} mm³ — <b>{pct}%</b> of the voxel volume</p>',
    'eje':
        'axis',
    'el mejor por morfometria NO es el mejor al sumar la mecanica. Sin esta etapa el ajuste habria elegido una estructura con la rigidez equivocada.':
        'the best by morphometry is NOT the best once the mechanics is added. Without this stage the fit would have chosen a structure with the wrong stiffness.',
    'el mejor por morfometria lo sigue siendo. La etapa no decidio nada, que tambien es informacion: la morfometria bastaba.':
        'the best by morphometry still is. The stage decided nothing, which is also information: morphometry was enough.',
    'elem.':
        'elem.',
    'ni':
        'nor',
    'no homogeneizo':
        'did not homogenise',
    'spin vs VOI':
        'spin vs VOI',
    '{etq}: E<sub>max</sub>/E<sub>min</sub> = {r} (mas rigido en <b>{eje}</b>)':
        '{etq}: E<sub>max</sub>/E<sub>min</sub> = {r} (stiffer in <b>{eje}</b>)',
    '{n} triángulos':
        '{n} triangles',

    # -- dialogos de figura, barras de progreso y avisos sueltos --------
    ' (con desempate mecanico a {r}³)':
        ' (with mechanical tie-break at {r}³)',
    ' (y el VOI)…':
        ' (and the VOI)…',
    '(DA2 {da2}: rigido en un PLANO, no en un eje)':
        '(DA2 {da2}: stiff in a PLANE, not along an axis)',
    '<b>Estas {n} filas no son {n} especimenes.</b> Son replicas de UN ajuste: no anaden grados de libertad a ninguna comparacion biologica. Para eso esta el N efectivo del panel de lote.':
        '<b>These {n} rows are not {n} specimens.</b> They are replicates of ONE fit: they add no degrees of freedom to any biological comparison. That is what the effective N of the batch panel is for.',
    '<b>Las replicas son tecnicas.</b> Reducen la incertidumbre de cada espécimen dividiendo la varianza DENTRO por K, y no tocan la varianza ENTRE, que es la biologica. Los grados de libertad para una afirmacion sobre la especie los ponen los animales, no las realizaciones.<br>La equivalencia es <b>TOST</b>, no una t de Student: «no se rechazo la nula» no demuestra que dos cosas sean iguales — con N pequeno nunca se rechaza y con N grande siempre. TOST declara un margen de antemano y demuestra que la diferencia cae dentro. El margen por defecto es el CV DENTRO observado, es decir el suelo de ruido del propio generador: pedir menos seria pedir que el sintetico se parezca al real mas de lo que el real se parece a si mismo.':
        '<b>Replicates are technical.</b> They reduce the uncertainty of each specimen by dividing the WITHIN variance by K, and leave the BETWEEN variance, the biological one, untouched. The degrees of freedom for a statement about the species come from the animals, not from the realisations.<br>Equivalence is <b>TOST</b>, not a Student t: «the null was not rejected» does not show that two things are equal — with small N it is never rejected and with large N always. TOST declares a margin in advance and shows the difference falls inside it. The default margin is the observed WITHIN CV, that is the noise floor of the generator itself: asking for less would be asking the synthetic to resemble the real more closely than the real resembles itself.',
    '<b>z</b> = (media de las realizaciones − VOI) / sd. Con |z| &lt; 2 la diferencia es indistinguible del ruido del propio generador: no se le puede pedir al ajuste que acierte mas de lo que el generador se repite a si mismo.<br>Un CV pequeno significa que el generador reproduce esa metrica de forma estable, <b>no</b> que la metrica sea exacta: BS es de las mas estables y arrastra el sesgo de marching cubes.':
        '<b>z</b> = (mean of the realisations − VOI) / sd. With |z| &lt; 2 the difference is indistinguishable from the noise of the generator itself: the fit cannot be asked to hit closer than the generator repeats itself.<br>A small CV means the generator reproduces that metric stably, <b>not</b> that the metric is accurate: BS is among the most stable and carries the marching cubes bias.',
    '<b>{n} replicas generadas</b> — {como}<br>':
        '<b>{n} replicates generated</b> — {como}<br>',
    "<p style='font-size:13px'><b>{nv} VOIs · {nr} replicas cada uno · {nf} filas en la tabla</b><br><span style='color:#b62324;font-size:15px'>N efectivo para inferencia biologica: {nef}</span></p>":
        "<p style='font-size:13px'><b>{nv} VOIs · {nr} replicates each · {nf} rows in the table</b><br><span style='color:#b62324;font-size:15px'>Effective N for biological inference: {nef}</span></p>",
    '<p>Orden observado <b>{orden}</b> · extrapolado <b>{ext} MPa</b> · error estimado del punto mas fino <b>{err}%</b></p>':
        '<p>Observed order <b>{orden}</b> · extrapolated <b>{ext} MPa</b> · estimated error of the finest point <b>{err}%</b></p>',
    '<p>{n} realizaciones de la MISMA parametrizacion, semillas {a}…{b}</p>':
        '<p>{n} realisations of the SAME parameters, seeds {a}…{b}</p>',
    '<tr><th></th><th>campo</th><th>mediana</th><th>p98</th><th>maximo</th><th>max/p98</th></tr>':
        '<tr><th></th><th>field</th><th>median</th><th>p98</th><th>maximum</th><th>max/p98</th></tr>',
    '<tr><th></th><th>grupos</th><th>filas</th><th>N efectivo</th><th>CV entre</th><th>CV dentro</th><th>ICC</th><th>dif. real−sint.</th><th>margen</th><th>equivalente</th></tr>':
        '<tr><th></th><th>groups</th><th>rows</th><th>effective N</th><th>CV between</th><th>CV within</th><th>ICC</th><th>diff. real−synth.</th><th>margin</th><th>equivalent</th></tr>',
    '<tr><th></th><th>media</th><th>mediana</th><th>sd</th><th>p05</th><th>p95</th><th>CV</th></tr>':
        '<tr><th></th><th>mean</th><th>median</th><th>sd</th><th>p05</th><th>p95</th><th>CV</th></tr>',
    '<tr><th></th><th>media</th><th>sd</th><th>CV</th><th>rango</th><th>VOI</th><th>dif.</th><th>z</th></tr>':
        '<tr><th></th><th>mean</th><th>sd</th><th>CV</th><th>range</th><th>VOI</th><th>diff.</th><th>z</th></tr>',
    '<tr><th>metrica</th><th>media</th><th>sd</th><th>CV</th></tr>':
        '<tr><th>metric</th><th>mean</th><th>sd</th><th>CV</th></tr>',
    'Ajustando con los tres metodos comparables…':
        'Fitting with the three comparable methods…',
    'Ajustando…':
        'Fitting…',
    'Calculando el espesor local…':
        'Computing the local thickness…',
    'Compresion en {eje}':
        'Compression along {eje}',
    'Convergencia de malla':
        'Mesh convergence',
    'Convergencia sobre el {cual}: {n} mallas hasta {res}³…':
        'Convergence on the {cual}: {n} meshes up to {res}³…',
    'Dispersion calculada':
        'Spread computed',
    'Dispersion de {n} realizaciones en {t} s — CV entre {a}% y {b}%':
        'Spread of {n} realisations in {t} s — CV between {a}% and {b}%',
    'Dispersion del generador':
        'Generator spread',
    'Distribucion del espesor local':
        'Local thickness distribution',
    'Distribuciones del ensayo de compresion':
        'Distributions from the compression test',
    "El VOI '{ruta}' no esta donde lo dejaste; cargalo a mano si lo necesitas.":
        "The VOI '{ruta}' is not where you left it; load it by hand if you need it.",
    'Generando a {res}³ y mallando…':
        'Generating at {res}³ and meshing…',
    'Homogeneizando a {n}³ (y el VOI si esta cargado)…':
        'Homogenising at {n}³ (and the VOI if loaded)…',
    'La linea de puntos es el <b>percentil 98</b>, que es el que decide el criterio de Pistoia; la verde es el 0.7% critico. La columna <b>max/p98</b> es el diagnostico: si es cercana a 1 la cola esta repartida y el criterio es estable; si es grande, unos pocos elementos —a menudo esquinas escalonadas de la malla, sin significado fisico— dominan el extremo, y conviene mirar la convergencia de malla antes de creerse la carga de fallo.<br>Escala logaritmica en x: los campos abarcan varios ordenes de magnitud entre el material descargado y los nudos.':
        'The dotted line is the <b>98th percentile</b>, which is what the Pistoia criterion decides on; the green one is the critical 0.7%. The <b>max/p98</b> column is the diagnostic: close to 1 the tail is spread out and the criterion is stable; if it is large, a few elements —often stair-stepped corners of the mesh, with no physical meaning— dominate the extreme, and it is worth looking at mesh convergence before believing the failure load.<br>Logarithmic x scale: the fields span several orders of magnitude between unloaded material and the nodes.',
    'Las mascaras no se guardan: se regeneran exactas desde la semilla y los parametros, que si estan en el CSV.<br>-> {archivo}':
        'The masks are not saved: they are regenerated exactly from the seed and the parameters, which are in the CSV.<br>-> {archivo}',
    'Los puntos que no resolvieron no son un fallo del estudio: a resolucion baja una estructura poco densa queda casi desconectada y el sistema mal condicionado. El solver lo detecta por el residuo y se descarta el punto en vez de devolver un campo equivocado.':
        'The points that did not solve are not a failure of the study: at low resolution a low-density structure is almost disconnected and the system ill-conditioned. The solver detects it from the residual and the point is discarded rather than returning a wrong field.',
    'Orientacion FIJADA por una sesion guardada (correccion F4). Mueve cualquier rotacion para liberarla.':
        'Orientation FIXED by a saved session (correction F4). Move any rotation to release it.',
    'Orientacion tomada del eje principal del VOI (Euler {e}°). Manda sobre los deslizadores de rotacion hasta que muevas uno.':
        'Orientation taken from the principal axis of the VOI (Euler {e}°). It overrides the rotation sliders until you move one.',
    'REORDENO':
        'RE-RANKED',
    'Resolviendo la compresion a {n}³ en {eje}…':
        'Solving compression at {n}³ along {eje}…',
    'Todo en mm. El <b>CV</b> (sd/media) es lo que el Tb.Th escalar no puede dar: mide cuanto varia el grosor dentro de la misma estructura.<br><b>La escala absoluta no es citable como Tb.Th</b>: a 3-4 voxeles de grosor el metodo subestima un ~30%. Usar como medida RELATIVA entre estructuras medidas al mismo spacing.':
        'All in mm. The <b>CV</b> (sd/mean) is what the scalar Tb.Th cannot give: it measures how much thickness varies inside the same structure.<br><b>The absolute scale is not quotable as Tb.Th</b>: at 3-4 voxels thick the method underestimates by ~30%. Use it as a RELATIVE measure between structures measured at the same spacing.',
    'Varianza del lote y equivalencia':
        'Batch variance and equivalence',
    'coeficiente de variacion [%]':
        'coefficient of variation [%]',
    'coeficiente de variacion entre realizaciones [%]':
        'coefficient of variation across realisations [%]',
    'confirmo la morfometria':
        'confirmed the morphometry',
    'densidad':
        'density',
    'densidad de probabilidad':
        'probability density',
    'dentro (ruido del generador)':
        'within (generator noise)',
    'desempate mecanico:':
        'mechanical tie-break:',
    'desempate por orientacion activo pero no aplicado ({n} candidato/s dentro del 5%)':
        'orientation tie-break active but not applied ({n} candidate(s) within 5%)',
    'eje MIL: {ang} grados respecto al VOI':
        'MIL axis: {ang} degrees from the VOI',
    'entre especimenes (biologica)':
        'between specimens (biological)',
    'espesor local [mm]  ·  la linea de trazos es la mediana':
        'local thickness [mm]  ·  the dashed line is the median',
    'fraccion':
        'fraction',
    'identicas (solo cambia la semilla)':
        'identical (only the seed changes)',
    'los tres ejes':
        'all three axes',
    'no':
        'no',
    'resolucion [voxeles por lado]':
        'resolution [voxels per side]',
    'rojo: el VOI queda a 2 sd o mas de la media de las realizaciones':
        'red: the VOI is 2 sd or more from the mean of the realisations',
    'si':
        'yes',
    'variando {mets} con amplitud {amp} %':
        'varying {mets} with amplitude {amp} %',
    '{n} hexaedros con volumen no positivo; la conectividad estaria mal orientada.':
        '{n} hexahedra with non-positive volume; the connectivity would be wrongly oriented.',
    '{n} metrica/s fuera de banda (|z|>=2)':
        '{n} metric(s) out of band (|z|>=2)',

    # -- ventana de comparacion de metodos (dialogo_metodos.py) ---------
    '<b>Ajustando el VOI con {n} metodos a la vez.</b>':
        '<b>Fitting the VOI with {n} methods at once.</b>',
    '<b>Elige con que ajuste seguir.</b> El resto se descarta.':
        '<b>Choose which fit to continue with.</b> The rest are discarded.',
    '<b>Ni las identicas ni las variadas son especimenes.</b> Las identicas miden el ruido del generador; las variadas son puntos de diseno de una familia sintetica. Tratar N replicas como N especimenes infla los grados de libertad y estrecha los intervalos de confianza: es pseudorreplicacion.<br><br>La amplitud se aplica al PARAMETRO que gobierna cada metrica, no a la metrica: el generador toma densidad, numero de onda y angulos, y las metricas salen de ahi. La dispersion realmente lograda se reporta al terminar, y no tiene por que coincidir con la pedida.':
        '<b>Neither the identical nor the varied ones are specimens.</b> The identical ones measure the noise of the generator; the varied ones are design points of a synthetic family. Treating N replicates as N specimens inflates the degrees of freedom and narrows the confidence intervals: it is pseudoreplication.<br><br>The amplitude applies to the PARAMETER that governs each metric, not to the metric: the generator takes density, wave number and angles, and the metrics follow from those. The spread actually achieved is reported at the end, and need not match the one requested.',
    '<b>{metodo}</b> — error {err}':
        '<b>{metodo}</b> — error {err}',
    'Ajustar con todos los metodos':
        'Fit with all methods',
    'Amplitud (desviacion relativa):':
        'Amplitude (relative deviation):',
    'Cancelar':
        'Cancel',
    'Con variacion en las metricas que elija':
        'Varying the metrics you choose',
    'Corren en hilos y comparten memoria, asi que el reparto de CPU no es perfecto: juntos tardan algo menos que en serie, no la cuarta parte. El desempate mecanico es el mas lento con diferencia, porque homogeneiza cada finalista.':
        'They run in threads and share memory, so the CPU split is not perfect: together they take somewhat less than in series, not a quarter of it. The mechanical tie-break is by far the slowest, because it homogenises every finalist.',
    'El <b>error no es comparable entre metodos que usan distinto numero de terminos</b>: el del desempate mecanico incluye E<sub>z</sub>/E<sub>s</sub> y E<sub>z</sub>/E<sub>x</sub> y los demas no. Para comparar de igual a igual, mira las diferencias por metrica.':
        "The <b>error is not comparable between methods that use a different number of terms</b>: the mechanical tie-break's includes E<sub>z</sub>/E<sub>s</sub> and E<sub>z</sub>/E<sub>x</sub> and the others do not. To compare like with like, look at the per-metric differences.",
    'Error al generar':
        'Error while generating',
    'Identicas (misma parametrizacion, distinta semilla)':
        'Identical (same parameters, different seed)',
    'Los {n} metodos fallaron. El primero dice:\n\n{msg}':
        'All {n} methods failed. The first says:\n\n{msg}',
    'Marca al menos una metrica que variar, o elige replicas identicas.':
        'Tick at least one metric to vary, or choose identical replicates.',
    'Marca con que ajuste quieres seguir.':
        'Tick which fit you want to continue with.',
    'Metricas que deben variar':
        'Metrics that must vary',
    'Ningun metodo termino':
        'No method finished',
    'Numero de replicas:':
        'Number of replicates:',
    'Sin metricas':
        'No metrics',
    'Sin seleccion':
        'Nothing selected',
    'Usar el seleccionado':
        'Use the selected one',
    'Usar y generar replicas…':
        'Use and generate replicates…',
    'en cola':
        'queued',
    'fase poro fragmentada: no es trabecula':
        'pore phase fragmented: not trabecular',
    '{t} s · {n} evaluaciones':
        '{t} s · {n} evaluations',

    # -- modo de medida: se MUESTRA traducido, se GUARDA en su forma interna
    # ("malla"/"voxel" viajan al JSON y al CSV; traducir el dato haria que dos
    # tablas del mismo estudio no se pudieran juntar).
    "malla": "mesh",
    "voxel": "voxel",

    "Superficie suavizada, reparada y cerrada: es la que se imprime.":
        "Smoothed, repaired and closed surface: this is the one that gets "
        "printed.",

    # -- menu Validacion y su ventana (dialogo_validacion.py) ----------
    '&Validacion':
        '&Validation',
    '(sin imagen)':
        '(no image)',
    '<b>Este método no es nuestro: es el de Kumar et al. (2020).</b> Lo que sigue es su figura publicada al lado de la misma figura regenerada con este programa, y las comprobaciones numéricas que fallarían si nuestra implementación se hubiera desviado del artículo.':
        '<b>This method is not ours: it is Kumar et al. (2020).</b> What follows is their published figure next to the same figure regenerated by this program, and the numerical checks that would fail if our implementation had drifted from the paper.',
    'Abrir la carpeta Test':
        'Open the Test folder',
    'Baja la malla de vista a 96³ y la de homogeneización a 16³. La réplica sale en algo más de un minuto en vez de cuatro, y las seis comprobaciones se siguen pasando; lo que se pierde es definición en la figura.':
        'Drops the view mesh to 96³ and the homogenisation mesh to 16³. The replication takes a little over a minute instead of four, and the six checks still hold; what is lost is definition in the figure.',
    'Código':
        'Code',
    'Ejecutar la réplica':
        'Run the replication',
    'La réplica falló':
        'The replication failed',
    'Las {n} comprobaciones se cumplen: la implementación reproduce el artículo.':
        'All {n} checks hold: the implementation reproduces the paper.',
    'Modo rápido (mallas menores)':
        'Fast mode (smaller meshes)',
    'NO cumple':
        'does NOT hold',
    'No se pudo abrir la validación':
        'Validation could not be opened',
    'Obtenido':
        'Obtained',
    'Original publicado (Fig. 2)':
        'Published original (Fig. 2)',
    'Predicción del artículo':
        'Prediction from the paper',
    'Regenerado con esta aplicación':
        'Regenerated with this application',
    'Replica de Kumar et al. (2020)…':
        'Replication of Kumar et al. (2020)…',
    'Validacion':
        'Validation',
    'Validación contra los resultados publicados':
        'Validation against the published results',
    'Veredicto':
        'Verdict',
    'cumple':
        'holds',
    '{n} comprobación/es no se cumplen.':
        '{n} check(s) do not hold.',

    # -- predicciones de la replica (Test/replicar_kumar2020.py) -------
    # Llegan a la interfaz desde un JSON, no como literales del codigo de
    # la ventana; el comprobador las saca del guion que las define.
    'Columnar (conos en e2 y e3): e1 es el eje RIGIDO':
        'Columnar (cones on e2 and e3): e1 is the STIFF axis',
    "Con conos DESIGUALES, 'equitativo' no es la ecuacion (2)":
        "With UNEQUAL cones, 'uniform' is not equation (2)",
    'Criterio: {c}':
        'Criterion: {c}',
    'Cubica (tres conos iguales): mucho menos anisotropa en los ejes que la lamelar y la columnar. Una realizacion suelta no es simetrica; ver la nota de C4 en la cabecera del guion.':
        'Cubic (three equal cones): far less anisotropic along the axes than the lamellar and columnar ones. A single realisation is not symmetric; see the C4 note in the script header.',
    'E1/E3 < 0.6':
        'E1/E3 < 0.6',
    'E1/E3 > 1.6':
        'E1/E3 > 1.6',
    'E_max/E_min < 1.25':
        'E_max/E_min < 1.25',
    'Isotropa (th = 90 grados): E1 = E2 = E3':
        'Isotropic (th = 90 degrees): E1 = E2 = E3',
    'Lamelar (cono solo en e1): e1 es el eje BLANDO':
        'Lamellar (cone on e1 only): e1 is the SOFT axis',
    'Lamelar y columnar ordenan al reves el eje e1':
        'Lamellar and columnar rank axis e1 in opposite order',
    'difieren mas del 5 %':
        'they differ by more than 5 %',
    'menos de la mitad de la mayor anisotropia':
        'less than half the largest anisotropy',
    'razon lamelar < razon columnar':
        'lamellar ratio < columnar ratio',

    "El proceso terminó con código {c}.":
        "The process exited with code {c}.",
    "No se pudo lanzar el proceso de réplica:\n{p}":
        "The replication process could not be started:\n{p}",

    # -- validacion: replicas de Zheng 2021 y Guo 2024 --
    # Cabeceras, criterios y predicciones de las dos replicas
    # anadidas en 2026-09-10. Las de Kumar siguen mas arriba.

    '<b>Aquí se comprueban COTAS, no proporciones.</b> Su Figura 1 dibuja cada superficie elástica dentro de la esfera de Voigt y, en la isótropa, de la de Hashin-Shtrikman. Son cotas con fórmula cerrada: un error de escala en nuestra homogeneización no las pasa. Se comprueba además la ortotropía del tensor y la razón por la que el artículo impone ρ ≥ 0.3.':
        '<b>What is checked here are BOUNDS, not proportions.</b> Their Figure 1 draws every elastic surface inside the Voigt sphere and, for the isotropic class, inside the Hashin-Shtrikman one. These are closed-form bounds: a scale error in our homogenisation does not get past them. The orthotropy of the tensor and the reason the article imposes ρ ≥ 0.3 are checked as well.',
    '<b>Aquí se replican las dos columnas «Target» de su Figura 7.</b> El spinodoide lleva parámetros publicados verbatim, y la superficie nodal periódica tiene curvaturas con fórmula cerrada: es la que valida nuestro estimador contra una respuesta exacta que no hemos calculado nosotros. No se replica su diseño inverso por redes neuronales, ni su muestra de hueso, que no tenemos.':
        '<b>What is replicated here are the two «Target» columns of their Figure 7.</b> The spinodoid uses verbatim published parameters, and the periodic nodal surface has closed-form curvatures: that is the one validating our estimator against an exact answer we did not compute ourselves. Their neural-network inverse design is not replicated, nor their bone sample, which we do not have.',
    'A rho = 0.5 solido y vacio son intercambiables: la curvatura media del ensemble es cero':
        'At rho = 0.5 solid and void are interchangeable: the ensemble mean curvature is zero',
    'Baja la malla de vista a 96³ y la de homogeneización a 24³. Las cotas se cumplen igual —una malla gruesa da rigideces más bajas, o sea más holgura—, pero las cifras que se reporten deberían salir de la ejecución completa, que homogeneiza a 40³ como pide el estudio de convergencia.':
        'Drops the view mesh to 96³ and the homogenisation mesh to 24³. The bounds hold either way — a coarse mesh gives lower stiffness, i.e. more slack — but any figure to be reported should come from the full run, which homogenises at 40³ as the convergence study requires.',
    'Baja las rejillas a 96³ y 120³. Sirve para comprobar que la tubería funciona: el error del estimador discreto sube al 10 %, que es justo el criterio, y pasar por una centésima es pasar por suerte. Para dar por buena la implementación, ejecución completa.':
        "Drops the grids to 96³ and 120³. Good for checking that the pipeline runs: the discrete estimator's error rises to 10 %, which is exactly the criterion, and passing by a hundredth is passing by luck. To call the implementation sound, use the full run.",
    'Calculando. Tarda unos minutos; la ventana responde igual.':
        'Computing. It takes a few minutes; the window stays responsive.',
    'Con SUS ternas —traspuestas respecto de las de Kumar 2020— el eje del cono activo es el blando en la lamelar (disco) y el rigido en la columnar (huso)':
        'With THEIR triples — transposed with respect to those of Kumar 2020 — the axis of the active cone is the soft one in the lamellar class (disc) and the stiff one in the columnar class (spindle)',
    'Doblar beta dobla todas las curvaturas: es el mismo campo con las longitudes a la mitad':
        'Doubling beta doubles every curvature: it is the same field with all lengths halved',
    'E(d) <= E_HS+ para la unica clase donde la cota aplica':
        'E(d) <= E_HS+ for the only class where the bound applies',
    'E(d) <= rho*E_s en las cuatro clases':
        'E(d) <= rho*E_s in all four classes',
    'El metodo en si: las cuatro clases de anisotropia y el muestreo por rechazo de su ecuacion (2).':
        'The method itself: the four anisotropy classes and the rejection sampling of their equation (2).',
    'El perfil del spinodoide cae en la mancha de su panel (b): silla, con el centro dentro de la caja leida de la figura':
        'The spinodoid profile lands on the blob of their panel (b): saddle, with its centre inside the box read off the published figure',
    'En el spinodoide, el estimador discreto coincide con la formula cerrada del propio campo':
        'On the spinodoid, the discrete estimator agrees with the closed-form answer from the field itself',
    'En la superficie nodal periodica —cuya formula NO es nuestra— el estimador discreto coincide con la respuesta cerrada':
        'On the periodic nodal surface — whose formula is NOT ours — the discrete estimator agrees with the closed-form answer',
    'La PNS no es una superficie minima; el articulo la elige por eso, porque k2 = -k1 seria un objetivo trivial':
        'The PNS is not a minimal surface; the article picks it for that reason, since k2 = -k1 would be a trivial target',
    'La clase isotropa no supera la cota superior de Hashin-Shtrikman, que es la esfera gris oscura de su panel (e)':
        'The isotropic class does not exceed the Hashin-Shtrikman upper bound, which is the dark grey sphere of their panel (e)',
    'Las cotas: cada superficie elastica dentro de Voigt y de Hashin-Shtrikman, la ortotropia del tensor, y por que el articulo impone rho >= 0.3.':
        'The bounds: every elastic surface inside Voigt and Hashin-Shtrikman, the orthotropy of the tensor, and why the article imposes rho >= 0.3.',
    'Las curvaturas: el perfil (k1, k2) de un spinodoide de parametros publicados y el de una superficie nodal periodica, cuya respuesta es cerrada.':
        'The curvatures: the (k1, k2) profile of a spinodoid with published parameters and that of a periodic nodal surface, whose answer is closed-form.',
    'Marching cubes sobre la mascara BINARIA sobreestima el area frente a la isosuperficie del campo continuo (correccion C4)':
        'Marching cubes on the BINARY mask overestimates the area relative to the isosurface of the continuous field (correction C4)',
    'Ninguna superficie elastica sale de la esfera de Voigt: es la esfera gris clara de su Fig. 1, de radio rho*E_s':
        'No elastic surface leaves the Voigt sphere: it is the light grey sphere of their Fig. 1, of radius rho*E_s',
    'No se encuentra la carpeta <b>Test</b>. La validación necesita los guiones de réplica y <code>Test/referencia/</code>.':
        'The <b>Test</b> folder cannot be found. Validation needs the replication scripts and <code>Test/referencia/</code>.',
    'Original publicado (Fig. 1)':
        'Published original (Fig. 1)',
    'Original publicado (Fig. 7)':
        'Published original (Fig. 7)',
    'Replica de Guo et al. (2024)…':
        'Guo et al. (2024) replication…',
    'Replica de Zheng et al. (2021)…':
        'Zheng et al. (2021) replication…',
    'Todavía no se ha ejecutado esta réplica en este equipo. Pulsa «Ejecutar la réplica».':
        'This replication has not been run on this machine yet. Press «Run the replication».',
    'Validación':
        'Validation',
    'Y la restriccion tiene sentido: por debajo de rho_min el material desconectado crece en todas las clases':
        'And the restriction makes sense: below rho_min the disconnected material grows in every class',
    'entre 1.02 y 1.20; medido +8.5 % sobre una esfera':
        'between 1.02 and 1.20; measured +8.5 % on a sphere',
    'error mediano <= 10 % de la escala de curvatura':
        'median error <= 10 % of the curvature scale',
    'la fraccion desconectada crece al bajar rho, en las cuatro':
        'the disconnected fraction grows as rho drops, in all four',
    'lamelar < 0.6 y columnar > 1.6':
        'lamellar < 0.6 and columnar > 1.6',
    'razon = 2.00 con 5 % de holgura':
        'ratio = 2.00 with 5 % slack',
    '|H medio| <= 15 % de la desviacion tipica de H':
        '|mean H| <= 15 % of the standard deviation of H',
    '|H medio| > 10 % de la escala de curvatura':
        '|mean H| > 10 % of the curvature scale',
    'ρ = {rho} · β = {beta_pi}π · θ = {thetas} · N = {num_waves} · spinodoide {resolucion_spinodoide}³ · PNS {resolucion_pns}³ · semilla {semilla} · {unidad}':
        'ρ = {rho} · β = {beta_pi}π · θ = {thetas} · N = {num_waves} · spinodoid {resolucion_spinodoide}³ · PNS {resolucion_pns}³ · seed {semilla} · {unidad}',
    'ρ = {rho} · β = {beta_pi}π · ν_s = {nu_s} · N = {num_waves} · vista {resolucion_vista}³ · homogeneización {resolucion_homogeneizacion}³ · semilla {semilla} · contorno: {condiciones_contorno}':
        'ρ = {rho} · β = {beta_pi}π · ν_s = {nu_s} · N = {num_waves} · view {resolucion_vista}³ · homogenisation {resolucion_homogeneizacion}³ · seed {semilla} · boundary: {condiciones_contorno}',

    'el centro cae dentro de la caja leida de su panel (b), en 1/(lado del cubo)':
        'the centre falls inside the box read off their panel (b), in 1/(cube side)',
    'ρ = {rho} · β = {beta_pi}π · N = {num_waves} · vista {resolucion_vista}³ · homogeneización {resolucion_homogeneizacion}³ · semilla {semilla} · muestreo: {esquema}':
        'ρ = {rho} · β = {beta_pi}π · N = {num_waves} · view {resolucion_vista}³ · homogenisation {resolucion_homogeneizacion}³ · seed {semilla} · sampling: {esquema}',

    'A rho = 0.3 —el rho_min que el articulo impone para evitar dominios disjuntos— queda poco hueso fuera de la mayor componente, con el theta_min = pi/6 que fija la MISMA frase':
        'At rho = 0.3 — the rho_min the article imposes to avoid disjoint domains — little bone is left outside the largest component, under the theta_min = pi/6 set by the SAME sentence',
    'El tensor homogeneizado es ORTOTROPO (p. 11): al promediar realizaciones, el acoplamiento normal-cortante baja hacia cero. Una realizacion suelta NO lo cumple, y el articulo tampoco lo afirma de ella':
        'The homogenised tensor is ORTHOTROPIC (p. 11): averaging realisations drives the normal-shear coupling towards zero. A single realisation does NOT satisfy it, and the article does not claim it of one either',
    'HALLAZGO, no una prediccion del articulo: la lamelar con los 15 grados del pie de su Fig. 1 NO cumple el proposito de rho_min; con los 30 grados del texto si, y la diferencia es de un orden de magnitud':
        'FINDING, not a prediction of the article: the lamellar class with the 15 degrees of their Fig. 1 caption does NOT fulfil the purpose of rho_min; with the 30 degrees of the text it does, and the difference is an order of magnitude',
    'el promedio baja respecto de las realizaciones sueltas y queda por debajo del 5 % de la media de C11, C22, C33':
        'the average drops below the single realisations and stays under 5 % of the mean of C11, C22, C33',

    'Anadir el tamano de poro (Po.Dm)':
        'Add the pore size (Po.Dm)',
    'Po.Dm es el espesor local de la fase PORO, medido con esferas\ninscritas (Hildebrand & Ruegsegger sobre el complemento). Es la\nseparacion MEDIDA, frente a Tb.Sp, que es la separacion que\nDEDUCE el modelo de placas de Parfitt. En una pila de placas los\ndos coinciden; en hueso trabecular la diferencia entre ambos dice\ncuanto se aparta la estructura de ese modelo.\n\nCasilla aparte porque es la medida mas cara de la tabla: una\ntransformada de distancia por radio sobre el poro, que a BV/TV\n0.3 es el 70 % del volumen. Medido: 0.7 s a 64^3, 3.2 s a 96^3 y\n12 s a 128^3.\n\nAviso de borde: los poros que tocan la cara del cubo estan\ncortados y su esfera inscrita no esta acotada por ese lado. No se\ncorrige -no hay forma honesta de saber cuanto seguia el poro\nfuera- y el informe lleva la fraccion afectada.':
        'Po.Dm is the local thickness of the PORE phase, measured with\ninscribed spheres (Hildebrand & Ruegsegger on the complement). It is\nthe MEASURED separation, as opposed to Tb.Sp, which is the separation\nthe Parfitt plate model DEDUCES. On a stack of plates the two agree;\nin trabecular bone the difference between them says how far the\nstructure departs from that model.\n\nIt has its own checkbox because it is the most expensive measurement\nin the table: one distance transform per radius over the pore phase,\nwhich at BV/TV 0.3 is 70 % of the volume. Measured: 0.7 s at 64^3,\n3.2 s at 96^3 and 12 s at 128^3.\n\nBoundary caveat: pores touching the face of the cube are cut, and\ntheir inscribed sphere is unbounded on that side. This is not\ncorrected -there is no honest way to know how far the pore went on-\nand the report carries the affected fraction.',

    '>= 95 % de material portante en todas las semillas de las clases con theta >= 30 grados':
        '>= 95 % load-bearing material in every seed of the classes with theta >= 30 degrees',
    'al menos dos realizaciones de cinco rotas a 15 grados, y ninguna a 30':
        'at least two realisations out of five broken at 15 degrees, and none at 30',

    # -- familias de microestructura (spinodoide / dual-lattice) ------------
    'Familia': 'Family',
    'Familia:': 'Family:',
    'Familia: {f}': 'Family: {f}',
    'Ambas': 'Both',
    'Dual-lattice': 'Dual-lattice',
    'Dual-lattice (vista)': 'Dual-lattice (view)',
    'Dual-lattice ajustado ({m})': 'Fitted dual-lattice ({m})',
    'dual-lattice vs VOI': 'dual-lattice vs VOI',
    'Microestructura': 'Microstructure',
    'Parametros del dual-lattice': 'Dual-lattice parameters',
    'Celdas por lado': 'Cells per side',
    'Estiramiento X': 'Stretch X',
    'Estiramiento Y': 'Stretch Y',
    'Estiramiento Z': 'Stretch Z',
    'Irregularidad': 'Irregularity',
    'Exportar el dual-lattice como solido': 'Export the dual-lattice as a solid',
    'Ajustando con los dos metodos comparables…':
        'Fitting with the two comparable methods…',
    'Spinodoide: conjunto de nivel de un campo aleatorio gaussiano\n(Kumar et al. 2020).\nDual-lattice: red dual de una teselacion de Delaunay, con cuatro\nbarras por nudo (Vafaeefar et al. 2022).\n\nAmbas: se generan, miden y ajustan las DOS a la vez, con un panel\n3D y una columna de la tabla por familia. Las acciones pesadas\n(tensor, ensayo, exportacion, lote) se hacen una familia detras\nde otra, con el mismo material y la misma resolucion.':
        'Spinodoid: level set of a Gaussian random field\n(Kumar et al. 2020).\nDual-lattice: dual network of a Delaunay tessellation, with four\nstruts per node (Vafaeefar et al. 2022).\n\nBoth: the TWO are generated, measured and fitted together, with one\n3D panel and one table column per family. Heavy actions (tensor,\ntest, export, batch) run one family after the other, with the same\nmaterial and the same resolution.',
    'Red dual de una teselacion de Delaunay (Vafaeefar et al. 2022):\nel centroide de cada tetraedro unido a los de sus cuatro caras,\nengrosado hasta la densidad pedida. Cuatro barras por nudo.':
        'Dual network of a Delaunay tessellation (Vafaeefar et al. 2022):\nthe centroid of each tetrahedron joined to those of its four faces,\nthickened up to the requested density. Four struts per node.',
    'El estiramiento alarga la MALLA antes de engrosarla: el eje estirado es el eje rigido. Con 2.5 en un eje el DA ronda 1.5, el del VOI proximal equino. El grosor de las barras no se elige: sale de la densidad.':
        'The stretch elongates the MESH before thickening it: the stretched axis is the stiff axis. With 2.5 on one axis DA is about 1.5, that of the proximal equine VOI. Strut thickness is not chosen: it follows from the density.',
    'Anadir el Ellipsoid Factor (EF)': 'Add the Ellipsoid Factor (EF)',
    'Ellipsoid Factor (Doube 2015): en cada voxel, el mayor elipsoide\nque lo contiene y cabe en el hueso. EF = a/b - b/c vale -1 en una\nplaca, 0 en una esfera y +1 en una barra.\n\nMide lo mismo que el SMI —placa o barra— sin su defecto: el SMI\nse confunde con la CONCAVIDAD (Salmon et al. 2015), y el hueso\ntiene mucha mas superficie concava que un spinodoide.\n\nAVISO DE COSTE: es la medida mas cara del programa, del orden de\nminutos por estructura a 64^3 y de media hora a 97^3. No se\nrecorta para ir mas rapido: recortarlo sesga hacia barra.':
        'Ellipsoid Factor (Doube 2015): at each voxel, the largest ellipsoid\nthat contains it and fits inside the bone. EF = a/b - b/c is -1 for a\nplate, 0 for a sphere and +1 for a rod.\n\nIt measures the same thing as the SMI —plate or rod— without its\nflaw: the SMI is confounded by CONCAVITY (Salmon et al. 2015), and\nbone has far more concave surface than a spinodoid.\n\nCOST WARNING: it is the most expensive measurement in the program,\nminutes per structure at 64^3 and half an hour at 97^3. It is not\ntrimmed to run faster: trimming it biases toward rod.',
    '<p><b>Orientacion del eje principal</b> (MIL, medido sobre la malla que se ensaya): VOI a {av}° del eje Z (DA {dav}), dual-lattice a {as_}° (DA {das}). <b>Angulo entre ambos: {ang}°.</b></p>':
        '<p><b>Orientation of the principal axis</b> (MIL, measured on the mesh being tested): VOI at {av}° from the Z axis (DA {dav}), dual-lattice at {as_}° (DA {das}). <b>Angle between them: {ang}°.</b></p>',
    '<b>Aviso:</b> el dual-lattice actual no procede de un ajuste al VOI. A densidad igual, el estiramiento de la malla decide en que eje es rigido: ajusta al VOI antes de comparar, o estaras midiendo la orientacion que dejaron los deslizadores.':
        '<b>Warning:</b> the current dual-lattice does not come from a fit to the VOI. At equal density, the mesh stretch decides which axis is stiff: fit to the VOI before comparing, or you will be measuring the orientation the sliders left behind.',

    # -- avisos de interpretabilidad (spinpy/avisos.py) --------------------
    "BV/TV {bv} ≥ {umbral}: a esta densidad el VOI no es trabecular (poros aislados). Un ajuste con error bajo no significa que se este imitando trabecula.":
        "BV/TV {bv} ≥ {umbral}: at this density the VOI is not trabecular (isolated pores). A low fit error does not mean trabecular bone is being reproduced.",
    "VOI no trabecular (BV/TV {bv})":
        "non-trabecular VOI (BV/TV {bv})",
    "eje a {ang}° del VOI: revisalo antes de ensayar":
        "axis {ang}° away from the VOI: check it before testing",
    "Candidato mal orientado":
        "Misoriented candidate",
    "El eje principal del {fam} esta a {ang}° del eje del VOI (umbral {umbral}°).\n\nEl ensayo cargaria el candidato en una direccion que no es la del hueso, y su rigidez no seria comparable con la del VOI. Suele pasar cuando la fabrica del candidato es casi plana y la correccion de orientacion no puede alinearla; en el banco de VOIs, el dual-lattice quedo alineado en todos.\n\n¿Ensayar de todos modos?":
        "The principal axis of the {fam} is {ang}° away from the VOI axis (threshold {umbral}°).\n\nThe test would load the candidate in a direction that is not the bone's, and its stiffness would not be comparable with the VOI's. This usually happens when the candidate's fabric is nearly planar and the orientation correction cannot align it; across the VOI bank, the dual-lattice stayed aligned in every case.\n\nTest anyway?",

    # -- simulaciones in silico (spinpy/simulacion.py) ---------------------
    '<p><b>Al final:</b> rigidez {E} y conectividad {c} de las iniciales, con BV/TV {b} del inicial.</p>':
        '<p><b>At the end:</b> stiffness {E} and connectivity {c} of the initial values, with BV/TV {b} of the initial.</p>',
    '<p><b>Carga de fallo estimada</b> {F} N en el paso {k}, con un {d} % del tejido danado.</p>':
        '<p><b>Estimated failure load</b> {F} N at step {k}, with {d} % of the tissue damaged.</p>',
    '<p>La estructura <b>colapsa en el paso {c}</b>: su rigidez cae por debajo del {u} % de la inicial. Cuanto antes colapsa con menos dano, mas fragil es.</p>':
        '<p>The structure <b>collapses at step {c}</b>: its stiffness drops below {u} % of the initial. The earlier it collapses with less damage, the more brittle it is.</p>',
    '<p>La estructura no colapsa en los pasos simulados: reparte el dano.</p>':
        '<p>The structure does not collapse within the simulated steps: it spreads the damage.</p>',
    'tras el colapso':
        'after collapse',
    '<p><b>Recuperacion:</b> con la masa osea de vuelta al {h} %, la rigidez vuelve a {E} de la inicial (su minimo fue {m}).</p>':
        '<p><b>Recovery:</b> with bone mass back at {h} %, stiffness returns to {E} of the initial (its minimum was {m}).</p>',
    '<p>A igual masa osea, falta rigidez: es arquitectura perdida. Engrosar las trabeculas que quedan no devuelve las que desaparecieron.</p>':
        '<p>At equal bone mass, stiffness is missing: this is lost architecture. Thickening the remaining trabeculae does not bring back the ones that disappeared.</p>',
    '<p>A igual masa osea, la rigidez supera la inicial en esta direccion: el hueso se redistribuyo de trabeculas finas a gruesas. Comprueba otra direccion antes de leerlo como mejora: puede haberse perdido rigidez transversal.</p>':
        '<p>At equal bone mass, stiffness exceeds the initial value in this direction: bone was redistributed from thin to thick trabeculae. Check another direction before reading it as an improvement: transverse stiffness may have been lost.</p>',
    'La malla mecanica tiene {v} elementos por trabecula (Tb.Th/h), por debajo de {u}: los cocientes de rigidez pueden ser artefacto de resolucion. Sube la resolucion del ensayo.':
        'The mechanical mesh has {v} elements per trabecula (Tb.Th/h), below {u}: the stiffness ratios may be a resolution artefact. Increase the test resolution.',
    'Adelgazamiento uniforme':
        'Uniform thinning',
    'BV/TV {b} · E/E0 {e}':
        'BV/TV {b} · E/E0 {e}',
    'CSV (*.csv)':
        'CSV (*.csv)',
    'Candidato activo':
        'Active candidate',
    'Carga un VOI de referencia o elige el candidato activo.':
        'Load a reference VOI or choose the active candidate.',
    'Carga, ablanda el tejido que el criterio de Pistoia da por roto y\nvuelve a cargar. La carga maxima de la serie estima la resistencia\nultima; como cae dice si la estructura es fragil o reparte el dano.':
        'Loads, softens the tissue the Pistoia criterion deems broken, and\nloads again. The maximum load of the series estimates the ultimate\nstrength; how it drops tells whether the structure is brittle or spreads the damage.',
    'Demasiada perdida':
        'Too much loss',
    'Desuso guiado por carga':
        'Load-driven disuse',
    'E / E inicial':
        'E / initial E',
    'El VOI es el gemelo digital del animal: la simulacion responde\nque le pasaria a ESE hueso. El candidato sirve para comparar\nsi una familia sintetica pierde rigidez como el hueso real.':
        'The VOI is the digital twin of the animal: the simulation answers\nwhat would happen to THAT bone. The candidate lets you compare\nwhether a synthetic family loses stiffness like real bone.',
    'En cada paso el tejido que el criterio de Pistoia da por roto conserva solo el 5 % de su rigidez, y se vuelve a cargar. Tras el colapso (zona gris) la carga calculada vuelve a subir porque la sostiene ese tejido ablandado, no resistencia real: por eso la carga de fallo se toma antes. Elastico lineal con dano en un escalon, sin plasticidad ni contacto. El 3D muestra el tejido intacto. Compara estructuras; no da una carga de rotura absoluta.':
        'At each step the tissue the Pistoia criterion deems broken keeps only 5 % of its stiffness, and the structure is loaded again. After collapse (grey area) the computed load rises again because that softened tissue carries it, not real strength: that is why the failure load is taken before collapse. Linear elastic with one-step damage, no plasticity and no contact. The 3D view shows the intact tissue. It compares structures; it does not give an absolute failure load.',
    'Estructura:':
        'Structure:',
    'Estructura: {e} · malla mecanica {n}³ · apoyo {a} · eje {x} · {t} s':
        'Structure: {e} · mechanical mesh {n}³ · support {a} · axis {x} · {t} s',
    'Exportar tabla':
        'Export table',
    'Exportar tabla CSV…':
        'Export CSV table…',
    'Fallo progresivo':
        'Progressive failure',
    'Fallo progresivo en {t} s: carga maxima {F} N en el paso {k}.':
        'Progressive failure in {t} s: maximum load {F} N at step {k}.',
    'Fallo progresivo sobre el {cual}: {n} pasos…':
        'Progressive failure on the {cual}: {n} steps…',
    'Fallo progresivo…':
        'Progressive failure…',
    'Fraccion del hueso INICIAL que se retira en cada paso, igual en\ntodos los protocolos. En el fallo progresivo no se usa: alli cada\npaso ablanda el tejido que rompe.':
        'Fraction of the INITIAL bone removed at each step, the same in\nall protocols. Not used in progressive failure: there each\nstep softens the tissue that breaks.',
    'Genera o ajusta primero un candidato.':
        'Generate or fit a candidate first.',
    'Hueso por paso:':
        'Bone per step:',
    'Librerias: numpy · scipy.ndimage (distancias, espesor local) · scipy.sparse.linalg (ensayo) · matplotlib · pyvistaqt':
        'Libraries: numpy · scipy.ndimage (distances, local thickness) · scipy.sparse.linalg (test) · matplotlib · pyvistaqt',
    'Paso:':
        'Step:',
    'Pasos:':
        'Steps:',
    'Perdida de trabeculas finas':
        'Loss of thin trabeculae',
    'Perdida y recuperacion':
        'Loss and recovery',
    'Protocolo:':
        'Protocol:',
    'Que hueso se pierde primero:\n\nAdelgazamiento: la superficie, capa a capa. Envejecimiento.\nTrabeculas finas: las de menor espesor, enteras. OVX,\n  remodelado acelerado.\nDesuso: la superficie menos deformada en el ensayo de cada paso\n  (mecanostato). Inmovilizacion, microgravedad.\nRecuperacion: perdida por trabeculas finas y despues\n  engrosamiento hasta la masa inicial. Farmaco anabolico.':
        'Which bone is lost first:\n\nThinning: the surface, layer by layer. Ageing.\nThin trabeculae: the thinnest ones, whole. OVX,\n  accelerated remodelling.\nDisuse: the least strained surface in each step\'s test\n  (mechanostat). Immobilisation, microgravity.\nRecovery: loss of thin trabeculae, then\n  thickening back to the initial mass. Anabolic drug.',
    'Retira hueso en pasos y, en cada uno, mide la morfometria y\nresuelve el ensayo de compresion. Devuelve la rigidez y la carga\nde fallo frente a BV/TV, y el 3D de cada paso.\n\nTarda: un ensayo y una morfometria completa por paso.':
        'Removes bone in steps and, at each one, measures morphometry and\nsolves the compression test. Returns stiffness and failure load\nversus BV/TV, and the 3D of every step.\n\nSlow: one test and a full morphometry per step.',
    'Simulacion de perdida osea':
        'Bone loss simulation',
    'Simulacion en {t} s: rigidez final {E} de la inicial.':
        'Simulation in {t} s: final stiffness {E} of the initial.',
    'Simulacion — {etapa} ({i} de {n})':
        'Simulation — {etapa} ({i} of {n})',
    'Simulaciones in silico':
        'In silico simulations',
    'Simulando perdida osea sobre el {cual}: {n} pasos…':
        'Simulating bone loss on the {cual}: {n} steps…',
    'Simular perdida osea…':
        'Simulate bone loss…',
    'Tabla guardada en {f}':
        'Table saved to {f}',
    'Todos los protocolos retiran la misma cantidad de hueso por paso: a igual paso, igual BV/TV, y las diferencias de rigidez entre protocolos son de arquitectura. Un paso es una cantidad de hueso, no un tiempo. E<sub>app</sub> absoluto a esta resolucion no es citable; los cocientes respecto al paso 0 si, porque todos los pasos comparten el sesgo de malla.':
        'All protocols remove the same amount of bone per step: at equal step, equal BV/TV, and stiffness differences between protocols are architectural. A step is an amount of bone, not a time. Absolute E<sub>app</sub> at this resolution is not citable; ratios to step 0 are, because every step shares the mesh bias.',
    'Todos los protocolos retiran la misma cantidad de hueso por paso: las diferencias de rigidez entre ellos son de arquitectura, no de masa. Un paso es una cantidad de hueso, no un tiempo. La mecanica usa la resolucion, el apoyo y la direccion del ensayo de compresion.':
        'All protocols remove the same amount of bone per step: stiffness differences between them are architectural, not due to mass. A step is an amount of bone, not a time. Mechanics use the resolution, support and direction of the compression test.',
    'carga de fallo':
        'failure load',
    'carga de fallo (Pistoia) [N]':
        'failure load (Pistoia) [N]',
    'inicial':
        'initial',
    'paso':
        'step',
    'perdida':
        'loss',
    'recuperacion':
        'recovery',
    'relativo al paso 0':
        'relative to step 0',
    'rigidez':
        'stiffness',
    'tejido danado acumulado [%]':
        'cumulative damaged tissue [%]',
    '{p} pasos de {q} % retiran el {t} % del hueso: no quedaria estructura que ensayar. Baja los pasos o el hueso por paso.':
        '{p} steps of {q} % remove {t} % of the bone: no structure would be left to test. Reduce the steps or the bone per step.',

    # -- informe para publicacion ------------------------------------------
    'Publicacion':
        'Publication',
    'Informe para publicación…':
        'Publication report…',
    'Párrafo de métodos (ES/EN) con los valores usados, citabilidad\nde cada resultado y paquete de reproducción verificable.':
        'Methods paragraph (ES/EN) with the values used, citability\nof every result and a verifiable reproduction package.',
    'Nada que informar':
        'Nothing to report',
    'Mide, ajusta u homogeneiza algo primero: el informe describe lo que ya se calculo.':
        'Measure, fit or homogenise something first: the report describes what has already been computed.',
    'Elige la carpeta del informe para publicacion':
        'Choose the folder for the publication report',
    'La carpeta ya tiene un informe':
        'The folder already contains a report',
    'Se sobrescribiran {n} archivo(s) de un informe anterior en esa carpeta. ¿Continuar?':
        '{n} file(s) from a previous report in that folder will be overwritten. Continue?',
    'Informe: regenerando las estructuras y dibujando las figuras…':
        'Report: regenerating the structures and drawing the figures…',
    'Informe: renderizando las estructuras en 3D…':
        'Report: rendering the structures in 3D…',
    'Informe: componiendo los PDF…':
        'Report: typesetting the PDFs…',
    '{n} archivos: informe en PDF y Markdown (ES/EN), {f} figuras y el paquete de reproduccion.':
        '{n} files: report as PDF and Markdown (ES/EN), {f} figures and the reproduction package.',
    'Carpeta: {c}':
        'Folder: {c}',
    'Abrir la carpeta':
        'Open the folder',
    '<b>Informe para publicacion</b>':
        '<b>Publication report</b>',
    '{c} citables, {r} con reservas, {n} no citables.':
        '{c} citable, {r} with caveats, {n} not citable.',
    '<b>No citables</b> (no deben aparecer como resultado):':
        '<b>Not citable</b> (must not appear as a result):',
    'Cada resultado, con su motivo, esta en el informe.':
        'Every result, with its reason, is in the report.',
    'Informe para publicacion':
        'Publication report',
    'Informe escrito en {d}':
        'Report written to {d}',

    # -- informe para publicacion automatico --------------------------------
    'Informe para publicación (auto)…':
        'Publication report (auto)…',
    'Corre sobre el VOI cargado todo el recorrido —mejor ajuste,\nmorfometría completa, tensor elástico, ensayo de compresión y\nanálisis comparado, más las etapas opcionales que marques— y\nescribe el informe al final. Puede tardar horas.':
        'Runs the whole workflow on the loaded VOI —best fit, full\nmorphometry, elastic tensor, compression test and compared\nanalysis, plus any optional stages you tick— and writes the\nreport at the end. It may take hours.',
    'Detener tras esta etapa':
        'Stop after this stage',
    'Un calculo en marcha no se puede interrumpir sin perderlo: la\netapa en curso termina y las siguientes ya no se lanzan. Lo\ncalculado queda en la sesión.':
        'A running computation cannot be interrupted without losing it:\nthe current stage finishes and the following ones are not\nlaunched. What was computed stays in the session.',
    'Informe para publicacion (automatico)':
        'Publication report (automatic)',
    'Corre sobre el VOI <b>{voi}</b> y la familia <b>{fam}</b> todo el recorrido marcado, con los mismos calculos que los botones del panel, y escribe el informe al final. Los dialogos de resultados no se abren: todo queda en la sesion y en el informe. Puede tardar <b>horas</b>; la ventana sigue respondiendo y se puede detener tras la etapa en curso.':
        'Runs the ticked workflow on the VOI <b>{voi}</b> and the family <b>{fam}</b>, with the same computations as the panel buttons, and writes the report at the end. Result dialogs are not opened: everything stays in the session and in the report. It may take <b>hours</b>; the window stays responsive and the run can be stopped after the current stage.',
    'Carpeta del informe:':
        'Report folder:',
    'Etapas que gradua el informe':
        'Stages graded by the report',
    'Mejor ajuste al VOI: el de menor error entre rapido, completo y equitativo':
        'Best fit to the VOI: the lowest error among fast, full and equitable',
    'El desempate mecanico queda fuera: su error suma dos terminos y no es comparable.':
        'The mechanical tie-break is left out: its error adds two terms and is not comparable.',
    '<b>Tiempo estimado: {t}</b> · terminaria hacia las {h}':
        '<b>Estimated time: {t}</b> · would finish around {h}',
    'Con esta resolucion del ensayo no salen tres mallas por encima de 12³.':
        'At this test resolution there are not three meshes above 12³.',
    'Modelo medido en el equipo de desarrollo y corregido con {n} ejecucion(es) en este equipo. Suele acertar dentro de un ±20 %; en un equipo nuevo, hasta un factor 2 hasta que corra una vez.':
        'Model measured on the development machine and corrected with {n} run(s) on this machine. It is usually within ±20 %; on a new machine, up to a factor of 2 until it has run once.',
    'Modelo medido en el equipo de desarrollo; aun no hay ejecuciones en este equipo para corregirlo, asi que puede errar hasta un factor 2. Se ajusta solo al terminar cada etapa.':
        'Model measured on the development machine; there are no runs on this machine yet to correct it, so it may be off by up to a factor of 2. It adjusts itself as each stage finishes.',
    'Ambas (VOI / Spinodoide / Dual-lattice)':
        'Both (VOI / Spinodoid / Dual-lattice)',
    'Con «Ambas» cada etapa se corre para las dos familias, las figuras las muestran lado a lado y el informe añade una tabla VOI / Spinodoide / Dual-lattice. El tiempo casi se duplica.':
        'With «Both» every stage runs for the two families, the figures show them side by side and the report adds a VOI / Spinodoid / Dual-lattice table. The time nearly doubles.',
    'Figuras 3D del informe':
        '3D figures of the report',
    'Estilos':
        'Styles',
    'Vistas':
        'Views',
    'Figura principal, estilo:':
        'Main figure, style:',
    'Figura principal, vista:':
        'Main figure, view:',
    'Suavizar la superficie (Taubin, solo para la figura)':
        'Smooth the surface (Taubin, figure only)',
    'Quita la escalera de los voxeles sin encoger la pieza ni adelgazar las trabeculas. Las medidas usan siempre la malla sin suavizar, y el pie de la figura lo dice.':
        'Removes the voxel staircase without shrinking the part or thinning the trabeculae. Measurements always use the unsmoothed mesh, and the figure caption says so.',
    'Todas las vistas usan proyeccion paralela y la misma escala. La figura 4 reune las vistas marcadas en el estilo principal, la 5 son cortes 2D, la 6 la anisotropia (fabrica y E direccional), y cada combinacion de estilo y vista queda en figuras/3d/ para elegir otra sin recalcular.':
        'All views use parallel projection and the same scale. Figure 4 gathers the ticked views in the main style, figure 5 shows 2D sections, figure 6 the anisotropy (fabric and directional E), and every style and view combination is saved in figuras/3d/ to pick another without recomputing.',
    'Figuras adicionales':
        'Additional figures',
    'Distribuciones de espesor trabecular y tamaño de poro':
        'Trabecular thickness and pore size distributions',
    'Espesor local por esferas inscritas en hueso y poro. Es la figura cara: del orden de un minuto para un VOI de 188³.':
        'Local thickness by inscribed spheres in bone and pore. This is the expensive figure: about a minute for a 188³ VOI.',
    'Mapa 3D de von Mises (del análisis comparado)':
        '3D von Mises map (from the compared analysis)',
    'Figura del método (del micro-CT al candidato)':
        'Method figure (from micro-CT to candidate)',
    "Figura 0, en 16:9: rebanada y segmentación, el cubo en la pila, el VOI "
    "y, por cada familia, su aleatoriedad, su campo, el umbral por densidad "
    "y el sólido, con la morfometría frente al VOI. La fila del micro-CT "
    "solo aparece si el VOI se recortó de una pila en esta sesión.":
        "Figure 0, in 16:9: slice and segmentation, the cube in the stack, "
        "the VOI and, for each family, its randomness, its field, the "
        "density threshold and the solid, with the morphometry against the "
        "VOI. The micro-CT row only appears if the VOI was cropped from a "
        "stack in this session.",
    'Informe y figuras':
        'Report and figures',
    'no se lanzaria':
        'would not run',
    'quedan {t}':
        '{t} left',
    'estimado {t}':
        'estimated {t}',
    'etapa':
        'stage',
    'estado':
        'status',
    'real':
        'actual',
    'estimado':
        'estimated',
    'Desmarcala para usar el ajuste que ya hay en la sesion.':
        'Untick it to use the fit already in the session.',
    'Morfometria completa de la estructura ajustada y del VOI':
        'Full morphometry of the fitted structure and of the VOI',
    'Conn.D y SMI':
        'Conn.D and SMI',
    'Po.Dm':
        'Po.Dm',
    'EF (muy caro)':
        'EF (very expensive)',
    'La medida mas cara del programa: minutos por estructura a 64³ y media hora a 97³. No se recorta para ir mas rapido.':
        'The most expensive measurement in the program: minutes per structure at 64³ and half an hour at 97³. It is not cut short to go faster.',
    'Tensor elastico (homogeneizacion periodica)':
        'Elastic tensor (periodic homogenisation)',
    'Ensayo de compresion y fallo (Pistoia)':
        'Compression test and failure (Pistoia)',
    'Analisis comparado (protocolo de Tapia et al. 2026), a la resolucion del ensayo':
        'Compared analysis (Tapia et al. 2026 protocol), at the test resolution',
    'El tensor elastico a {n}³ queda por debajo de {m}³: el informe lo dara con reservas (la interfaz lo limita a {max}³).':
        'The elastic tensor at {n}³ is below {m}³: the report will give it with caveats (the interface caps it at {max}³).',
    'El ensayo a {n}³ queda por debajo de {m}³: el informe dara la mecanica con reservas.':
        'The test at {n}³ is below {m}³: the report will give the mechanics with caveats.',
    'Etapas opcionales (el informe las gradua y las dibuja si se corren)':
        'Optional stages (the report grades and plots them when they are run)',
    'Dispersion entre semillas':
        'Scatter across seeds',
    'Convergencia de malla del ensayo sobre el VOI':
        'Mesh convergence of the test on the VOI',
    'Simulacion de perdida osea ({est}, {prot}, {n} pasos)':
        'Bone loss simulation ({est}, {prot}, {n} steps)',
    'Fallo progresivo ({est}, {n} pasos)':
        'Progressive failure ({est}, {n} steps)',
    'Las simulaciones usan la estructura, el protocolo y los pasos de su seccion del panel; cambialos alli antes si hace falta.':
        'The simulations use the structure, protocol and steps set in their panel section; change them there first if needed.',
    'Empezar':
        'Start',
    'El informe automatico parte del VOI de referencia. Cargalo primero.':
        'The automatic report starts from the reference VOI. Load it first.',
    'Informe automatico: empezando…':
        'Automatic report: starting…',
    'Mejor ajuste al VOI':
        'Best fit to the VOI',
    'Morfometria completa':
        'Full morphometry',
    'Ensayo de compresion y fallo':
        'Compression test and failure',
    'Convergencia de malla (VOI)':
        'Mesh convergence (VOI)',
    'Analisis comparado':
        'Compared analysis',
    'Informe automatico {i}/{n}: {etapa}':
        'Automatic report {i}/{n}: {etapa}',
    'Informe automatico: escribiendo el informe…':
        'Automatic report: writing the report…',
    'se detiene al acabar esta etapa':
        'stops when this stage ends',
    'Informe automatico detenido. Lo calculado sigue en la sesion.':
        'Automatic report stopped. What was computed stays in the session.',
    'Informe automatico detenido':
        'Automatic report stopped',
    'No se ha escrito el informe. Lo calculado hasta aqui sigue en la sesion: «Informe para publicación…» lo describe tal como esta.':
        'The report was not written. What was computed so far stays in the session: «Publication report…» describes it as it is.',
    'Error en la etapa; se sigue con la siguiente':
        'Error in the stage; moving on to the next one',
    '{fam}: eje a {ang}° del VOI; ensayado igualmente, el informe lo marca como desalineado':
        '{fam}: axis at {ang}° from the VOI; tested anyway, the report flags it as misaligned',
    'Recorrido automatico':
        'Automatic run',
    '{t} en total':
        '{t} in total',
    'hecha':
        'done',
    'omitida':
        'skipped',
    'con error':
        'failed',
    'no lanzada':
        'not launched',
    'Ajustando con los metodos comparables…':
        'Fitting with the comparable methods…',

    # --- Pantalla de bienvenida -------------------------------------------
    'spinpy — visor de spinodoides y VOIs':
        'spinpy — spinodoid and VOI viewer',
    'Bienvenido a spinpy':
        'Welcome to spinpy',
    'Acerca de spinpy…':
        'About spinpy…',
    'Ajuste de microestructuras espinodales a VOIs de hueso trabecular':
        'Fitting spinodoid microstructures to trabecular bone VOIs',
    'versión':
        'version',
    'AUTORES':
        'AUTHORS',
    'Dos familias':
        'Two families',
    'Spinodoide (campo aleatorio gaussiano) y dual-lattice, sobre la misma '
    'rejilla de vóxeles. Se ajustan y se comparan a la vez.':
        'Spinodoid (Gaussian random field) and dual-lattice, on the same voxel '
        'grid. They are fitted and compared at once.',
    'Morfometría':
        'Morphometry',
    'BV/TV, BS/BV, Tb.Th, Tb.Sp, Tb.N, anisotropía por tensor MIL, Conn.D, '
    'SMI, tamaño de poro y factor de elipsoide.':
        'BV/TV, BS/BV, Tb.Th, Tb.Sp, Tb.N, anisotropy from the MIL tensor, '
        'Conn.D, SMI, pore size and ellipsoid factor.',
    'Mecánica':
        'Mechanics',
    'Homogeneización elástica periódica, ensayo de compresión y carga de '
    'fallo por el criterio de Pistoia.':
        'Periodic elastic homogenization, compression test and failure load from '
        'the Pistoia criterion.',
    'In silico':
        'In silico',
    'Simulaciones de pérdida ósea y de fallo progresivo sobre el gemelo '
    'digital del VOI.':
        'Simulations of bone loss and of progressive failure on the digital twin '
        'of the VOI.',
    'Informe':
        'Report',
    'Dice de cada número si es citable, con qué reservas y con qué parámetros '
    'se obtuvo.':
        'For every number it states whether it is quotable, with what caveats, '
        'and with which parameters it was obtained.',
    'Exportación':
        'Export',
    'STL, VTU, Abaqus, ANSYS y FEBio.':
        'STL, VTU, Abaqus, ANSYS and FEBio.',
    'No volver a mostrar':
        'Do not show again',
    'Entrar':
        'Enter',
    'Código:':
        'Code:',
    'Código: el repositorio público aún no está publicado; el enlace '
    'aparecerá aquí cuando lo esté.':
        'Code: the public repository is not published yet; the link will appear '
        'here once it is.',
}
