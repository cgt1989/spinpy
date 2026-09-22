"""voi.py — De la pila del escaner al VOI cubico, sin salir de Python.

QUE CIERRA ESTE MODULO
-----------------------
Hasta aqui el port empezaba en «cargame un .vtk que alguien ya extrajo». Todo
el camino anterior —umbralizar las rebanadas, apilarlas, orientar el hueso por
PCA y recortar los cubos— vivia solo en `segmentacion_por_lote.m`,
`Segmentacion_a_3D.m` y `VOIsTIFF.m`. Eso obligaba a tener MATLAB para poder
usar la herramienta, que es tanto como decir que otra persona no podia usarla.

Y hay una razon mas fuerte que la comodidad. Medida la sensibilidad al umbral
sobre H4, mover la segmentacion un +-15 % mueve Tb.Th un 62 % y BS/BV un 58 %,
mientras que la diferencia entre spinpy y BoneJ es del 6 %. **La segmentacion
es la mayor fuente de incertidumbre de todo el proceso**, y era justo el unico
paso que no se podia tocar desde aqui. Con esto se puede barrer el umbral y
reextraer los VOIs REALES —alineados por PCA, no cubos de conveniencia— para
poner barras de error a cada metrica.

CONVENCION DE EJES
------------------
Las mascaras salen con **dimension 1 = X**, igual que `readVTKVOI` y que todo
el resto de spinpy. `imread` devuelve (fila, columna) = (y, x), asi que cada
rebanada se traspone al apilar. Mezclar las dos convenciones es el bug B1/B2
del historial del proyecto.

EL SIGNO DE LOS EJES PRINCIPALES ES ARBITRARIO
-----------------------------------------------
La PCA determina las DIRECCIONES, no su sentido: -v es tan valido como v.
MATLAB y numpy pueden elegir signos distintos, y eso no es un error de ninguno
de los dos. Consecuencias practicas:

  * el cubo extraido puede salir reflejado respecto al de MATLAB. BV/TV, BS,
    Tb.Th, Tb.Sp, Tb.N, Conn.D y el DA son invariantes ante reflexion, asi que
    las metricas coinciden igual.
  * lo que SI puede intercambiarse es la etiqueta proximal/distal, porque los
    tercios se definen a lo largo de PC1. Por eso `vois_por_tercios` devuelve
    los tercios ORDENADOS a lo largo del eje y avisa de que quien nombra los
    extremos es la anatomia, no el algoritmo.

Aqui se fija ademas un convenio de signo determinista (la componente de mayor
magnitud de cada eje se hace positiva) para que dos ejecuciones de este codigo
den siempre lo mismo.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Lectura y segmentacion de la pila
# ---------------------------------------------------------------------------

def _orden_natural(rutas):
    """Ordena por el numero del nombre, no alfabeticamente.

    `..._rec00000999.tif` y `..._rec00001000.tif` se ordenan bien alfabetica-
    mente solo mientras el relleno de ceros sea constante. Si alguna vez no lo
    es, el volumen sale con las rebanadas barajadas y nada avisa.
    """
    def clave(p):
        n = re.findall(r"\d+", Path(p).stem)
        return (int(n[-1]) if n else 0, str(p))
    return sorted(rutas, key=clave)


def listar_rebanadas(carpeta, patron="*.tif"):
    rutas = _orden_natural(str(p) for p in Path(carpeta).glob(patron))
    if not rutas:
        raise FileNotFoundError(f"{carpeta}: ningun archivo coincide con {patron}")
    return rutas


def _leer_imagen(ruta):
    """Lee una rebanada, con Pillow de reserva.

    Las `BW_*.tif` que escribe MATLAB con `imwrite` sobre una imagen logica son
    de 1 bit con compresion CCITT, y `tifffile` —que es el lector que usa
    scikit-image para TIFF— no la soporta: lanza
    `<COMPRESSION.CCITTRLE: 2> not supported`.

    Pillow si la lee. Se intenta primero scikit-image porque maneja mejor los
    16 bits de las reconstrucciones de NRecon, que es el otro caso de uso, y se
    cae a Pillow solo cuando hace falta.
    """
    from skimage import io as skio
    try:
        im = skio.imread(ruta)
    except Exception:
        from PIL import Image
        im = np.array(Image.open(ruta))
    if im.ndim == 3:
        im = im[..., 0]
    return im


def apilar(carpeta, patron="*.tif", umbral=None, progreso=None):
    """Apila las rebanadas en una mascara booleana (nx, ny, nz).

    `umbral` en las unidades del propio archivo (16 bits para las
    reconstrucciones de NRecon). Si es None se supone que las imagenes ya son
    binarias y se toma > 0, que es el caso de las `BW_*.tif` que produce
    `segmentacion_por_lote.m`.

    Se lee rebanada a rebanada y se escribe directamente en el volumen ya
    reservado: apilar una lista de 438 imagenes de 740x656 duplicaria el pico
    de memoria justo cuando el volumen ya ocupa lo suyo.
    """
    rutas = listar_rebanadas(carpeta, patron)
    im0 = _leer_imagen(rutas[0])
    ny, nx = im0.shape
    vol = np.zeros((nx, ny, len(rutas)), dtype=bool)

    for k, r in enumerate(rutas):
        if progreso and k % 25 == 0:
            progreso(k, len(rutas), Path(r).name)
        im = _leer_imagen(r)
        if im.shape != (ny, nx):
            raise ValueError(
                f"{Path(r).name}: forma {im.shape} distinta de {(ny, nx)}. "
                f"MATLAB reescalaba en silencio; aqui se para, porque "
                f"reescalar una rebanada cambia la morfometria sin avisar.")
        # .T lleva (fila, columna) = (y, x) a (x, y): dimension 1 = X.
        vol[:, :, k] = (im > 0).T if umbral is None else (im > umbral).T
    if progreso:
        progreso(len(rutas), len(rutas), "listo")
    return vol


# ---------------------------------------------------------------------------
# Pila TIFF -> volumen con su tamano fisico real
# ---------------------------------------------------------------------------

# En milimetros, que es la unidad interna de TODO spinpy: los spacing van en
# mm/voxel, Tb.Th sale en mm, TV en mm^3 y Conn.D en 1/mm^3. Se acepta la
# entrada en la unidad que traiga el escaner y se convierte aqui, una sola vez.
# Cambiar la unidad interna en vez de convertir a la entrada obligaria a tocar
# todas las metricas y todos los resultados publicados.
UNIDADES_MM = {"um": 1e-3, "µm": 1e-3, "micras": 1e-3, "micron": 1e-3,
               "mm": 1.0, "cm": 10.0, "m": 1000.0, "nm": 1e-6}

# Claves de tamano de voxel en el `*_rec.log` de SkyScan/NRecon, de mas a menos
# fiable. "Image Pixel Size" es el del volumen reconstruido; "Camera Pixel Size"
# es el del detector y NO sirve — son 75.0 um frente a los 51.489 reales en H4.
_CLAVES_LOG = ("Image Pixel Size (um)", "Scaled Image Pixel Size (um)",
               "Pixel Size (um)")


def tam_voxel_desde_log(origen):
    """Busca el tamano de voxel en un `*_rec.log` de SkyScan. Devuelve mm.

    `origen` puede ser el propio .log, la carpeta de rebanadas o su carpeta
    padre: las reconstrucciones de NRecon dejan el log junto a las imagenes o
    un nivel por encima, y en este proyecto esta arriba (`H4/H 1-4_rec.log`
    con las rebanadas en `H4/` y las binarizadas en `H4/Segmentadas/`).

    Devuelve `(mm, info)` o `(None, info)` si no se encuentra. Nunca inventa un
    valor por defecto: un tamano de voxel equivocado escala TODAS las metricas
    en longitud y no da ningun sintoma visible.
    """
    p = Path(origen)
    info = {"log": None, "clave": None, "linea": None, "buscado": []}

    candidatos = []
    if p.is_file() and p.suffix.lower() == ".log":
        candidatos = [p]
    else:
        base = p if p.is_dir() else p.parent
        for d in (base, base.parent):
            info["buscado"].append(str(d))
            candidatos += sorted(d.glob("*_rec.log")) + sorted(d.glob("*.log"))

    vistos = set()
    for log in candidatos:
        if log in vistos or not log.is_file():
            continue
        vistos.add(log)
        try:
            texto = log.read_text(encoding="latin-1", errors="ignore")
        except Exception:
            continue
        for clave in _CLAVES_LOG:
            for linea in texto.splitlines():
                if linea.strip().startswith(clave):
                    try:
                        um = float(linea.split("=", 1)[1].strip())
                    except (IndexError, ValueError):
                        continue
                    if um > 0:
                        info.update({"log": str(log), "clave": clave,
                                     "linea": linea.strip()})
                        return um * UNIDADES_MM["um"], info
    return None, info


def _forma_imagen(ruta):
    """Forma (filas, columnas) sin leer los pixeles.

    Pillow abre de forma perezosa: `.size` sale de la cabecera y no cuesta
    nada. Hacen falta las formas de TODOS los candidatos antes de decidir
    cuales son rebanadas, y leerlos enteros para eso seria absurdo.
    """
    try:
        from PIL import Image
        with Image.open(str(ruta)) as im:
            w, h = im.size          # Pillow da (ancho, alto)
        return (h, w)
    except Exception:
        try:
            return _leer_imagen(ruta).shape
        except Exception:
            return None


def _puede_leer(ruta):
    try:
        _leer_imagen(ruta)
        return True
    except Exception:
        return False


def _separar_rebanadas(rutas):
    """Se queda con la serie de rebanadas y aparta lo que no lo es.

    POR QUE HACE FALTA. Los escaneres no dejan solo rebanadas en la carpeta de
    reconstruccion. En `H4/` conviven 438 rebanadas de 664x628 con ocho
    archivos que NO lo son: dos proyecciones `_pp1/_pp2` de 200x936, cinco
    previsualizaciones `_prev_*` de 936x936 y la hoja `_rec_spr` de 1944x936.

    Filtrar por nombre no sirve: `_prev_00000235` tambien lleva numero de
    rebanada, y cada fabricante inventa sus sufijos. Lo que si distingue a la
    serie real es que es la MAYORIA y comparte forma. Asi que se agrupa por
    forma y se toma el grupo mas numeroso.

    Lo apartado se DEVUELVE con nombre y forma, no se descarta en silencio:
    quitar rebanadas cambia la morfometria igual que reescalarlas, y quien
    carga tiene que poder ver que se dejo fuera.
    """
    # La cabecera es barata pero no siempre es de fiar: Pillow interpreta el
    # ultimo eje de tamano 3 o 4 como canales RGB/RGBA, de modo que una imagen
    # de 4 columnas le sale como (ancho, 1) en modo RGBA. Con rebanadas reales
    # de micro-CT no pasa —nunca miden 3 o 4 pixeles de ancho—, pero basta un
    # caso para leer mal toda una carpeta. Se contrasta la cabecera del primer
    # archivo contra su lectura completa y, si discrepan, se leen todas
    # enteras: mas lento y correcto antes que rapido y equivocado.
    sonda = _forma_imagen(rutas[0])
    try:
        real = _leer_imagen(rutas[0]).shape
    except Exception:
        real = None
    fiable = real is not None and sonda == real
    medir = _forma_imagen if fiable else (
        lambda r: (_leer_imagen(r).shape if _puede_leer(r) else None))

    formas = {r: medir(r) for r in rutas}
    validas = {r: f for r, f in formas.items() if f is not None}
    if not validas:
        raise ValueError("Ninguno de los archivos se pudo leer como imagen.")

    conteo = {}
    for f in validas.values():
        conteo[f] = conteo.get(f, 0) + 1
    forma_serie = max(conteo, key=conteo.get)

    serie = [r for r in rutas if validas.get(r) == forma_serie]
    apartados = [{"nombre": Path(r).name,
                  "forma": (list(formas[r]) if formas[r] else None)}
                 for r in rutas if validas.get(r) != forma_serie]

    # Si la "mayoria" no es tal, la carpeta no contiene una serie clara y
    # adivinar seria peor que parar.
    if len(serie) < 0.5 * len(rutas):
        detalle = ", ".join(f"{f[0]}x{f[1]}: {n}" for f, n in
                            sorted(conteo.items(), key=lambda kv: -kv[1]))
        raise ValueError(
            f"No hay una serie de rebanadas dominante en la carpeta "
            f"({len(rutas)} archivos, formas -> {detalle}). Afina el patron "
            f"para seleccionar una sola serie.")
    return serie, list(forma_serie), apartados


def _leer_pila_multipagina(ruta):
    """Lee un TIFF multipagina como (nx, ny, nz), dimension 1 = X."""
    from skimage import io as skio
    try:
        pila = skio.imread(str(ruta))
    except Exception:
        from PIL import Image, ImageSequence
        im = Image.open(str(ruta))
        pila = np.stack([np.array(p) for p in ImageSequence.Iterator(im)])
    pila = np.asarray(pila)
    if pila.ndim == 4:
        pila = pila[..., 0]
    if pila.ndim == 2:
        pila = pila[None, ...]
    # skimage entrega (z, fila, columna) = (z, y, x); a (x, y, z)
    return np.transpose(pila, (2, 1, 0))


def _umbral_otsu(vol, clases=3):
    """Otsu sobre un submuestreo, para no duplicar el pico de memoria.

    `clases=3` es el valor por defecto y NO es un capricho. En una
    reconstruccion de micro-CT el campo de vision es aire en su mayor parte, y
    hay tres poblaciones, no dos: aire, medula/tejido blando y hueso. Otsu de
    dos clases separa el aire de TODO LO DEMAS, de modo que la medula acaba
    contada como hueso.

    Medido sobre la pila cruda de H4, comparando contra la segmentacion manual
    del proyecto (BV/TV global 0.0657; tercios 0.282 / 0.546 / 0.766):

        2 clases (umbral 12358)  ->  0.0889 ;  0.612 / 0.791 / 0.938
        3 clases (umbral 18983)  ->  0.0717 ;  0.351 / 0.588 / 0.806

    Ninguno reproduce la segmentacion manual —con 3 clases el tercio menos
    denso todavia sale un 25 % alto—, y eso NO es un defecto de la
    implementacion: es la razon por la que la segmentacion es la mayor fuente
    de incertidumbre de todo el proceso. Un +-15 % de umbral mueve Tb.Th un
    62 %. Por eso el valor se reporta siempre y se puede fijar a mano.

    Con `clases=2` se recupera el Otsu clasico, que es el correcto cuando el
    volumen ya viene recortado al hueso y solo hay dos poblaciones.
    """
    v = np.asarray(vol[::2, ::2, ::2] if vol.size > 8_000_000 else vol)
    if int(clases) <= 2:
        from skimage.filters import threshold_otsu
        return float(threshold_otsu(v))

    # El casteo a float NO es cosmetico. `threshold_multiotsu` histograma con
    # `source_range='image'`, y para imagenes ENTERAS scikit-image ignora
    # `nbins` y usa el rango completo del tipo: 65536 bins para uint16. La
    # busqueda es O(bins^2), o sea 4.3e9 combinaciones — medido, 393 s sobre
    # un volumen de 14 400 voxeles. En float si respeta `nbins`, y con 256 el
    # umbral sale igual (18983 sobre la pila cruda de H4) en menos de un
    # segundo.
    from skimage.filters import threshold_multiotsu
    umbrales = threshold_multiotsu(v.astype(np.float64), classes=int(clases),
                                   nbins=256)
    return float(np.asarray(umbrales)[-1])


def leer_pila_tiff(origen, patron="*.tif", umbral=None, tam_voxel=None,
                   unidad="um", progreso=None, clases_otsu=3):
    """Carga una pila TIFF de micro-CT como mascara booleana con su escala real.

    Es la puerta que faltaba: hasta ahora solo se podia entrar por un .vtk o un
    .mat que alguien ya habia extraido. Devuelve exactamente lo mismo que
    `leer_vtk_voi` —`(BW, spacing)` con dimension 1 = X y `spacing` en mm— de
    modo que TODO lo que hay aguas abajo (morfometria, ajuste, homogeneizacion,
    exportacion) funciona sin cambiar una linea.

    origen : carpeta con las rebanadas, o un unico TIFF multipagina.
    umbral : None intenta deducirlo. Si el volumen solo tiene dos valores ya es
        binario y se toma > minimo; si no, se aplica Otsu y se DECLARA en
        `info`, porque el umbral es la mayor fuente de incertidumbre de todo el
        proceso (mover la segmentacion un +-15 % mueve Tb.Th un 62 %).
    tam_voxel : tamano de voxel en `unidad`. Si es None se busca en el
        `*_rec.log` de SkyScan. Si tampoco aparece, se LANZA excepcion en vez
        de suponer 1: un tamano inventado escala en silencio toda la
        morfometria y la rigidez.
    unidad : 'um', 'mm', 'cm', 'nm' o 'm'. Se convierte a mm internamente.

    Devuelve (BW, spacing_mm, info).
    """
    p = Path(origen)
    info = {"origen": str(p), "multipagina": False, "n_rebanadas": 0,
            "umbral": None, "umbral_metodo": None, "dtype": None,
            "tam_voxel_mm": None, "tam_voxel_origen": None, "log": None}

    # --- volumen crudo ----------------------------------------------------
    if p.is_file():
        info["multipagina"] = True
        vol = _leer_pila_multipagina(p)
        info["n_rebanadas"] = int(vol.shape[2])
    else:
        candidatas = listar_rebanadas(p, patron)
        rutas, forma, apartados = _separar_rebanadas(candidatas)
        info["n_rebanadas"] = len(rutas)
        info["apartados"] = apartados
        info["forma_rebanada"] = forma

        ny, nx = forma
        im0 = _leer_imagen(rutas[0])
        vol = np.zeros((nx, ny, len(rutas)), dtype=im0.dtype)
        for k, r in enumerate(rutas):
            if progreso and k % 25 == 0:
                progreso(k, len(rutas), Path(r).name)
            im = _leer_imagen(r) if k else im0
            if im.shape != (ny, nx):
                # No deberia ocurrir: la forma ya se comprobo por cabecera. Si
                # pasa, la cabecera y el dato no coinciden y el archivo esta
                # corrupto — parar sigue siendo lo correcto.
                raise ValueError(
                    f"{Path(r).name}: forma {im.shape} distinta de {(ny, nx)} "
                    f"pese a que la cabecera decia lo contrario. El archivo "
                    f"parece corrupto.")
            vol[:, :, k] = im.T
        if progreso:
            progreso(len(rutas), len(rutas), "listo")
    info["dtype"] = str(vol.dtype)

    # --- binarizacion -----------------------------------------------------
    if umbral is None:
        unicos = np.unique(vol[::4, ::4, ::4])
        if unicos.size <= 2 or vol.dtype == bool:
            umbral = float(unicos.min()) if unicos.size else 0.0
            info["umbral_metodo"] = "ya binaria (> minimo)"
        else:
            umbral = _umbral_otsu(vol, clases_otsu)
            info["umbral_metodo"] = (
                f"Otsu automatico, {int(clases_otsu)} clases"
                + (" (aire/medula/hueso)" if int(clases_otsu) >= 3 else ""))
    else:
        info["umbral_metodo"] = "indicado por el usuario"
    info["umbral"] = float(umbral)
    BW = vol > umbral
    del vol

    # --- escala fisica ----------------------------------------------------
    if tam_voxel is None:
        mm, info_log = tam_voxel_desde_log(p)
        info["log"] = info_log
        if mm is None:
            raise ValueError(
                "No se encontro el tamano de voxel. Indica `tam_voxel` con su "
                "unidad, o deja junto a las imagenes el `*_rec.log` del "
                f"escaner. Se busco en: {', '.join(info_log['buscado'])}")
        info["tam_voxel_origen"] = f"log de SkyScan ({info_log['clave']})"
    else:
        u = str(unidad).strip().lower()
        if u not in UNIDADES_MM:
            raise ValueError(f"Unidad no reconocida: {unidad!r}. "
                             f"Use una de {sorted(set(UNIDADES_MM))}.")
        mm = float(tam_voxel) * UNIDADES_MM[u]
        if not np.isfinite(mm) or mm <= 0:
            raise ValueError(f"Tamano de voxel invalido: {tam_voxel} {unidad}")
        info["tam_voxel_origen"] = f"indicado por el usuario ({tam_voxel} {unidad})"

    info["tam_voxel_mm"] = float(mm)
    spacing = np.array([mm, mm, mm], dtype=float)
    info["tam_fisico_mm"] = (np.array(BW.shape, float) * spacing).tolist()
    info["BVTV"] = float(BW.mean())
    return BW, spacing, info


# ---------------------------------------------------------------------------
# Marco anatomico por PCA
# ---------------------------------------------------------------------------

def marco_pca(mask, spacing):
    """Ejes principales de la nube de voxeles oseos, en mm.

    Devuelve (ejes, centro, proyecciones) con `ejes` de 3x3 y los vectores
    propios EN COLUMNAS, ordenados de mayor a menor varianza. `proyecciones`
    son las coordenadas de cada voxel oseo en ese marco.

    OJO — esto NO es el tensor MIL ni sustituye al DA. Es la forma EXTERNA del
    hueso (donde esta el eje largo de la pieza), no la orientacion de las
    trabeculas. Confundir las dos es exactamente el error que la correccion del
    proyecto prohibe: el DA sale del MIL, nunca de la covarianza de la nube.
    Aqui la PCA se usa solo para orientar el recorte del VOI.
    """
    mask = np.asarray(mask, bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    idx = np.argwhere(mask)              # (x, y, z), dimension 1 = X
    if idx.shape[0] < 10:
        raise ValueError("la mascara no tiene hueso suficiente para una PCA")
    pts = idx * spacing[None, :]
    centro = pts.mean(axis=0)
    _, _, Vt = np.linalg.svd(pts - centro, full_matrices=False)
    ejes = Vt.T                          # columnas = componentes principales

    # Convenio de signo determinista: la componente de mayor magnitud de cada
    # eje se hace positiva. La PCA no fija el sentido, y sin esto dos corridas
    # sobre datos casi iguales pueden devolver ejes opuestos.
    for j in range(3):
        k = int(np.argmax(np.abs(ejes[:, j])))
        if ejes[k, j] < 0:
            ejes[:, j] *= -1

    return ejes, centro, (pts - centro) @ ejes


def centros_por_tercios(proyecciones, n=3):
    """Centro de cada tercio a lo largo del eje principal, en el marco PCA.

    Los limites se ponen sobre el RANGO de la proyeccion, no sobre percentiles.
    Es lo que hace `VOIsTIFF.m` y hay que conservarlo: con percentiles los tres
    tercios tendrian el mismo numero de voxeles y sus centros caerian en otro
    sitio, de modo que los VOIs no serian los mismos.
    """
    p = np.asarray(proyecciones, float)
    eje = p[:, 0]
    lo, hi = float(eje.min()), float(eje.max())
    cortes = [lo + (hi - lo) * i / n for i in range(n + 1)]
    centros = []
    for i in range(n):
        m = (eje > cortes[i]) if i else (eje <= cortes[1])
        if i:
            m &= eje <= cortes[i + 1]
        if not m.any():
            raise ValueError(f"el tercio {i} quedo vacio")
        centros.append(p[m].mean(axis=0))
    return np.array(centros)


def extraer_cubo(mask, spacing, centro_pca, ejes, centro, lado_mm=5.0,
                 lado_vox=None):
    """Recorta un cubo orientado segun los ejes principales.

    El cubo se define en el marco de la PCA y se lleva de vuelta a indices de
    voxel por vecino mas proximo, igual que `VOIsTIFF.m`. El lado se fuerza a
    IMPAR para que exista un voxel central exacto.

    Los puntos que caen fuera del volumen se toman como fondo. No es lo mismo
    que recortarlos: si un VOI asoma por el borde, quedaria hueso artificial-
    mente ausente y el BV/TV bajaria sin que nada lo indique. Por eso se
    devuelve tambien `fuera`, la fraccion del cubo que cayo fuera, y quien
    llama decide si ese VOI es utilizable.
    """
    mask = np.asarray(mask, bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)

    if lado_vox is None:
        lado_vox = int(round(lado_mm / spacing[0]))
    if lado_vox % 2 == 0:
        lado_vox += 1
    h = (lado_vox - 1) // 2

    r = (np.arange(lado_vox) - h)
    dx, dy, dz = np.meshgrid(r * spacing[0], r * spacing[1], r * spacing[2],
                             indexing="ij")
    delta = np.stack([dx.ravel(), dy.ravel(), dz.ravel()], axis=1)

    # marco PCA -> mm -> indices de voxel
    pts_mm = (delta + np.asarray(centro_pca, float)) @ ejes.T + centro
    idx = np.rint(pts_mm / spacing[None, :]).astype(np.int64)

    dentro = np.ones(idx.shape[0], bool)
    for k in range(3):
        dentro &= (idx[:, k] >= 0) & (idx[:, k] < mask.shape[k])

    vals = np.zeros(idx.shape[0], bool)
    if dentro.any():
        i = idx[dentro]
        vals[dentro] = mask[i[:, 0], i[:, 1], i[:, 2]]

    cubo = vals.reshape(lado_vox, lado_vox, lado_vox)

    # ORDEN DE LOS EJES DEL CUBO — determinado empiricamente, no supuesto.
    # Comparando contra los VOI_*_PCAaligned.vtk de VOIsTIFF.m sobre H4, los
    # cubos coinciden en 912 673 de 912 673 voxeles en los tres VOIs, con
    # permutacion (1, 0, 2) y SIN reflexion: lo que cambia es que eje principal
    # cae en la primera dimension del array, no el contenido.
    #
    # Se transpone para reproducir el orden de MATLAB bit a bit. Da igual para
    # BV/TV y BS, que son invariantes, pero no para el DA: el MIL lanza rayos
    # sobre un conjunto fijo de direcciones y permutar los ejes cambia que
    # direccion atraviesa que, lo que movia el DA un 1-2 %.
    cubo = np.transpose(cubo, (1, 0, 2))
    return cubo, float(1.0 - dentro.mean())


def rango_pca(proyecciones):
    """Extremos del hueso en el marco PCA: (lo, hi) por eje, en mm.

    Es lo que hace falta para colocar un cubo a mano: `centro_en` mueve el
    recorte dentro de esta caja, y la interfaz la usa para acotar los mandos.
    """
    p = np.asarray(proyecciones, float)
    return p.min(axis=0), p.max(axis=0)


def centro_en(proyecciones, t, off2=0.0, off3=0.0, ventana=0.15):
    """Centro de recorte en el marco PCA, a la fraccion `t` del eje principal.

    `t` va de 0 a 1 sobre el RANGO de la proyeccion, el mismo criterio que
    `centros_por_tercios` (rango, no percentiles).

    EN PC1 manda `t` exactamente: es la posicion que ha pedido quien recorta.
    EN PC2 y PC3 el centro SIGUE AL HUESO — sale del centroide de los voxeles
    oseos que caen en una ventana de anchura `ventana` (fraccion del eje)
    alrededor de esa posicion, mas el desplazamiento que se pida. Sin eso, un
    hueso curvado dejaria el cubo fuera del material en los extremos, porque el
    eje PC1 pasa por el centroide global y no por el centro de cada seccion.

    NO reproduce exactamente `centros_por_tercios`, y no debe: aquella toma el
    centroide del tercio entero en las TRES coordenadas, asi que su PC1 tampoco
    es el punto medio del tercio. Aqui PC1 lo fija el usuario.
    """
    p = np.asarray(proyecciones, float)
    lo, hi = rango_pca(p)
    L = hi[0] - lo[0]
    c1 = lo[0] + float(np.clip(t, 0.0, 1.0)) * L

    # Si la ventana queda vacia se ensancha antes que devolver algo arbitrario.
    for w in (ventana, 2 * ventana, 4 * ventana):
        sel = np.abs(p[:, 0] - c1) <= 0.5 * w * L
        if sel.any():
            loc = p[sel].mean(axis=0)
            break
    else:
        loc = p.mean(axis=0)

    return np.array([c1, loc[1] + float(off2), loc[2] + float(off3)])


def vois_por_tercios(mask, spacing, lado_mm=5.0, n=3):
    """Los `n` VOIs cubicos alineados por PCA, uno por tercio del hueso.

    Port de `VOIsTIFF.m`. Devuelve una lista de dicts con el cubo, el centro,
    la fraccion fuera del volumen y el BV/TV, **ordenados a lo largo del eje
    principal**.

    NO se les pone nombre anatomico. El sentido de PC1 es arbitrario (ver la
    cabecera del modulo), asi que llamar «proximal» al primero seria acertar la
    mitad de las veces. Quien conoce la pieza decide, y para ayudar se devuelve
    `bvtv`: en estos sesamoideos la densidad crece de proximal a distal.
    """
    ejes, centro, proy = marco_pca(mask, spacing)
    centros = centros_por_tercios(proy, n)
    salida = []
    for i, c in enumerate(centros):
        cubo, fuera = extraer_cubo(mask, spacing, c, ejes, centro, lado_mm)
        salida.append({
            "indice": i,
            "cubo": cubo,
            "centro_pca": c,
            "fuera": fuera,
            "bvtv": float(cubo.mean()),
            "lado_vox": int(cubo.shape[0]),
            "spacing": np.atleast_1d(np.asarray(spacing, float)).ravel(),
        })
    return salida, {"ejes": ejes, "centro": centro}


# ---------------------------------------------------------------------------
# Escritura en VTK legacy, el formato que ya lee `leer_vtk_voi`
# ---------------------------------------------------------------------------

def escribir_vtk_voi(mask, spacing, ruta, comentario="VOI cubico"):
    """Escribe la mascara como STRUCTURED_POINTS ASCII, 0/255.

    Mismo formato que producen los scripts de MATLAB, para que los VOIs de las
    dos rutas sean intercambiables y `leer_vtk_voi` los lea sin distinguirlos.
    """
    mask = np.asarray(mask, bool)
    spacing = np.atleast_1d(np.asarray(spacing, float)).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    nx, ny, nz = mask.shape

    # order='F' porque VTK recorre X como indice mas rapido y la mascara tiene
    # dimension 1 = X.
    vals = np.where(mask.ravel(order="F"), 255, 0).astype(np.uint8)
    with open(ruta, "w", encoding="ascii", newline="\n") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(comentario[:200] + "\n")
        f.write("ASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        f.write("ORIGIN 0 0 0\n")
        f.write("SPACING %.6f %.6f %.6f\n" % tuple(spacing))
        f.write(f"POINT_DATA {vals.size}\n")
        f.write("SCALARS scalars unsigned_char\nLOOKUP_TABLE default\n")
        f.write("\n".join(str(int(v)) for v in vals))
        f.write("\n")
    return Path(ruta)
