"""
visor.py — Interfaz grafica del port: VOI real y spinodoide candidato, lado a lado.

ALCANCE — menor que el de AppFinal_V2, y solo lo que esta validado
------------------------------------------------------------------
Esta ventana expone lo que el nucleo portado respalda: generacion del
spinodoide, carga de un VOI, morfometria comparada, ajuste por busqueda
escalonada, tensor elastico homogeneizado, ensayo de compresion en uno o en
los tres ejes con campos de deformacion efectiva y de von Mises, y exportacion
como SOLIDO (malla hexaedrica o TET10) para Abaqus y ANSYS.

NO hay informes ni los otros cuatro optimizadores (Pareto, bayesiano, MOBO,
best-fit rapido): esa logica sigue viviendo solo en MATLAB. Un boton apagado
seria peor que un boton ausente.

SESION REPRODUCIBLE
-------------------
El menu Sesion guarda en JSON todo lo que determina el estado: parametros,
semilla, esquema de muestreo, resoluciones y la ruta del VOI. Como el generador
es estocastico pero esta sembrado, ese archivo regenera la MISMA realizacion bit
a bit. Sin el, un caso interesante que se pierde al mover un deslizador no se
recupera.

EXPORTAR: HEXAEDRICA O TET10
-----------------------------
La hexaedrica pone un elemento por voxel. Es EXACTA —reproduce el 100% del
volumen—, instantanea y no puede fallar al mallar, pero el elemento es lineal
y la superficie queda escalonada.

La TET10 suaviza la superficie y la tetraedraliza: elemento cuadratico y piel
lisa, a cambio de tiempo y de perder algo de material. Cuanto se pierde depende
de la RESOLUCION: medido sobre el mismo spinodoide, 88.8% del volumen a 40^3 y
93.8% a 48^3. Con trabeculas de uno o dos voxeles la superficie no tiene con
que trabajar, asi que conviene exportar a la resolucion mas alta que se pueda
pagar.

ORIENTACION IMPUESTA POR EL AJUSTE (correccion F4)
--------------------------------------------------
Los deslizadores de rotacion solo cubren +-90 grados y no pueden representar
todas las orientaciones de un VOI. Cuando el ajuste deriva una R del eje
principal, esa R manda sobre los deslizadores y se avisa en el panel; mover
cualquier rotacion la libera, porque ese gesto expresa que el usuario quiere
controlar la orientacion a mano.

DOS SUPERFICIES DE MEDIDA, Y HAY QUE DECLARAR CUAL
--------------------------------------------------
La seccion Morfometria tiene un desplegable "Superficie":

  * voxeles — marching cubes crudo sobre la mascara. Es el modo de la app de
    MATLAB, el que usa el ajuste, y el que sobreestima el area un ~8.5 %
    (medido sobre una esfera). Como VOI y candidato lo sufren igual, el
    sesgo se cancela al comparar.
  * malla — superficie cerrada, suavizada con Taubin y reparada con
    PyMeshFix, recortada medio voxel por dentro y sin las seis tapas
    (`spinpy.morphometry.morfometria_malla`). Quita ese sesgo: sobre el VOI
    proximal de H4, BS pasa de 412 a 377 mm2 (-8.7 %) y Tb.Th sube un 4.7 %.
    BV/TV sale ~1-3 % menor, que es la diferencia entre la superficie a
    nivel 0.5 y el conteo de cubos enteros, no una perdida: la perdida real
    de suavizar y reparar se mide (-0.4 % en ese VOI) y se ensena en la
    barra de estado.

Cambiar de modo vuelve a medir el VOI: una tabla con el VOI en un modo y el
candidato en otro mezcla dos definiciones de BS (correccion C4, otra vez).
El modo viaja en la sesion y en `_res["morfometria"]["modo"]`, para que
cualquier tabla exportada diga con que superficie se midio.

EL PANEL VA POR ETAPAS DEL TRABAJO
-----------------------------------
El panel izquierdo esta dividido en secciones con cabecera propia, en el orden
en que se usan: VOI de referencia, Spinodoide, Visualizacion, Morfometria,
Ajuste al VOI, Analisis mecanico, Exportacion, y Lote y varianza. Antes todo
colgaba del mismo nivel: los grupos existian, pero no se veia de un vistazo en
que etapa estaba cada control ni donde correspondia anadir uno nuevo.

Cada seccion se abre con `_seccion`, que dibuja la raya separadora, el titulo y
la raya bajo el titulo, y devuelve el layout donde va el contenido. Anadir una
seccion es una linea; anadir un control es meterlo en el layout que devuelve.

DOS RESOLUCIONES, QUE ES LO QUE HACE USABLE LA INTERFAZ
-------------------------------------------------------
Generar y medir a 64^3 cuesta ~9 s en esta maquina, lo que vuelve inservible un
deslizador. Aqui se separan:

  * Resolucion de VISTA   — se recalcula al mover un deslizador. A 32^3 con
    pocas ondas baja a decimas de segundo y sirve para explorar la forma.
  * Resolucion de MEDIDA  — solo se usa al pulsar "Medir". Es la que produce
    numeros comparables con el VOI, y es la que hay que citar.

Ambas usan la MISMA semilla, asi que muestrean el mismo campo aleatorio a
distinta finura: la vista previa es una version gruesa de lo que se medira, no
otra estructura.

El esquema de muestreo de ondas es un desplegable y no tiene valor por defecto
silencioso: 'rechazo' (GIBBON/AppFinal_V2) y 'equitativo' (TPMS-Scaffolds-
generator) producen anisotropias distintas para los mismos thetas cuando los
conos tienen angulos desiguales. Ver README.

DOS IDIOMAS: ESPANOL E INGLES
------------------------------
El menu **Idioma** cambia toda la interfaz sin reconstruir la ventana ni perder
nada: el VOI cargado, las metricas medidas y los campos de deformacion siguen
donde estaban. Son minutos u horas de calculo y no se pueden tirar por cambiar
de idioma.

COMO. `spinpy.idioma.capturar()` recorre el arbol de widgets UNA vez al
arrancar y anota el texto ORIGINAL de cada uno; `aplicar()` los reescribe
traducidos. Como el original queda guardado, ir y volver entre idiomas es
idempotente. Lo que ese recorrido no ve -mensajes de la barra de estado,
cuadros de dialogo, informes en HTML- lleva `_()` explicito donde se construye.

Los textos cuyo valor cambia en marcha -la etiqueta de un deslizador dice
"Densidad relativa: 35 %"- se marcan con la propiedad `i18n_dinamico` para que
el recorrido los salte, y se retraducen ellos solos.

EL DICCIONARIO VA INDEXADO POR LA CADENA ESPANOLA (`spinpy/idioma_textos.py`).
Eso mantiene el codigo legible y hace que un texto sin traducir salga en
espanol en vez de romperse, pero tiene un filo: cambiar una tilde en un texto
espanol y no tocar el diccionario pierde la traduccion inglesa SIN NINGUN
ERROR. Para eso esta `idioma_revisar.py`, que construye la ventana de verdad,
recorre sus widgets, saca del codigo las cadenas envueltas en `_()` y lista lo
que falta, lo que sobra y los marcadores {} que no cuadran entre las dos
versiones. Hay que ejecutarlo despues de tocar cualquier texto:

    python idioma_revisar.py              (resumen; devuelve 1 si falta algo)
    python idioma_revisar.py --esqueleto  (deja las claves listas para pegar)

QUE NO SE TRADUCE, A PROPOSITO:
  * los nombres de metrica -BV/TV, Tb.Th, Tb.Sp, Tb.N, BS/BV, DA, SMI,
    Conn.D-, que son notacion internacional (Parfitt / Bouxsein et al. 2010);
  * las unidades y los simbolos (mm, MPa, GPa, sigma, epsilon, E_app);
  * las citas bibliograficas;
  * el informe de `--autocomprobacion`, que es una herramienta de diagnostico
    para enviar por correo, no parte de la interfaz.

La eleccion se guarda en QSettings("spinpy", "visor") y se restaura al abrir.

LA FABRICA SE DIBUJA, NO SOLO SE TABULA
----------------------------------------
La casilla "Mostrar fabrica MIL" pinta en cada panel el elipsoide del tensor
MIL en alambre y, en magenta, su eje mayor: la direccion en la que la
estructura es mas continua y por tanto mas rigida. En el panel del spinodoide
se rotula ademas el angulo entre ese eje y el del VOI.

POR QUE HACIA FALTA. El tensor MIL se calcula desde el principio, pero de el
solo se ensenaban dos escalares: DA y DA2. Dos estructuras pueden tener el
MISMO DA y estar orientadas a 90 grados una de otra, y la tabla no las
distingue. En el ajuste del VOI proximal de H4 el candidato ganador resulto
tener su eje rigido a 87 grados del hueso, con lo que su modulo aparente caia
de 633 a 11 MPa; eso se supo despues de 477 s de ajuste y de leer el dialogo
del analisis comparado. Dibujado, se ve antes de ajustar nada.

El rotulo pasa a rojo si el angulo supera 30 grados, y avisa aparte cuando
DA2 < 1.06: por debajo de ese umbral la fabrica es un PLANO y no un eje, la
flecha apunta a una direccion cualquiera de ese plano, y ninguna rotacion
convierte esa estructura en axial. Es el mismo criterio con el que la
correccion L1 (`spinpy/fit.py`) decide abstenerse de alinear.

CONVENCION, que es donde se equivoca uno. El ajuste del tensor es
1/MIL(n)^2 = n' M n, de modo que el MIL en la direccion del autovector i vale
1/sqrt(lambda_i): el autovalor MENOR corresponde al MIL MAYOR. El elipsoide se
dibuja con semiejes proporcionales a 1/sqrt(lambda), no a lambda, que lo
pintaria aplastado justo en la direccion en la que la estructura es continua.
El semieje mayor se lleva al 62 % del lado y las flechas al 85 %, para que
asomen por las caras: inscritos quedarian enterrados dentro de una estructura
opaca y no se veria nada.

VIGENCIA. El elipsoide del spinodoide describe la ULTIMA MEDIDA, no la posicion
de los deslizadores. Al regenerar la vista desaparece hasta que se vuelva a
medir, en vez de quedarse pegado: una flecha de la estructura anterior sobre la
geometria nueva es indistinguible a la vista de un resultado valido.

RENDER: SOMBRAS DE CONTACTO Y FONDO
------------------------------------
Una trabecula es una marana de barras cruzadas. Sin sombra donde una pasa por
detras de otra, la profundidad se lee solo por el escorzo y el conjunto parece
plano. La casilla "Sombras de contacto (SSAO)" activa oclusion ambiental de
pantalla; el radio NO es fijo, sino una fraccion de la diagonal de la caja,
porque en milimetros absolutos un VOI de 5 mm y un cubo unidad de 1 mm darian
resultados opuestos. Se anade tambien suavizado SSAA -las aristas de la
isosuperficie salian dentadas en las capturas- y un punto de especularidad en
el material.

Todo esto es APARIENCIA: no toca ninguna medida, y las tres casillas se pueden
desmarcar si el driver de video lo dibuja mal o si el equipo va lento.

El fondo por defecto es un degradado suave, que da a la silueta algo contra lo
que recortarse. "Fondo blanco (para figuras)" devuelve el blanco puro, que es
el correcto para una figura de tesis.

Uso:  python visor.py
"""

from __future__ import annotations

import html
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pyvista as pv
from PyQt5 import QtCore, QtGui, QtWidgets
from pyvistaqt import QtInteractor

sys.path.insert(0, str(Path(__file__).parent))
from spinpy import procedencia  # noqa: E402
from spinpy import (ajustar_dual_lattice, ajustar_spinodoide,  # noqa: E402
                    constantes_ingenieria, euler_R, generar_dual_lattice,
                    generar_mascara, homogeneizar, leer_voi, morfometria)
from spinpy.elastic import remuestrear_bw   # noqa: E402
from spinpy.estadistica import informe      # noqa: E402
from spinpy.lote import (METRICAS, correr_lote,  # noqa: E402
                         dispersion_semillas)
from spinpy.metodos import COMPARADAS as COMPARADAS_MET  # noqa: E402
from spinpy.metodos import REGISTRO                       # noqa: E402
from spinpy.metodos import REGISTROS                      # noqa: E402
from spinpy.metodos import comparar_resultados            # noqa: E402
# `_` es la funcion de traduccion en todo este modulo. NO usarla como
# variable descartable: hacerlo la convierte en local de esa funcion y
# cualquier `_("...")` anterior revienta con UnboundLocalError. Para los
# descartes se usa `_x`.
from spinpy.idioma import (IDIOMAS, _, aplicar, capturar,  # noqa: E402
                           fijar_idioma)
from spinpy.idioma import idioma as idioma_actual             # noqa: E402
from spinpy.espesor import (espesor_local, muestrear_en_puntos,  # noqa: E402
                            estadisticas as estadisticas_esp)
from spinpy.resistencia import (APOYOS, EJES, EPS_CRITICA,  # noqa: E402
                                _solo_portante,
                                criterio_pistoia, ensayo_compresion,
                                cuantiles_vm_superficie,
                                ensayo_compresion_eje, estadisticos_vm,
                                estudio_convergencia)
from spinpy.escribe import (escribir_abaqus, escribir_apdl,  # noqa: E402
                            escribir_febio, escribir_stl, escribir_vtu)
from spinpy.morphometry import morfometria_malla, tensor_mil  # noqa: E402
from spinpy.avisos import desalineacion, voi_no_trabecular   # noqa: E402
from spinpy import informe as informe_pub                   # noqa: E402
from spinpy.informe import (etiqueta_item, sha256_archivo,   # noqa: E402
                            texto_aviso, texto_motivos)
from spinpy.simulacion import PROTOCOLOS as PROTOCOLOS_SIM   # noqa: E402
from spinpy.simulacion import TBTH_H_MIN                     # noqa: E402
from spinpy.simulacion import fallo_progresivo, simular_perdida  # noqa: E402
from spinpy.solido import (malla_hex, malla_tet10,           # noqa: E402
                           volumen_hex)
from dialogo_febio import (DialogoFEBio, DialogoFEMAuto,     # noqa: E402
                           DialogoResultadosFEBio, escribir_csv,
                           figura_validacion)

RAIZ = Path(__file__).parent


def _carpeta_documentos():
    """Ruta real de "Documentos", preguntandosela a Windows.

    No vale `Path.home() / "Documents"`. En un Windows en espanol la carpeta
    del disco sigue llamandose `Documents` -"Documentos" es solo el nombre que
    ensena el explorador-, pero OneDrive la REDIRIGE con frecuencia a
    `~/OneDrive/Documentos`, y entonces la ruta ingenua apunta a una carpeta
    que no existe o, peor, que existe y nadie mira. Se pregunta por el
    identificador de carpeta conocida, que es lo unico que sobrevive a la
    redireccion, y solo si eso falla se prueban los nombres a mano.
    """
    try:
        import ctypes
        import ctypes.wintypes as wt
        from uuid import UUID

        class GUID(ctypes.Structure):
            _fields_ = [("d1", wt.DWORD), ("d2", wt.WORD), ("d3", wt.WORD),
                        ("d4", ctypes.c_ubyte * 8)]

            def __init__(self, texto):
                super().__init__()
                u = UUID(texto)
                self.d1 = u.time_low
                self.d2 = u.time_mid
                self.d3 = u.time_hi_version
                for i, b in enumerate(u.bytes[8:]):   # los 8 ultimos bytes
                    self.d4[i] = b

        FOLDERID_Documents = GUID("FDD39AD0-238F-46AF-ADB4-6C85480369C7")
        p = ctypes.c_wchar_p()
        if ctypes.windll.shell32.SHGetKnownFolderPath(
                ctypes.byref(FOLDERID_Documents), 0, None,
                ctypes.byref(p)) == 0:
            d = Path(p.value)
            ctypes.windll.ole32.CoTaskMemFree(p)
            if d.exists():
                return d
    except Exception:
        pass
    for cand in (Path.home() / "Documents", Path.home() / "Documentos",
                 Path.home() / "OneDrive" / "Documentos",
                 Path.home() / "OneDrive" / "Documents"):
        if cand.exists():
            return cand
    return Path.home()


def _carpeta_datos():
    """Donde la aplicacion guarda y donde abre los dialogos por defecto.

    Sin empaquetar es la carpeta del codigo, que es lo que se ha usado siempre
    y lo que espera quien trabaja desde el repositorio.

    EMPAQUETADA es otra cosa, y no es un detalle: `__file__` apunta entonces a
    la carpeta de instalacion -normalmente bajo Archivos de programa, de solo
    lectura para un usuario sin privilegios- o, si se empaqueta en un solo
    archivo, al directorio temporal que Windows BORRA al cerrar la aplicacion.
    Guardar ahi no da error visible: el usuario exporta sus resultados, la
    aplicacion dice que los ha guardado, y no estan. Por eso todo lo que se
    escribe va a Documentos/spinpy.
    """
    if getattr(sys, "frozen", False):
        d = _carpeta_documentos() / "spinpy"
        try:
            d.mkdir(parents=True, exist_ok=True)
            return d
        except Exception:
            return Path.home()
    return RAIZ


def _ruta_icono():
    """El .ico de la aplicacion, este el codigo suelto o empaquetado.

    Empaquetado, PyInstaller deja los datos en `sys._MEIPASS` (la carpeta
    `_internal` en modo carpeta). Suelto, vive en `instalador/`. Se prueban
    los dos y se devuelve None si no esta: un icono que falta no puede
    impedir que la aplicacion abra.
    """
    cand = []
    empaquetado = getattr(sys, "_MEIPASS", None)
    if empaquetado:
        cand.append(Path(empaquetado) / "spinpy.ico")
    cand += [RAIZ / "instalador" / "spinpy.ico", RAIZ / "spinpy.ico"]
    for c in cand:
        try:
            if c.is_file():
                return c
        except Exception:
            pass
    return None


def _icono_de_clase(w):
    """Pone el .ico tambien como icono de la CLASE de ventana (solo Windows).

    `setWindowIcon` deja el icono en la VENTANA, y con eso la barra de titulo
    ya sale bien. La barra de tareas, en cambio, unas veces lee ese icono y
    otras el de la clase de ventana que Qt registro, que es `IDI_APPLICATION`
    -el rectangulo generico de Windows- porque el ejecutable que arranca es
    `python.exe`. De ahi que el fallo fuera intermitente: el shell cachea el
    icono POR AppUserModelID, asi que un arranque que leyo el generico deja
    generico tambien al siguiente.

    Medido en esta maquina, con un AppUserModelID nuevo en cada prueba para
    saltarse esa cache: sin esto, 2 de 4 arranques dan el icono generico; con
    esto, 9 de 9 dan el de spinpy. El icono de clase es del PROCESO y lo
    comparten todas las ventanas de Qt, pero se llama desde la portada y
    desde el visor porque la portada se puede desactivar.

    Silencioso a proposito: un icono es cosmetico y no puede impedir que la
    aplicacion abra.
    """
    if sys.platform != "win32":
        return
    ruta = _ruta_icono()
    if ruta is None:
        return
    try:
        import ctypes
        u = ctypes.windll.user32
        u.LoadImageW.restype = ctypes.c_void_p
        poner = getattr(u, "SetClassLongPtrW", None) or u.SetClassLongW
        poner.restype = ctypes.c_void_p
        poner.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
        # `winId()` fuerza la creacion de la ventana nativa: la clase tiene
        # que estar corregida ANTES de `show()`, porque es entonces cuando
        # Windows crea el boton de la barra y le busca icono.
        h = ctypes.c_void_p(int(w.winId()))
        IMAGE_ICON, LR_LOADFROMFILE = 1, 0x00000010
        GCLP_HICON, GCLP_HICONSM = -14, -34
        SM_CXICON, SM_CXSMICON = 11, 49
        for indice, metrica in ((GCLP_HICON, SM_CXICON),
                                (GCLP_HICONSM, SM_CXSMICON)):
            tam = u.GetSystemMetrics(metrica)
            ico = u.LoadImageW(None, str(ruta), IMAGE_ICON, tam, tam,
                               LR_LOADFROMFILE)
            if ico:
                poner(h, indice, ctypes.c_void_p(ico))
    except Exception:
        pass


# DATOS: donde se escribe.  RAIZ: donde vive el codigo.  Coinciden mientras se
# trabaja desde el repositorio y se separan al empaquetar.
DATOS = _carpeta_datos()

# Carpeta inicial del dialogo de VOIs. Fuera del repositorio no existe, y
# entonces el dialogo se abre en DATOS.
VOIDIR = RAIZ.parent / "H4" / "Segmentadas"

# Limites de la homogeneizacion. Los valores salen de tiempos MEDIDOS con el
# solver multigrid sobre un spinodoide, con la maquina en reposo (ver la tabla
# completa en spinpy/elastic.py):
#     16^3 -> 5.7 s     24^3 -> 34 s     32^3 -> 99 s
# El maximo es 32 porque a partir de ahi la espera dentro de una interfaz deja
# de ser razonable, no porque el metodo falle.
RES_HOMOG_MAX = 32
RES_HOMOG_DEF = 16

# Material base de la app (localBaseMaterial): matriz mineralizada.
E_S_PA, NU_S = 20e9, 0.30

def _fmt_num(v, spec, sufijo=""):
    """Numero formateado, o «—» si no hay numero (None, NaN, texto)."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return "—"
    return format(x, spec) + sufijo if np.isfinite(x) else "—"


def _n_efectivo_min(inf):
    """El menor N efectivo del informe de varianza, o NaN si no lo hay.

    Con un solo VOI no hay varianza entre especimenes que estimar y
    `estadistica.informe` deja `n_efectivo` en None. `np.nanmin` no sabe
    comparar None con None y reventaba justo al terminar el lote, despues de
    minutos de CPU y antes de mostrar nada: un lote de un VOI -el caso tipico
    para probar- no llegaba nunca a su dialogo.
    """
    if inf is None or not len(inf) or "n_efectivo" not in inf:
        return float("nan")
    v = np.array([np.nan if x is None else float(x)
                  for x in inf["n_efectivo"]], float)
    return float(np.nanmin(v)) if np.isfinite(v).any() else float("nan")


def mejor_ajuste(VOI, spacing, registro, claves, m_voi=None,
                 semilla=20260720, num_waves=700, informar=None):
    """El ajuste de MENOR ERROR entre los metodos `claves` de `registro`.

    Una sola definicion de «el mejor ajuste» para los dos sitios que la usan:
    el analisis comparado con «Ajustar antes» y el informe automatico. Si cada
    uno eligiera a su manera, el informe podria describir un ganador distinto
    del que el analisis comparado ensayo con el mismo VOI.

    Solo tiene sentido con metodos que comparten la definicion de error
    (`METODOS_COMPARABLES`): el desempate mecanico suma dos terminos y su
    cifra no se puede poner en la misma cola.

    Corre en el hilo que llama, en serie. `informar(i, n, etapa)` recibe el
    avance sobre el total de evaluaciones de todos los metodos. Devuelve
    `(ganador, hechos)`: el resultado ganador, con `_todos` resumiendo el
    error y el tiempo de cada metodo, y el dict metodo -> resultado entero.
    """
    mv = m_voi if m_voi is not None else morfometria(VOI, spacing)
    total = sum(registro[c]["n_eval"] for c in claves)
    hechos, base = {}, 0
    for clave in claves:
        met = registro[clave]

        def pr(i, n_tot, etapa, _d=base, _e=met["etiqueta"]):
            if informar is not None:
                informar(_d + int(i), total, f"{_e} — etapa {etapa}")

        hechos[clave] = met["correr"](VOI, spacing, m_voi=mv,
                                      semilla=semilla, num_waves=num_waves,
                                      progreso=pr)
        base += met["n_eval"]
    ganador = min(hechos.values(), key=lambda r: r["error"])
    ganador["_todos"] = {k: {"error": float(r["error"]),
                             "tiempo_s": float(r["tiempo_s"]),
                             "etiqueta": r["etiqueta"]}
                         for k, r in hechos.items()}
    return ganador, hechos


def opciones_figuras():
    """Estilos, vistas y suavizado de las figuras 3D del informe.

    Se guardan en QSettings al aceptar la ventana del informe automatico y los
    lee tambien el informe manual. Lo guardado se valida contra
    `figuras.opciones_3d`: una clave vieja o un estilo que ya no exista no
    rompen el informe, toman el valor por omision.
    """
    from spinpy import figuras as F
    try:
        crudo = QtCore.QSettings("spinpy", "visor").value("figuras_3d", "")
        guardadas = json.loads(crudo) if crudo else None
    except (TypeError, ValueError):
        guardadas = None
    return F.opciones_3d(guardadas)


def guardar_opciones_figuras(o):
    QtCore.QSettings("spinpy", "visor").setValue("figuras_3d", json.dumps(o))


def factores_tiempo():
    """({etapa: real/estimado}, numero de ejecuciones) de ESTE equipo.

    `spinpy.tiempos` esta medido en el equipo de desarrollo; cada etapa que
    termina en el informe automatico corrige su factor (`guardar_factor_
    tiempo`), asi que la estimacion se acerca a este equipo con el uso.
    """
    ajustes = QtCore.QSettings("spinpy", "visor")
    try:
        crudo = ajustes.value("tiempos/factores", "")
        factores = json.loads(crudo) if crudo else {}
    except (TypeError, ValueError):
        factores = {}
    try:
        n = int(ajustes.value("tiempos/ejecuciones", 0))
    except (TypeError, ValueError):
        n = 0
    return factores, n


def guardar_factor_tiempo(clave, real, estimado):
    """Corrige el factor de `clave` con una etapa terminada.

    Se descartan las etapas de menos de 20 s: ahi manda el coste fijo y el
    cociente es ruido. `estimado` ya incluia el factor anterior, asi que el
    nuevo factor es el anterior por el cociente de esta ejecucion.
    """
    from spinpy import tiempos
    if not estimado or real is None or real < 20:
        return
    factores, _n = factores_tiempo()
    anterior = factores.get(clave)
    base = estimado / (anterior or 1.0)
    factores[clave] = tiempos.actualizar_factor(anterior, real, base)
    QtCore.QSettings("spinpy", "visor").setValue(
        "tiempos/factores", json.dumps(factores))


# Metricas mostradas en la tabla, con su formato y unidad
TABLA = [
    ("BVTV",     "BV/TV",     "{:.4f}",  ""),
    ("PoTot",    "Po.tot",    "{:.2f}",  "%"),
    ("BSBV",     "BS/BV",     "{:.3f}",  "1/mm"),
    # BS/PV es la MISMA interfaz dividida por el volumen de poro en vez de por
    # el de hueso. En hueso se cita BS/BV (Parfitt); la literatura de
    # materiales porosos para implantes usa la razon referida al poro, que es
    # la que gobierna transporte y adhesion celular. Se dan las dos y no se
    # elige por el lector.
    ("BSPV",     "BS/PV",     "{:.3f}",  "1/mm"),
    ("BS",       "BS",        "{:.3f}",  "mm2"),
    # La misma area, separada. La interna es la interfaz de verdad -y es la
    # que entra en BS/BV y Tb.Th-; la externa son las seis caras por donde el
    # cubo corta trabeculas, que no son interfaz de nada. `%BS.ext` dice cuanto
    # del area medida es huella del corte: si sale grande, el VOI es pequeno
    # para la estructura que contiene.
    ("BS_interna", "BS int.", "{:.3f}",  "mm2"),
    ("BS_externa", "BS ext.", "{:.3f}",  "mm2"),
    ("BS_frac_externa", "%BS.ext", "{:.1%}", ""),
    ("TbTh",     "Tb.Th",     "{:.5f}",  "mm"),
    ("TbSp",     "Tb.Sp",     "{:.5f}",  "mm"),
    ("TbN",      "Tb.N",      "{:.4f}",  "1/mm"),
    ("DA",       "DA",        "{:.4f}",  ""),
    ("DA2",      "DA2",       "{:.4f}",  ""),
    ("FracPort", "Frac.port", "{:.4f}",  ""),
    # Conn.D y SMI solo se rellenan si se marca la casilla: cada una cuesta una
    # marching cubes o un etiquetado completo de mas, y ninguna entra en la
    # funcion de error del ajuste. Se muestran siempre como fila para que se
    # vea que existen y estan sin calcular, no para que parezcan ausentes.
    ("ConnD",    "Conn.D",    "{:.1f}",  "1/mm3"),
    ("SMI",      "SMI",       "{:.3f}",  ""),
    # Po.Dm tiene casilla propia porque es la mas cara de todas: una
    # transformada de distancia por radio sobre la fase PORO, que a BV/TV 0.3
    # es el 70 % del volumen. Medido: 0.7 s a 64^3, 3.2 s a 96^3, 12 s a 128^3.
    ("PoDm",     "Po.Dm",     "{:.5f}",  "mm"),
    ("PoDm_p95", "Po.Dm p95", "{:.5f}",  "mm"),
    # Ellipsoid Factor (Doube 2015): placa -1, barra +1. Casilla propia porque
    # es, con diferencia, la medida mas cara del programa (minutos, no
    # segundos), y se da la fraccion de hueso que es placa y la que es barra
    # porque la media sola esconde justo la mezcla que interesa.
    ("EF",       "EF",        "{:+.3f}", ""),
    ("EF_frac_placa", "EF placa", "{:.1%}", ""),
    ("EF_frac_barra", "EF barra", "{:.1%}", ""),
]

# Titulo de la barra de color, por indice de `cmb_color`. Va en una tabla y no
# en un condicional porque cada campo nuevo tiene que obligar a declarar su
# unidad: un mapa de tensiones rotulado "deformacion efectiva [-]" se lee mal y
# no da ninguna senal de estar mal.
TITULO_BARRA = {1: "espesor local [mm]",
                2: "deformacion efectiva [-]",
                3: "von Mises [Pa]",
                4: "deformacion total [mm]"}

# --- Protocolo comparativo microCT-FEA -------------------------------------
# Tapia D, Gonzalez A, Vidal F, Salinas P. A Specimen-Based Comparative
# MicroCT-FEA Analysis of Vertebral Trabecular Bone Microarchitecture and
# Mechanical Response in Two South American Cervids. Biology 2026;15:722.
#
# Del articulo se toman los cuatro valores que definen el ensayo. NO se toma
# el mallado: alli son tetraedros SOLID187 de 0.05 mm generados en ANSYS sobre
# el STL, y aqui son hexaedros de un voxel. Los campos son comparables en
# PATRON y en orden de magnitud, no cifra a cifra.
PAPER_E_S = 18e9          # Pa, modulo del tejido
PAPER_NU = 0.30
PAPER_CARGA_N = 100.0     # N, carga axial de compresion
PAPER_APOYO = "empotrado"  # "fixed support" en la cara opuesta


def apoyo_de_combo(cmb):
    """Clave CANONICA del apoyo ("deslizante"/"empotrado") de un combo de apoyo.

    Se lee por INDICE, nunca por texto: con la interfaz en ingles el texto es
    "fixed"/"sliding", y `ensayo_compresion` no lo entiende (antes lo tomaba en
    silencio por deslizante; ahora lanza ValueError). Lo que se pasa al ensayo
    y lo que se guarda en los registros es siempre esta clave; para mostrarla
    en pantalla, `_(clave)`.
    """
    return APOYOS[max(0, int(cmb.currentIndex()))]
PAPER_REF = ("Tapia, Gonzalez, Vidal &amp; Salinas, <i>Biology</i> 2026;15:722 "
             "— 18 GPa, &nu; = 0.30, 100 N axiales, apoyo empotrado")

# Metodos que compiten por el menor error antes del analisis comparado.
#
# NO estan los cuatro: el error del desempate mecanico incluye Ez_rel y Ez_Ex
# y el de los demas no, asi que su numero NO es comparable con el de estos
# tres —tiene mas sumandos y otra normalizacion—. Elegir el minimo entre los
# cuatro seria elegir por una cifra que no mide lo mismo en cada fila. Los
# tres de aqui comparten exactamente la misma definicion de error, asi que el
# minimo entre ellos si significa algo.
#
# Medido sobre el VOI proximal de H4: rapido 0.01449 (152 s), completo
# 0.00283 (491 s), equitativo 0.00542 (493 s). El completo suele ganar, pero
# no siempre: el equitativo acierta mejor la densidad y en un VOI distinto
# puede quedar por delante. Por eso se corren los tres y se elige, en vez de
# dar por hecho el ganador.
METODOS_COMPARABLES = ("rapido", "completo", "equitativo")

# Familias de microestructura y modos del selector de la barra de herramientas.
FAMILIAS = ("spinodoide", "dual-lattice")
MODOS_FAMILIA = ("spinodoide", "dual-lattice", "ambas")

# Ancho del texto de las filas que necesitan aviso al leerse
AVISO_SMI = ("El SMI esta confundido por la CONCAVIDAD (Salmon et al. 2015):\n"
             "a BV/TV alto sale negativo y pierde la interpretacion\n"
             "placa-barra. Ademas marching cubes sobre voxeles lo sesga\n"
             "un ~-15%. Reportarlo con esa salvedad, no como un numero limpio.")

# El coste crece con el CUBO del lado, y la memoria tambien. MEDIDO en esta
# maquina con 700 ondas, generar + morfometria completa:
#
#       lado    generar   morfometria   total    pico
#        48^3     1.5 s       1.2 s      2.7 s   137 MB
#        64^3     3.5 s       1.1 s      4.7 s   149 MB
#        96^3    11.9 s       1.7 s     13.7 s   199 MB
#       128^3    28.1 s       2.6 s     30.7 s   296 MB
#       192^3   100   s       5   s    105   s   694 MB
#       256^3   257   s       6   s    263   s  1470 MB
#
# Lo que domina es GENERAR, no medir: la morfometria apenas crece porque
# marching cubes trabaja sobre la superficie, que va con el cuadrado del lado.
#
# El tope es 256 y no 512 a proposito. Extrapolando desde el 256 medido,
# 512^3 son 134 millones de voxeles: del orden de 34 minutos y ~12 GB de pico.
# Un deslizador que ofrece una opcion capaz de agotar la memoria de la maquina
# no es una opcion, es una trampa.
AVISO_RES = ("El coste crece con el CUBO del lado. MEDIDO aqui, generar y\n"
             "medir tarda 4.7 s a 64 vox, 31 s a 128, 105 s a 192 y 263 s a\n"
             "256, con 1.5 GB de pico. Lo caro es generar, no medir.\n\n"
             "Por encima de ~96 vox conviene DESMARCAR 'Actualizar vista\n"
             "automaticamente': si no, cada deslizador que muevas lanza un\n"
             "calculo de minutos.\n\n"
             "La resolucion de MEDIDA por defecto es 64 porque es la que usa\n"
             "la app de MATLAB para comparar. Cambiarla cambia la\n"
             "discretizacion del candidato, asi que hay que declararla.")


# ---------------------------------------------------------------------------
# Hilo de calculo
# ---------------------------------------------------------------------------

class Trabajador(QtCore.QThread):
    """Ejecuta una funcion lenta fuera del hilo de la interfaz.

    Sin esto la ventana se congela durante la generacion y Windows la marca
    como 'no responde', que es justo lo que ocurre en la app de MATLAB con
    `drawnow` cuando la resolucion sube.
    """
    listo = QtCore.pyqtSignal(object)
    fallo = QtCore.pyqtSignal(str)
    avance = QtCore.pyqtSignal(int, int, str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self._fn, self._a, self._kw = fn, args, kwargs

    def informar(self, i, n, etapa):
        """Callback de progreso para tareas largas.

        Se invoca DESDE EL HILO DE TRABAJO. Emitir una senal es la forma
        correcta de cruzar al hilo de la interfaz: Qt la encola y la entrega
        en el bucle de eventos. Tocar los widgets directamente desde aqui
        provoca caidas intermitentes muy dificiles de reproducir.
        """
        self.avance.emit(int(i), int(n), str(etapa))

    def run(self):
        try:
            self.listo.emit(self._fn(*self._a, **self._kw))
        except Exception:
            self.fallo.emit(traceback.format_exc())


class DialogoPilaTiff(QtWidgets.QDialog):
    """Carga de una pila TIFF de micro-CT, con su escala fisica.

    POR QUE UN DIALOGO Y NO UN `getOpenFileName` A SECAS. Un .vtk trae dentro
    el tamano de voxel; una pila de imagenes no trae NADA. Si se abriera sin
    preguntar habria que suponer una escala, y una escala equivocada multiplica
    en silencio Tb.Th, Tb.Sp, TV, Conn.D y la rigidez aparente sin producir un
    solo sintoma visible: las figuras salen igual de bonitas y todos los
    numeros estan mal. Por eso la escala es obligatoria y se ensena de donde
    sale.

    Lo mismo con el umbral. Medido sobre H4, mover la segmentacion un +-15 %
    mueve Tb.Th un 62 %: es la mayor fuente de incertidumbre de todo el
    proceso, asi que se declara cual se aplico en vez de esconderlo.
    """

    UNIDADES = [("um", "µm  (micras — lo habitual en micro-CT)"),
                ("mm", "mm"), ("cm", "cm"), ("nm", "nm")]

    def __init__(self, padre, dir_inicial=""):
        super().__init__(padre)
        self.setWindowTitle(_("Cargar pila de imagenes (micro-CT)"))
        self.setMinimumWidth(560)
        lay = QtWidgets.QVBoxLayout(self)

        aviso = QtWidgets.QLabel(_(
            "<b>La imagen debe estar en formato TIFF (.tif / .tiff)</b>, "
            "que es lo que exportan los micro-CT.<br><br>"
            "Se admiten dos formas:<br>"
            "&nbsp;&nbsp;• una <b>carpeta</b> con una rebanada por archivo "
            "(se ordenan por el numero del nombre), o<br>"
            "&nbsp;&nbsp;• un unico <b>TIFF multipagina</b>.<br><br>"
            "Todas las rebanadas deben tener el <b>mismo tamano</b>: reescalar "
            "una cambiaria la morfometria sin avisar, asi que la carga se "
            "detiene si alguna difiere. Las imagenes pueden estar ya "
            "binarizadas o en escala de grises."))
        aviso.setWordWrap(True)
        aviso.setStyleSheet(
            "background:#fbf7e8; border:1px solid #e0d5a8; padding:9px;")
        lay.addWidget(aviso)

        form = QtWidgets.QFormLayout()
        lay.addLayout(form)

        fila = QtWidgets.QHBoxLayout()
        self.ed_origen = QtWidgets.QLineEdit()
        self.ed_origen.setPlaceholderText(_("carpeta de rebanadas o TIFF multipagina"))
        b_dir = QtWidgets.QPushButton(_("Carpeta…"))
        b_fic = QtWidgets.QPushButton(_("Archivo…"))
        b_dir.clicked.connect(self._elegir_carpeta)
        b_fic.clicked.connect(self._elegir_archivo)
        for w in (self.ed_origen, b_dir, b_fic):
            fila.addWidget(w)
        form.addRow(_("Origen"), fila)

        self.ed_patron = QtWidgets.QLineEdit("*.tif")
        self.ed_patron.setToolTip(_(
            "Solo para carpetas. Util cuando conviven varias series en el mismo "
            "sitio: en H4/Segmentadas hay BW_*.tif y SEG_*.tif mezcladas."))
        form.addRow(_("Patron de archivo"), self.ed_patron)

        fila = QtWidgets.QHBoxLayout()
        self.sp_tam = QtWidgets.QDoubleSpinBox()
        self.sp_tam.setDecimals(4)
        self.sp_tam.setRange(1e-4, 1e6)
        self.sp_tam.setValue(51.489)
        self.cmb_uni = QtWidgets.QComboBox()
        for clave, texto in self.UNIDADES:
            self.cmb_uni.addItem(texto, clave)
        fila.addWidget(self.sp_tam, 1)
        fila.addWidget(self.cmb_uni, 2)
        form.addRow(_("Tamano de voxel"), fila)

        self.lab_escala = QtWidgets.QLabel(_("sin comprobar"))
        self.lab_escala.setWordWrap(True)
        self.lab_escala.setStyleSheet("color:#555; font-size:11px;")
        form.addRow("", self.lab_escala)

        fila = QtWidgets.QHBoxLayout()
        self.chk_auto = QtWidgets.QCheckBox(_("automatico"))
        self.chk_auto.setChecked(True)
        self.chk_auto.setToolTip(_(
            "Si la pila ya es binaria toma > minimo. Si viene en escala de "
            "grises aplica Otsu y lo deja anotado."))
        self.cmb_clases = QtWidgets.QComboBox()
        self.cmb_clases.addItem(_("3 clases: aire / medula / hueso"), 3)
        self.cmb_clases.addItem(_("2 clases: Otsu clasico"), 2)
        self.sp_umbral = QtWidgets.QDoubleSpinBox()
        self.sp_umbral.setDecimals(1)
        self.sp_umbral.setRange(0, 1e9)
        self.sp_umbral.setEnabled(False)
        self.chk_auto.toggled.connect(
            lambda v: (self.sp_umbral.setEnabled(not v),
                       self.cmb_clases.setEnabled(v)))
        fila.addWidget(self.chk_auto)
        fila.addWidget(self.cmb_clases, 2)
        fila.addWidget(self.sp_umbral, 1)
        form.addRow(_("Umbral"), fila)

        aviso_u = QtWidgets.QLabel(_(
            "En una reconstruccion cruda el campo es aire en su mayor parte y "
            "hay <b>tres</b> poblaciones: con 2 clases la medula acaba contada "
            "como hueso. Medido en H4 frente a la segmentacion manual "
            "(tercios 0.28 / 0.55 / 0.77): con 2 clases salen 0.61 / 0.79 / "
            "0.94, con 3 clases 0.35 / 0.59 / 0.81. <b>Ningun automatico "
            "reproduce una segmentacion manual</b> — comprueba el valor."))
        aviso_u.setWordWrap(True)
        aviso_u.setStyleSheet("color:#555; font-size:11px;")
        form.addRow("", aviso_u)

        bb = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self._aceptar)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

        self._dir_inicial = dir_inicial

    # -- seleccion -------------------------------------------------------
    def _elegir_carpeta(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, _("Carpeta con las rebanadas"),
            self.ed_origen.text() or self._dir_inicial)
        if d:
            self.ed_origen.setText(d)
            self._buscar_escala(d)

    def _elegir_archivo(self):
        f, _x = QtWidgets.QFileDialog.getOpenFileName(
            self, _("TIFF multipagina"),
            self.ed_origen.text() or self._dir_inicial,
            _("TIFF (*.tif *.tiff);;Todos (*)"))
        if f:
            self.ed_origen.setText(f)
            self._buscar_escala(f)

    def _buscar_escala(self, origen):
        """Rellena el tamano de voxel desde el `*_rec.log` de SkyScan, si lo hay."""
        from spinpy.voi import tam_voxel_desde_log
        try:
            mm, info = tam_voxel_desde_log(origen)
        except Exception:
            mm, info = None, {}
        if mm is None:
            self.lab_escala.setText(_(
                "No se encontro ningun <tt>*_rec.log</tt> del escaner. "
                "<b>Escribe el tamano de voxel a mano</b> — no hay valor por "
                "defecto que sea seguro."))
            return
        self.sp_tam.setValue(mm / 1e-3)
        self.cmb_uni.setCurrentIndex(0)          # um
        self.lab_escala.setText(_(
            "Leido de <tt>{log}</tt> — <tt>{linea}</tt>").format(
                log=Path(info["log"]).name, linea=info["linea"]))

    def _aceptar(self):
        if not self.ed_origen.text().strip():
            QtWidgets.QMessageBox.warning(
                self, _("Falta el origen"),
                _("Elige la carpeta de rebanadas o el TIFF multipagina."))
            return
        self.accept()

    def valores(self):
        return {
            "origen": self.ed_origen.text().strip(),
            "patron": self.ed_patron.text().strip() or "*.tif",
            "tam_voxel": float(self.sp_tam.value()),
            "unidad": self.cmb_uni.currentData(),
            "umbral": None if self.chk_auto.isChecked()
                      else float(self.sp_umbral.value()),
            "clases_otsu": int(self.cmb_clases.currentData()),
        }


def _pixmap_secciones(cubo, alto=104):
    """Las tres secciones centrales de un cubo, una al lado de otra.

    Se eligen cortes ORTOGONALES por el centro en vez de una vista 3D porque
    para decidir un VOI hay que ver el interior: una isosuperficie solo enseña
    la piel, y dos cubos con interiores muy distintos se ven casi iguales por
    fuera. Ademas hay que pintar tres de estos por cada movimiento del mando y
    una vista 3D no daria esa cadencia.
    """
    from PyQt5 import QtGui
    c = np.asarray(cubo, bool)
    n = c.shape[0] // 2
    cortes = [c[n, :, :], c[:, n, :], c[:, :, n]]

    sep = 4
    lado = cortes[0].shape[0]
    ancho = 3 * lado + 2 * sep
    lienzo = np.full((lado, ancho), 235, np.uint8)
    for i, s in enumerate(cortes):
        # .T lleva (x, y) a (fila, columna) para que se vea con X horizontal
        img = np.where(np.asarray(s, bool).T, 30, 250).astype(np.uint8)
        x0 = i * (lado + sep)
        lienzo[:, x0:x0 + lado] = img

    qi = QtGui.QImage(lienzo.tobytes(), ancho, lado, ancho,
                      QtGui.QImage.Format_Grayscale8).copy()
    return QtGui.QPixmap.fromImage(qi).scaledToHeight(
        alto, QtCore.Qt.SmoothTransformation)


class DialogoRecorte(QtWidgets.QDialog):
    """Recorte del VOI cubico sobre el volumen cargado, con vista previa.

    POR QUE ESTE PASO NECESITA INTERFAZ. Una pila de escaner no es un VOI:
    628x664x438 son 1.8e8 voxeles, y ni la homogeneizacion ni el ajuste son
    viables ahi. El recorte ya existia (`vois_por_tercios`, `extraer_cubo`)
    pero solo desde un script, que es tanto como decir que no existia para
    quien no programa.

    SOBRE LOS NOMBRES ANATOMICOS. El signo de PC1 es arbitrario: la PCA
    determina las direcciones, no su sentido, y dos ejecuciones sobre datos
    casi iguales pueden devolver ejes opuestos. Por eso `vois_por_tercios`
    devuelve los tercios ORDENADOS a lo largo del eje y sin nombre, y por eso
    aqui la asignacion proximal/distal la hace quien conoce la pieza, con el
    BV/TV como pista: en estos sesamoideos la densidad crece de proximal a
    distal. Poner el nombre por nuestra cuenta seria acertar la mitad de las
    veces y no avisar de la otra mitad.
    """

    def __init__(self, padre, mask, spacing, marco, tercios, lado_mm):
        super().__init__(padre)
        self.setWindowTitle(_("Recortar VOI cubico"))
        self.setMinimumWidth(760)
        self.mask, self.spacing = mask, spacing
        self.ejes = marco["ejes"]
        self.centro = marco["centro"]
        self.proy = marco["proy"]
        self.tercios = tercios
        self._cache_libre = None

        lay = QtWidgets.QVBoxLayout(self)

        cab = QtWidgets.QLabel(_(
            "El hueso se orienta por PCA y se corta en tres tramos a lo largo "
            "de su eje mayor. <b>El sentido del eje es arbitrario</b>: los "
            "tramos van ordenados, pero cual es proximal lo decides tu. Pista: "
            "en estos sesamoideos la densidad crece de proximal a distal."))
        cab.setWordWrap(True)
        cab.setStyleSheet(
            "background:#fbf7e8; border:1px solid #e0d5a8; padding:8px;")
        lay.addWidget(cab)

        fila = QtWidgets.QHBoxLayout()
        fila.addWidget(QtWidgets.QLabel(_("Lado del cubo (mm)")))
        self.sp_lado = QtWidgets.QDoubleSpinBox()
        self.sp_lado.setDecimals(2)
        self.sp_lado.setRange(0.2, 50.0)
        self.sp_lado.setValue(lado_mm)
        self.sp_lado.setSingleStep(0.5)
        fila.addWidget(self.sp_lado)
        self.b_recalc = QtWidgets.QPushButton(_("Recalcular tramos"))
        fila.addWidget(self.b_recalc)
        fila.addStretch(1)
        fila.addWidget(QtWidgets.QLabel(_("El tramo 1 es")))
        self.cmb_orden = QtWidgets.QComboBox()
        self.cmb_orden.addItem(_("proximal"), "prox")
        self.cmb_orden.addItem(_("distal"), "dist")
        fila.addWidget(self.cmb_orden)
        lay.addLayout(fila)

        self.grupo = QtWidgets.QButtonGroup(self)
        self.filas = []
        self.caja_tramos = QtWidgets.QVBoxLayout()
        lay.addLayout(self.caja_tramos)

        # --- modo libre ---------------------------------------------------
        gb = QtWidgets.QGroupBox(_("Elegir libremente"))
        gl = QtWidgets.QGridLayout(gb)
        self.rb_libre = QtWidgets.QRadioButton(_("usar esta posicion"))
        self.grupo.addButton(self.rb_libre, 99)
        gl.addWidget(self.rb_libre, 0, 0)

        self.lab_libre = QtWidgets.QLabel("—")
        self.lab_libre.setStyleSheet("color:#555;")
        gl.addWidget(self.lab_libre, 0, 1, 1, 3)

        self.vista_libre = QtWidgets.QLabel()
        self.vista_libre.setFixedHeight(104)
        gl.addWidget(self.vista_libre, 1, 0, 4, 1)

        self.mandos = {}
        for r, (clave, texto, lo, hi, val) in enumerate([
                ("t",  _("posicion en el eje (%)"), 0, 100, 50),
                ("o2", _("desplazar eje 2 (mm)"), -100, 100, 0),
                ("o3", _("desplazar eje 3 (mm)"), -100, 100, 0)]):
            gl.addWidget(QtWidgets.QLabel(texto), r + 1, 1)
            s = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            s.setRange(lo, hi)
            s.setValue(val)
            s.sliderReleased.connect(self._preview_libre)
            gl.addWidget(s, r + 1, 2)
            e = QtWidgets.QLabel(str(val))
            e.setFixedWidth(46)
            s.valueChanged.connect(lambda v, w=e: w.setText(str(v)))
            gl.addWidget(e, r + 1, 3)
            self.mandos[clave] = s
        lay.addWidget(gb)

        bb = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self._aceptar)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

        self.b_recalc.clicked.connect(self._recalcular)
        self.cmb_orden.currentIndexChanged.connect(self._reetiquetar)
        self._pintar_tramos()
        self._preview_libre()

    # -- tramos ----------------------------------------------------------
    def _nombres(self):
        base = [_("proximal"), _("centro"), _("distal")]
        return base if self.cmb_orden.currentData() == "prox" else base[::-1]

    def _pintar_tramos(self):
        while self.caja_tramos.count():
            w = self.caja_tramos.takeAt(0).widget()
            if w:
                w.deleteLater()
        self.filas = []
        for i, t in enumerate(self.tercios):
            w = QtWidgets.QWidget()
            h = QtWidgets.QHBoxLayout(w)
            rb = QtWidgets.QRadioButton()
            self.grupo.addButton(rb, i)
            if i == 1:
                rb.setChecked(True)
            h.addWidget(rb)
            v = QtWidgets.QLabel()
            v.setPixmap(_pixmap_secciones(t["cubo"]))
            h.addWidget(v)
            txt = QtWidgets.QLabel()
            txt.setWordWrap(True)
            h.addWidget(txt, 1)
            self.caja_tramos.addWidget(w)
            self.filas.append((rb, txt, t))
        self._reetiquetar()

    def _reetiquetar(self):
        nom = self._nombres()
        for i, (rb, txt, t) in enumerate(self.filas):
            rb.setText(f"{i + 1}")
            aviso = ""
            if t["fuera"] > 0.005:
                aviso = _("  ·  <b style='color:#c0392b;'>{p:.1f} % del cubo "
                          "cae fuera del volumen</b>").format(p=100 * t["fuera"])
            txt.setText(
                f"<b>Tramo {i + 1} — {nom[i]}</b><br>"
                + _("BV/TV {b:.4f} · lado {lv} vox").format(
                    b=t["bvtv"], lv=t["lado_vox"]) + aviso)

    def _recalcular(self):
        from spinpy.voi import extraer_cubo, centros_por_tercios
        lado = float(self.sp_lado.value())
        centros = centros_por_tercios(self.proy, 3)
        self.tercios = []
        for i, c in enumerate(centros):
            cubo, fuera = extraer_cubo(self.mask, self.spacing, c, self.ejes,
                                       self.centro, lado)
            self.tercios.append({"indice": i, "cubo": cubo, "centro_pca": c,
                                 "fuera": fuera, "bvtv": float(cubo.mean()),
                                 "lado_vox": int(cubo.shape[0])})
        self._pintar_tramos()
        self._preview_libre()

    # -- libre -----------------------------------------------------------
    def _cubo_libre(self):
        from spinpy.voi import centro_en, extraer_cubo
        t = self.mandos["t"].value() / 100.0
        c = centro_en(self.proy, t,
                      off2=self.mandos["o2"].value() / 10.0,
                      off3=self.mandos["o3"].value() / 10.0)
        cubo, fuera = extraer_cubo(self.mask, self.spacing, c, self.ejes,
                                   self.centro, float(self.sp_lado.value()))
        return {"cubo": cubo, "centro_pca": c, "fuera": fuera,
                "bvtv": float(cubo.mean()), "lado_vox": int(cubo.shape[0])}

    def _preview_libre(self):
        try:
            r = self._cubo_libre()
        except Exception as e:
            self.lab_libre.setText(f"<span style='color:#c0392b;'>{e}</span>")
            return
        self._cache_libre = r
        self.vista_libre.setPixmap(_pixmap_secciones(r["cubo"]))
        aviso = ""
        if r["fuera"] > 0.005:
            aviso = _("  ·  <b style='color:#c0392b;'>{p:.1f} % fuera</b>"
                      ).format(p=100 * r["fuera"])
        self.lab_libre.setText(
            _("BV/TV {b:.4f} · lado {lv} vox").format(
                b=r["bvtv"], lv=r["lado_vox"]) + aviso)

    def _aceptar(self):
        i = self.grupo.checkedId()
        if i == 99:
            self.elegido = self._cache_libre
            self.etiqueta = _("libre")
        elif 0 <= i < len(self.filas):
            self.elegido = self.filas[i][2]
            self.etiqueta = self._nombres()[i]
        else:
            QtWidgets.QMessageBox.warning(self, _("Nada seleccionado"),
                                          _("Elige un tramo o el modo libre."))
            return
        if self.elegido is None:
            QtWidgets.QMessageBox.warning(
                self, _("Sin vista previa"),
                _("Mueve algun mando para calcular el cubo libre."))
            return
        self.accept()


# ---------------------------------------------------------------------------
# Utilidades de render
# ---------------------------------------------------------------------------

def malla_de_mascara(BW, spacing):
    """Isosuperficie 0.5 de la mascara como PolyData.

    Se construye un ImageData con la mascara en los NODOS y se contornea: es
    marching cubes de VTK, el mismo algoritmo que mide `area_superficie`, de
    modo que lo que se ve es lo que se mide.

    `flatten(order='F')` porque VTK recorre X como indice mas rapido y BW tiene
    dimension 1 = X - la misma convencion que `readVTKVOI`.

    Las normales se calculan EXPLICITAMENTE: `contour()` devuelve la malla sin
    ellas y, sin normales, `smooth_shading` no tiene con que sombrear. La
    estructura se renderiza entonces como una silueta plana de un solo color y
    parece un bloque macizo aunque tenga BV/TV de 0.28.
    """
    if BW is None or not BW.any():
        return pv.PolyData()
    grid = pv.ImageData(dimensions=BW.shape,
                        spacing=tuple(float(s) for s in spacing))
    grid.point_data["v"] = BW.astype(np.float32).flatten(order="F")
    malla = grid.contour([0.5], scalars="v")
    if malla.n_points:
        malla = malla.compute_normals(auto_orient_normals=True,
                                      consistent_normals=False)
    return malla


# ---------------------------------------------------------------------------
# Dialogos de figura: espesor, convergencia de malla y campos del ensayo
# ---------------------------------------------------------------------------

class DialogoFigura(QtWidgets.QDialog):
    """Base de los dialogos que muestran una figura con su lectura al pie.

    Todos comparten la misma estructura —lienzo de matplotlib, tabla o nota en
    HTML debajo, y un boton para guardar la figura— y la misma regla: **ninguna
    figura se muestra sin el texto que dice como leerla**. Un histograma de
    espesor sin la advertencia del sesgo del 30%, o una curva de convergencia
    sin el veredicto, invitan a citar un numero que no es citable.
    """

    def __init__(self, titulo, parent=None, tam=(780, 560)):
        super().__init__(parent)
        self.setWindowTitle(titulo)
        self.resize(*tam)
        self._fig = None

    @staticmethod
    def _lienzo(figsize=(7.4, 4.4)):
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        fig = Figure(figsize=figsize, tight_layout=True)
        return fig, FigureCanvasQTAgg(fig)

    def _montar(self, fig, lienzo, html):
        self._fig = fig
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(lienzo, 1)

        txt = QtWidgets.QLabel(html)
        txt.setTextFormat(QtCore.Qt.RichText)
        txt.setWordWrap(True)
        area = QtWidgets.QScrollArea()
        area.setWidget(txt)
        area.setWidgetResizable(True)
        area.setFrameShape(QtWidgets.QFrame.NoFrame)
        area.setMaximumHeight(220)
        lay.addWidget(area)

        botones = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton(_("Guardar figura…"))
        b.clicked.connect(self._guardar)
        botones.addWidget(b)
        botones.addStretch(1)
        c = QtWidgets.QPushButton(_("Cerrar"))
        c.clicked.connect(self.accept)
        botones.addWidget(c)
        lay.addLayout(botones)

    def _guardar(self):
        f, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Guardar la figura"), str(DATOS / "figura.png"),
            _("PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"))
        if f and self._fig is not None:
            self._fig.savefig(f, dpi=200)
            p = self.parent()
            if p is not None:
                p.statusBar().showMessage(_("Figura guardada en {f}").format(f=f))


class DialogoEspesor(DialogoFigura):
    """Distribucion del espesor local, VOI y spinodoide superpuestos.

    POR QUE UN HISTOGRAMA Y NO SOLO Tb.Th. El Tb.Th de la tabla es un escalar
    derivado de Parfitt (2*BV/BS): supone un modelo de placas y colapsa toda la
    estructura en un numero. Dos VOIs con el mismo Tb.Th pueden tener uno las
    trabeculas homogeneas y el otro una mezcla de trabeculas gruesas y finas
    —que es justo la firma de la perdida osea— y el escalar no los distingue.
    La distribucion si.

    QUE NO SE PUEDE HACER CON ESTOS NUMEROS. La escala absoluta NO es citable
    como Tb.Th: a 3-4 voxeles de grosor el metodo de esferas inscritas
    subestima un ~30%. Sirve para comparar dos estructuras medidas con el MISMO
    spacing y para leer la FORMA de la distribucion (moda, cola, bimodalidad).
    Los dos paneles comparten eje x por esa razon; con ejes distintos la
    superposicion no significaria nada.
    """

    def __init__(self, esp, parent=None, etq_spin=None):
        super().__init__(_("Distribucion del espesor local"), parent)
        fig, lienzo = self._lienzo()
        ax = fig.add_subplot(111)

        colores = {"voi": "#b62324", "spin": "#1f4e9c"}
        etiquetas = {"voi": _("VOI real"), "spin": etq_spin or _("Spinodoide")}
        series, resumen = {}, []
        for cual in ("voi", "spin"):
            if cual not in esp:
                continue
            v = np.asarray(esp[cual][0], float)
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                series[cual] = v

        if series:
            # Bordes COMUNES a las dos series. Con bordes propios cada
            # histograma elegiria su ancho de barra y las alturas dejarian de
            # ser comparables aunque compartan eje.
            todo = np.concatenate(list(series.values()))
            bordes = np.linspace(0.0, float(np.percentile(todo, 99.5)), 60)
            for cual, v in series.items():
                # densidad: los dos volumenes tienen distinto numero de
                # voxeles solidos, asi que la frecuencia bruta compararia
                # tamanos y no formas.
                ax.hist(v, bins=bordes, density=True, histtype="stepfilled",
                        alpha=0.40, color=colores[cual], label=etiquetas[cual])
                ax.hist(v, bins=bordes, density=True, histtype="step",
                        lw=1.6, color=colores[cual])
                med = float(np.median(v))
                ax.axvline(med, color=colores[cual], ls="--", lw=1.2)
                s = estadisticas_esp(v)
                resumen.append(
                    f"<tr><td><b>{etiquetas[cual]}</b></td>"
                    f"<td>{s['media']:.4f}</td><td>{s['mediana']:.4f}</td>"
                    f"<td>{s['sd']:.4f}</td><td>{s['p05']:.4f}</td>"
                    f"<td>{s['p95']:.4f}</td>"
                    f"<td>{s['sd']/max(s['media'],1e-30):.3f}</td></tr>")
            ax.legend(frameon=False)

        ax.set_xlabel(_("espesor local [mm]  ·  la linea de trazos es la "
                        "mediana"))
        ax.set_ylabel(_("densidad de probabilidad"))
        ax.spines[["top", "right"]].set_visible(False)

        self._montar(fig, lienzo,
            "<table cellpadding=4 cellspacing=0 border=1 "
            "style='border-collapse:collapse'>"
            + _("<tr><th></th><th>media</th><th>mediana</th><th>sd</th>"
                "<th>p05</th><th>p95</th><th>CV</th></tr>")
            + "".join(resumen) + "</table>"
            + "<p style='color:#666;font-size:10px'>"
            + _("Todo en mm. El <b>CV</b> (sd/media) es lo que el Tb.Th "
                "escalar no puede dar: mide cuanto varia el grosor dentro de "
                "la misma estructura.<br>"
                "<b>La escala absoluta no es citable como Tb.Th</b>: a 3-4 "
                "voxeles de grosor el metodo subestima un ~30%. Usar como "
                "medida RELATIVA entre estructuras medidas al mismo spacing.")
            + "</p>")


class DialogoConvergencia(DialogoFigura):
    """E_app frente al tamano de elemento, con el veredicto al pie.

    La app ya advertia de que `E_app` no converge sin mas en estructuras poco
    densas, pero no ofrecia forma de comprobarlo. Esta figura es esa forma.

    SE TRAZAN DOS COSAS EN EL MISMO EJE. Arriba `E_app`, que es lo que se
    quiere; abajo la **fraccion portante**, que es el diagnostico. Si la
    fraccion portante sube al refinar, la estructura conectada esta cambiando,
    y entonces la variacion de E_app no es error de discretizacion: son
    geometrias distintas. Sin esa segunda curva, una serie que sube y baja
    parece ruido numerico y se descarta por el motivo equivocado.
    """

    def __init__(self, res, parent=None):
        super().__init__(_("Convergencia de malla"), parent,
                         tam=(820, 620))
        fig, lienzo = self._lienzo(figsize=(7.6, 4.8))
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax)

        val = [p for p in res["puntos"] if p["ok"]]
        malos = [p for p in res["puntos"] if not p["ok"]]

        if val:
            n = [p["n"] for p in val]
            E = [p["E_app"] / 1e6 for p in val]
            ax.plot(n, E, "o-", color="#1f4e9c", lw=1.8, ms=7)
            for x, y in zip(n, E):
                ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points",
                            xytext=(0, 8), ha="center", fontsize=8,
                            color="#1f4e9c")
            if "E_extrapolado" in res:
                ax.axhline(res["E_extrapolado"] / 1e6, color="#1a7f37",
                           ls="--", lw=1.4,
                           label=f"extrapolado {res['E_extrapolado']/1e6:.0f} MPa")
                ax.legend(frameon=False, fontsize=9)
            ax2.plot(n, [p["frac_portante"] for p in val], "s-",
                     color="#b06000", lw=1.6, ms=6, label="fraccion portante")
            ax2.plot(n, [p["rho"] for p in val], "^--", color="#888", lw=1.2,
                     ms=5, label="BV/TV")
            ax2.legend(frameon=False, fontsize=9)

        ax.set_ylabel("E$_{app}$ [MPa]")
        ax.set_title(_("Compresion en {eje}").format(
            eje=res.get("eje", "Z")), fontsize=10)
        ax2.set_xlabel(_("resolucion [voxeles por lado]"))
        ax2.set_ylabel(_("fraccion"))
        for a in (ax, ax2):
            a.spines[["top", "right"]].set_visible(False)
            a.grid(alpha=0.25)

        filas = []
        for p in res["puntos"]:
            if p["ok"]:
                filas.append(
                    f"<tr><td><b>{p['n']}³</b></td><td>{p['h']*1000:.1f} µm</td>"
                    f"<td>{p['E_app']/1e6:.1f} MPa</td>"
                    f"<td>{p['rho']:.4f}</td>"
                    f"<td>{p['frac_portante']:.3f}</td>"
                    f"<td>{p['n_elem']:,}</td>"
                    f"<td>{p['residuo']:.1e}</td></tr>")
            else:
                filas.append(f"<tr><td><b>{p['n']}³</b></td>"
                             f"<td colspan=6 style='color:#b62324'>"
                             f"{p['msg']}</td></tr>")

        extra = ""
        if "orden" in res:
            extra = _("<p>Orden observado <b>{orden}</b> · extrapolado "
                      "<b>{ext} MPa</b> · error estimado del punto mas fino "
                      "<b>{err}%</b></p>").format(
                          orden=f"{res['orden']:.2f}",
                          ext=f"{res['E_extrapolado']/1e6:.1f}",
                          err=f"{100*res['error_estimado_rel']:.1f}")
        if malos:
            extra += (
                "<p style='color:#666;font-size:10px'>"
                + _("Los puntos que no resolvieron no son un fallo del "
                    "estudio: a resolucion baja una estructura poco densa "
                    "queda casi desconectada y el sistema mal condicionado. El "
                    "solver lo detecta por el residuo y se descarta el punto "
                    "en vez de devolver un campo equivocado.")
                + "</p>")

        self._montar(fig, lienzo,
            "<table cellpadding=4 cellspacing=0 border=1 "
            "style='border-collapse:collapse'>"
            "<tr><th>malla</th><th>h</th><th>E<sub>app</sub></th><th>BV/TV</th>"
            "<th>frac.port.</th><th>elem.</th><th>residuo</th></tr>"
            + "".join(filas) + "</table>"
            + extra
            + f"<p><b>Veredicto:</b> {res.get('veredicto','')}</p>")


class DialogoDistribucionFE(DialogoFigura):
    """Distribucion de la deformacion efectiva y de la tension de von Mises.

    POR QUE LA DISTRIBUCION Y NO EL MAXIMO. El criterio de Pistoia no mira el
    maximo —seria un solo elemento, y en una malla de voxeles ese elemento suele
    ser una esquina escalonada sin significado fisico—: mira el **percentil 98**,
    el valor que supera el 2% del tejido. Eso es un punto de la cola, y para
    saber si es robusto hay que ver la cola entera. Si esta concentrada en unos
    pocos nudos, el criterio depende de esos nudos; si esta repartida, el
    resultado es estable.

    Se traza en escala logaritmica en x: los dos campos abarcan varios ordenes
    de magnitud entre el material apenas cargado y los nudos, y en escala
    lineal todo se amontona contra el cero.
    """

    def __init__(self, datos, parent=None, etq_spin=None):
        super().__init__(_("Distribuciones del ensayo de compresion"),
                         parent, tam=(860, 620))
        fig, lienzo = self._lienzo(figsize=(8.0, 4.6))
        ejes = [fig.add_subplot(121), fig.add_subplot(122)]
        colores = {"voi": "#b62324", "spin": "#1f4e9c"}
        etiquetas = {"voi": _("VOI real"), "spin": etq_spin or _("Spinodoide")}

        campos = [("eps", "deformacion efectiva [-]", EPS_CRITICA,
                   "0.7% critico"),
                  ("vm", "von Mises [MPa]", None, None)]
        filas = []
        for ax, (clave, etiqueta, marca, marca_txt) in zip(ejes, campos):
            series = {}
            for cual in ("voi", "spin"):
                v = datos.get(cual, {}).get(clave)
                if v is None:
                    continue
                v = np.asarray(v, float)
                v = v[np.isfinite(v) & (v > 0)]
                if v.size:
                    series[cual] = v / 1e6 if clave == "vm" else v
            if series:
                todo = np.concatenate(list(series.values()))
                bordes = np.logspace(np.log10(max(todo.min(), todo.max() * 1e-5)),
                                     np.log10(todo.max()), 60)
                for cual, v in series.items():
                    ax.hist(v, bins=bordes, density=True, histtype="stepfilled",
                            alpha=0.38, color=colores[cual],
                            label=etiquetas[cual])
                    ax.hist(v, bins=bordes, density=True, histtype="step",
                            lw=1.5, color=colores[cual])
                    p98 = float(np.percentile(v, 98.0))
                    ax.axvline(p98, color=colores[cual], ls=":", lw=1.5)
                    filas.append(
                        f"<tr><td><b>{etiquetas[cual]}</b></td>"
                        f"<td>{etiqueta.split(' [')[0]}</td>"
                        f"<td>{np.median(v):.4g}</td>"
                        f"<td>{p98:.4g}</td>"
                        f"<td>{v.max():.4g}</td>"
                        f"<td>{v.max()/max(p98,1e-30):.2f}</td></tr>")
                ax.set_xscale("log")
                ax.legend(frameon=False, fontsize=8)
            if marca:
                ax.axvline(marca, color="#1a7f37", lw=1.6)
                ax.annotate(marca_txt, (marca, 0), rotation=90, fontsize=8,
                            color="#1a7f37", va="bottom", ha="right",
                            xycoords=("data", "axes fraction"))
            ax.set_xlabel(etiqueta)
            ax.set_ylabel(_("densidad"))
            ax.spines[["top", "right"]].set_visible(False)

        self._montar(fig, lienzo,
            "<table cellpadding=4 cellspacing=0 border=1 "
            "style='border-collapse:collapse'>"
            + _("<tr><th></th><th>campo</th><th>mediana</th><th>p98</th>"
                "<th>maximo</th><th>max/p98</th></tr>")
            + "".join(filas) + "</table>"
            + "<p style='color:#666;font-size:10px'>"
            + _("La linea de puntos es el <b>percentil 98</b>, que es el que "
                "decide el criterio de Pistoia; la verde es el 0.7% critico. "
                "La columna <b>max/p98</b> es el diagnostico: si es cercana a "
                "1 la cola esta repartida y el criterio es estable; si es "
                "grande, unos pocos elementos —a menudo esquinas escalonadas "
                "de la malla, sin significado fisico— dominan el extremo, y "
                "conviene mirar la convergencia de malla antes de creerse la "
                "carga de fallo.<br>"
                "Escala logaritmica en x: los campos abarcan varios ordenes "
                "de magnitud entre el material descargado y los nudos.")
            + "</p>")


class DialogoDispersion(DialogoFigura):
    """Suelo de ruido del generador: K realizaciones de los mismos parametros.

    LA PREGUNTA QUE CONTESTA. Cuando el candidato difiere del VOI en un 3% en
    Tb.Sp, ¿es que el ajuste falla o es que el generador no reproduce Tb.Sp
    mejor que eso? Sin esta medida no hay forma de saberlo, y la diferencia
    porcentual sola invita a interpretar ruido como senal.

    La columna **z** es el criterio con el que se valido el port: cuantas
    desviaciones tipicas del propio generador separan al VOI de la media de las
    realizaciones. |z| < 2 significa indistinguible del ruido, y no se puede
    pedir a un ajuste que acierte mas de lo que el generador se repite a si
    mismo.
    """

    def __init__(self, disp, parent=None):
        super().__init__(_("Dispersion del generador"), parent,
                         tam=(860, 620))
        fig, lienzo = self._lienzo(figsize=(8.0, 4.4))
        ax = fig.add_subplot(111)

        r = disp["resumen"]
        # Se grafican los CV, no los valores: cada metrica tiene sus unidades y
        # su orden de magnitud, y ponerlas en un mismo eje sin normalizar solo
        # mostraria cual es mas grande.
        mets = [m for m in r if np.isfinite(r[m].get("cv_pct", np.nan))]
        mets.sort(key=lambda m: r[m]["cv_pct"])
        if mets:
            y = np.arange(len(mets))
            cv = [r[m]["cv_pct"] for m in mets]
            hay_z = any("z" in r[m] for m in mets)
            colores = ["#b62324" if abs(r[m].get("z", 0)) >= 2 else "#1f4e9c"
                       for m in mets]
            ax.barh(y, cv, color=colores, height=0.62)
            ax.set_yticks(y)
            ax.set_yticklabels(mets, fontsize=9)
            for i, (m, c) in enumerate(zip(mets, cv)):
                etq = f"{c:.2f} %"
                if "z" in r[m] and np.isfinite(r[m]["z"]):
                    etq += f"   z = {r[m]['z']:+.1f}"
                ax.text(c, i, "  " + etq, va="center", fontsize=8,
                        color="#333")
            ax.set_xlabel(_("coeficiente de variacion entre "
                            "realizaciones [%]"))
            if hay_z:
                ax.set_title(_("rojo: el VOI queda a 2 sd o mas de la media "
                               "de las realizaciones"),
                             fontsize=9, color="#555")
            ax.margins(x=0.22)
            ax.spines[["top", "right"]].set_visible(False)

        filas = []
        for m in mets:
            d = r[m]
            z = d.get("z")
            col = ("#b62324" if z is not None and np.isfinite(z)
                   and abs(z) >= 2 else "#1a7f37")
            filas.append(
                f"<tr><td><b>{m}</b></td>"
                f"<td>{d['media']:.5g}</td><td>{d['sd']:.3g}</td>"
                f"<td>{d['cv_pct']:.3f} %</td>"
                f"<td>{d['min']:.5g} – {d['max']:.5g}</td>"
                + (f"<td>{d['ref']:.5g}</td>"
                   f"<td>{d.get('dif_rel_pct', float('nan')):+.2f} %</td>"
                   f"<td style='color:{col}'><b>"
                   f"{z:+.2f}</b></td>" if z is not None and np.isfinite(z)
                   else "<td>—</td><td>—</td><td>—</td>")
                + "</tr>")

        self._montar(fig, lienzo,
            _("<p>{n} realizaciones de la MISMA parametrizacion, semillas "
              "{a}…{b}</p>").format(
                n=disp["n_semillas"], a=disp["semilla_base"],
                b=disp["semilla_base"] + disp["n_semillas"] - 1)
            + "<table cellpadding=4 cellspacing=0 border=1 "
              "style='border-collapse:collapse'>"
            + _("<tr><th></th><th>media</th><th>sd</th><th>CV</th>"
                "<th>rango</th><th>VOI</th><th>dif.</th><th>z</th></tr>")
            + "".join(filas) + "</table>"
            + "<p style='color:#666;font-size:10px'>"
            + _("<b>z</b> = (media de las realizaciones − VOI) / sd. Con "
                "|z| &lt; 2 la diferencia es indistinguible del ruido del "
                "propio generador: no se le puede pedir al ajuste que acierte "
                "mas de lo que el generador se repite a si mismo.<br>"
                "Un CV pequeno significa que el generador reproduce esa "
                "metrica de forma estable, <b>no</b> que la metrica sea "
                "exacta: BS es de las mas estables y arrastra el sesgo de "
                "marching cubes.")
            + "</p>")


class DialogoVarianza(DialogoFigura):
    """Descomposicion de varianza del lote, con el N efectivo en primer plano.

    EL ERROR QUE ESTE DIALOGO EXISTE PARA IMPEDIR. Generar K replicas de cada
    VOI y tratarlas despues como K especimenes independientes es
    pseudorreplicacion: infla los grados de libertad, estrecha los intervalos
    de confianza y convierte cualquier diferencia en significativa. El numero
    de filas de la tabla NO es el N del estudio.

    Por eso el N efectivo va arriba y en grande, y no escondido en una columna.
    """

    def __init__(self, inf, meta, parent=None):
        super().__init__(_("Varianza del lote y equivalencia"), parent,
                         tam=(920, 660))
        fig, lienzo = self._lienzo(figsize=(8.4, 4.2))
        ax = fig.add_subplot(111)

        met = list(inf["metrica"])
        entre = np.nan_to_num(np.asarray(inf["cv_entre_pct"], float))
        dentro = np.nan_to_num(np.asarray(inf["cv_dentro_pct"], float))
        y = np.arange(len(met))
        ax.barh(y - 0.19, entre, height=0.36, color="#1f4e9c",
                label=_("entre especimenes (biologica)"))
        ax.barh(y + 0.19, dentro, height=0.36, color="#b06000",
                label=_("dentro (ruido del generador)"))
        ax.set_yticks(y)
        ax.set_yticklabels(met, fontsize=9)
        ax.set_xlabel(_("coeficiente de variacion [%]"))
        ax.legend(frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(x=0.08)

        filas = []
        for _x, f in inf.iterrows():
            eq = f.get("eq_equivalente")
            marca = ("<span style='color:#1a7f37'>" + _("si") + "</span>"
                     if eq is True
                     else "<span style='color:#b62324'>" + _("no") + "</span>"
                     if eq is False else "—")
            # Con un solo VOI no hay varianza entre especimenes, y el informe
            # deja varias columnas en None: se muestran como «—» en vez de
            # reventar al formatearlas.
            filas.append(
                f"<tr><td><b>{f['metrica']}</b></td>"
                f"<td>{_fmt_num(f.get('n_grupos'), '.0f')}</td>"
                f"<td>{_fmt_num(f.get('n_total'), '.0f')}</td>"
                f"<td><b>{_fmt_num(f.get('n_efectivo'), '.1f')}</b></td>"
                f"<td>{_fmt_num(f.get('cv_entre_pct'), '.2f', ' %')}</td>"
                f"<td>{_fmt_num(f.get('cv_dentro_pct'), '.2f', ' %')}</td>"
                f"<td>{_fmt_num(f.get('ICC'), '.3f')}</td>"
                f"<td>{_fmt_num(f.get('eq_dif_rel_pct'), '+.2f', ' %')}</td>"
                f"<td>{_fmt_num(f.get('eq_margen_rel_pct'), '.2f', ' %')}</td>"
                f"<td>{marca}</td></tr>")

        n_ef = _n_efectivo_min(inf)
        self._montar(fig, lienzo,
            _("<p style='font-size:13px'><b>{nv} VOIs · {nr} replicas cada "
              "uno · {nf} filas en la tabla</b><br>"
              "<span style='color:#b62324;font-size:15px'>N efectivo para "
              "inferencia biologica: {nef}</span></p>").format(
                nv=meta["n_vois"], nr=meta["n_replicas"],
                nf=meta.get("n_filas", "?"), nef=f"{n_ef:.1f}")
            + "<table cellpadding=4 cellspacing=0 border=1 "
              "style='border-collapse:collapse'>"
            + _("<tr><th></th><th>grupos</th><th>filas</th>"
                "<th>N efectivo</th><th>CV entre</th><th>CV dentro</th>"
                "<th>ICC</th><th>dif. real−sint.</th><th>margen</th>"
                "<th>equivalente</th></tr>")
            + "".join(filas) + "</table>"
            + "<p style='color:#666;font-size:10px'>"
            + _("<b>Las replicas son tecnicas.</b> Reducen la incertidumbre "
                "de cada espécimen dividiendo la varianza DENTRO por K, y no "
                "tocan la varianza ENTRE, que es la biologica. Los grados de "
                "libertad para una afirmacion sobre la especie los ponen los "
                "animales, no las realizaciones.<br>"
                "La equivalencia es <b>TOST</b>, no una t de Student: «no se "
                "rechazo la nula» no demuestra que dos cosas sean iguales — "
                "con N pequeno nunca se rechaza y con N grande siempre. TOST "
                "declara un margen de antemano y demuestra que la diferencia "
                "cae dentro. El margen por defecto es el CV DENTRO observado, "
                "es decir el suelo de ruido del propio generador: pedir menos "
                "seria pedir que el sintetico se parezca al real mas de lo "
                "que el real se parece a si mismo.")
            + "</p>")


# ---------------------------------------------------------------------------
# Controles
# ---------------------------------------------------------------------------

class DialogoSimulacion(DialogoFigura):
    """Resultado de una simulacion in silico: curvas, tabla y el 3D de cada paso.

    Es el unico dialogo de figura con visor 3D, y no por adorno: una curva de
    rigidez que cae en picado no dice POR QUE. Ver que en ese paso desaparecio
    la ultima trabecula vertical, o que el dano se concentro en un cuello, es
    lo que convierte la curva en una explicacion.

    La tabla lleva el residuo y la fraccion portante de cada paso por la misma
    razon que `DialogoConvergencia`: un E que cae porque la estructura se
    desconecto y un E que cae porque el solver no convergio no son el mismo
    resultado, y solo esas dos columnas los separan.
    """

    COLOR = {"inicial": "#1f4e9c", "perdida": "#1f4e9c",
             "recuperacion": "#1a7f37"}

    def __init__(self, res, parent=None):
        self._fallo = res.get("tipo") == "fallo_progresivo"
        super().__init__(_("Fallo progresivo") if self._fallo
                         else _("Simulacion de perdida osea"), parent,
                         tam=(1240, 760))
        self._res = res
        fig, lienzo = self._lienzo(figsize=(6.6, 5.2))
        if self._fallo:
            self._figura_fallo(fig, res["pasos"])
        else:
            self._figura_perdida(fig, res["pasos"])
        self._montar(fig, lienzo, self._html(res))

        # El lienzo que `_montar` puso arriba pasa a la izquierda de un
        # divisor, con el visor 3D del paso elegido a la derecha.
        lay = self.layout()
        lay.removeWidget(lienzo)
        div = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        div.addWidget(lienzo)
        caja = QtWidgets.QWidget()
        cl = QtWidgets.QVBoxLayout(caja)
        cl.setContentsMargins(0, 0, 0, 0)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Paso:")))
        self._spin = QtWidgets.QSpinBox()
        self._spin.setRange(0, max(0, len(res.get("mascaras", [])) - 1))
        f.addWidget(self._spin)
        self._lab = QtWidgets.QLabel("")
        f.addWidget(self._lab, 1)
        cl.addLayout(f)
        self._vista = QtInteractor(caja)
        self._vista.set_background("white")
        # La misma iluminacion que los paneles principales: sin ella la malla
        # sale como una silueta plana y no se ve ninguna trabecula.
        self._vista.enable_lightkit()
        cl.addWidget(self._vista.interactor, 1)
        b = QtWidgets.QPushButton(_("Exportar tabla CSV…"))
        b.clicked.connect(self._csv)
        cl.addWidget(b)
        div.addWidget(caja)
        div.setSizes([640, 560])
        lay.insertWidget(0, div, 1)

        self._primera = True
        self._spin.valueChanged.connect(self._pintar)
        if res.get("mascaras"):
            self._spin.setValue(len(res["mascaras"]) - 1)
            self._pintar(self._spin.value())

    # -- figuras -----------------------------------------------------------

    def _figura_perdida(self, fig, filas):
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for fase in ("perdida", "recuperacion"):
            serie = [f for f in filas if f["fase"] in (fase, "inicial")]
            if fase == "recuperacion":
                ult = [f for f in filas if f["fase"] == "perdida"]
                serie = ult[-1:] + [f for f in filas if f["fase"] == fase]
            if len(serie) < 2:
                continue
            x = [f["BVTV"] for f in serie]
            y = [f["E_rel"] for f in serie]
            ax.plot(x, y, "o-" if fase == "perdida" else "s--",
                    color=self.COLOR[fase], lw=1.8, ms=6,
                    label=_("perdida") if fase == "perdida"
                    else _("recuperacion"))
            for f in serie:
                ax.annotate(str(f["paso"]), (f["BVTV"], f["E_rel"]),
                            textcoords="offset points", xytext=(4, 5),
                            fontsize=7, color=self.COLOR[fase])
        ax.axhline(1.0, color="#999", lw=0.8)
        ax.set_xlabel("BV/TV")
        ax.set_ylabel(_("E / E inicial"))
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, fontsize=8)

        pasos = [f["paso"] for f in filas]
        f0 = filas[0]

        def rel(k):
            v0 = f0.get(k)
            if v0 is None or not np.isfinite(v0) or v0 == 0:
                return None
            return [f.get(k, np.nan) / v0 for f in filas]

        for clave, etq, color, ls in (
                ("E_rel", _("rigidez"), "#1f4e9c", "-"),
                ("sigma_fallo_rel", _("carga de fallo"), "#b62324", "-"),
                ("ConnD", "Conn.D", "#b06000", "--"),
                ("TbN", "Tb.N", "#6b3fa0", ":"),
                ("BVTV", "BV/TV", "#888888", "-.")):
            y = ([f.get(clave, np.nan) for f in filas]
                 if clave.endswith("_rel") else rel(clave))
            if y is not None:
                ax2.plot(pasos, y, ls, color=color, lw=1.6, label=etq)
        ax2.set_xlabel(_("paso"))
        ax2.set_ylabel(_("relativo al paso 0"))
        ax2.set_ylim(bottom=0)
        ax2.legend(frameon=False, fontsize=8, ncol=3)
        for a in (ax, ax2):
            a.spines[["top", "right"]].set_visible(False)
            a.grid(alpha=0.25)

    def _figura_fallo(self, fig, filas):
        ax = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ok = [f for f in filas if f["ok"]]
        if ok:
            ax.plot([f["paso"] for f in ok], [f["F_fallo"] for f in ok],
                    "o-", color="#b62324", lw=1.8, ms=6)
            antes = [f for f in ok if not f.get("tras_colapso")]
            if antes:
                pico = max(antes, key=lambda f: f["F_fallo"])
                ax.plot([pico["paso"]], [pico["F_fallo"]], "*",
                        color="#b62324", ms=14)
            tras = [f["paso"] for f in filas if f.get("tras_colapso")]
            if tras:
                # Tras el colapso la carga que sube no es resistencia: la
                # sostiene el tejido ablandado. Se sombrea para que no se lea.
                ax.axvspan(min(tras) - 0.5, max(f["paso"] for f in filas) + 0.5,
                           color="#999999", alpha=0.18, lw=0,
                           label=_("tras el colapso"))
                ax.legend(frameon=False, fontsize=8)
            ax2.plot([100 * f["dano_acumulado"] for f in ok],
                     [f["E_rel"] for f in ok], "s-", color="#1f4e9c", lw=1.6,
                     label=_("rigidez"))
            ax2.plot([100 * f["dano_acumulado"] for f in ok],
                     [f["F_rel"] for f in ok], "o--", color="#b62324", lw=1.4,
                     label=_("carga de fallo"))
            ax2.legend(frameon=False, fontsize=8)
        ax.set_xlabel(_("paso"))
        ax.set_ylabel(_("carga de fallo (Pistoia) [N]"))
        ax2.set_xlabel(_("tejido danado acumulado [%]"))
        ax2.set_ylabel(_("relativo al paso 0"))
        ax2.set_ylim(bottom=0)
        for a in (ax, ax2):
            a.spines[["top", "right"]].set_visible(False)
            a.grid(alpha=0.25)

    # -- texto -------------------------------------------------------------

    def _html(self, res):
        p = res.get("parametros", {})
        s = res.get("resumen", {})
        filas = []
        if self._fallo:
            cab = ("<tr><th>paso</th><th>dano</th><th>BV/TV</th>"
                   "<th>E<sub>app</sub></th><th>E/E0</th><th>F fallo</th>"
                   "<th>rotos</th><th>frac.port.</th><th>residuo</th></tr>")
            for f in res["pasos"]:
                if not f["ok"]:
                    filas.append(f"<tr><td>{f['paso']}</td><td colspan=8 "
                                 f"style='color:#b62324'>{f['msg']}</td></tr>")
                    continue
                filas.append(
                    f"<tr><td>{f['paso']}</td>"
                    f"<td>{100*f['dano_acumulado']:.1f} %</td>"
                    f"<td>{f['BVTV']:.3f}</td><td>{f['E_app']/1e6:.1f} MPa</td>"
                    f"<td>{f['E_rel']:.3f}</td><td>{f['F_fallo']:.1f} N</td>"
                    f"<td>{f['rotos']}</td>"
                    f"<td>{(f['frac_portante'] or 0):.3f}</td>"
                    f"<td>{(f['residuo'] or 0):.1e}</td></tr>")
            cabecera = _("<p><b>Carga de fallo estimada</b> {F} N en el paso "
                         "{k}, con un {d} % del tejido danado.</p>").format(
                F=f"{s.get('F_max', float('nan')):.1f}", k=s.get("paso_F_max"),
                d=f"{100*s.get('dano_en_F_max', float('nan')):.1f}")
            c = s.get("colapso_en_paso")
            if c is not None:
                cabecera += _("<p>La estructura <b>colapsa en el paso {c}</b>: "
                              "su rigidez cae por debajo del {u} % de la "
                              "inicial. Cuanto antes colapsa con menos dano, "
                              "mas fragil es.</p>").format(
                    c=c, u=int(round(100 * s.get("E_rel_colapso", 0.5))))
            else:
                cabecera += _("<p>La estructura no colapsa en los pasos "
                              "simulados: reparte el dano.</p>")
            tbh = s.get("TbTh_h_mec", float("nan"))
            if np.isfinite(tbh) and tbh < TBTH_H_MIN:
                cabecera += (
                    "<p style='color:#b62324'>"
                    + _("La malla mecanica tiene {v} elementos por trabecula "
                        "(Tb.Th/h), por debajo de {u}: los cocientes de "
                        "rigidez pueden ser artefacto de resolucion. Sube la "
                        "resolucion del ensayo.").format(
                        v=f"{tbh:.1f}", u=f"{TBTH_H_MIN:.1f}")
                    + "</p>")
            nota = _("En cada paso el tejido que el criterio de Pistoia da "
                     "por roto conserva solo el 5 % de su rigidez, y se vuelve "
                     "a cargar. Tras el colapso (zona gris) la carga calculada "
                     "vuelve a subir porque la sostiene ese tejido ablandado, "
                     "no resistencia real: por eso la carga de fallo se toma "
                     "antes. Elastico lineal con dano en un escalon, sin "
                     "plasticidad ni contacto. El 3D muestra el tejido "
                     "intacto. Compara estructuras; no da una carga de rotura "
                     "absoluta.")
        else:
            cab = ("<tr><th>paso</th><th>fase</th><th>hueso</th><th>BV/TV</th>"
                   "<th>Tb.Th</th><th>Tb.N</th><th>Conn.D</th>"
                   "<th>E<sub>app</sub></th><th>E/E0</th><th>&sigma; fallo</th>"
                   "<th>frac.port.</th><th>residuo</th></tr>")
            for f in res["pasos"]:
                mec = (f"<td>{f['E_app']/1e6:.1f} MPa</td><td>{f['E_rel']:.3f}</td>"
                       f"<td>{f['sigma_fallo']/1e6:.2f} MPa</td>"
                       f"<td>{(f['frac_portante'] or 0):.3f}</td>"
                       f"<td>{(f['residuo'] or 0):.1e}</td>"
                       if f["ok"] else
                       f"<td colspan=5 style='color:#b62324'>{f['msg']}</td>")
                filas.append(
                    f"<tr><td>{f['paso']}</td><td>{_(f['fase'])}</td>"
                    f"<td>{100*f['hueso_rel']:.0f} %</td><td>{f['BVTV']:.3f}</td>"
                    f"<td>{f['TbTh']:.3f}</td><td>{f['TbN']:.2f}</td>"
                    f"<td>{f['ConnD']:.2f}</td>" + mec + "</tr>")
            cabecera = _("<p><b>Al final:</b> rigidez {E} y conectividad {c} "
                         "de las iniciales, con BV/TV {b} del inicial.</p>"
                         ).format(E=f"{s.get('E_rel_final', float('nan')):.2f}",
                                  c=f"{s.get('ConnD_rel_final', float('nan')):.2f}",
                                  b=f"{s.get('BVTV_rel_final', float('nan')):.2f}")
            if "E_rel_recuperado" in s:
                cabecera += _("<p><b>Recuperacion:</b> con la masa osea de "
                              "vuelta al {h} %, la rigidez vuelve a {E} de la "
                              "inicial (su minimo fue {m}).</p>").format(
                    h=f"{100*s['hueso_rel_recuperado']:.0f}",
                    E=f"{s['E_rel_recuperado']:.2f}",
                    m=f"{s.get('E_rel_minimo', float('nan')):.2f}")
                # La lectura depende del resultado: no se afirma de antemano
                # que engrosar no recupera. En la primera prueba del dialogo
                # la rigidez volvio POR ENCIMA de la inicial.
                e = s["E_rel_recuperado"]
                if np.isfinite(e) and e < 0.95:
                    cabecera += _("<p>A igual masa osea, falta rigidez: es "
                                  "arquitectura perdida. Engrosar las "
                                  "trabeculas que quedan no devuelve las que "
                                  "desaparecieron.</p>")
                elif np.isfinite(e) and e > 1.05:
                    cabecera += _("<p>A igual masa osea, la rigidez supera la "
                                  "inicial en esta direccion: el hueso se "
                                  "redistribuyo de trabeculas finas a "
                                  "gruesas. Comprueba otra direccion antes de "
                                  "leerlo como mejora: puede haberse perdido "
                                  "rigidez transversal.</p>")
            tbh = s.get("TbTh_h_mec", float("nan"))
            if np.isfinite(tbh) and tbh < TBTH_H_MIN:
                cabecera += (
                    "<p style='color:#b62324'>"
                    + _("La malla mecanica tiene {v} elementos por trabecula "
                        "(Tb.Th/h), por debajo de {u}: los cocientes de "
                        "rigidez pueden ser artefacto de resolucion. Sube la "
                        "resolucion del ensayo.").format(
                        v=f"{tbh:.1f}", u=f"{TBTH_H_MIN:.1f}")
                    + "</p>")
            nota = _("Todos los protocolos retiran la misma cantidad de hueso "
                     "por paso: a igual paso, igual BV/TV, y las diferencias "
                     "de rigidez entre protocolos son de arquitectura. Un paso "
                     "es una cantidad de hueso, no un tiempo. E<sub>app</sub> "
                     "absoluto a esta resolucion no es citable; los cocientes "
                     "respecto al paso 0 si, porque todos los pasos comparten "
                     "el sesgo de malla.")
        pie = (f"<p style='color:#666;font-size:10px'>{nota}<br>"
               + _("Estructura: {e} · malla mecanica {n}³ · apoyo {a} · "
                   "eje {x} · {t} s").format(
                   e=res.get("estructura", ""), n=p.get("n_mec"),
                   a=_(p.get("apoyo") or ""), x="XYZ"[int(p.get("eje", 2))],
                   t=res.get("tiempo_s")) + "</p>")
        return (cabecera + "<table cellpadding=3 cellspacing=0 border=1 "
                "style='border-collapse:collapse;font-size:10px'>" + cab
                + "".join(filas) + "</table>" + pie)

    # -- 3D y exportacion --------------------------------------------------

    def _pintar(self, i):
        masc = self._res.get("mascaras") or []
        if not (0 <= i < len(masc)):
            return
        f = self._res["pasos"][i]
        self._vista.clear()
        # `clear()` se lleva tambien las luces, igual que en los paneles
        # principales: sin volver a ponerlas la malla es una silueta plana.
        self._vista.enable_lightkit()
        malla = malla_de_mascara(masc[i], np.asarray(self._res["spacing"]))
        if malla.n_points:
            self._vista.add_mesh(malla, color="#d8d2c4", smooth_shading=True,
                                 specular=0.25)
        if self._primera:
            self._vista.camera_position = "iso"
            self._vista.reset_camera()
            self._primera = False
        self._vista.render()
        e = f.get("E_rel", float("nan"))
        self._lab.setText(_("BV/TV {b} · E/E0 {e}").format(
            b=f"{f['BVTV']:.3f}",
            e=f"{e:.3f}" if np.isfinite(e) else "—"))

    def _csv(self):
        ruta, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Exportar tabla"), str(DATOS / "simulacion.csv"),
            _("CSV (*.csv)"))
        if not ruta:
            return
        import csv
        filas = self._res["pasos"]
        claves = [k for k in filas[0] if np.ndim(filas[0][k]) == 0]
        with open(ruta, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=claves, extrasaction="ignore")
            w.writeheader()
            for f in filas:
                w.writerow({k: f.get(k) for k in claves})
        p = self.parent()
        if p is not None:
            p.statusBar().showMessage(_("Tabla guardada en {f}").format(f=ruta))

    def done(self, r):
        # El QtInteractor retiene un contexto de VTK: cerrarlo a mano evita
        # que se quede vivo al cerrar el dialogo y reviente al salir.
        try:
            self._vista.close()
        except Exception:
            pass
        super().done(r)


class Deslizador(QtWidgets.QWidget):
    """Deslizador con etiqueta y valor, en unidades reales (no pasos enteros).

    ALTURA MINIMA FIJA. Sin ella, cuando el panel no cabe a lo alto Qt reparte
    el deficit comprimiendo cada hijo, y la etiqueta acaba solapada sobre la
    barra: ilegible. Con un minimo, el panel prefiere desbordar —y entonces la
    QScrollArea que lo envuelve muestra una barra de desplazamiento— antes que
    aplastar los controles.
    """
    cambiado = QtCore.pyqtSignal()
    ALTO_MIN = 46

    def __init__(self, texto, minimo, maximo, valor, decimales=0, sufijo=""):
        super().__init__()
        self._esc = 10 ** decimales
        self._dec = decimales
        self._suf = sufijo

        self.setMinimumHeight(self.ALTO_MIN)
        self.setSizePolicy(QtWidgets.QSizePolicy.Preferred,
                           QtWidgets.QSizePolicy.Fixed)

        self.lab = QtWidgets.QLabel()
        # Su texto lleva el VALOR dentro, asi que no puede capturarse como
        # cadena fija. Se marca para que el recorrido de traduccion lo salte
        # y se retraduce por su cuenta (`retraducir`).
        self.lab.setProperty("i18n_dinamico", True)
        self.sld = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sld.setMinimumHeight(20)
        self.sld.setMinimum(int(round(minimo * self._esc)))
        self.sld.setMaximum(int(round(maximo * self._esc)))
        self.sld.setValue(int(round(valor * self._esc)))
        self.sld.valueChanged.connect(self._actualiza)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(1)
        lay.addWidget(self.lab)
        lay.addWidget(self.sld)
        self._texto = texto
        self._actualiza()

    def _actualiza(self):
        self._pinta()
        self.cambiado.emit()

    def _pinta(self):
        self.lab.setText(
            f"{_(self._texto)}: {self.valor():.{self._dec}f}{self._suf}")

    def retraducir(self):
        """Reescribe la etiqueta SIN emitir `cambiado`.

        La senal esta conectada a `_pedir_vista` y a `_soltar_R`: emitirla al
        cambiar de idioma relanzaria el calculo de la vista y, peor, liberaria
        la orientacion que impuso el ajuste. Cambiar de idioma no puede tener
        efectos secundarios.
        """
        self._pinta()

    def valor(self):
        return self.sld.value() / self._esc

    def fijar(self, v):
        """Coloca un valor, acotandolo al rango del control.

        Se acota en vez de fallar porque una sesion guardada con otra version
        puede traer un valor fuera de los limites de ahora; quedarse en el
        extremo y seguir es preferible a no cargar nada. Qt no emite
        `valueChanged` si el valor no cambia, asi que el aviso al usuario no
        puede colgar de la senal.
        """
        self.sld.setValue(int(round(
            min(max(float(v), self.sld.minimum() / self._esc),
                self.sld.maximum() / self._esc) * self._esc)))


# ---------------------------------------------------------------------------
# Pantalla de bienvenida
# ---------------------------------------------------------------------------

# Nombres de quienes han contribuido sin figurar como autores del software.
# Vacio mientras no haya una lista acordada: esta pantalla NO es sitio para
# poner a nadie de memoria. Los AUTORES no se escriben aqui -salen de
# `informe.AUTORES`, que `tests/test_19` cuadra con `CITATION.cff`-, para que
# no puedan discrepar de la cita del software.
COLABORADORES = ()


def _version_spinpy():
    try:
        import spinpy
        return getattr(spinpy, "__version__", "?")
    except Exception:
        return "?"


class DialogoBienvenida(QtWidgets.QDialog):
    """Portada: que es esto, quien lo firma y por donde se entra.

    Se abre ANTES de construir la ventana principal, no despues. Construir el
    `Visor` entero -con su render de VTK- tarda unos segundos, y durante ese
    rato la pantalla se queda vacia; asi la primera respuesta es inmediata.

    Las caracteristicas que se enumeran son las que la aplicacion tiene DE
    VERDAD y estan verificadas: esto es la portada de una herramienta de
    medida, no un folleto. El enlace al repositorio se dibuja a partir de
    `informe.REPOSITORIO`; mientras siga siendo `[PENDIENTE]` sale como texto
    apagado, y el dia que se rellene esa constante se convierte en enlace solo,
    sin tocar esta clase. Los autores salen de `informe.AUTORES`, que
    `tests/test_19` cuadra con `CITATION.cff`, para que la portada y la cita
    del software no puedan discrepar.
    """

    CLAVE = "mostrar_bienvenida"
    ANCHO_IZQ = 300
    ANCHO = 860

    def __init__(self, padre=None):
        super().__init__(padre)
        self.setModal(True)
        self.setStyleSheet("QDialog { background: white; }")
        self._ico = _ruta_icono()
        if self._ico is not None:
            self.setWindowIcon(QtGui.QIcon(str(self._ico)))
        _icono_de_clase(self)

        self._h = QtWidgets.QHBoxLayout(self)
        self._h.setContentsMargins(0, 0, 0, 0)
        self._h.setSpacing(0)
        self.setFixedWidth(self.ANCHO)
        self._construir()

    def _construir(self):
        """Pone las dos columnas en el idioma activo.

        Se llama al crear el dialogo y otra vez al cambiar de idioma desde el
        propio selector. Rehacer las columnas enteras es mas corto y mas
        seguro que capturar y retraducir: aqui casi todo son QLabel con HTML
        compuesto a mano (clave en negrita + texto), que `idioma.capturar` no
        sabria devolver al original.
        """
        marcado = getattr(self, "chk_no", None) is not None \
            and self.chk_no.isChecked()
        while self._h.count():
            it = self._h.takeAt(0)
            if it.widget() is not None:
                it.widget().deleteLater()
        self.setWindowTitle(_("Bienvenido a spinpy"))
        self._h.addWidget(self._identidad(self._ico))
        self._h.addWidget(self._contenido(), 1)
        self.chk_no.setChecked(marcado)

    # -- columnas ----------------------------------------------------------
    def _identidad(self, ico):
        """Columna oscura: logo, nombre, version, autores y repositorio."""
        w = QtWidgets.QFrame()
        w.setFixedWidth(self.ANCHO_IZQ)
        w.setStyleSheet("QFrame { background: #16324a; }")
        v = QtWidgets.QVBoxLayout(w)
        v.setContentsMargins(26, 30, 26, 22)
        v.setSpacing(8)

        if ico is not None:
            lg = QtWidgets.QLabel()
            # El .ico ya viaja dentro del paquete y ya lo localiza
            # `_ruta_icono`; pedirle a QIcon un pixmap evita anadir un PNG
            # suelto que habria que acordarse de meter en `spinpy.spec`.
            lg.setPixmap(QtGui.QIcon(str(ico)).pixmap(96, 96))
            lg.setFixedSize(96, 96)
            v.addWidget(lg, 0, QtCore.Qt.AlignLeft)
            v.addSpacing(6)

        v.addWidget(self._txt("spinpy", "color:white; font-size:28pt;"
                              " font-weight:bold;", False))
        v.addWidget(self._txt(_("versión") + " " + _version_spinpy(),
                              "color:#8fb3cc; font-size:10pt;", False))
        v.addSpacing(18)
        v.addWidget(self._txt(_("AUTORES"),
                              "color:#5f93bb; font-size:9pt; font-weight:bold;"
                              " letter-spacing:1px;", False))
        v.addWidget(self._txt(
            "<br>".join(a["nombre"] + " " + a["apellidos"]
                        for a in informe_pub.AUTORES),
            "color:#dbe8f2; font-size:10pt;"))
        v.addSpacing(10)
        # Nombre propio de la institucion: no se traduce.
        v.addWidget(self._txt("Universidad de Valparaíso",
                              "color:#8fb3cc; font-size:9pt;"))
        v.addStretch(1)
        v.addWidget(self._repositorio())
        return w

    def _contenido(self):
        """Columna clara: lema, lo que hace y la puerta de entrada."""
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        v.setContentsMargins(30, 30, 30, 20)
        v.setSpacing(10)

        v.addWidget(self._txt(
            _("Ajuste de microestructuras espinodales a VOIs de hueso "
              "trabecular"),
            "color:#16324a; font-size:13pt; font-weight:bold;"))
        ln = QtWidgets.QFrame()
        ln.setFixedHeight(1)
        ln.setStyleSheet("background: #d5dde4;")
        v.addWidget(ln)
        v.addSpacing(4)

        # Clave y texto van como literales sueltos dentro de `_()` a proposito:
        # `idioma_revisar.py` recoge las cadenas de marcha leyendo el CODIGO
        # (AST) y solo ve las que se le pasan literalmente. Sacarlas a una
        # constante las volveria invisibles para el comprobador y sus
        # traducciones se darian por sobrantes.
        for clave, texto in (
            (_("Dos familias"),
             _("Spinodoide (campo aleatorio gaussiano) y dual-lattice, sobre "
               "la misma rejilla de vóxeles. Se ajustan y se comparan a la "
               "vez.")),
            (_("Morfometría"),
             _("BV/TV, BS/BV, Tb.Th, Tb.Sp, Tb.N, anisotropía por tensor MIL, "
               "Conn.D, SMI, tamaño de poro y factor de elipsoide.")),
            (_("Mecánica"),
             _("Homogeneización elástica periódica, ensayo de compresión y "
               "carga de fallo por el criterio de Pistoia.")),
            (_("In silico"),
             _("Simulaciones de pérdida ósea y de fallo progresivo sobre el "
               "gemelo digital del VOI.")),
            (_("Informe"),
             _("Dice de cada número si es citable, con qué reservas y con qué "
               "parámetros se obtuvo.")),
            (_("Exportación"),
             _("STL, VTU, Abaqus, ANSYS y FEBio.")),
        ):
            v.addWidget(self._txt(
                "<span style='color:#1f5c8b; font-weight:bold;'>" + clave
                + "</span> &nbsp;" + texto, "font-size:10pt;"))

        v.addStretch(1)
        fila = QtWidgets.QHBoxLayout()
        self.chk_no = QtWidgets.QCheckBox(_("No volver a mostrar"))
        # Fija: sin esto el layout la estrecha antes que al boton y el texto
        # sale cortado ("No volver a mostra").
        self.chk_no.setSizePolicy(QtWidgets.QSizePolicy.Fixed,
                                  QtWidgets.QSizePolicy.Fixed)
        fila.addWidget(self.chk_no)
        fila.addSpacing(14)
        # El idioma se elige AQUI y no solo en el menu de la ventana principal:
        # la portada es lo primero que se ve, antes de que ese menu exista.
        fila.addWidget(QtWidgets.QLabel(_("Idioma")))
        self.cmb_idioma = QtWidgets.QComboBox()
        for cod, nombre in IDIOMAS.items():
            self.cmb_idioma.addItem(nombre, cod)
        self.cmb_idioma.setCurrentIndex(
            max(0, self.cmb_idioma.findData(idioma_actual())))
        self.cmb_idioma.currentIndexChanged.connect(self._cambiar_idioma)
        fila.addWidget(self.cmb_idioma)
        fila.addStretch(1)
        b = QtWidgets.QPushButton(_("Entrar") + "  →")
        b.setMinimumHeight(34)
        b.setMinimumWidth(120)
        b.setDefault(True)
        b.clicked.connect(self.accept)
        fila.addWidget(b)
        v.addLayout(fila)
        return w

    # -- piezas menores ----------------------------------------------------
    @staticmethod
    def _txt(t, css="", wrap=True):
        lb = QtWidgets.QLabel(t)
        lb.setTextFormat(QtCore.Qt.RichText)
        lb.setWordWrap(wrap)
        if css:
            lb.setStyleSheet(css)
        return lb

    def _repositorio(self):
        repo = informe_pub.REPOSITORIO
        if repo and repo != informe_pub.PENDIENTE:
            lb = self._txt(_("Código:") + ' <a href="' + repo + '">' + repo
                           + "</a>", "color:#7fa6c4; font-size:8pt;")
            lb.setOpenExternalLinks(True)
            return lb
        return self._txt(
            _("Código: el repositorio público aún no está publicado; el "
              "enlace aparecerá aquí cuando lo esté."),
            "color:#7fa6c4; font-size:8pt;")

    def _cambiar_idioma(self, _indice):
        """Retraduce la portada y deja el idioma elegido como el guardado.

        Con la ventana principal ya abierta (Sesion -> Acerca de) se traduce
        tambien esa, por su camino de siempre: si solo se cambiara el idioma
        global, la ventana quedaria a medias entre los dos. Al arrancar no hay
        ventana todavia; basta con guardar el ajuste, que `Visor.__init__` lee.
        """
        cod = self.cmb_idioma.itemData(self.cmb_idioma.currentIndex())
        if cod not in IDIOMAS or cod == idioma_actual():
            return
        padre = self.parent()
        if isinstance(padre, Visor):
            padre.cambiar_idioma(cod)
        else:
            fijar_idioma(cod)
            QtCore.QSettings("spinpy", "visor").setValue("idioma", cod)
        # Se difiere: este metodo lo llama una senal del combo que
        # `_construir` va a destruir.
        QtCore.QTimer.singleShot(0, self._reconstruir)

    def _reconstruir(self):
        self._construir()
        QtCore.QTimer.singleShot(0, self._medir)

    # -- entrada y salida --------------------------------------------------
    def showEvent(self, ev):
        super().showEvent(ev)
        if not getattr(self, "_medido", False):
            self._medido = True
            self._medir()

    def _medir(self):
        # Un QLabel con ajuste de linea no sabe su altura hasta que conoce su
        # ancho, y el layout ya ha decidido la del dialogo: sin esta segunda
        # pasada las ultimas lineas de cada parrafo salen cortadas. Se repite
        # tras cambiar de idioma, porque el ingles corta las lineas en otro
        # sitio.
        self.layout().activate()
        for lb in self.findChildren(QtWidgets.QLabel):
            if lb.wordWrap() and lb.width() > 0:
                lb.setMinimumHeight(lb.heightForWidth(lb.width()))
        self.adjustSize()

    def closeEvent(self, ev):
        # Cerrar con la cruz entra igual que pulsar el boton: la portada
        # informa, no guarda la puerta. El estado se guarda AQUI y no en
        # `accept`, porque este es el unico camino por el que pasan los dos.
        QtCore.QSettings("spinpy", "visor").setValue(
            self.CLAVE, "no" if self.chk_no.isChecked() else "si")
        super().closeEvent(ev)

    @classmethod
    def mostrar_al_arrancar(cls):
        return str(QtCore.QSettings("spinpy", "visor").value(
            cls.CLAVE, "si")).lower() != "no"


class DialogoInformeAuto(QtWidgets.QDialog):
    """Que etapas corre el informe automatico, cuanto tardara y que figuras.

    Las resoluciones y casillas NO son copias propias: salen de los controles
    del panel y, al aceptar, se escriben de vuelta en ellos. Asi las etapas
    leen los mismos controles que al pulsar sus botones, y la sesion
    guardada o exportada declara los valores con los que se calculo.

    TIEMPO ESTIMADO. Cada etapa lleva al lado lo que se espera que tarde, y el
    total se rehace con cada cambio (`spinpy.tiempos`). El modelo esta medido
    en el equipo de desarrollo y se corrige con lo que tardaron las etapas en
    ESTE equipo en ejecuciones anteriores; la nota de abajo dice cuantas hay.

    La carpeta se pide AQUI, al principio, y no al final: la cadena puede
    tardar horas y tiene que poder terminar sin nadie delante.
    """

    def __init__(self, padre):
        super().__init__(padre)
        from spinpy import figuras as F
        self.v = padre
        self.F = F
        self.setWindowTitle(_("Informe para publicacion (automatico)"))
        self._etq = {}
        raiz = QtWidgets.QVBoxLayout(self)

        self._ocupado_act = False
        self.intro = QtWidgets.QLabel()
        self.intro.setWordWrap(True)
        self.intro.setTextFormat(QtCore.Qt.RichText)
        raiz.addWidget(self.intro)

        # Familia: la misma eleccion que el selector de la barra, y al aceptar
        # se escribe en el. «Ambas» corre cada etapa para las dos familias y
        # el informe las compara con el VOI en una sola tabla.
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Familia:")))
        self.cmb_fam = QtWidgets.QComboBox()
        for modo, etq in zip(MODOS_FAMILIA, (_("Spinodoide"),
                                             _("Dual-lattice"),
                                             _("Ambas (VOI / Spinodoide / "
                                               "Dual-lattice)"))):
            self.cmb_fam.addItem(etq, modo)
        self.cmb_fam.setCurrentIndex(MODOS_FAMILIA.index(padre._modo_fam))
        self.cmb_fam.setToolTip(_(
            "Con «Ambas» cada etapa se corre para las dos familias, las "
            "figuras las muestran lado a lado y el informe añade una tabla "
            "VOI / Spinodoide / Dual-lattice. El tiempo casi se duplica."))
        f.addWidget(self.cmb_fam, 1)
        f.addStretch(2)
        raiz.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Carpeta del informe:")))
        self.ed_carpeta = QtWidgets.QLineEdit()
        self.ed_carpeta.setReadOnly(True)
        f.addWidget(self.ed_carpeta, 1)
        b = QtWidgets.QPushButton(_("Carpeta…"))
        b.clicked.connect(self._elegir_carpeta)
        f.addWidget(b)
        raiz.addLayout(f)

        columnas = QtWidgets.QHBoxLayout()
        izq = QtWidgets.QVBoxLayout()
        der = QtWidgets.QVBoxLayout()
        columnas.addLayout(izq, 3)
        columnas.addLayout(der, 2)
        raiz.addLayout(columnas)

        # -- etapas que gradua el informe --
        g = QtWidgets.QGroupBox(_("Etapas que gradua el informe"))
        gl = QtWidgets.QGridLayout(g)
        gl.setColumnStretch(0, 1)
        fila = 0

        # Sin un ajuste previo de TODAS las familias elegidas no hay candidato
        # que proceda del VOI: la casilla solo se puede quitar si ya lo hay
        # (`_actualizar`, porque depende de la familia).
        self.chk_ajuste = QtWidgets.QCheckBox(_(
            "Mejor ajuste al VOI: el de menor error entre rapido, completo y "
            "equitativo"))
        self.chk_ajuste.setChecked(True)
        gl.addWidget(self.chk_ajuste, fila, 0, 1, 3)
        self._tiempo(gl, fila, "ajuste")
        fila += 1

        self.chk_morfo = QtWidgets.QCheckBox(_(
            "Morfometria completa de la estructura ajustada y del VOI"))
        self.chk_morfo.setChecked(True)
        gl.addWidget(self.chk_morfo, fila, 0, 1, 3)
        self._tiempo(gl, fila, "morfometria")
        fila += 1
        sub = QtWidgets.QHBoxLayout()
        sub.setContentsMargins(24, 0, 0, 0)
        self.chk_extra = QtWidgets.QCheckBox(_("Conn.D y SMI"))
        self.chk_extra.setChecked(True)
        self.chk_poro = QtWidgets.QCheckBox(_("Po.Dm"))
        self.chk_poro.setChecked(True)
        self.chk_ef = QtWidgets.QCheckBox(_("EF (muy caro)"))
        self.chk_ef.setChecked(padre.chk_ef.isChecked())
        self.chk_ef.setToolTip(_(
            "La medida mas cara del programa: minutos por estructura a 64³ y "
            "media hora a 97³. No se recorta para ir mas rapido."))
        for c in (self.chk_extra, self.chk_poro, self.chk_ef):
            sub.addWidget(c)
            self.chk_morfo.toggled.connect(c.setEnabled)
        sub.addStretch(1)
        gl.addLayout(sub, fila, 0, 1, 3)
        fila += 1

        self.chk_ela = QtWidgets.QCheckBox(_("Tensor elastico (homogeneizacion "
                                             "periodica)"))
        self.chk_ela.setChecked(True)
        gl.addWidget(self.chk_ela, fila, 0)
        gl.addWidget(QtWidgets.QLabel(_("Resolucion:")), fila, 1)
        self.spin_homog = QtWidgets.QSpinBox()
        self.spin_homog.setRange(padre.spin_homog.minimum(),
                                 padre.spin_homog.maximum())
        self.spin_homog.setValue(padre.spin_homog.maximum())
        self.spin_homog.setSuffix(" vox")
        gl.addWidget(self.spin_homog, fila, 2)
        self._tiempo(gl, fila, "elastico")
        fila += 1

        self.chk_fe = QtWidgets.QCheckBox(_("Ensayo de compresion y fallo "
                                            "(Pistoia)"))
        self.chk_fe.setChecked(True)
        gl.addWidget(self.chk_fe, fila, 0)
        gl.addWidget(QtWidgets.QLabel(_("Resolucion:")), fila, 1)
        self.spin_fe = QtWidgets.QSpinBox()
        self.spin_fe.setRange(padre.spin_res_fe.minimum(),
                              padre.spin_res_fe.maximum())
        self.spin_fe.setValue(max(padre.spin_res_fe.value(),
                                  informe_pub.N_MECANICA_MIN))
        self.spin_fe.setSuffix(" vox")
        gl.addWidget(self.spin_fe, fila, 2)
        self._tiempo(gl, fila, "ensayo")
        fila += 1
        sub = QtWidgets.QHBoxLayout()
        sub.setContentsMargins(24, 0, 0, 0)
        sub.addWidget(QtWidgets.QLabel(_("Direccion:")))
        self.cmb_eje = QtWidgets.QComboBox()
        self.cmb_eje.addItems([padre.cmb_eje_fe.itemText(i)
                               for i in range(padre.cmb_eje_fe.count())])
        self.cmb_eje.setCurrentIndex(padre.cmb_eje_fe.currentIndex())
        sub.addWidget(self.cmb_eje, 1)
        sub.addWidget(QtWidgets.QLabel(_("Apoyo:")))
        self.cmb_apoyo = QtWidgets.QComboBox()
        self.cmb_apoyo.addItems([padre.cmb_apoyo.itemText(i)
                                 for i in range(padre.cmb_apoyo.count())])
        self.cmb_apoyo.setCurrentIndex(padre.cmb_apoyo.currentIndex())
        sub.addWidget(self.cmb_apoyo, 1)
        gl.addLayout(sub, fila, 0, 1, 3)
        fila += 1

        self.chk_comp = QtWidgets.QCheckBox(_(
            "Analisis comparado (protocolo de Tapia et al. 2026), a la "
            "resolucion del ensayo"))
        self.chk_comp.setChecked(True)
        gl.addWidget(self.chk_comp, fila, 0, 1, 3)
        self._tiempo(gl, fila, "comparado")
        fila += 1

        self.lab_res = QtWidgets.QLabel()
        self.lab_res.setWordWrap(True)
        self.lab_res.setStyleSheet("color:#9a6700; font-size:10px;")
        gl.addWidget(self.lab_res, fila, 0, 1, 4)
        izq.addWidget(g)

        # -- etapas que el informe aun no gradua --
        g = QtWidgets.QGroupBox(_("Etapas opcionales (el informe las "
                                  "gradua y las dibuja si se corren)"))
        gl = QtWidgets.QGridLayout(g)
        gl.setColumnStretch(0, 1)
        self.chk_disp = QtWidgets.QCheckBox(_("Dispersion entre semillas"))
        gl.addWidget(self.chk_disp, 0, 0)
        gl.addWidget(QtWidgets.QLabel("K:"), 0, 1)
        self.spin_k = QtWidgets.QSpinBox()
        self.spin_k.setRange(padre.spin_k.minimum(), padre.spin_k.maximum())
        self.spin_k.setValue(padre.spin_k.value())
        gl.addWidget(self.spin_k, 0, 2)
        self._tiempo(gl, 0, "dispersion")
        self.chk_conv = QtWidgets.QCheckBox(_(
            "Convergencia de malla del ensayo sobre el VOI"))
        gl.addWidget(self.chk_conv, 1, 0, 1, 3)
        self._tiempo(gl, 1, "convergencia")
        est = padre.cmb_sim_estructura.currentText()
        self.chk_perdida = QtWidgets.QCheckBox(_(
            "Simulacion de perdida osea ({est}, {prot}, {n} pasos)").format(
                est=est, prot=padre.cmb_sim_protocolo.currentText(),
                n=padre.spin_sim_pasos.value()))
        gl.addWidget(self.chk_perdida, 2, 0, 1, 3)
        self._tiempo(gl, 2, "perdida")
        self.chk_fallo = QtWidgets.QCheckBox(_(
            "Fallo progresivo ({est}, {n} pasos)").format(
                est=est, n=padre.spin_sim_pasos.value()))
        gl.addWidget(self.chk_fallo, 3, 0, 1, 3)
        self._tiempo(gl, 3, "fallo")
        nota = QtWidgets.QLabel(_(
            "Las simulaciones usan la estructura, el protocolo y los pasos "
            "de su seccion del panel; cambialos alli antes si hace falta."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota, 4, 0, 1, 4)
        izq.addWidget(g)
        izq.addStretch(1)

        # -- figuras --
        op = opciones_figuras()
        li = 1 if idioma_actual() == "en" else 0     # nombres de figuras.py
        g = QtWidgets.QGroupBox(_("Figuras 3D del informe"))
        gl = QtWidgets.QVBoxLayout(g)
        gl.addWidget(QtWidgets.QLabel("<b>" + _("Estilos") + "</b>"))
        self.chk_estilos = {}
        for clave, e in F.ESTILOS_3D.items():
            c = QtWidgets.QCheckBox(e["nombre"][li])
            c.setChecked(clave in op["estilos"])
            self.chk_estilos[clave] = c
            gl.addWidget(c)
        gl.addWidget(QtWidgets.QLabel("<b>" + _("Vistas") + "</b>"))
        rej = QtWidgets.QGridLayout()
        self.chk_vistas = {}
        for k, (clave, vv) in enumerate(F.VISTAS_3D.items()):
            c = QtWidgets.QCheckBox(vv[2][li])
            c.setChecked(clave in op["vistas"])
            self.chk_vistas[clave] = c
            rej.addWidget(c, k // 2, k % 2)
        gl.addLayout(rej)
        f = QtWidgets.QFormLayout()
        self.cmb_estilo_ppal = QtWidgets.QComboBox()
        for clave, e in F.ESTILOS_3D.items():
            self.cmb_estilo_ppal.addItem(e["nombre"][li], clave)
        self.cmb_estilo_ppal.setCurrentIndex(
            list(F.ESTILOS_3D).index(op["estilo_principal"]))
        self.cmb_vista_ppal = QtWidgets.QComboBox()
        for clave, vv in F.VISTAS_3D.items():
            self.cmb_vista_ppal.addItem(vv[2][li], clave)
        self.cmb_vista_ppal.setCurrentIndex(
            list(F.VISTAS_3D).index(op["vista_principal"]))
        f.addRow(_("Figura principal, estilo:"), self.cmb_estilo_ppal)
        f.addRow(_("Figura principal, vista:"), self.cmb_vista_ppal)
        gl.addLayout(f)
        self.chk_suave = QtWidgets.QCheckBox(_(
            "Suavizar la superficie (Taubin, solo para la figura)"))
        self.chk_suave.setChecked(op["suavizar"])
        self.chk_suave.setToolTip(_(
            "Quita la escalera de los voxeles sin encoger la pieza ni "
            "adelgazar las trabeculas. Las medidas usan siempre la malla sin "
            "suavizar, y el pie de la figura lo dice."))
        gl.addWidget(self.chk_suave)
        gl.addWidget(QtWidgets.QLabel("<b>" + _("Figuras adicionales")
                                      + "</b>"))
        self.chk_dist = QtWidgets.QCheckBox(_(
            "Distribuciones de espesor trabecular y tamaño de poro"))
        self.chk_dist.setChecked(op["distribuciones"])
        self.chk_dist.setToolTip(_(
            "Espesor local por esferas inscritas en hueso y poro. Es la "
            "figura cara: del orden de un minuto para un VOI de 188³."))
        gl.addWidget(self.chk_dist)
        self.chk_vm = QtWidgets.QCheckBox(_(
            "Mapa 3D de von Mises (del análisis comparado)"))
        self.chk_vm.setChecked(op["von_mises"])
        gl.addWidget(self.chk_vm)
        self.chk_metodo = QtWidgets.QCheckBox(_(
            "Figura del método (del micro-CT al candidato)"))
        self.chk_metodo.setChecked(op["metodo"])
        self.chk_metodo.setToolTip(_(
            "Figura 0, en 16:9: rebanada y segmentación, el cubo en la pila, "
            "el VOI y, por cada familia, su aleatoriedad, su campo, el umbral "
            "por densidad y el sólido, con la morfometría frente al VOI. La "
            "fila del micro-CT solo aparece si el VOI se recortó de una pila "
            "en esta sesión."))
        gl.addWidget(self.chk_metodo)
        nota = QtWidgets.QLabel(_(
            "Todas las vistas usan proyeccion paralela y la misma escala. La "
            "figura 4 reune las vistas marcadas en el estilo principal, la 5 "
            "son cortes 2D, la 6 la anisotropia (fabrica y E direccional), y "
            "cada combinacion de estilo y vista queda en figuras/3d/ para "
            "elegir otra sin recalcular."))
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel(_("Informe y figuras")), 1)
        self._etq["informe"] = QtWidgets.QLabel()
        self._etq["informe"].setAlignment(QtCore.Qt.AlignRight)
        f.addWidget(self._etq["informe"])
        gl.addLayout(f)
        der.addWidget(g)
        der.addStretch(1)

        # -- total --
        marco = QtWidgets.QFrame()
        marco.setStyleSheet("QFrame { background:#f3f2ee; border-radius:4px; }")
        fl = QtWidgets.QVBoxLayout(marco)
        self.lab_total = QtWidgets.QLabel()
        self.lab_total.setTextFormat(QtCore.Qt.RichText)
        fl.addWidget(self.lab_total)
        self.lab_nota_tiempo = QtWidgets.QLabel()
        self.lab_nota_tiempo.setWordWrap(True)
        self.lab_nota_tiempo.setStyleSheet("color:#666; font-size:10px;")
        fl.addWidget(self.lab_nota_tiempo)
        raiz.addWidget(marco)

        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel)
        self.b_ok = bb.addButton(_("Empezar"),
                                 QtWidgets.QDialogButtonBox.AcceptRole)
        self.b_ok.setEnabled(False)
        bb.accepted.connect(self._aceptar)
        bb.rejected.connect(self.reject)
        raiz.addWidget(bb)

        # Cualquier cambio rehace el aviso de resolucion y la estimacion.
        for c in self.findChildren(QtWidgets.QCheckBox):
            c.toggled.connect(self._actualizar)
        for s in self.findChildren(QtWidgets.QSpinBox):
            s.valueChanged.connect(self._actualizar)
        for c in self.findChildren(QtWidgets.QComboBox):
            c.currentIndexChanged.connect(self._actualizar)
        self._actualizar()

    def _tiempo(self, rejilla, fila, clave):
        lab = QtWidgets.QLabel()
        lab.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        lab.setMinimumWidth(80)
        lab.setStyleSheet("color:#52514e;")
        rejilla.addWidget(lab, fila, 3)
        self._etq[clave] = lab

    # -- estimacion --

    def familias(self):
        modo = self.cmb_fam.currentData()
        return list(FAMILIAS) if modo == "ambas" else [modo]

    def _hay_ajuste(self):
        return all((self.v._res.get("ajuste" if f == "spinodoide"
                                    else "ajuste_dual") or {})
                   .get("parametros") for f in self.familias())

    def plan(self):
        """Lo que se va a correr, en los terminos de `spinpy.tiempos`."""
        from spinpy.fit import resolucion_comparacion
        v = self.v
        fams = self.familias()
        n_fit = resolucion_comparacion(v.VOI.shape)
        mascaras = [v._de(f, "BW") for f in fams]
        if (self.chk_ajuste.isChecked()
                or any(m is None for m in mascaras)):
            n_cand = n_fit
        else:
            n_cand = max(int(m.shape[0]) for m in mascaras)
        extra, poro, ef = (self.chk_extra.isChecked(),
                           self.chk_poro.isChecked(), self.chk_ef.isChecked())
        mv = v.m_voi or {}
        # Tras ajustar, el VOI solo trae la morfometria basica: con cualquier
        # extra se vuelve a medir (la misma regla que `Visor.medir`).
        medir_voi = (self.chk_ajuste.isChecked() or not mv
                     or (extra and "ConnD" not in mv)
                     or (poro and "PoDm" not in mv)
                     or (ef and "EF" not in mv))
        sim_voi = v.cmb_sim_estructura.currentIndex() == 0
        estilos = [k for k, c in self.chk_estilos.items() if c.isChecked()]
        vistas = [k for k, c in self.chk_vistas.items() if c.isChecked()]
        o = self.F.opciones_3d({
            "estilos": estilos, "vistas": vistas,
            "estilo_principal": self.cmb_estilo_ppal.currentData(),
            "vista_principal": self.cmb_vista_ppal.currentData(),
            "suavizar": self.chk_suave.isChecked(),
            "distribuciones": self.chk_dist.isChecked(),
            "von_mises": self.chk_vm.isChecked(),
            "metodo": self.chk_metodo.isChecked()})
        return {"familias": fams, "vox_voi": int(v.VOI.size),
                "n_fit": n_fit, "ajuste": self.chk_ajuste.isChecked(),
                "morfometria": self.chk_morfo.isChecked(), "extra": extra,
                "poro": poro, "ef": ef, "n_cand": n_cand,
                "medir_voi": medir_voi,
                "elastico": self.chk_ela.isChecked(),
                "n_homog": self.spin_homog.value(),
                "ensayo": self.chk_fe.isChecked(),
                "n_fe": self.spin_fe.value(),
                "n_ejes": 3 if self.cmb_eje.currentIndex() == 3 else 1,
                "comparado": self.chk_comp.isChecked(),
                "dispersion": self.chk_disp.isChecked(),
                "K": self.spin_k.value(),
                "n_resm": int(v.sl["resm"].valor()),
                "convergencia": self.chk_conv.isChecked(),
                "perdida": self.chk_perdida.isChecked(),
                "fallo": self.chk_fallo.isChecked(),
                "pasos": int(v.spin_sim_pasos.value()),
                "protocolo": PROTOCOLOS_SIM[
                    v.cmb_sim_protocolo.currentIndex()],
                "n_sim": (int(v.VOI.shape[0]) if sim_voi else n_cand),
                "n_estilos": len(o["estilos"]),
                "n_vistas": len(o["vistas"]),
                "suavizar": o["suavizar"],
                "distribuciones": o["distribuciones"],
                "von_mises": o["von_mises"] and self.chk_comp.isChecked(),
                "metodo": o["metodo"],
                "pila": v.VOI_contexto is not None,
                "figuras": o}

    def _actualizar(self, *_a):
        # `setChecked` dentro vuelve a emitir y entraria aqui otra vez.
        if self._ocupado_act:
            return
        self._ocupado_act = True
        try:
            self._actualizar_una()
        finally:
            self._ocupado_act = False

    def _actualizar_una(self):
        from spinpy import tiempos
        v = self.v
        fams = self.familias()
        ya = self._hay_ajuste()
        self.chk_ajuste.setEnabled(ya)
        if not ya:
            self.chk_ajuste.setChecked(True)
        self.chk_ajuste.setToolTip(
            _("El desempate mecanico queda fuera: su error suma dos terminos "
              "y no es comparable.")
            + ("\n" + _("Desmarcala para usar el ajuste que ya hay en la "
                        "sesion.") if ya else ""))
        self.intro.setText(_(
            "Corre sobre el VOI <b>{voi}</b> y la familia <b>{fam}</b> todo "
            "el recorrido marcado, con los mismos calculos que los botones "
            "del panel, y escribe el informe al final. Los dialogos de "
            "resultados no se abren: todo queda en la sesion y en el "
            "informe. Puede tardar <b>horas</b>; la ventana sigue respondiendo "
            "y se puede detener tras la etapa en curso.").format(
                voi=html.escape(str(v.VOI_nombre)),
                fam=", ".join(v._etq_fam(f) for f in fams)))
        # Sin el analisis comparado no hay campos de von Mises que pintar.
        self.chk_vm.setEnabled(self.chk_comp.isChecked())
        self._avisar_resolucion()
        plan = self.plan()
        factores, n_ejec = factores_tiempo()
        self.estimacion = tiempos.estimar(plan, factores)
        activas = {"ajuste": self.chk_ajuste, "morfometria": self.chk_morfo,
                   "elastico": self.chk_ela, "ensayo": self.chk_fe,
                   "comparado": self.chk_comp, "dispersion": self.chk_disp,
                   "convergencia": self.chk_conv,
                   "perdida": self.chk_perdida, "fallo": self.chk_fallo}
        for clave, lab in self._etq.items():
            s = sum((self.estimacion.get(clave) or {}).values())
            marcada = clave == "informe" or activas[clave].isChecked()
            if not marcada:
                lab.setText("—")
            elif clave == "convergencia" and s == 0:
                lab.setText(_("no se lanzaria"))
                lab.setToolTip(_("Con esta resolucion del ensayo no salen "
                                 "tres mallas por encima de 12³."))
            else:
                lab.setText(tiempos.texto(s))
        t = tiempos.total(self.estimacion)
        fin = QtCore.QDateTime.currentDateTime().addSecs(int(t))
        hoy = fin.date() == QtCore.QDate.currentDate()
        self.lab_total.setText(
            _("<b>Tiempo estimado: {t}</b> · terminaria hacia las {h}").format(
                t=tiempos.texto(t),
                h=fin.toString("HH:mm") if hoy
                else fin.toString("dd/MM HH:mm")))
        self.lab_nota_tiempo.setText(
            _("Modelo medido en el equipo de desarrollo y corregido con {n} "
              "ejecucion(es) en este equipo. Suele acertar dentro de un "
              "±20 %; en un equipo nuevo, hasta un factor 2 hasta que corra "
              "una vez.").format(n=n_ejec) if n_ejec else
            _("Modelo medido en el equipo de desarrollo; aun no hay "
              "ejecuciones en este equipo para corregirlo, asi que puede "
              "errar hasta un factor 2. Se ajusta solo al terminar cada "
              "etapa."))

    def _avisar_resolucion(self):
        n_min = informe_pub.N_MECANICA_MIN
        t = []
        if self.chk_ela.isChecked() and self.spin_homog.value() < n_min:
            t.append(_("El tensor elastico a {n}³ queda por debajo de {m}³: "
                       "el informe lo dara con reservas (la interfaz lo "
                       "limita a {max}³).").format(
                n=self.spin_homog.value(), m=n_min,
                max=self.v.spin_homog.maximum()))
        if ((self.chk_fe.isChecked() or self.chk_comp.isChecked())
                and self.spin_fe.value() < n_min):
            t.append(_("El ensayo a {n}³ queda por debajo de {m}³: el informe "
                       "dara la mecanica con reservas.").format(
                n=self.spin_fe.value(), m=n_min))
        self.lab_res.setText("\n".join(t))
        self.lab_res.setVisible(bool(t))

    def _elegir_carpeta(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, _("Elige la carpeta del informe para publicacion"),
            self.ed_carpeta.text() or str(DATOS))
        if d:
            self.ed_carpeta.setText(d)
            self.b_ok.setEnabled(True)

    def _aceptar(self):
        v = self.v
        # La familia primero, sin regenerar vistas: la cadena empieza en
        # cuanto se cierra esta ventana, y un hilo de vista en marcha haria
        # que las primeras etapas se saltaran por encontrar el hilo ocupado.
        i = MODOS_FAMILIA.index(self.cmb_fam.currentData())
        if v.cmb_familia.currentIndex() != i:
            previo = v._aplicando
            v._aplicando = True
            try:
                v.cmb_familia.setCurrentIndex(i)
            finally:
                v._aplicando = previo
        # De vuelta al panel: son los controles que leen las etapas.
        v.spin_homog.setValue(self.spin_homog.value())
        v.spin_res_fe.setValue(self.spin_fe.value())
        v.cmb_eje_fe.setCurrentIndex(self.cmb_eje.currentIndex())
        v.cmb_apoyo.setCurrentIndex(self.cmb_apoyo.currentIndex())
        v.spin_k.setValue(self.spin_k.value())
        if self.chk_morfo.isChecked():
            v.chk_extra.setChecked(self.chk_extra.isChecked())
            v.chk_poro.setChecked(self.chk_poro.isChecked())
            v.chk_ef.setChecked(self.chk_ef.isChecked())
        guardar_opciones_figuras(self.plan()["figuras"])
        self.accept()

    def opciones(self):
        return {"carpeta": self.ed_carpeta.text(),
                "familias": self.familias(),
                "ajuste": self.chk_ajuste.isChecked(),
                "morfometria": self.chk_morfo.isChecked(),
                "extra": self.chk_extra.isChecked(),
                "poro": self.chk_poro.isChecked(),
                "ef": self.chk_ef.isChecked(),
                "elastico": self.chk_ela.isChecked(),
                "res_homog": self.spin_homog.value(),
                "ensayo": self.chk_fe.isChecked(),
                "res_fe": self.spin_fe.value(),
                "eje_fe": self.cmb_eje.currentText(),
                "apoyo": apoyo_de_combo(self.cmb_apoyo),
                "comparado": self.chk_comp.isChecked(),
                "dispersion": self.chk_disp.isChecked(),
                "k": self.spin_k.value(),
                "convergencia": self.chk_conv.isChecked(),
                "perdida": self.chk_perdida.isChecked(),
                "fallo": self.chk_fallo.isChecked(),
                "figuras": self.plan()["figuras"],
                "estimacion": {k: {("" if f is None else f): s
                                   for f, s in d.items()}
                               for k, d in self.estimacion.items()}}


# ---------------------------------------------------------------------------
# Ventana principal
# ---------------------------------------------------------------------------

class Visor(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("spinpy — visor de spinodoides y VOIs")
        self.resize(1500, 900)
        _ico = _ruta_icono()
        if _ico is not None:
            self.setWindowIcon(QtGui.QIcon(str(_ico)))
        _icono_de_clase(self)

        self.VOI = None
        self.VOI_spacing = None
        self.VOI_nombre = "—"
        self.m_voi = None
        self.m_spin = None
        self._m_spin_vigente = False   # ¿`m_spin` corresponde a `BW_vista`?
        self.BW_vista = None
        self.malla_voi = None
        self.malla_spin = None
        self.hilo = None
        self.R_forzada = None      # R impuesta por el ajuste (F4)
        self._aplicando = False    # evita regenerar mientras se colocan valores
        self._obs = []             # observadores de sincronizacion de camaras
        self._esp = {}             # campos de espesor local cacheados
        self._fe = {}              # campos de deformacion efectiva
        self._vm = {}              # campos de tension de von Mises
        # FEBio: evento para detener la corrida en curso y mapas de la ultima
        # (malla + von Mises por elemento; en memoria, no en `_res`).
        self._febio_evento = None
        self._febio_mapas = []
        self._desp = {}            # campos de deformacion total (mm)
        self._clim = None          # escala de color COMUN a los dos paneles
        self._esp_luego = None     # que hacer cuando termine el espesor
        self.VOI_ruta = None       # para poder guardarla en la sesion
        # Figura 0 del informe: de donde salio el VOI. `_pila_param` es la
        # carga de la pila en curso; `VOI_contexto` (rebanada, pila reducida,
        # esquinas del cubo) se toma al recortar, porque despues la pila se
        # suelta; `VOI_recorte` es lo que el JSON guarda para rehacerlo.
        self._pila_param = None
        self.VOI_contexto = None
        self.VOI_recorte = None
        # Ultimo resultado de cada calculo, para poder exportarlos. Se guardan
        # segun se producen y no se recalculan al exportar: exportar tiene que
        # volcar EXACTAMENTE los numeros que el usuario vio, y como el
        # generador es estocastico, recalcular daria otros.
        self._res = {}
        # Estado del informe automatico mientras recorre sus etapas; None
        # fuera de el. Mientras no es None, los dialogos de resultados no se
        # abren y los errores se anotan en vez de interrumpir.
        self._auto = None
        # Campos de von Mises del ANALISIS COMPARADO por familia, para la
        # figura 8 del informe: {familia: {"voi": (campo, sp), "spin": ...}}.
        # Aparte de `_vm`, que tambien escribe el ensayo propio: la figura
        # declara el protocolo de 100 N y no puede pintar el de otro ensayo.
        self._vm_comparado = {}

        # FAMILIAS. Todo el codigo de los botones se escribio para UN
        # candidato: `BW_vista`, `m_spin`, `malla_spin`, `vis_spin`, `R_forzada`
        # y la entrada "spin" de los campos cacheados. En vez de reescribirlo,
        # esos atributos pasan a ser los de la familia ACTIVA, y `_activar`
        # guarda y carga cada familia en `_cand`. Asi el mismo boton sirve para
        # el spinodoide, para el dual-lattice o para los dos, y la ruta de
        # calculo de cada familia es LITERALMENTE la misma (correccion C4).
        self._modo_fam = "spinodoide"
        self._fam = "spinodoide"
        self._cand = {f: self._cand_vacio() for f in FAMILIAS}

        cen = QtWidgets.QWidget()
        self.setCentralWidget(cen)
        lay = QtWidgets.QHBoxLayout(cen)

        # El panel va dentro de una QScrollArea: si la pantalla es baja, se
        # desplaza en vez de comprimir los controles hasta hacerlos ilegibles.
        panel = self._panel_control()
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(panel)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        scroll.setFixedWidth(panel.width() + 22)   # + hueco de la barra
        lay.addWidget(scroll, 0)

        div = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        vistas = QtWidgets.QSplitter(QtCore.Qt.Horizontal)

        self.vis_voi = QtInteractor(self)
        # Un panel por familia. Solo se ven los de las familias activas: con
        # «Ambas» son tres paneles, VOI y los dos candidatos.
        self._vis = {f: QtInteractor(self) for f in FAMILIAS}
        self.vis_spin = self._vis["spinodoide"]
        self._caja = {}
        for clave, v, t in (("voi", self.vis_voi, _("VOI real")),
                            ("spinodoide", self._vis["spinodoide"],
                             _("Spinodoide")),
                            ("dual-lattice", self._vis["dual-lattice"],
                             _("Dual-lattice"))):
            self._estilo_render(v)
            v.add_text(t, font_size=11, color="black")
            # Triada de orientacion: sin ella no se sabe que eje es cual, y la
            # direccion principal del MIL solo tiene sentido respecto a los ejes.
            try:
                v.add_axes(line_width=2, labels_off=False)
            except Exception:
                pass
            caja = QtWidgets.QWidget()
            cl = QtWidgets.QVBoxLayout(caja)
            cl.setContentsMargins(0, 0, 0, 0)
            cl.addWidget(v.interactor)
            vistas.addWidget(caja)
            self._caja[clave] = caja
        self._caja["dual-lattice"].setVisible(False)
        vistas.setSizes([600, 600, 600])

        div.addWidget(vistas)
        div.addWidget(self._tabla())
        div.setSizes([620, 260])
        lay.addWidget(div, 1)

        self._menu()
        self._barra_familia()
        self._barra_publicacion()
        self.statusBar().showMessage(_("Listo. Carga un VOI o pulsa Generar."))
        self._sincronizar(True)
        self._sync_formatos()      # el .stl arranca deshabilitado con la hex

        # La captura va AQUI: con la interfaz entera construida y todavia en
        # espanol. Guarda el texto original de cada widget para poder ir y
        # volver entre idiomas sin degradar nada. Si se hiciera despues de
        # traducir, el ingles quedaria grabado como original.
        self._textos = capturar(self)
        cod = QtCore.QSettings("spinpy", "visor").value("idioma", "es")
        if cod in IDIOMAS and cod != "es":
            self.cambiar_idioma(cod, guardar=False)

        self._generar_vista()

    def abrir_bienvenida(self):
        """La portada, ya con la ventana abierta y en el idioma en curso."""
        DialogoBienvenida(self).exec_()

    def _menu(self):
        m = self.menuBar().addMenu("&Sesion")
        a = m.addAction("Guardar sesion…")
        a.setShortcut("Ctrl+S")
        a.triggered.connect(self.guardar_sesion)
        a = m.addAction("Cargar sesion…")
        a.setShortcut("Ctrl+O")
        a.triggered.connect(self.cargar_sesion)
        m.addSeparator()
        a = m.addAction("Exportar resultados (CSV + JSON)…")
        a.setShortcut("Ctrl+E")
        a.triggered.connect(self.exportar_resultados)
        m.addSeparator()
        a = m.addAction("Cargar VOI…")
        a.triggered.connect(self.cargar_voi)
        a = m.addAction("Cargar pila TIFF (micro-CT)…")
        a.triggered.connect(self.cargar_pila_tiff)
        m.addSeparator()
        # La portada se puede desactivar con su propia casilla, asi que tiene
        # que haber una forma de volver a ella: si no, quien la apague pierde
        # los creditos y el enlace al repositorio para siempre.
        a = m.addAction("Acerca de spinpy…")
        a.triggered.connect(self.abrir_bienvenida)

        # --- Validacion ---
        # Un menu propio y no un boton del panel: no es una operacion sobre
        # los datos cargados, es documentacion ejecutable sobre el metodo. El
        # panel esta ordenado por etapas del trabajo y esto no es una etapa.
        mv = self.menuBar().addMenu("&Validacion")
        # Una entrada por articulo, y no una sola que abra "la validacion":
        # cada replica comprueba una cosa distinta y el menu es donde se ve
        # que son tres cosas. Todas abren la misma ventana, en su pestana.
        for clave, texto, ayuda in (
                ("kumar2020", "Replica de Kumar et al. (2020)…",
                 "El metodo en si: las cuatro clases de anisotropia y el "
                 "muestreo por rechazo de su ecuacion (2)."),
                ("zheng2021", "Replica de Zheng et al. (2021)…",
                 "Las cotas: cada superficie elastica dentro de Voigt y de "
                 "Hashin-Shtrikman, la ortotropia del tensor, y por que el "
                 "articulo impone rho >= 0.3."),
                ("guo2024", "Replica de Guo et al. (2024)…",
                 "Las curvaturas: el perfil (k1, k2) de un spinodoide de "
                 "parametros publicados y el de una superficie nodal "
                 "periodica, cuya respuesta es cerrada.")):
            a = mv.addAction(texto)
            a.setToolTip(ayuda)
            a.triggered.connect(lambda _c, k=clave: self.abrir_validacion(k))

        # --- Idioma ---
        # Un grupo excluyente y marcable, no dos acciones sueltas: asi se ve
        # cual esta activo sin tener que cambiarlo para averiguarlo.
        mi = self.menuBar().addMenu("&Idioma")
        self._grupo_idioma = QtWidgets.QActionGroup(self)
        self._grupo_idioma.setExclusive(True)
        actual = QtCore.QSettings("spinpy", "visor").value("idioma", "es")
        for cod, nombre in IDIOMAS.items():
            a = mi.addAction(nombre)
            a.setCheckable(True)
            a.setChecked(cod == actual)
            a.setData(cod)
            self._grupo_idioma.addAction(a)
            a.triggered.connect(lambda _c, k=cod: self.cambiar_idioma(k))

    def abrir_validacion(self, clave=None):
        """Abre la ventana de replicas, en la pestana `clave` si se pide.

        El dialogo vive en `dialogo_validacion.py` y cada replica en su
        `Test/replicar_*.py`, que se pueden correr sin interfaz. Aqui solo se
        abre.
        """
        try:
            from dialogo_validacion import DialogoValidacion
        except Exception:
            QtWidgets.QMessageBox.critical(
                self, _("No se pudo abrir la validación"),
                traceback.format_exc())
            return
        DialogoValidacion(self, clave=clave).exec_()

    def cambiar_idioma(self, cod, guardar=True):
        """Reescribe la interfaz en `cod` SIN reconstruirla ni perder estado.

        Reconstruir la ventana seria mas simple de programar y peor de usar:
        se perderian el VOI cargado, las metricas medidas y los campos de
        deformacion, que son minutos u horas de calculo. En vez de eso se
        vuelven a escribir los textos sobre los mismos widgets.

        Lo que hay que refrescar a mano es lo que NO es un widget: los
        rotulos dibujados dentro de los paneles 3D y las celdas de la tabla,
        que se escriben en cada medida y no en la construccion.
        """
        if fijar_idioma(cod) != cod:
            return
        if guardar:
            QtCore.QSettings("spinpy", "visor").setValue("idioma", cod)
        aplicar(self._textos)
        for d in self.findChildren(Deslizador):
            d.retraducir()
        if self.VOI is None:
            self.lab_voi.setText(_("ninguno"))
        for a in self._grupo_idioma.actions():
            a.setChecked(a.data() == cod)
        # Las cabeceras de la tabla dependen de la familia y no se pueden
        # restaurar desde la captura del arranque, que se hizo con 4 columnas.
        self._reconstruir_tabla()
        if self.m_voi or any(self._de(f, "m") for f in FAMILIAS):
            self._rellenar_tabla()
        self._redibujar()          # los rotulos de los paneles 3D
        self.statusBar().showMessage(
            _("Idioma cambiado. Los resultados ya calculados se conservan."))

    # -- panel izquierdo ---------------------------------------------------

    def _seccion(self, lay, titulo, primera=False):
        """Cabecera de una seccion del panel: raya, titulo y raya.

        El panel es largo y antes todo colgaba del mismo nivel: los grupos
        estaban ahi, pero no habia forma de ver de un vistazo en que etapa
        del trabajo estaba cada control. Las secciones agrupan por ETAPA
        -cargar, generar, mirar, medir, ajustar, ensayar, exportar- y
        dejan sitio evidente donde colgar lo que se anada despues.
        """
        if not primera:
            lay.addSpacing(12)
            raya = QtWidgets.QFrame()
            raya.setFrameShape(QtWidgets.QFrame.HLine)
            raya.setStyleSheet("color:#c0c0c0;")
            lay.addWidget(raya)
        t = QtWidgets.QLabel(titulo)
        t.setStyleSheet("font-size:13px; font-weight:bold; color:#1b4f72;"
                        "padding:6px 0 2px 0;")
        lay.addWidget(t)
        raya = QtWidgets.QFrame()
        raya.setFrameShape(QtWidgets.QFrame.HLine)
        raya.setStyleSheet("color:#1b4f72;")
        lay.addWidget(raya)
        cont = QtWidgets.QVBoxLayout()
        cont.setContentsMargins(0, 4, 0, 0)
        lay.addLayout(cont)
        return cont

    def _libs(self, sec, texto):
        """Pie de seccion con las librerias de Python que sostienen el bloque.

        En gris y a 9 px a proposito. Quien viene a USAR la aplicacion no
        necesita saber que la morfometria sale de scikit-image; quien viene a
        reproducir el calculo, a portarlo o a discutir un numero lo encuentra
        sin abrir el codigo. Sigue siendo texto legible, no texto oculto:
        simplemente no compite por la atencion con los controles.
        """
        l = QtWidgets.QLabel(texto)
        l.setWordWrap(True)
        l.setStyleSheet("color:#9a9a9a; font-size:9px; padding:3px 0 0 0;")
        l.setToolTip("Librerias de Python que hacen el calculo de esta seccion.")
        sec.addWidget(l)

    def _panel_control(self):
        w = QtWidgets.QWidget()
        w.setFixedWidth(330)
        lay = QtWidgets.QVBoxLayout(w)

        sec = self._seccion(lay, "VOI de referencia", primera=True)
        gl = sec
        b = QtWidgets.QPushButton("Cargar VOI (.vtk / .mat)…")
        b.clicked.connect(self.cargar_voi)
        gl.addWidget(b)
        # Segunda puerta de entrada: las imagenes tal como salen del escaner.
        # Antes habia que pasar por MATLAB para llegar a un .vtk.
        b = QtWidgets.QPushButton("Cargar pila TIFF (micro-CT)…")
        b.setToolTip("Carpeta de rebanadas o TIFF multipagina, con su tamano "
                     "de voxel real")
        b.clicked.connect(self.cargar_pila_tiff)
        gl.addWidget(b)
        # Una pila de escaner no es un VOI. Este es el paso que faltaba para
        # cerrar el camino desde el escaner sin salir de la interfaz.
        self.b_recorte = QtWidgets.QPushButton("Recortar VOI cubico…")
        self.b_recorte.setToolTip("Orienta el hueso por PCA y recorta un cubo: "
                                  "tres tramos sugeridos o posicion libre")
        self.b_recorte.clicked.connect(self.recortar_voi)
        self.b_recorte.setEnabled(False)
        gl.addWidget(self.b_recorte)
        self.lab_voi = QtWidgets.QLabel(_("ninguno"))
        self.lab_voi.setProperty("i18n_dinamico", True)
        self.lab_voi.setWordWrap(True)
        self.lab_voi.setStyleSheet("color:#555;")
        gl.addWidget(self.lab_voi)

        self._libs(sec, "Librerias: numpy · scipy.io para los VOI en .mat · "
                        "scikit-image y Pillow para las pilas TIFF")

        sec = self._seccion(lay, "Microestructura")
        g = QtWidgets.QGroupBox("Parametros del campo")
        self.grp_spin = g
        gl = QtWidgets.QVBoxLayout(g)
        self.sl = {}
        for clave, texto, lo, hi, val, dec, suf in [
            ("dens",  "Densidad relativa", 5, 95, 35, 0, " %"),
            ("wave",  "Numero de onda",    4, 30, 15, 1, " pi"),
            ("nw",    "Numero de ondas", 100, 1500, 700, 0, ""),
            ("thx",   "Theta X",           0, 90, 15, 0, " deg"),
            ("thy",   "Theta Y",           0, 90, 15, 0, " deg"),
            ("thz",   "Theta Z",           0, 90, 45, 0, " deg"),
            ("rotx",  "Rotacion X",      -90, 90, 0, 0, " deg"),
            ("roty",  "Rotacion Y",      -90, 90, 0, 0, " deg"),
            ("rotz",  "Rotacion Z",      -90, 90, 0, 0, " deg"),
        ]:
            d = Deslizador(texto, lo, hi, val, dec, suf)
            # Las rotaciones ademas liberan la R impuesta por el ajuste (F4),
            # y SOLO la de esta familia.
            d.cambiado.connect(
                (lambda: self._soltar_R_fam("spinodoide"))
                if clave.startswith("rot") else self._pedir_vista)
            self.sl[clave] = d
            gl.addWidget(d)

        self.lab_R = QtWidgets.QLabel("")
        self.lab_R.setProperty("i18n_dinamico", True)
        self.lab_R.setWordWrap(True)
        self.lab_R.setStyleSheet("color:#b06000; font-size:10px;")
        gl.addWidget(self.lab_R)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Muestreo de ondas:"))
        self.cmb_esq = QtWidgets.QComboBox()
        self.cmb_esq.addItems(["rechazo (GIBBON)", "equitativo (TPMS-Scaffolds)"])
        self.cmb_esq.setToolTip(
            "Con conos de angulos DESIGUALES los dos esquemas producen\n"
            "anisotropias distintas para los mismos thetas. No es un ajuste\n"
            "cosmetico: hay que declarar cual se uso.")
        self.cmb_esq.currentIndexChanged.connect(self._pedir_vista)
        f.addWidget(self.cmb_esq, 1)
        gl.addLayout(f)

        sec.addWidget(g)

        # --- dual-lattice ---
        g = QtWidgets.QGroupBox("Parametros del dual-lattice")
        g.setToolTip(
            "Red dual de una teselacion de Delaunay (Vafaeefar et al. 2022):\n"
            "el centroide de cada tetraedro unido a los de sus cuatro caras,\n"
            "engrosado hasta la densidad pedida. Cuatro barras por nudo.")
        self.grp_dual = g
        gl = QtWidgets.QVBoxLayout(g)
        for clave, texto, lo, hi, val, dec, suf in [
            ("d_dens",   "Densidad relativa", 5, 95, 30, 0, " %"),
            ("d_celdas", "Celdas por lado",   2, 20, 5, 2, ""),
            ("d_ex",     "Estiramiento X",    0.5, 4, 1, 2, ""),
            ("d_ey",     "Estiramiento Y",    0.5, 4, 1, 2, ""),
            ("d_ez",     "Estiramiento Z",    0.5, 4, 1, 2, ""),
            ("d_irr",    "Irregularidad",     0.05, 1, 0.5, 2, ""),
            ("d_rotx",   "Rotacion X",      -90, 90, 0, 0, " deg"),
            ("d_roty",   "Rotacion Y",      -90, 90, 0, 0, " deg"),
            ("d_rotz",   "Rotacion Z",      -90, 90, 0, 0, " deg"),
        ]:
            d = Deslizador(texto, lo, hi, val, dec, suf)
            d.cambiado.connect(
                (lambda: self._soltar_R_fam("dual-lattice"))
                if clave.startswith("d_rot") else self._pedir_vista)
            self.sl[clave] = d
            gl.addWidget(d)
        self.lab_R_dual = QtWidgets.QLabel("")
        self.lab_R_dual.setProperty("i18n_dinamico", True)
        self.lab_R_dual.setWordWrap(True)
        self.lab_R_dual.setStyleSheet("color:#b06000; font-size:10px;")
        gl.addWidget(self.lab_R_dual)
        self._lab_R = {"spinodoide": self.lab_R,
                       "dual-lattice": self.lab_R_dual}
        nota = QtWidgets.QLabel(
            "El estiramiento alarga la MALLA antes de engrosarla: el eje "
            "estirado es el eje rigido. Con 2.5 en un eje el DA ronda 1.5, el "
            "del VOI proximal equino. El grosor de las barras no se elige: "
            "sale de la densidad.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        g.setVisible(False)
        sec.addWidget(g)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Semilla:"))
        self.spin_semilla = QtWidgets.QSpinBox()
        self.spin_semilla.setRange(0, 10**8)
        self.spin_semilla.setValue(20260720)
        self.spin_semilla.valueChanged.connect(self._pedir_vista)
        f.addWidget(self.spin_semilla, 1)
        b = QtWidgets.QPushButton("Otra")
        b.setToolTip("El generador es estocastico: la misma parametrizacion da\n"
                     "realizaciones distintas. Cambiar la semilla muestra esa\n"
                     "dispersion, que es el suelo de ruido de cualquier ajuste.")
        b.clicked.connect(lambda: self.spin_semilla.setValue(
            int(np.random.default_rng().integers(0, 10**8))))
        f.addWidget(b)
        sec.addLayout(f)

        g = QtWidgets.QGroupBox("Resolucion")
        gl = QtWidgets.QVBoxLayout(g)
        self.sl["resv"] = Deslizador("Vista (interactiva)", 16, 256, 32, 0, " vox")
        self.sl["resv"].cambiado.connect(self._pedir_vista)
        self.sl["resv"].setToolTip(AVISO_RES)
        gl.addWidget(self.sl["resv"])
        self.sl["resm"] = Deslizador("Medida", 32, 256, 64, 0, " vox")
        self.sl["resm"].setToolTip(AVISO_RES)
        gl.addWidget(self.sl["resm"])
        nota = QtWidgets.QLabel(
            "La vista se recalcula al mover un deslizador; la medida solo al "
            "pulsar Medir. Misma semilla, asi que la vista previa es una "
            "version gruesa de lo que se medira.\n"
            "El coste crece con el CUBO del lado: medido aqui, 64³ tarda 5 s, "
            "128³ 31 s y 256³ 263 s con 1.5 GB. Por encima de ~96 vox conviene "
            "quitar la actualizacion automatica.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        self.chk_auto = QtWidgets.QCheckBox("Actualizar vista automaticamente")
        self.chk_auto.setChecked(True)
        gl.addWidget(self.chk_auto)
        sec.addWidget(g)

        self.btn_gen = QtWidgets.QPushButton("Generar vista")
        self.btn_gen.clicked.connect(self._generar_vista)
        sec.addWidget(self.btn_gen)

        self._libs(sec, "Librerias: numpy · scipy.special (el erf⁻¹ del umbral) · pyvista/VTK para la isosuperficie de la vista")

        sec = self._seccion(lay, "Visualizacion")
        g = QtWidgets.QGroupBox("Camara y color")
        gl = QtWidgets.QVBoxLayout(g)
        self.chk_sync = QtWidgets.QCheckBox("Enlazar las camaras de los dos paneles")
        self.chk_sync.setChecked(True)
        self.chk_sync.setToolTip(
            "Con encuadres distintos, comparar grosores a ojo entre los dos\n"
            "paneles no significa nada: una estructura puede parecer mas\n"
            "gruesa solo por estar mas cerca de la camara.")
        self.chk_sync.toggled.connect(self._sincronizar)
        gl.addWidget(self.chk_sync)
        self.chk_ejes = QtWidgets.QCheckBox("Mostrar ejes y cotas en mm")
        self.chk_ejes.setChecked(True)
        self.chk_ejes.toggled.connect(lambda _: self._redibujar())
        gl.addWidget(self.chk_ejes)

        self.chk_fab = QtWidgets.QCheckBox("Mostrar fabrica MIL (eje y elipsoide)")
        self.chk_fab.setChecked(True)
        self.chk_fab.setToolTip(
            "Dibuja el elipsoide de fabrica del tensor MIL y, en magenta, su\n"
            "eje mayor: la direccion en la que la estructura es mas continua\n"
            "y, por tanto, mas rigida.\n\n"
            "Es el mismo dato que da el DA de la tabla, pero orientado. Dos\n"
            "estructuras pueden tener el MISMO DA y estar orientadas a 90\n"
            "grados una de otra: eso la tabla no lo distingue y esto si.\n\n"
            "Necesita haber medido. El elipsoide del spinodoide corresponde a\n"
            "la ultima medida, no a los deslizadores: si cambias parametros,\n"
            "desaparece hasta que vuelvas a medir.")
        self.chk_fab.toggled.connect(lambda _: self._redibujar())
        gl.addWidget(self.chk_fab)

        self.chk_ssao = QtWidgets.QCheckBox("Sombras de contacto (SSAO)")
        self.chk_ssao.setChecked(True)
        self.chk_ssao.setToolTip(
            "Oclusion ambiental: oscurece los huecos donde una trabecula pasa\n"
            "por detras de otra. Sin ella la marana se lee plana y no se sabe\n"
            "que esta delante.\n\n"
            "Es un efecto de RENDER. No toca ninguna medida. Desmarcala si el\n"
            "equipo va lento o si el driver de video la dibuja mal.")
        self.chk_ssao.toggled.connect(lambda _: self._redibujar())
        gl.addWidget(self.chk_ssao)

        self.chk_fondo = QtWidgets.QCheckBox("Fondo blanco (para figuras)")
        self.chk_fondo.setToolTip(
            "En pantalla el degradado da a la silueta algo contra lo que\n"
            "recortarse y el volumen se lee mejor. Para una figura de tesis o\n"
            "una captura de un informe conviene el blanco puro.")
        self.chk_fondo.toggled.connect(lambda _: self._aplicar_fondo())
        gl.addWidget(self.chk_fondo)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Colorear por:"))
        self.cmb_color = QtWidgets.QComboBox()
        self.cmb_color.addItems(["Material (liso)", "Espesor local",
                                 "Deformacion efectiva",
                                 "Tension de von Mises",
                                 "Deformacion total (mm)"])
        self.cmb_color.setToolTip(
            "El espesor local muestra DONDE la estructura es gruesa o fina.\n"
            "Los dos paneles comparten escala de color, sin lo cual la\n"
            "comparacion no significaria nada.\n\n"
            "AVISO: a 3-4 voxeles de grosor el metodo subestima un 30%, asi\n"
            "que los valores absolutos NO son citables como Tb.Th. Para eso\n"
            "esta el Tb.Th de Parfitt de la tabla.\n\n"
            "La deformacion efectiva y von Mises salen del mismo ensayo de\n"
            "compresion, pero NO son el mismo mapa reescalado: la primera\n"
            "viene de la energia total y la segunda solo de la parte\n"
            "desviadora. Se separan en los nudos, donde el estado es triaxial.")
        self.cmb_color.currentIndexChanged.connect(self._cambiar_color)
        f.addWidget(self.cmb_color, 1)
        gl.addLayout(f)

        self.btn_hist = QtWidgets.QPushButton("Histograma de espesor…")
        self.btn_hist.setToolTip(
            "Distribucion completa del espesor local, VOI y spinodoide\n"
            "superpuestos. Dos estructuras con el mismo Tb.Th pueden tener\n"
            "distribuciones muy distintas, y eso el escalar no lo distingue.")
        self.btn_hist.clicked.connect(self.histograma_espesor)
        gl.addWidget(self.btn_hist)
        sec.addWidget(g)

        # --- recorte ---
        # Una estructura con BV/TV alto se ve desde fuera como un bloque macizo:
        # la superficie exterior tapa el interior. Recortar es la unica forma de
        # inspeccionar la trabecula. No afecta a NINGUNA medida, solo al render.
        g = QtWidgets.QGroupBox("Recorte (solo visual)")
        gl = QtWidgets.QVBoxLayout(g)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Eje:"))
        self.cmb_eje = QtWidgets.QComboBox()
        self.cmb_eje.addItems(["X", "Y", "Z"])
        self.cmb_eje.currentIndexChanged.connect(self._redibujar)
        f.addWidget(self.cmb_eje, 1)
        gl.addLayout(f)
        self.sl["clip"] = Deslizador("Fraccion recortada", 0, 90, 0, 0, " %")
        self.sl["clip"].cambiado.connect(self._redibujar)
        gl.addWidget(self.sl["clip"])
        sec.addWidget(g)

        self._libs(sec, "Librerias: pyvista y pyvistaqt (los dos paneles) · matplotlib (histograma) · scipy.ndimage (espesor local por transformada de distancia)")

        sec = self._seccion(lay, "Morfometria")
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Superficie:"))
        self.cmb_modo_med = QtWidgets.QComboBox()
        self.cmb_modo_med.addItems(["voxeles (marching cubes)",
                                    "malla suavizada + PyMeshFix"])
        self.cmb_modo_med.setToolTip(
            "voxeles: marching cubes crudo sobre la mascara. Es el modo de la\n"
            "app de MATLAB y el que usa el ajuste. Sobreestima el area un\n"
            "~8.5 % (medido sobre una esfera), igual en el VOI y en el\n"
            "candidato, asi que el sesgo se cancela al comparar.\n\n"
            "malla: superficie cerrada, suavizada (Taubin) y reparada con\n"
            "PyMeshFix, recortada medio voxel por dentro (caja de centros de\n"
            "voxel, la misma de la superficie abierta) y sin las seis tapas.\n"
            "Baja ese sesgo a <1 %. El BV/TV sale ~3-5 % menor que el de\n"
            "voxeles: es la diferencia entre la superficie a nivel 0.5 y el\n"
            "conteo de cubos enteros, no una perdida. DA y fraccion portante\n"
            "siguen saliendo de los voxeles.\n\n"
            "VOI y candidato se miden SIEMPRE con el mismo modo. Cambiarlo\n"
            "vuelve a medir el VOI. Cada tabla debe declarar cual se uso.")
        f.addWidget(self.cmb_modo_med, 1)
        sec.addLayout(f)
        self.btn_med = QtWidgets.QPushButton("Medir (morfometria completa)")
        self.btn_med.setStyleSheet("font-weight:bold; padding:6px;")
        self.btn_med.clicked.connect(self.medir)
        sec.addWidget(self.btn_med)

        f = QtWidgets.QHBoxLayout()
        self.btn_disp = QtWidgets.QPushButton("Dispersion…")
        self.btn_disp.setToolTip(
            "Repite la MISMA parametrizacion con K semillas distintas y da\n"
            "media +- sd de cada metrica. Es el suelo de ruido por debajo del\n"
            "cual ninguna diferencia significa nada.\n\n"
            "Con un VOI cargado anade la z: cuantas sd separan al VOI de la\n"
            "media de las realizaciones. |z| < 2 es indistinguible del ruido.")
        self.btn_disp.clicked.connect(self.dispersion)
        f.addWidget(self.btn_disp, 1)
        f.addWidget(QtWidgets.QLabel("K:"))
        self.spin_k = QtWidgets.QSpinBox()
        self.spin_k.setRange(3, 40)
        self.spin_k.setValue(8)
        f.addWidget(self.spin_k)
        sec.addLayout(f)

        self.chk_extra = QtWidgets.QCheckBox("Anadir Conn.D y SMI a la medida")
        self.chk_extra.setToolTip(
            "Conn.D (Odgaard & Gundersen) mide las conexiones REDUNDANTES de\n"
            "la red: separa 'hueso mas fino' de 'hueso roto', que es algo que\n"
            "BV/TV y Tb.Th no distinguen.\n\n"
            "SMI (Hildebrand & Ruegsegger) mide si la estructura tiende a\n"
            "placa o a barra. Leerlo con la salvedad de Salmon et al. 2015:\n"
            "a BV/TV alto la concavidad lo hace negativo y pierde sentido.\n\n"
            "Estan fuera del bucle de ajuste a proposito: cuestan una\n"
            "marching cubes de mas cada una y no entran en el error.")
        sec.addWidget(self.chk_extra)

        self.chk_poro = QtWidgets.QCheckBox("Anadir el tamano de poro (Po.Dm)")
        self.chk_poro.setToolTip(
            "Po.Dm es el espesor local de la fase PORO, medido con esferas\n"
            "inscritas (Hildebrand & Ruegsegger sobre el complemento). Es la\n"
            "separacion MEDIDA, frente a Tb.Sp, que es la separacion que\n"
            "DEDUCE el modelo de placas de Parfitt. En una pila de placas los\n"
            "dos coinciden; en hueso trabecular la diferencia entre ambos dice\n"
            "cuanto se aparta la estructura de ese modelo.\n\n"
            "Casilla aparte porque es la medida mas cara de la tabla: una\n"
            "transformada de distancia por radio sobre el poro, que a BV/TV\n"
            "0.3 es el 70 % del volumen. Medido: 0.7 s a 64^3, 3.2 s a 96^3 y\n"
            "12 s a 128^3.\n\n"
            "Aviso de borde: los poros que tocan la cara del cubo estan\n"
            "cortados y su esfera inscrita no esta acotada por ese lado. No se\n"
            "corrige -no hay forma honesta de saber cuanto seguia el poro\n"
            "fuera- y el informe lleva la fraccion afectada.")
        sec.addWidget(self.chk_poro)

        self.chk_ef = QtWidgets.QCheckBox("Anadir el Ellipsoid Factor (EF)")
        self.chk_ef.setToolTip(
            "Ellipsoid Factor (Doube 2015): en cada voxel, el mayor elipsoide\n"
            "que lo contiene y cabe en el hueso. EF = a/b - b/c vale -1 en una\n"
            "placa, 0 en una esfera y +1 en una barra.\n\n"
            "Mide lo mismo que el SMI —placa o barra— sin su defecto: el SMI\n"
            "se confunde con la CONCAVIDAD (Salmon et al. 2015), y el hueso\n"
            "tiene mucha mas superficie concava que un spinodoide.\n\n"
            "AVISO DE COSTE: es la medida mas cara del programa, del orden de\n"
            "minutos por estructura a 64^3 y de media hora a 97^3. No se\n"
            "recorta para ir mas rapido: recortarlo sesga hacia barra.")
        sec.addWidget(self.chk_ef)

        self._libs(sec, "Librerias: numpy · scipy.ndimage (componentes conexas y distancias) · scikit-image (marching cubes y numero de Euler) · pyvista y pymeshfix en el modo malla")

        sec = self._seccion(lay, "Ajuste al VOI")
        self.btn_fit = QtWidgets.QPushButton("Ajustar al VOI")
        self.btn_fit.setStyleSheet("font-weight:bold; padding:6px;")
        self.btn_fit.setToolTip(
            "Busqueda escalonada de fitSpinodoidToVOI: densidad y numero de\n"
            "onda, despues angulos conicos, despues refinado. Son 53\n"
            "evaluaciones y tarda varios minutos a la resolucion del VOI.")
        self.btn_fit.clicked.connect(self.ajustar)
        sec.addWidget(self.btn_fit)

        self.btn_todos = QtWidgets.QPushButton("Ajustar con TODOS los metodos…")
        self.btn_todos.setToolTip(
            "Lanza los cuatro metodos disponibles a la vez, con una barra de\n"
            "progreso por metodo, y al terminar los compara lado a lado con\n"
            "sus metricas frente al VOI para que elijas con cual seguir.\n\n"
            "Desde ahi se pueden generar ademas N replicas del elegido,\n"
            "identicas o con variacion en las metricas que marques.")
        self.btn_todos.clicked.connect(self.ajustar_todos)
        sec.addWidget(self.btn_todos)

        f = QtWidgets.QHBoxLayout()
        self.chk_mec = QtWidgets.QCheckBox("Desempate mecanico")
        self.chk_mec.setToolTip(
            "Anade una etapa D: reordena los mejores candidatos incluyendo\n"
            "Ez/Es y Ez/Ex, obtenidos por homogeneizacion.\n\n"
            "NO es una busqueda con el termino mecanico dentro: homogeneizar\n"
            "las ~53 evaluaciones seria de tres a cuatro ordenes de magnitud\n"
            "mas caro. Solo se paga entre candidatos que ya son buenos\n"
            "morfometricamente, asi que si el optimo mecanico esta en una\n"
            "region que la busqueda descarto pronto, esto no lo encuentra.")
        self.chk_mec.toggled.connect(
            lambda v: [w.setEnabled(v) for w in (self.spin_peso_mec,
                                                 self.spin_res_mec)])
        f.addWidget(self.chk_mec, 1)
        sec.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Peso:"))
        self.spin_peso_mec = QtWidgets.QDoubleSpinBox()
        self.spin_peso_mec.setRange(0.1, 5.0)
        self.spin_peso_mec.setSingleStep(0.5)
        self.spin_peso_mec.setValue(1.0)
        self.spin_peso_mec.setEnabled(False)
        f.addWidget(self.spin_peso_mec)
        f.addWidget(QtWidgets.QLabel("Malla:"))
        self.spin_res_mec = QtWidgets.QSpinBox()
        self.spin_res_mec.setRange(8, 24)
        self.spin_res_mec.setValue(16)
        self.spin_res_mec.setSuffix(" vox")
        self.spin_res_mec.setEnabled(False)
        self.spin_res_mec.setToolTip(
            "Ez/Es depende MUCHO de esta malla: sobre el mismo VOI salio\n"
            "0.108 a 12 voxeles y 0.064 a 16. No es citable como rigidez.\n"
            "Sirve porque VOI y candidatos se miden con la MISMA malla y el\n"
            "sesgo se cancela al comparar (mismo argumento que C4).")
        f.addWidget(self.spin_res_mec)
        sec.addLayout(f)

        self._libs(sec, "Librerias: numpy · scipy.ndimage · scipy.sparse solo si se activa el desempate mecanico")

        sec = self._seccion(lay, "Analisis mecanico")
        self.btn_ela = QtWidgets.QPushButton("Tensor elastico")
        self.btn_ela.setToolTip(
            "Homogeneizacion periodica sobre la rejilla de voxeles.\n"
            "El coste crece con el CUBO del lado, asi que la estructura se\n"
            "remuestrea antes al limite indicado.")
        self.btn_ela.clicked.connect(self.calcular_elastico)
        sec.addWidget(self.btn_ela)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Resolucion del tensor:"))
        self.spin_homog = QtWidgets.QSpinBox()
        self.spin_homog.setRange(8, RES_HOMOG_MAX)
        self.spin_homog.setValue(RES_HOMOG_DEF)
        self.spin_homog.setSuffix(" vox")
        f.addWidget(self.spin_homog, 1)
        sec.addLayout(f)

        g = QtWidgets.QGroupBox("Ensayo de compresion")
        gl = QtWidgets.QVBoxLayout(g)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Apoyo:"))
        self.cmb_apoyo = QtWidgets.QComboBox()
        # Mismo orden que resistencia.APOYOS: `apoyo_de_combo` lee el INDICE.
        self.cmb_apoyo.addItems(["deslizante", "empotrado"])
        self.cmb_apoyo.setToolTip(
            "deslizante: uz=0 en la base, sin friccion. Estado uniaxial.\n"
            "empotrado: base totalmente fija, como un ensayo con friccion\n"
            "entre probeta y platos. Coacciona lateralmente y rigidiza.\n\n"
            "La eleccion cambia E_app y la carga de fallo: hay que declararla.")
        f.addWidget(self.cmb_apoyo, 1)
        gl.addLayout(f)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Resolucion:"))
        self.spin_res_fe = QtWidgets.QSpinBox()
        self.spin_res_fe.setRange(8, 64)
        self.spin_res_fe.setValue(32)
        self.spin_res_fe.setSuffix(" vox")
        f.addWidget(self.spin_res_fe, 1)
        gl.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Direccion:"))
        self.cmb_eje_fe = QtWidgets.QComboBox()
        self.cmb_eje_fe.addItems(["Z (axial)", "X", "Y", "Los tres ejes"])
        self.cmb_eje_fe.setToolTip(
            "Con los tres ejes sale E_max/E_min, una anisotropia MECANICA\n"
            "independiente del DA del tensor MIL, que es puramente geometrico,\n"
            "y contrastable con Ex/Ey/Ez de la homogeneizacion periodica.\n\n"
            "Cuesta el triple de tiempo: son tres sistemas completos.")
        f.addWidget(self.cmb_eje_fe, 1)
        gl.addLayout(f)

        self.btn_fe = QtWidgets.QPushButton("Resolver y estimar el fallo")
        self.btn_fe.clicked.connect(self.ensayo_fe)
        gl.addWidget(self.btn_fe)

        self.btn_febio = QtWidgets.QPushButton("Analizar con FEBio…")
        self.btn_febio.setToolTip(
            "El mismo ensayo resuelto en FEBio, con la malla de ladrillos de\n"
            "la app y/o una malla suave de tetraedros cuadráticos, en lineal\n"
            "y no lineal. Usa la estructura activa, el VOI y la resolución,\n"
            "dirección y apoyo de este panel.")
        self.btn_febio.clicked.connect(self.analizar_febio)
        gl.addWidget(self.btn_febio)

        f = QtWidgets.QHBoxLayout()
        self.btn_dist = QtWidgets.QPushButton("Distribuciones…")
        self.btn_dist.setToolTip(
            "Histograma de la deformacion efectiva y de von Mises.\n"
            "Pistoia no mira el maximo sino el percentil 98, que es un punto\n"
            "de la cola: para saber si es robusto hay que ver la cola entera.")
        self.btn_dist.clicked.connect(self.distribuciones_fe)
        f.addWidget(self.btn_dist)

        self.btn_conv = QtWidgets.QPushButton("Convergencia…")
        self.btn_conv.setToolTip(
            "Resuelve a varias resoluciones y traza E_app frente al tamano\n"
            "de elemento. Es lo que convierte un E_app en citable o no.\n\n"
            "Tarda: son varios ensayos completos, el mas fino el mas caro.")
        self.btn_conv.clicked.connect(self.convergencia_fe)
        f.addWidget(self.btn_conv)
        gl.addLayout(f)

        # El rotulo cabe en los 330 px del panel; el protocolo va en el
        # tooltip y en la cabecera del dialogo, no en el boton.
        self.btn_paper = QtWidgets.QPushButton(
            "Analisis comparado VOI / spinodoide…")
        self.btn_paper.setStyleSheet("font-weight:bold; padding:6px;")
        self.btn_paper.setToolTip(
            "Aplica a las DOS estructuras el mismo ensayo del articulo de\n"
            "Tapia, Gonzalez, Vidal y Salinas (Biology 2026;15:722):\n"
            "tejido lineal elastico isotropo de 18 GPa y nu = 0.30, carga\n"
            "axial de compresion de 100 N y apoyo empotrado en la cara\n"
            "opuesta.\n\n"
            "Devuelve las dos variables que ese trabajo reporta:\n"
            "  · distribucion espacial de la tension equivalente de von Mises\n"
            "  · distribucion espacial de la deformacion total (mm)\n\n"
            "Los dos paneles quedan coloreados con la MISMA escala, que es lo\n"
            "que permite compararlos a ojo. La malla aqui es hexaedrica de un\n"
            "voxel, no TET10 de ANSYS: los patrones son comparables, las\n"
            "cifras absolutas no.")
        self.btn_paper.clicked.connect(self.analisis_comparado)
        gl.addWidget(self.btn_paper)

        self.chk_ajustar_antes = QtWidgets.QCheckBox(
            "Ajustar antes al VOI (elige el de menor error)")
        self.chk_ajustar_antes.setChecked(True)
        self.chk_ajustar_antes.setToolTip(
            "Corre los tres metodos que comparten la misma definicion de\n"
            "error —rapido, completo y equitativo—, se queda con el de menor\n"
            "error y ensaya ESE. El desempate mecanico se excluye a proposito:\n"
            "su error incluye dos terminos mas y no es comparable con el de\n"
            "los otros.\n\n"
            "Sin esto se ensaya el spinodoide que haya en pantalla, que puede\n"
            "no tener nada que ver con el VOI: a densidad igual, los angulos\n"
            "de cono cambian E_app en un factor de 15.\n\n"
            "TARDA: del orden de veinte minutos a la resolucion del VOI, mas\n"
            "el ensayo. Desmarcalo si ya has ajustado.")
        gl.addWidget(self.chk_ajustar_antes)
        nota = QtWidgets.QLabel(
            "Criterio de Pistoia: la carga a la que el 2% del tejido supera el "
            "0.7% de deformacion efectiva. Los dos valores son convenciones "
            "calibradas en radio distal humano, no constantes fisicas.\n"
            "El ensayo deja ademas los campos de deformacion efectiva y de "
            "von Mises listos para colorear.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)
        sec.addWidget(g)

        self._libs(sec, "Librerias: numpy · scipy.sparse y scipy.sparse.linalg (cg, splu) · matplotlib (distribuciones y convergencia)")

        sec = self._seccion(lay, "Simulaciones in silico")
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Estructura:"))
        self.cmb_sim_estructura = QtWidgets.QComboBox()
        self.cmb_sim_estructura.addItems(["VOI de referencia",
                                          "Candidato activo"])
        self.cmb_sim_estructura.setToolTip(
            "El VOI es el gemelo digital del animal: la simulacion responde\n"
            "que le pasaria a ESE hueso. El candidato sirve para comparar\n"
            "si una familia sintetica pierde rigidez como el hueso real.")
        f.addWidget(self.cmb_sim_estructura, 1)
        sec.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Protocolo:"))
        self.cmb_sim_protocolo = QtWidgets.QComboBox()
        # El orden es el de `spinpy.simulacion.PROTOCOLOS`.
        self.cmb_sim_protocolo.addItems([
            "Adelgazamiento uniforme",
            "Perdida de trabeculas finas",
            "Desuso guiado por carga",
            "Perdida y recuperacion"])
        self.cmb_sim_protocolo.setToolTip(
            "Que hueso se pierde primero:\n\n"
            "Adelgazamiento: la superficie, capa a capa. Envejecimiento.\n"
            "Trabeculas finas: las de menor espesor, enteras. OVX,\n"
            "  remodelado acelerado.\n"
            "Desuso: la superficie menos deformada en el ensayo de cada paso\n"
            "  (mecanostato). Inmovilizacion, microgravedad.\n"
            "Recuperacion: perdida por trabeculas finas y despues\n"
            "  engrosamiento hasta la masa inicial. Farmaco anabolico.")
        f.addWidget(self.cmb_sim_protocolo, 1)
        sec.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Pasos:"))
        self.spin_sim_pasos = QtWidgets.QSpinBox()
        self.spin_sim_pasos.setRange(1, 20)
        self.spin_sim_pasos.setValue(8)
        f.addWidget(self.spin_sim_pasos)
        f.addWidget(QtWidgets.QLabel("Hueso por paso:"))
        self.spin_sim_paso = QtWidgets.QSpinBox()
        self.spin_sim_paso.setRange(1, 25)
        self.spin_sim_paso.setValue(5)
        self.spin_sim_paso.setSuffix(" %")
        self.spin_sim_paso.setToolTip(
            "Fraccion del hueso INICIAL que se retira en cada paso, igual en\n"
            "todos los protocolos. En el fallo progresivo no se usa: alli cada\n"
            "paso ablanda el tejido que rompe.")
        f.addWidget(self.spin_sim_paso)
        sec.addLayout(f)

        f = QtWidgets.QHBoxLayout()
        self.btn_sim = QtWidgets.QPushButton("Simular perdida osea…")
        self.btn_sim.setToolTip(
            "Retira hueso en pasos y, en cada uno, mide la morfometria y\n"
            "resuelve el ensayo de compresion. Devuelve la rigidez y la carga\n"
            "de fallo frente a BV/TV, y el 3D de cada paso.\n\n"
            "Tarda: un ensayo y una morfometria completa por paso.")
        self.btn_sim.clicked.connect(self.simular_perdida_gui)
        f.addWidget(self.btn_sim)
        self.btn_fallo = QtWidgets.QPushButton("Fallo progresivo…")
        self.btn_fallo.setToolTip(
            "Carga, ablanda el tejido que el criterio de Pistoia da por roto y\n"
            "vuelve a cargar. La carga maxima de la serie estima la resistencia\n"
            "ultima; como cae dice si la estructura es fragil o reparte el dano.")
        self.btn_fallo.clicked.connect(self.fallo_progresivo_gui)
        f.addWidget(self.btn_fallo)
        sec.addLayout(f)

        nota = QtWidgets.QLabel(
            "Todos los protocolos retiran la misma cantidad de hueso por paso: "
            "las diferencias de rigidez entre ellos son de arquitectura, no de "
            "masa. Un paso es una cantidad de hueso, no un tiempo. La mecanica "
            "usa la resolucion, el apoyo y la direccion del ensayo de "
            "compresion.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        # En el panel de ancho fijo el layout le daba una linea de menos y se
        # comia "compresion." (visto en la captura del informe). Se reserva la
        # altura que pide el texto a ese ancho, con margen para el ingles.
        nota.setMinimumHeight(nota.heightForWidth(300) + 16)
        nota.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        sec.addWidget(nota)
        self._libs(sec, "Librerias: numpy · scipy.ndimage (distancias, espesor local) · scipy.sparse.linalg (ensayo) · matplotlib · pyvistaqt")

        sec = self._seccion(lay, "Exportacion")
        gl = sec
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Malla:"))
        self.cmb_solido = QtWidgets.QComboBox()
        self.cmb_solido.addItems(["Hexaedrica (voxeles)", "Tetraedrica TET10"])
        self.cmb_solido.setToolTip(
            "Hexaedrica: un elemento por voxel. Exacta, instantanea y no puede\n"
            "fallar al mallar, pero el elemento es lineal y la superficie queda\n"
            "escalonada. Es el metodo estandar de micro-elementos finitos.\n\n"
            "TET10: superficie suavizada y tetraedralizada. Elemento cuadratico\n"
            "y superficie lisa, pero mas lenta y pierde algo de volumen.")
        f.addWidget(self.cmb_solido, 1)
        self.cmb_solido.currentIndexChanged.connect(self._sync_formatos)
        gl.addLayout(f)

        # Casillas de formato. Antes se escribian los cuatro archivos siempre,
        # y con TET10 a resolucion alta el .inp solo son cientos de MB que
        # muchas veces no se quieren: para mirar la malla basta el .vtu y para
        # imprimir basta el .stl.
        gl.addWidget(QtWidgets.QLabel("Formatos:"))
        rej = QtWidgets.QGridLayout()
        self.chk_vtu = QtWidgets.QCheckBox(".vtu")
        self.chk_vtu.setToolTip(
            "Malla volumetrica para ParaView o pyvista. Lleva la calidad y el\n"
            "volumen de cada elemento, para localizar los degenerados antes de\n"
            "mandar la malla a un solver.")
        self.chk_inp = QtWidgets.QCheckBox(".inp  (Abaqus)")
        self.chk_inp.setToolTip(
            "Abaqus / CalculiX, con material y los conjuntos de nodos BASE y\n"
            "TECHO. Es el archivo mas pesado con diferencia.")
        self.chk_apdl = QtWidgets.QCheckBox(".apdl  (ANSYS)")
        self.chk_stl = QtWidgets.QCheckBox(".stl")
        self.chk_feb = QtWidgets.QCheckBox(".feb  (FEBio)")
        # Columna izquierda los formatos para mirar e imprimir, derecha los de
        # solver. Los dos rotulos largos van juntos en la derecha: en la
        # izquierda invadian la casilla vecina, que solo son 165 px.
        for c, fila, col in ((self.chk_vtu, 0, 0), (self.chk_inp, 0, 1),
                             (self.chk_stl, 1, 0), (self.chk_apdl, 1, 1),
                             (self.chk_feb, 2, 1)):
            c.setChecked(True)
            rej.addWidget(c, fila, col)
        rej.setColumnStretch(0, 1)
        rej.setColumnStretch(1, 1)
        gl.addLayout(rej)

        self.btn_exp = QtWidgets.QPushButton("Exportar solido…")
        self.btn_exp.clicked.connect(self.exportar_solido)
        gl.addWidget(self.btn_exp)
        nota = QtWidgets.QLabel(
            "Se exporta a la resolucion de MEDIDA, no a la de vista, "
            "regenerando con la misma semilla. El .stl solo existe por la via "
            "TET10: es la superficie suavizada y reparada, cerrada y por tanto "
            "imprimible. La hexaedrica no la produce porque seria la piel "
            "escalonada de los voxeles.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#666; font-size:10px;")
        gl.addWidget(nota)

        self._libs(sec, "Librerias: pyvista (.vtu y .stl) · tetgen (TET10) · pymeshfix (estanqueidad de la superficie)")

        sec = self._seccion(lay, "Lote y varianza")
        gl = sec
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Replicas por VOI:"))
        self.spin_repl = QtWidgets.QSpinBox()
        self.spin_repl.setRange(2, 60)
        self.spin_repl.setValue(10)
        f.addWidget(self.spin_repl, 1)
        gl.addLayout(f)
        f = QtWidgets.QHBoxLayout()
        f.addWidget(QtWidgets.QLabel("Modo:"))
        self.cmb_modo_lote = QtWidgets.QComboBox()
        self.cmb_modo_lote.addItems(["rapido (A+B)", "completo (A+B+C)"])
        f.addWidget(self.cmb_modo_lote, 1)
        gl.addLayout(f)
        self.btn_lote = QtWidgets.QPushButton("Elegir VOIs y correr lote…")
        self.btn_lote.setToolTip(
            "Ajusta cada VOI, genera las replicas y descompone la varianza\n"
            "en entre-especimenes y dentro, con ICC, N efectivo y TOST.\n\n"
            "TARDA MUCHO: es un ajuste completo por VOI mas K generaciones.")
        self.btn_lote.clicked.connect(self.correr_lote_gui)
        gl.addWidget(self.btn_lote)
        nota = QtWidgets.QLabel(
            "El numero de filas NO es el N del estudio. Las replicas son "
            "tecnicas: reducen la incertidumbre de cada espécimen pero no "
            "anaden grados de libertad para comparar entre animales.")
        nota.setWordWrap(True)
        nota.setStyleSheet("color:#b62324; font-size:10px;")
        gl.addWidget(nota)

        self._libs(sec, "Librerias: numpy · pandas · joblib (paralelo) · scipy.stats (TOST)")

        self.barra = QtWidgets.QProgressBar()
        self.barra.setRange(0, 0)
        self.barra.hide()
        lay.addWidget(self.barra)

        lay.addStretch(1)
        return w
    def _tabla(self):
        self.tabla = QtWidgets.QTableWidget(len(TABLA), 4)
        self.tabla.setHorizontalHeaderLabels(
            ["Metrica", "VOI", "Spinodoide", "Dif. relativa"])
        self.tabla.verticalHeader().setVisible(False)
        self.tabla.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.Stretch)
        for i, (_x, etq, _y, uni) in enumerate(TABLA):
            it = QtWidgets.QTableWidgetItem(f"{etq}  [{uni}]" if uni else etq)
            it.setFlags(QtCore.Qt.ItemIsEnabled)
            self.tabla.setItem(i, 0, it)
            for j in (1, 2, 3):
                self.tabla.setItem(i, j, QtWidgets.QTableWidgetItem("—"))
        return self.tabla

    # -- parametros --------------------------------------------------------

    def _params(self, resolucion):
        # Correccion F4: los deslizadores de rotacion solo cubren +-90 grados y
        # no pueden representar todas las orientaciones de un VOI. Cuando el
        # ajuste deriva una R del eje principal, esa R MANDA sobre los
        # deslizadores hasta que el usuario mueva uno, que es el momento en que
        # expresa que quiere controlar la orientacion a mano.
        dual = self._fam == "dual-lattice"
        if self.R_forzada is not None:
            R = self.R_forzada
        else:
            pre = "d_rot" if dual else "rot"
            R = euler_R(self.sl[pre + "x"].valor(), self.sl[pre + "y"].valor(),
                        self.sl[pre + "z"].valor())
        if dual:
            return dict(
                familia="dual-lattice",
                resolution=int(resolucion),
                celdas=float(self.sl["d_celdas"].valor()),
                rho=self.sl["d_dens"].valor() / 100.0,
                estiramiento=(self.sl["d_ex"].valor(), self.sl["d_ey"].valor(),
                              self.sl["d_ez"].valor()),
                irregularidad=float(self.sl["d_irr"].valor()),
                R=R,
                seed=int(self.spin_semilla.value()),
            )
        return dict(
            familia="spinodoide",
            resolution=int(resolucion),
            wave_number=self.sl["wave"].valor() * np.pi,
            num_waves=int(self.sl["nw"].valor()),
            thetas=[self.sl["thx"].valor(), self.sl["thy"].valor(),
                    self.sl["thz"].valor()],
            rho=self.sl["dens"].valor() / 100.0,
            R=R,
            esquema="rechazo" if self.cmb_esq.currentIndex() == 0 else "equitativo",
            seed=int(self.spin_semilla.value()),
        )

    def _soltar_R(self):
        """El usuario toca una rotacion: deja de mandar la R del ajuste."""
        self._soltar_R_fam(self._fam)

    def _soltar_R_fam(self, fam):
        """Libera la R del ajuste de UNA familia: la de la rotacion tocada.

        Cada familia tiene sus deslizadores de rotacion y su propia R forzada.
        Tocar la rotacion del dual-lattice no puede soltar la orientacion que
        el ajuste impuso al spinodoide.
        """
        if fam == self._fam:
            self.R_forzada = None
        else:
            self._cand[fam]["R"] = None
        self._lab_R[fam].setText("")
        self._pedir_vista()

    # -- familias de microestructura ---------------------------------------

    # Campos cacheados cuya entrada "spin" pertenece a la familia activa.
    _CACHES = (("_esp", "esp"), ("_fe", "fe"), ("_vm", "vm"),
               ("_desp", "desp"))

    @staticmethod
    def _cand_vacio():
        return {"BW": None, "m": None, "vigente": False, "malla": None,
                "R": None, "esp": None, "fe": None, "vm": None, "desp": None}

    def _guardar_fam(self):
        c = self._cand[self._fam]
        c.update(BW=self.BW_vista, m=self.m_spin,
                 vigente=self._m_spin_vigente, malla=self.malla_spin,
                 R=self.R_forzada)
        for attr, k in self._CACHES:
            c[k] = getattr(self, attr).get("spin")

    def _cargar_fam(self, f):
        c = self._cand[f]
        self.BW_vista, self.m_spin = c["BW"], c["m"]
        self._m_spin_vigente, self.malla_spin = c["vigente"], c["malla"]
        self.R_forzada = c["R"]
        for attr, k in self._CACHES:
            cache = getattr(self, attr)
            if c[k] is None:
                cache.pop("spin", None)
            else:
                cache["spin"] = c[k]
        self.vis_spin = self._vis[f]
        self.lab_R = self._lab_R[f]
        self._fam = f

    def _activar(self, f):
        """Deja la familia `f` en los atributos de candidato."""
        if f == self._fam:
            return
        self._guardar_fam()
        self._cargar_fam(f)

    def _de(self, f, clave):
        """Un dato de la familia `f`, este activa o no."""
        if f == self._fam:
            self._guardar_fam()
        return self._cand[f][clave]

    def _familias(self):
        return list(FAMILIAS) if self._modo_fam == "ambas" else [self._modo_fam]

    def _etq_fam(self, f=None):
        f = f or self._fam
        return _("Dual-lattice") if f == "dual-lattice" else _("Spinodoide")

    def _clave_res(self, nombre):
        """Clave de `_res` de la familia activa. La del spinodoide conserva el
        nombre de siempre, para que un JSON exportado antes se siga leyendo."""
        return nombre if self._fam == "spinodoide" else f"{nombre}_dual"

    def _densidad_actual(self):
        return self.sl["d_dens" if self._fam == "dual-lattice"
                       else "dens"].valor()

    def _para_cada(self, fn):
        """Ejecuta `fn` con cada familia activa, sin hilos de por medio."""
        previa = self._fam
        fams = self._familias()
        for f in fams:
            self._activar(f)
            fn()
        self._activar(previa if previa in fams else fams[0])

    def _en_cada_familia(self, accion):
        """Lanza `accion` para cada familia activa, UNA DETRAS DE OTRA.

        Las acciones pesadas corren en un hilo y terminan en un metodo
        `*_listo` que lee el estado de la familia activa. Por eso no se pueden
        lanzar las dos a la vez: la siguiente familia se activa cuando el hilo
        de la anterior ha terminado, no antes. Si la accion no lanza hilo
        -porque el usuario cancelo un dialogo, o porque lo pedido ya estaba
        calculado- se pasa a la siguiente sin esperar.
        """
        fams = self._familias()
        cola = list(fams)

        def siguiente():
            if not cola:
                self._activar(fams[0])
                self._redibujar()
                return
            self._activar(cola.pop(0))
            previo = self.hilo
            accion()
            h = self.hilo
            if h is not None and h is not previo and h.isRunning():
                h.finished.connect(
                    lambda: QtCore.QTimer.singleShot(0, siguiente))
            else:
                siguiente()

        siguiente()

    def _confirmar_orientacion(self):
        """Pregunta antes de ensayar un candidato girado respecto al VOI.

        El error del ajuste mira el DA como escalar y no ve la orientacion: en
        10 de los 12 VOIs equinos el spinodoide ganador quedo a mas de 30° del
        eje del hueso (Estudio_Familias, 4.5). La morfometria no lo nota; un
        ensayo si, porque carga el candidato en otra direccion que el hueso.

        Devuelve True si se puede seguir. Sin VOI o sin candidato no hay nada
        que comparar, y con estructuras casi isotropas el angulo es ruido:
        en esos casos no se pregunta (ver `spinpy.avisos.desalineacion`).
        """
        if self.VOI is None or self.BW_vista is None:
            return True

        def con_eje(m, BW, sp):
            # La fabrica ya medida sirve si esta al dia; si no, basta el MIL,
            # que es mucho mas barato que la morfometria entera.
            if (isinstance(m, dict) and m.get("dir_principal") is not None
                    and np.isfinite(m.get("DA", np.nan))):
                return m
            return tensor_mil(BW, sp)

        m_c = self.m_spin if getattr(self, "_m_spin_vigente", False) else None
        d = desalineacion(
            con_eje(self.m_voi, self.VOI, self.VOI_spacing),
            con_eje(m_c, self.BW_vista, self._spacing(self.BW_vista.shape[0])))
        if not d["aviso"]:
            return True
        if self._auto is not None:
            # Sin nadie delante a quien preguntar: se ensaya igual y se anota.
            # El informe marca ese resultado como desalineado por su cuenta
            # (`informe.comprobar`), asi que no se pierde la advertencia.
            self._auto["avisos"].append(
                _("{fam}: eje a {ang}° del VOI; ensayado igualmente, el "
                  "informe lo marca como desalineado").format(
                    fam=self._etq_fam(), ang=f"{d['angulo_deg']:.0f}"))
            return True
        r = QtWidgets.QMessageBox.question(
            self, _("Candidato mal orientado"),
            _("El eje principal del {fam} esta a {ang}° del eje del VOI "
              "(umbral {umbral}°).\n\nEl ensayo cargaria el candidato en una "
              "direccion que no es la del hueso, y su rigidez no seria "
              "comparable con la del VOI. Suele pasar cuando la fabrica del "
              "candidato es casi plana y la correccion de orientacion no "
              "puede alinearla; en el banco de VOIs, el dual-lattice quedo "
              "alineado en todos.\n\n¿Ensayar de todos modos?").format(
                fam=self._etq_fam(), ang=f"{d['angulo_deg']:.0f}",
                umbral=f"{d['umbral_deg']:.0f}"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No)
        return r == QtWidgets.QMessageBox.Yes

    def _campos(self, attr):
        """Campos 3D cacheados en `attr` de TODAS las familias activas y del
        VOI, para fijar una escala de color comun a todos los paneles."""
        k = dict(self._CACHES)[attr]
        vals = [c for c, _sp in getattr(self, attr).values()]
        for f in self._familias():
            if f != self._fam and self._cand[f][k] is not None:
                vals.append(self._cand[f][k][0])
        return vals

    def _fijar_R(self, fam, R):
        arr = None if R is None else np.asarray(R, float)
        if fam == self._fam:
            self.R_forzada = arr
        else:
            self._cand[fam]["R"] = arr

    @staticmethod
    def _generar(p):
        """Mascara de cualquier familia a partir de `_params`. Sin Qt: corre
        dentro de los hilos de trabajo."""
        p = dict(p)
        if p.pop("familia", "spinodoide") == "dual-lattice":
            BW, campo, info = generar_dual_lattice(**p)
            info.setdefault("esquema", "dual-lattice")
            return BW, campo, info
        return generar_mascara(**p)

    def _barra_familia(self):
        """Selector de familia de microestructura.

        Va en una barra de herramientas y no en el panel: no es un parametro
        de una etapa sino el objeto sobre el que trabajan TODAS las etapas, y
        tiene que estar a la vista mientras se recorre el panel entero.
        """
        tb = self.addToolBar("Familia")
        tb.setObjectName("barra_familia")
        tb.setMovable(False)
        lab = QtWidgets.QLabel("Familia:")
        lab.setContentsMargins(6, 0, 6, 0)
        tb.addWidget(lab)
        self.cmb_familia = QtWidgets.QComboBox()
        self.cmb_familia.addItems(["Spinodoide", "Dual-lattice", "Ambas"])
        self.cmb_familia.setToolTip(
            "Spinodoide: conjunto de nivel de un campo aleatorio gaussiano\n"
            "(Kumar et al. 2020).\n"
            "Dual-lattice: red dual de una teselacion de Delaunay, con cuatro\n"
            "barras por nudo (Vafaeefar et al. 2022).\n\n"
            "Ambas: se generan, miden y ajustan las DOS a la vez, con un panel\n"
            "3D y una columna de la tabla por familia. Las acciones pesadas\n"
            "(tensor, ensayo, exportacion, lote) se hacen una familia detras\n"
            "de otra, con el mismo material y la misma resolucion.")
        self.cmb_familia.currentIndexChanged.connect(self._cambiar_familia)
        tb.addWidget(self.cmb_familia)

    def _barra_publicacion(self):
        """Boton del informe para publicacion.

        En la barra de herramientas y no en el panel de Exportacion: no
        exporta una etapa, describe la sesion ENTERA —VOI, ajuste, mecanica— y
        se usa al final de cualquier recorrido, este donde este el panel.
        """
        tb = self.addToolBar("Publicacion")
        tb.setObjectName("barra_publicacion")
        tb.setMovable(False)
        self.act_informe = QtWidgets.QAction("Informe para publicación…", self)
        self.act_informe.setToolTip(
            "Párrafo de métodos (ES/EN) con los valores usados, citabilidad\n"
            "de cada resultado y paquete de reproducción verificable.")
        self.act_informe.triggered.connect(self.informe_publicacion)
        tb.addAction(self.act_informe)

        # Al lado y no dentro del mismo dialogo: el de la izquierda describe lo
        # que YA se calculo; este calcula antes todo el recorrido. Son dos
        # decisiones distintas —una tarda segundos, la otra horas— y un boton
        # que a veces tarda horas segun una casilla sorprende.
        self.act_informe_auto = QtWidgets.QAction(
            "Informe para publicación (auto)…", self)
        self.act_informe_auto.setToolTip(
            "Corre sobre el VOI cargado todo el recorrido —mejor ajuste,\n"
            "morfometría completa, tensor elástico, ensayo de compresión y\n"
            "análisis comparado, más las etapas opcionales que marques— y\n"
            "escribe el informe al final. Puede tardar horas.")
        self.act_informe_auto.triggered.connect(self.informe_auto)
        tb.addAction(self.act_informe_auto)

        # Al lado del informe automatico y con su misma logica: se elige todo
        # en una ventana y se espera. Encadena las MISMAS `_febio_una` que el
        # boton «Analizar con FEBio…» del panel.
        self.act_fem_auto = QtWidgets.QAction("FEM automático (FEBio)…", self)
        self.act_fem_auto.setToolTip(
            "Resuelve en FEBio las estructuras, protocolos (ensayo de la app,\n"
            "Tapia et al. 2026, homogeneización) y mallas (ladrillos y/o\n"
            "tetraedros suaves) marcados, y los compara con la app.\n"
            "Puede tardar horas.")
        self.act_fem_auto.triggered.connect(self.fem_auto)
        tb.addAction(self.act_fem_auto)

        # Avance de la cadena y boton para pararla. En la barra de estado y no
        # en un dialogo: la ventana tiene que seguir a la vista y usable para
        # mirar resultados intermedios durante horas.
        self.lab_auto = QtWidgets.QLabel()
        self.lab_auto.setVisible(False)
        self.btn_auto_detener = QtWidgets.QPushButton("Detener tras esta etapa")
        self.btn_auto_detener.setToolTip(
            "Un calculo en marcha no se puede interrumpir sin perderlo: la\n"
            "etapa en curso termina y las siguientes ya no se lanzan. Lo\n"
            "calculado queda en la sesión.")
        self.btn_auto_detener.clicked.connect(self._auto_detener)
        self.btn_auto_detener.setVisible(False)
        self.statusBar().addPermanentWidget(self.lab_auto)
        self.statusBar().addPermanentWidget(self.btn_auto_detener)

    def _cambiar_familia(self, i):
        modo = MODOS_FAMILIA[int(i)]
        if modo == self._modo_fam:
            return
        self._modo_fam = modo
        fams = self._familias()
        self._activar(fams[0])
        self.grp_spin.setVisible("spinodoide" in fams)
        self.grp_dual.setVisible("dual-lattice" in fams)
        for f in FAMILIAS:
            self._caja[f].setVisible(f in fams)
        self._reconstruir_tabla()
        self._rellenar_tabla()
        faltan = [f for f in fams if self._de(f, "BW") is None]
        if faltan and not self._aplicando:
            self._generar_vista()
        else:
            self._redibujar()
        self.statusBar().showMessage(_("Familia: {f}").format(
            f=", ".join(self._etq_fam(f) for f in fams)))

    def _spacing(self, resolucion):
        """mm/voxel del candidato, escalado al tamano fisico del VOI.

        Correccion F1: la malla se genera en un cubo unidad adimensional. Para
        que Tb.Th, Tb.Sp y BS/BV salgan en las mismas unidades que las del VOI
        hay que escalar al tamano fisico de este antes de medir. Sin VOI
        cargado se usa 1 mm de lado y las longitudes son relativas.
        """
        if self.VOI is not None:
            lado = float(self.VOI.shape[0] * self.VOI_spacing[0])
        else:
            lado = 1.0
        return np.full(3, lado / float(resolucion))

    # -- acciones ----------------------------------------------------------

    def _ocupado(self, si, msg="", determinada=False, total=0):
        # Durante el informe automatico los botones siguen bloqueados ENTRE
        # etapas: cada `*_listo` libera la interfaz al terminar la suya, y un
        # clic en ese instante lanzaria un calculo en medio de la cadena.
        bloquear = si or self._auto is not None
        for b in (self.btn_gen, self.btn_med, self.btn_fit, self.btn_ela,
                  self.btn_exp, self.btn_fe, self.btn_hist, self.btn_dist,
                  self.btn_conv, self.btn_disp, self.btn_lote, self.btn_todos,
                  self.btn_paper, self.cmb_familia, self.btn_sim,
                  self.btn_fallo, self.act_informe, self.act_informe_auto,
                  self.btn_febio, self.act_fem_auto):
            b.setEnabled(not bloquear)
        if si and determinada:
            self.barra.setRange(0, max(1, int(total)))
            self.barra.setValue(0)
        elif si:
            self.barra.setRange(0, 0)      # indeterminada
        self.barra.setVisible(si)
        if msg:
            self.statusBar().showMessage(msg)

    def _avance(self, i, n, etapa):
        if self.barra.maximum() != n:
            self.barra.setRange(0, max(1, n))
        self.barra.setValue(i)
        self.statusBar().showMessage(
            _("Ajustando — etapa {etapa}: {i} de {n}").format(
                etapa=etapa, i=i, n=n))

    def _pedir_vista(self):
        if self._aplicando:
            return
        if self.chk_auto.isChecked():
            # Se reinicia un temporizador corto para no lanzar un calculo por
            # cada paso intermedio mientras se arrastra el deslizador.
            if not hasattr(self, "_temp"):
                self._temp = QtCore.QTimer(self)
                self._temp.setSingleShot(True)
                self._temp.timeout.connect(self._generar_vista)
            self._temp.start(250)

    def _generar_vista(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        res = int(self.sl["resv"].valor())
        fams = self._familias()
        ps = {}
        for f in fams:
            self._activar(f)
            ps[f] = self._params(res)
        self._activar(fams[0])
        self._ocupado(True, f"Generando vista a {res}^3…")
        self._t0 = time.time()

        def tarea():
            out = {}
            for f, p in ps.items():
                BW, _x, info = Visor._generar(p)
                out[f] = (BW, info)
            return out

        self.hilo = Trabajador(tarea)
        self.hilo.listo.connect(self._vista_lista)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    def _recortar(self, malla):
        """Aplica el plano de recorte visual. Devuelve la malla tal cual si no
        hay recorte pedido o si el corte dejaria la vista vacia."""
        frac = self.sl["clip"].valor() / 100.0 if "clip" in self.sl else 0.0
        if malla is None or not malla.n_points or frac <= 0:
            return malla
        eje = "xyz"[self.cmb_eje.currentIndex()]
        i = "xyz".index(eje)
        lo, hi = malla.bounds[2 * i], malla.bounds[2 * i + 1]
        origen = [0.0, 0.0, 0.0]
        origen[i] = lo + frac * (hi - lo)
        try:
            cortada = malla.clip(normal=eje, origin=origen, invert=True)
            return cortada if cortada.n_points else malla
        except Exception:
            return malla

    @staticmethod
    def _vaciar(vis):
        """Quita los actores CONSERVANDO las luces.

        `Plotter.clear()` destruye tambien el light kit: tras llamarlo el
        renderer se queda con 0 luces y todo lo que se anada despues se dibuja
        SIN ILUMINAR — una silueta plana de un solo color que parece un bloque
        macizo aunque la estructura tenga BV/TV de 0.28. `clear_actors()` no
        toca las luces. El rearme es por si acaso: una vez destruidas,
        `clear_actors()` no las recupera.
        """
        vis.clear_actors()
        if not vis.renderer.lights:
            vis.enable_lightkit()

    # Mismo material en los dos paneles. Colorear uno con un mapa de color y el
    # otro liso ROMPE la comparacion visual: el ojo compara colores antes que
    # formas, y la altura —que es lo que codificaba el viridis— ya se ve en la
    # geometria, asi que ese color no aportaba informacion.
    COLOR_HUESO = "#d9c7a0"

    # Fondo de TRABAJO frente a fondo de PUBLICACION. El blanco puro es el
    # correcto para una figura de tesis, pero en pantalla deja la silueta de la
    # estructura sin nada contra lo que recortarse y el volumen se lee mal. El
    # degradado da esa referencia sin meter color en el material.
    FONDO_PANTALLA = ("#dfe5ec", "#f7f9fb")
    FONDO_FIGURA = ("white", None)

    # Fabrica MIL. Un solo color en los DOS paneles a proposito: lo que se
    # compara es la direccion, y colores distintos invitarian a leer el color
    # como si significase algo.
    COLOR_EJE = "#c2185b"
    COLOR_ELIPSOIDE = "#1b4f72"

    # Umbral de DA2 por debajo del cual la fabrica es un PLANO y no un eje. Es
    # el mismo `da2_max` con el que la correccion L1 se abstiene de rotar
    # (spinpy/fit.py): si el candidato es rigido en un plano, ninguna rotacion
    # lo convierte en axial y su "eje mayor" es uno cualquiera de ese plano.
    DA2_DEGENERADO = 1.06

    def _estilo_render(self, vis):
        """Fondo y suavizado de un panel.

        Se aplica al crear el panel y cada vez que se cambia el fondo. Va en
        try/except porque son pasos de renderizado que dependen del driver de
        video: si el equipo no los soporta, la aplicacion tiene que seguir
        dibujando la geometria correcta aunque salga mas fea.
        """
        a, b = (self.FONDO_FIGURA if self.chk_fondo.isChecked()
                else self.FONDO_PANTALLA)
        try:
            vis.set_background(a, top=b)
        except Exception:
            vis.set_background(a)
        # SSAA: la isosuperficie de una trabecula es todo aristas finas y sin
        # suavizado salen dentadas, muy visible en las capturas del informe.
        try:
            vis.enable_anti_aliasing("ssaa")
        except Exception:
            pass

    def _aplicar_fondo(self):
        for v in (self.vis_voi, *self._vis.values()):
            self._estilo_render(v)
        self._redibujar()

    def _ssao(self, vis, malla):
        """Oclusion ambiental de pantalla, con el radio a escala de la escena.

        Es la mejora de render que mas cambia la lectura: una trabecula es una
        marana de barras cruzadas y, sin sombra de contacto donde una pasa por
        detras de otra, no hay forma de saber cual esta delante; la
        profundidad se veia solo por el escorzo.

        El radio NO puede ser fijo. SSAO oscurece los puntos que tienen
        geometria dentro de una esfera de ese radio, asi que en milimetros
        absolutos un VOI de 5 mm de lado y un cubo unidad de 1 mm darian
        resultados opuestos: uno casi sin sombra, el otro entero en sombra. Se
        fija como fraccion de la diagonal de la caja.
        """
        if not self.chk_ssao.isChecked():
            try:
                vis.disable_ssao()
            except Exception:
                pass
            return
        try:
            if malla is None or not malla.n_points:
                return
            b = np.asarray(malla.bounds, float)
            diag = float(np.linalg.norm([b[1] - b[0], b[3] - b[2], b[5] - b[4]]))
            if not np.isfinite(diag) or diag <= 0:
                return
            vis.enable_ssao(radius=0.045 * diag, bias=0.0009 * diag,
                            kernel_size=128, blur=True)
        except Exception:
            pass

    def _fabrica(self, vis, malla, m):
        """Dibuja el eje principal del MIL y el elipsoide de fabrica.

        POR QUE ESTA VISTA EXISTE
          El tensor MIL se calcula desde el principio y su eje mayor decide si
          una estructura es rigida en Z o en Y. Pero ese dato solo aparecia
          como un escalar (DA) en la tabla y como una direccion enterrada en un
          informe. En el ajuste del VOI proximal de H4 el candidato ganador
          tenia su eje rigido a 89.8 grados del hueso, y eso se supo despues de
          477 s de ajuste y de abrir un dialogo. Dibujado, se ve de golpe.

        CONVENCION, que es la parte que se presta a equivocarse
          El ajuste del tensor es 1/MIL(n)^2 = n' M n, de modo que la longitud
          del MIL en la direccion del autovector i vale 1/sqrt(lambda_i): el
          autovalor MENOR corresponde al MIL MAYOR. Los autovalores vienen
          ordenados de menor a mayor, asi que el semieje LARGO del elipsoide es
          el primero, el mismo que `dir_principal`. Dibujarlo con semiejes
          proporcionales a lambda -y no a 1/sqrt(lambda)- lo pintaria aplastado
          justo en la direccion en la que la estructura es mas continua.

        El elipsoide va en alambre y no en superficie translucida: relleno
        taparia la estructura, que es lo que se ha venido a mirar.
        """
        if malla is None or not malla.n_points or not m:
            return
        if not m.get("MIL_valid"):
            return
        lam = np.asarray(m.get("eigenvalues", []), float)
        V = np.asarray(m.get("eigenvectors", []), float)
        if lam.shape != (3,) or V.shape != (3, 3) or np.any(lam <= 0):
            return

        b = np.asarray(malla.bounds, float)
        ctr = np.array([(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2])
        lado = float(np.mean([b[1] - b[0], b[3] - b[2], b[5] - b[4]]))
        if not np.isfinite(lado) or lado <= 0:
            return

        # El semieje mayor se lleva al 62% del lado, es decir MAS ALLA de la
        # media caja: si se dibujase inscrito quedaria enterrado dentro de la
        # estructura -que es opaca- y no se veria nada. Asomando por las caras
        # se lee desde fuera, y la parte que queda dentro la tapa la trabecula,
        # que es lo correcto.
        r = 1.0 / np.sqrt(lam)
        r = r / r.max() * (0.62 * lado)

        R = V.copy()
        if np.linalg.det(R) < 0:             # eigh puede devolver terna levogira
            R[:, 2] *= -1.0

        T = np.eye(4)
        T[:3, :3] = R @ np.diag(r)
        T[:3, 3] = ctr
        try:
            elip = pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=10)
            elip.transform(T, inplace=True)
            vis.add_mesh(elip, style="wireframe", line_width=2,
                         color=self.COLOR_ELIPSOIDE, opacity=0.55,
                         lighting=False, pickable=False)
        except Exception:
            pass

        # Eje mayor: dos flechas opuestas desde el centro. Una sola sugeriria
        # un sentido, y un eje de fabrica no lo tiene.
        d = np.asarray(m.get("dir_principal", R[:, 0]), float)
        n = float(np.linalg.norm(d))
        d = d / n if n > 0 else np.array([0.0, 0.0, 1.0])
        L = 0.85 * lado          # sale por las dos caras: siempre visible
        for signo in (1.0, -1.0):
            try:
                vis.add_mesh(pv.Arrow(start=ctr, direction=signo * d,
                                      tip_length=0.22, tip_radius=0.055,
                                      shaft_radius=0.018, scale=L),
                             color=self.COLOR_EJE, pickable=False)
            except Exception:
                pass

    def _texto_fabrica(self):
        """Rotulo del angulo entre los dos ejes principales, para el panel spin.

        Devuelve (texto, color) o None. El angulo se toma entre EJES: el signo
        de un autovector es arbitrario, asi que se usa el valor absoluto del
        producto escalar y el resultado vive en [0, 90].
        """
        a, b = self.m_voi, self.m_spin
        if not (a and b and a.get("MIL_valid") and b.get("MIL_valid")):
            return None
        try:
            da = np.asarray(a["dir_principal"], float)
            db = np.asarray(b["dir_principal"], float)
            ca = abs(float(np.dot(da / np.linalg.norm(da),
                                  db / np.linalg.norm(db))))
            ang = float(np.degrees(np.arccos(min(1.0, ca))))
        except Exception:
            return None
        txt = _("eje MIL: {ang} grados respecto al VOI").format(
            ang="%.0f" % ang)
        col = "#b62324" if ang > 30.0 else "#1a7f37"
        da2 = float(b.get("DA2", np.nan))
        if np.isfinite(da2) and da2 < self.DA2_DEGENERADO:
            txt += "\n" + _("(DA2 {da2}: rigido en un PLANO, no en un "
                             "eje)").format(da2="%.3f" % da2)
            col = "#b62324"
        return txt, col

    def _ejes(self, vis, malla):
        """Caja delimitadora con cotas en mm.

        Sin una referencia de escala los dos paneles pueden estar a zooms
        distintos y no habria forma de notarlo: se concluiria que una
        estructura tiene trabeculas mas gruesas cuando es un artefacto de la
        vista.
        """
        if malla is None or not malla.n_points:
            return
        try:
            vis.show_bounds(mesh=malla, grid=False, location="outer",
                            ticks="outside", font_size=9, color="#555555",
                            xtitle="x [mm]", ytitle="y [mm]", ztitle="z [mm]",
                            use_3d_text=False)
        except Exception:
            pass

    def _pintar(self, vis, malla, campo, clim, barra):
        """Anade la malla, coloreada por `campo` si lo hay.

        `clim` es COMUN a los dos paneles. Con escalas independientes el azul
        de un panel y el azul del otro significarian espesores distintos, y la
        comparacion visual —que es todo el proposito de tener dos paneles—
        seria enganosa.
        """
        # Un punto de especularidad: con material puramente difuso las barras
        # de la trabecula se funden unas con otras y el relieve se pierde. Es
        # apariencia, no dato: no toca ninguna medida.
        mat = dict(smooth_shading=True, ambient=0.18, diffuse=0.82,
                   specular=0.22, specular_power=18)
        if campo is None:
            vis.add_mesh(malla, color=self.COLOR_HUESO, **mat)
        else:
            vis.add_mesh(malla, scalars=campo, cmap="turbo", clim=clim,
                         show_scalar_bar=barra, **mat,
                         scalar_bar_args={"title": TITULO_BARRA.get(
                                              self.cmb_color.currentIndex(), ""),
                                          "n_labels": 5, "color": "black",
                                          "vertical": False, "height": 0.06,
                                          "position_y": 0.02})

    def _redibujar_spin(self, reset=False):
        self._vaciar(self.vis_spin)
        self.vis_spin.add_text(
            _("Dual-lattice (vista)") if self._fam == "dual-lattice"
            else _("Spinodoide (vista)"), font_size=11, color="black")
        m = self._recortar(self.malla_spin)
        if m is not None and m.n_points:
            self._pintar(self.vis_spin, m, self._campo(m, "spin"),
                         self._clim, True)
            self._ssao(self.vis_spin, m)
            if self.chk_ejes.isChecked():
                self._ejes(self.vis_spin, m)
            # La fabrica se dibuja sobre la malla SIN recortar: el recorte es
            # solo visual y su centro no es el del volumen, asi que el
            # elipsoide saldria descolocado respecto a la estructura que
            # describe.
            if self.chk_fab.isChecked() and self._m_spin_vigente:
                self._fabrica(self.vis_spin, self.malla_spin, self.m_spin)
                t = self._texto_fabrica()
                if t is not None:
                    # Arriba a la derecha: abajo a la izquierda esta la
                    # triada de orientacion y abajo del todo la barra de color.
                    self.vis_spin.add_text(t[0], position="upper_right",
                                           font_size=9, color=t[1])
        if reset:
            self.vis_spin.reset_camera()

    def _redibujar_voi(self, reset=False):
        self._vaciar(self.vis_voi)
        self.vis_voi.add_text(f'{_("VOI")}: {self.VOI_nombre}',
                              font_size=11, color="black")
        m = self._recortar(self.malla_voi)
        if m is not None and m.n_points:
            self._pintar(self.vis_voi, m, self._campo(m, "voi"),
                         self._clim, False)
            self._ssao(self.vis_voi, m)
            if self.chk_ejes.isChecked():
                self._ejes(self.vis_voi, m)
            if self.chk_fab.isChecked():
                self._fabrica(self.vis_voi, self.malla_voi, self.m_voi)
        if reset:
            self.vis_voi.reset_camera()

    # -- coloreado por espesor local ---------------------------------------

    def _campo(self, malla, cual):
        """Valores del campo activo en los vertices, o None si no procede."""
        i = self.cmb_color.currentIndex()
        if i == 0:
            return None
        dat = {1: self._esp, 2: self._fe, 3: self._vm,
               4: self._desp}[i].get(cual)
        if dat is None:
            return None
        campo, sp = dat
        # El campo de deformacion lleva NaN donde no hay hueso portante: ese
        # material no se resolvio, y pintarlo como cero diria que esta
        # descargado cuando en realidad no se sabe. `muestrear_en_puntos` usa
        # un maximo 3x3x3, que ignora los NaN de alrededor.
        return muestrear_en_puntos(np.nan_to_num(campo, nan=0.0), sp,
                                   np.asarray(malla.points, float))

    def _cambiar_color(self):
        i = self.cmb_color.currentIndex()
        if i == 0:
            self._clim = None
            self._redibujar()
            self.statusBar().showMessage(_("Coloreado por material."))
            return
        if i in (2, 3, 4):
            # Los tres campos los produce el ensayo de compresion; no se lanza
            # solo porque puede tardar minutos y el usuario no lo espera al
            # tocar un desplegable de color.
            attr = {2: "_fe", 3: "_vm", 4: "_desp"}[i]
            cache = getattr(self, attr)
            if not cache:
                QtWidgets.QMessageBox.information(
                    self, _("Falta el ensayo"),
                    _("Pulsa antes «Resolver y estimar el fallo»: los campos "
                      "de deformacion efectiva y de von Mises salen de ese "
                      "calculo."))
                self.cmb_color.setCurrentIndex(0)
                return
            vals = [c[np.isfinite(c) & (c > 0)] for c in self._campos(attr)]
            vals = [v for v in vals if v.size]
            if vals:
                todos = np.concatenate(vals)
                self._clim = (float(np.percentile(todos, 1)),
                              float(np.percentile(todos, 99)))
            self._redibujar()
            if i == 2:
                self.statusBar().showMessage(
                    _("Coloreado por deformacion efectiva [escala comun "
                      "{a}–{b}; el criterio de fallo usa 7.0e-03]").format(
                        a=f"{self._clim[0]:.2e}", b=f"{self._clim[1]:.2e}"))
            elif i == 3:
                self.statusBar().showMessage(
                    _("Coloreado por von Mises [escala comun {a}–{b} MPa a la "
                      "carga de referencia; escala con la carga porque el "
                      "problema es lineal]").format(
                        a=f"{self._clim[0]/1e6:.3f}",
                        b=f"{self._clim[1]/1e6:.3f}"))
            else:
                # El campo de desplazamiento esta en las unidades de `spacing`,
                # es decir en milimetros: no lleva el /1e6 de las tensiones.
                self.statusBar().showMessage(
                    _("Coloreado por deformacion total [escala comun {a}–{b} "
                      "mm a la carga de referencia; el maximo esta en la cara "
                      "cargada porque acumula todo el desplazamiento]").format(
                        a=f"{self._clim[0]:.4f}", b=f"{self._clim[1]:.4f}"))
            return
        # Todas las familias activas: con dos paneles de candidato, colorear
        # solo uno dejaria el otro liso y la comparacion sin sentido.
        self._en_cada_familia(self._calcular_espesor)

    def _calcular_espesor(self, luego=None):
        """Lanza el campo de espesor local si no esta ya calculado.

        Devuelve True si el campo YA estaba listo —el que llama sigue sin
        esperar— y False si se ha lanzado el calculo, en cuyo caso `luego` se
        invoca al terminar. Lo usan el coloreado y el histograma: el campo es
        el mismo y calcularlo dos veces seria pagar dos veces la parte cara
        (una transformada de distancia por cada radio de prueba).
        """
        # Solo se calcula lo que falta: con varias familias el VOI ya puede
        # tener su campo y la familia activa no, y recalcularlo todo porque
        # "algo" falte pagaria dos veces la parte cara.
        falta_spin = self.BW_vista is not None and "spin" not in self._esp
        falta_voi = self.VOI is not None and "voi" not in self._esp
        if not (falta_spin or falta_voi):
            return bool(self._esp)
        if self.hilo is not None and self.hilo.isRunning():
            return False

        BW_s = self.BW_vista if falta_spin else None
        sp_s = self._spacing(BW_s.shape[0]) if BW_s is not None else None
        VOI = self.VOI if falta_voi else None
        sp_v = self.VOI_spacing

        self._esp_luego = luego
        self._ocupado(True, _("Calculando el espesor local…"))
        self._t0 = time.time()

        def tarea():
            out = {}
            if BW_s is not None:
                out["spin"] = (espesor_local(BW_s, sp_s), sp_s)
            if VOI is not None:
                out["voi"] = (espesor_local(VOI, sp_v), sp_v)
            return out

        self.hilo = Trabajador(tarea)
        self.hilo.listo.connect(self._espesor_listo)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()
        return False

    def histograma_espesor(self):
        self._en_cada_familia(self._histograma_una)

    def _histograma_una(self):
        if not self._calcular_espesor(luego=self._mostrar_histograma):
            return
        self._mostrar_histograma()

    def _mostrar_histograma(self):
        if not self._esp:
            return
        DialogoEspesor(self._esp, self, etq_spin=self._etq_fam()).exec_()

    def _espesor_listo(self, out):
        self._esp.update(out)
        # Escala COMUN a todos los paneles —el VOI y cada familia activa—, del
        # percentil 1 al 99 para que unos pocos voxeles extremos no aplasten
        # todo el resto del rango.
        vals = [e[e > 0] for e in self._campos("_esp") if (e > 0).any()]
        if vals:
            todos = np.concatenate(vals)
            self._clim = (float(np.percentile(todos, 1)),
                          float(np.percentile(todos, 99)))
        self._ocupado(False)
        self._redibujar()
        partes = []
        for cual, etq in (("voi", "VOI"), ("spin", "spin")):
            if cual in out:
                s = estadisticas_esp(out[cual][0])
                partes.append(f"{etq}: mediana {s['mediana']:.4f} mm")
        self.statusBar().showMessage(
            _("Espesor local — {partes}   [escala comun {a}–{b} mm; sesgo "
              "~30% a 3 voxeles: usar como medida RELATIVA]").format(
                partes="  ·  ".join(partes),
                a=f"{self._clim[0]:.3f}", b=f"{self._clim[1]:.3f}"))

        # Lo que pidio el campo (por ejemplo el histograma) se atiende ahora,
        # y se limpia antes de invocarlo para que un fallo dentro no deje la
        # peticion pegada y se repita en el siguiente calculo.
        luego, self._esp_luego = self._esp_luego, None
        if luego:
            luego()

    def _redibujar(self):
        self._redibujar_voi()
        self._para_cada(self._redibujar_spin)

    # -- sincronizacion de camaras -----------------------------------------

    def _sincronizar(self, activar):
        """Enlaza o separa las camaras de los dos paneles.

        Se comparte el MISMO objeto vtkCamera. Es correcto porque el candidato
        se escala al tamano fisico del VOI (correccion F1), de modo que las dos
        escenas ocupan el mismo volumen y una camara comun las encuadra igual.

        Sin esto cada panel tiene su zoom, y comparar grosores a ojo entre dos
        vistas con encuadres distintos es sencillamente invalido.
        """
        # Todos los paneles de candidato, se vean o no: una familia oculta que
        # se vuelve a mostrar tiene que aparecer con el mismo encuadre.
        cands = list(self._vis.values())
        todos = [self.vis_voi] + cands
        if activar:
            for v in cands:
                v.camera = self.vis_voi.camera
            for v in todos:
                otros = [o for o in todos if o is not v]
                obs = v.iren.add_observer(
                    "InteractionEvent",
                    lambda *a, os_=otros: [o.render() for o in os_])
                self._obs.append((v, obs))
            for v in todos:
                v.render()
        else:
            for v, obs in self._obs:
                try:
                    v.iren.remove_observer(obs)
                except Exception:
                    pass
            self._obs = []
            for v in cands:
                v.camera = pv.Camera()
                v.reset_camera()
        self.statusBar().showMessage(
            _("Camaras enlazadas: las dos vistas comparten encuadre.")
            if activar else _("Camaras independientes."))

    def _vista_lista(self, r):
        partes = []
        for f, (BW, info) in r.items():
            self._activar(f)
            self.BW_vista = BW
            self._vm_comparado.pop(f, None)
            # La estructura cambio: todos los campos calculados sobre la
            # anterior caducan. Dejarlos pegados no daria error —se muestrean
            # por posicion— sino algo peor: pintaria sobre la geometria nueva
            # la respuesta de la vieja, que es indistinguible a la vista de un
            # resultado valido.
            for cache in (self._esp, self._fe, self._vm, self._desp):
                cache.pop("spin", None)
            # ...y tambien la fabrica: la flecha del MIL describe la estructura
            # medida, no la que se acaba de generar. Se oculta hasta volver a
            # medir en vez de quedarse pegada.
            self._m_spin_vigente = False
            spacing = self._spacing(BW.shape[0])
            self.malla_spin = malla_de_mascara(BW, spacing)
            self._redibujar_spin(reset=True)
            partes.append(
                _("Vista a {n}^3 en {t} s — densidad pedida {ped}%, obtenida "
                  "{obt}%   [{esq}]").format(
                    n=BW.shape[0], t=f"{time.time()-self._t0:.1f}",
                    ped=f"{100*info['rho_objetivo']:.0f}",
                    obt=f"{100*info['rho_obtenida']:.1f}",
                    esq=info["esquema"]))
            # Region degenerada de los conos (`grf.region_degenerada`): mover
            # los angulos ahi no cambia la estructura, y sin decirlo parece
            # que el deslizador no funciona.
            if info.get("thetas_degenerados"):
                partes.append(_("thetas en la region degenerada: la "
                                "estructura es la isotropa (90, 90, 90) y "
                                "mover estos angulos no cambia nada"))
            elif (info.get("prob_realizacion_isotropa") or 0.0) >= 0.5:
                partes.append(_("thetas junto a la region degenerada: {p}% "
                                "de probabilidad de que esta realizacion sea "
                                "exactamente la isotropa").format(
                    p=f"{100 * info['prob_realizacion_isotropa']:.0f}"))
        self._activar(self._familias()[0])
        self._ocupado(False)
        self.statusBar().showMessage("   |   ".join(partes))

    def cargar_voi(self, ruta=None):
        if not ruta:
            ini = str(VOIDIR if VOIDIR.exists() else DATOS)
            # Los dos formatos, no solo VTK: los VOIs equinos son .vtk pero los
            # porcinos son .mat, y `spinpy.lote` ya procesaba ambos mientras la
            # interfaz solo sabia abrir uno.
            ruta, _x = QtWidgets.QFileDialog.getOpenFileName(
                self, _("Abrir VOI"), ini,
                _("VOI (*.vtk *.mat);;VTK legacy (*.vtk);;MATLAB (*.mat);;"
                  "Todos (*)"))
        if not ruta:
            return
        f = str(ruta)
        try:
            VOI, spacing = leer_voi(f)
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self, _("No se pudo leer el VOI"), str(e))
            return

        self._instalar_voi(VOI, spacing, Path(f).name, f)

    def cargar_pila_tiff(self):
        """Carga un volumen directamente desde las imagenes del micro-CT.

        Entra por el mismo `_instalar_voi` que un .vtk: en cuanto hay mascara y
        spacing en mm, todo lo de aguas abajo —morfometria, ajuste,
        homogeneizacion, von Mises, exportacion— funciona sin saber de donde
        vino.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        dlg = DialogoPilaTiff(self, str(VOIDIR if VOIDIR.exists() else DATOS))
        if dlg.exec_() != QtWidgets.QDialog.Accepted:
            return
        v = dlg.valores()

        from spinpy.voi import leer_pila_tiff

        prog = QtWidgets.QProgressDialog(
            _("Leyendo rebanadas…"), None, 0, 100, self)
        prog.setWindowTitle(_("Cargar pila TIFF"))
        prog.setWindowModality(QtCore.Qt.WindowModal)
        prog.setMinimumDuration(0)
        prog.setValue(0)

        def tarea(informar=None):
            return leer_pila_tiff(
                v["origen"], patron=v["patron"], umbral=v["umbral"],
                tam_voxel=v["tam_voxel"], unidad=v["unidad"],
                clases_otsu=v["clases_otsu"],
                progreso=(lambda i, n, et: informar(i, n, et)) if informar else None)

        hilo = Trabajador(lambda: tarea(hilo.informar))
        self.hilo = hilo

        def avance(i, n, etapa):
            prog.setValue(int(100 * i / max(n, 1)))
            prog.setLabelText(_("Leyendo rebanada {i} de {n}…").format(i=i, n=n))

        def listo(r):
            prog.close()
            BW, spacing, info = r
            self._instalar_pila(BW, spacing, info, v)

        def fallo(tb):
            prog.close()
            # Los fallos de carga son ESPERABLES y accionables: falta la
            # escala, la carpeta no tiene una serie clara, el patron no
            # selecciona nada. Enseñar el traceback de Python en esos casos
            # esconde la frase util debajo de quince lineas de ruido. El
            # traceback completo sigue disponible en «Ver detalle».
            ultima = tb.strip().splitlines()[-1] if tb.strip() else ""
            for marca in ("ValueError:", "FileNotFoundError:"):
                if marca in ultima:
                    m = QtWidgets.QMessageBox(self)
                    m.setIcon(QtWidgets.QMessageBox.Warning)
                    m.setWindowTitle(_("No se pudo cargar la pila"))
                    m.setText(ultima.split(marca, 1)[1].strip())
                    m.setDetailedText(tb)
                    m.exec_()
                    return
            self._error(tb)

        hilo.avance.connect(avance)
        hilo.listo.connect(listo)
        hilo.fallo.connect(fallo)
        hilo.start()

    def recortar_voi(self):
        """Orienta el volumen cargado por PCA y abre el recorte del cubo."""
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None:
            return

        from spinpy.voi import marco_pca, centros_por_tercios, extraer_cubo

        mask, spc = self.VOI, self.VOI_spacing
        # Lado por defecto: 5 mm es el del banco de este proyecto, pero si el
        # volumen cargado es mas pequeño que eso hay que bajarlo o el cubo
        # saldria casi entero fuera.
        menor = float(min(np.array(mask.shape) * np.asarray(spc).ravel()))
        lado = min(5.0, 0.6 * menor)

        prog = QtWidgets.QProgressDialog(
            _("Orientando el hueso por PCA…"), None, 0, 0, self)
        prog.setWindowTitle(_("Recortar VOI cubico"))
        prog.setWindowModality(QtCore.Qt.WindowModal)
        prog.setMinimumDuration(0)

        from spinpy import figura_metodo

        def tarea():
            ejes, centro, proy = marco_pca(mask, spc)
            centros = centros_por_tercios(proy, 3)
            tercios = []
            for i, c in enumerate(centros):
                cubo, fuera = extraer_cubo(mask, spc, c, ejes, centro, lado)
                tercios.append({"indice": i, "cubo": cubo, "centro_pca": c,
                                "fuera": fuera, "bvtv": float(cubo.mean()),
                                "lado_vox": int(cubo.shape[0])})
            # Pila reducida para la figura 0: aqui, en el hilo de trabajo,
            # porque tras el recorte la pila completa se suelta.
            reducido = figura_metodo.reducir(mask)
            return ({"ejes": ejes, "centro": centro, "proy": proy}, tercios,
                    reducido)

        def listo(r):
            prog.close()
            marco, tercios, reducido = r
            dlg = DialogoRecorte(self, mask, spc, marco, tercios, lado)
            if dlg.exec_() != QtWidgets.QDialog.Accepted:
                return
            e = dlg.elegido
            nota = ""
            if e["fuera"] > 0.005:
                nota = _("{p:.1f} % del cubo cae fuera del volumen: el BV/TV "
                         "esta subestimado.").format(p=100 * e["fuera"])
            nombre = f"{self.VOI_nombre} · {dlg.etiqueta} {e['lado_vox']}³"
            lado_elegido = float(dlg.sp_lado.value())
            pila = self._pila_param
            try:
                ctx = figura_metodo.contexto_pila(
                    mask, spc, marco["ejes"], marco["centro"],
                    e["centro_pca"], lado_elegido,
                    origen=(pila or {}).get("origen"),
                    patron=(pila or {}).get("patron") or "*.tif",
                    reducido=reducido)
            except Exception:
                ctx = None
            rec = (figura_metodo.registro_recorte(
                pila, marco["ejes"], marco["centro"], e["centro_pca"],
                lado_elegido) if pila else None)
            self._instalar_voi(e["cubo"], np.asarray(spc).ravel(), nombre,
                               self.VOI_ruta, nota=nota)
            self.VOI_contexto, self.VOI_recorte = ctx, rec
            self.statusBar().showMessage(
                _("VOI recortado: {et}, {lv}³ vox, BV/TV {b:.4f}. El volumen "
                  "completo ya no esta cargado.").format(
                      et=dlg.etiqueta, lv=e["lado_vox"], b=e["bvtv"]))

        def fallo(tb):
            prog.close()
            self._error(tb)

        hilo = Trabajador(tarea)
        self.hilo = hilo
        hilo.listo.connect(listo)
        hilo.fallo.connect(fallo)
        hilo.start()

    def _instalar_pila(self, BW, spacing, info, v):
        """Avisa de lo que hay que saber y deja la pila como VOI activo."""
        n_vox = int(np.prod(BW.shape))
        notas = []
        if info["umbral_metodo"].startswith("Otsu"):
            notas.append(_("Umbral por {m} = {u:.1f}. La segmentacion es la "
                           "mayor fuente de incertidumbre: comprueba el "
                           "valor.").format(m=info["umbral_metodo"],
                                            u=info["umbral"]))

        # Apartar archivos NO es inocuo —quitar rebanadas cambia la morfometria
        # igual que reescalarlas—, asi que se dice cuantos y cuales, no se
        # entierra en el log.
        apartados = info.get("apartados") or []
        if apartados:
            nombres = ", ".join(a["nombre"] for a in apartados[:6])
            if len(apartados) > 6:
                nombres += f", … (+{len(apartados) - 6})"
            notas.append(_("{n} archivo(s) apartados por no ser rebanadas de "
                           "la serie: {lista}").format(n=len(apartados),
                                                       lista=nombres))
            QtWidgets.QMessageBox.information(
                self, _("Archivos apartados"),
                _("Se cargaron {ok} rebanadas de {tot}x{ancho} px.\n\n"
                  "Estos {n} archivos se apartaron por tener otra forma — son "
                  "las proyecciones, previsualizaciones y hojas que deja el "
                  "escaner junto a la reconstruccion:\n\n{lista}").format(
                      ok=info["n_rebanadas"],
                      tot=info["forma_rebanada"][0],
                      ancho=info["forma_rebanada"][1],
                      n=len(apartados),
                      lista="\n".join(
                          f"  · {a['nombre']}   {a['forma'][0]}x{a['forma'][1]}"
                          if a["forma"] else f"  · {a['nombre']}"
                          for a in apartados)))
        nota = "  ".join(notas)

        # El coste de la homogeneizacion crece con el CUBO del lado. Una pila
        # entera de escaner (740x656x438 = 2.1e8 voxeles) no es un VOI: no se
        # bloquea la carga —ver y medir la pila completa es legitimo— pero se
        # dice con todas las letras antes de que alguien pulse «Analisis
        # mecanico» y se pregunte por que no vuelve.
        if n_vox > 20_000_000:
            QtWidgets.QMessageBox.warning(
                self, _("Volumen muy grande"),
                _("La pila tiene {nv:,} voxeles ({f}).\n\n"
                  "Se puede visualizar y medir la morfometria, pero el analisis "
                  "mecanico y el ajuste NO son viables a este tamano: el coste "
                  "va con el cubo del lado.\n\n"
                  "Para eso hay que recortar un VOI cubico "
                  "(spinpy.voi.vois_por_tercios / extraer_cubo).").format(
                      nv=n_vox,
                      f="x".join(str(s) for s in BW.shape)))

        nombre = Path(info["origen"]).name or info["origen"]
        self._instalar_voi(BW, spacing, nombre, info["origen"], nota=nota)
        self._pila_param = {"origen": info["origen"], "patron": v["patron"],
                            "umbral": info["umbral"],
                            "tam_voxel_mm": info["tam_voxel_mm"]}
        self.statusBar().showMessage(
            _("Pila cargada: {n} rebanadas · {mm:.6f} mm/vox ({orig}) · "
              "umbral {u} [{met}]").format(
                n=info["n_rebanadas"], mm=info["tam_voxel_mm"],
                orig=info["tam_voxel_origen"], u=f"{info['umbral']:g}",
                met=info["umbral_metodo"]))

    def _instalar_voi(self, VOI, spacing, nombre, ruta, nota=""):
        """Deja un VOI ya leido como VOI activo y rehace todo lo que depende.

        Se separo de `cargar_voi` para que la carga desde una pila TIFF entre
        exactamente por el mismo sitio: en cuanto hay `(BW, spacing)` en mm, el
        origen deja de importar y no debe haber dos caminos que puedan
        divergir.
        """
        self.VOI, self.VOI_spacing = VOI, spacing
        self.VOI_nombre = nombre
        self.VOI_ruta = ruta
        # Un VOI nuevo no hereda el origen del anterior: quien instala desde
        # una pila o un recorte lo vuelve a poner DESPUES de esta llamada.
        self._pila_param = None
        self.VOI_contexto = None
        self.VOI_recorte = None
        # El recorte tiene sentido sobre cualquier volumen cargado, no solo
        # sobre una pila: un .vtk demasiado grande tambien se puede recortar.
        try:
            self.b_recorte.setEnabled(True)
        except AttributeError:
            pass
        lado = VOI.shape[0] * spacing[0]
        # Aviso de densidad (Estudio_Familias, 4.6): a BV/TV >= 0.83 ninguna
        # familia sale bicontinua y el VOI ya no es trabecula. Se avisa al
        # cargar, con el BV/TV crudo, porque es cuando se decide si seguir.
        dens = voi_no_trabecular(float(VOI.mean()))
        if dens["aviso"]:
            txt = _("BV/TV {bv} ≥ {umbral}: a esta densidad el VOI no es "
                    "trabecular (poros aislados). Un ajuste con error bajo no "
                    "significa que se este imitando trabecula.").format(
                bv=f"{dens['BVTV']:.3f}", umbral=f"{dens['umbral']:.2f}")
            nota = f"{nota}<br>{txt}" if nota else txt
        tam = " x ".join(f"{s * spacing[i]:.3f}" for i, s in enumerate(VOI.shape))
        self.lab_voi.setText(
            f"<b>{self.VOI_nombre}</b><br>{VOI.shape[0]}x{VOI.shape[1]}x"
            f"{VOI.shape[2]} vox · {spacing[0]:.6f} mm/vox<br>"
            f"{tam} mm · BV/TV crudo {VOI.mean():.4f}"
            + (f"<br><span style='color:#a33;'>{nota}</span>" if nota else ""))

        self.malla_voi = malla_de_mascara(VOI, spacing)

        # El candidato se escala al tamano fisico del VOI (correccion F1), asi
        # que al cargar uno CAMBIA la escala del spinodoide ya generado: antes
        # media 1 mm de lado y ahora debe medir lo que el VOI. Si no se
        # reconstruye aqui, los dos paneles quedan a escalas distintas y, con
        # las camaras enlazadas, uno aparece diminuto mientras del otro se ve
        # el interior. No hace falta regenerar el campo —la mascara no cambia—,
        # solo rehacer la malla con el nuevo spacing.
        # Las dos familias, se vean o no: una familia oculta que se vuelva a
        # mostrar tiene que aparecer ya a la escala del VOI nuevo.
        previa = self._fam
        for f in FAMILIAS:
            self._activar(f)
            if self.BW_vista is not None:
                self.malla_spin = malla_de_mascara(
                    self.BW_vista, self._spacing(self.BW_vista.shape[0]))
                self._esp.pop("spin", None)      # el campo de espesor caduca
        self._activar(previa)
        for cache in (self._esp, self._fe, self._vm):
            cache.pop("voi", None)
        self._vm_comparado = {}
        self._clim = None
        if self.cmb_color.currentIndex() != 0:
            self.cmb_color.setCurrentIndex(0)   # se recalcula si se vuelve a pedir

        self._redibujar_voi(reset=True)
        self.vis_voi.reset_camera()

        def _panel():
            self._redibujar_spin(reset=True)
            self.vis_spin.reset_camera()
        self._para_cada(_panel)

        self.m_voi = None
        self.statusBar().showMessage(
            _("VOI cargado. Pulsa Medir para la morfometria (el candidato se "
              "escalara a {lado} mm de lado).").format(lado=f"{lado:.3f}"))

    def medir(self, sobre_mascara=False):
        """Morfometria de cada familia activa y, si hace falta, del VOI.

        `sobre_mascara=True` mide la mascara que ya tiene cada familia en vez
        de regenerarla a la resolucion de medida. Es lo que usa el informe
        automatico despues de ajustar: los deslizadores redondean densidad,
        numero de onda y angulos a su paso, asi que regenerar desde ellos
        mediria una estructura vecina de la ganadora, no la ganadora, y su
        huella no cuadraria con la del paquete de reproduccion.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        res = int(self.sl["resm"].valor())
        fams = self._familias()
        ps, mascaras = {}, {}
        for f in fams:
            self._activar(f)
            if sobre_mascara and self.BW_vista is not None:
                mascaras[f] = self.BW_vista
            else:
                ps[f] = self._params(res)
        self._activar(fams[0])
        if mascaras:
            res = int(next(iter(mascaras.values())).shape[0])
        esquemas = {f: ((self._res.get("ajuste" if f == "spinodoide"
                                       else "ajuste_dual") or {})
                        .get("parametros") or {}).get("esquema")
                    for f in mascaras}
        sps = {f: self._spacing(BW.shape[0]) for f, BW in mascaras.items()}
        sp_spin = self._spacing(res)
        VOI, sp_voi = self.VOI, self.VOI_spacing
        extra = self.chk_extra.isChecked()
        poro = self.chk_poro.isChecked()
        ef = self.chk_ef.isChecked()
        # Si se pide Conn.D/SMI y el VOI ya se midio SIN ellos, hay que
        # volver a medirlo: si no, la columna del VOI saldria vacia junto a la
        # del spinodoide llena, que se lee como "el VOI no tiene conectividad"
        # en vez de "no se calculo".
        # El VOI se vuelve a medir si nunca se midio, si ahora se piden las
        # metricas extra y no las tiene, o si el MODO cambio: una tabla con
        # el VOI en voxeles y el candidato en malla mezcla dos definiciones de
        # BS que difieren un 9 % (correccion C4, otra vez).
        modo = self._modo_medida()
        medir_voi = VOI is not None and (
            self.m_voi is None
            or (extra and "ConnD" not in self.m_voi)
            or (poro and not np.isfinite(self.m_voi.get("PoDm", np.nan)))
            or (ef and not np.isfinite(self.m_voi.get("EF", np.nan)))
            or self.m_voi.get("modo", "voxel") != modo)
        self._ocupado(True, f"Midiendo a {res}^3 [{modo}]"
                            + (_(" (y el VOI)…") if medir_voi else "…"))
        self._t0 = time.time()

        def medir_una(BW, sp):
            if modo == "malla":
                m = morfometria_malla(BW, sp)
                if extra or poro or ef:
                    # Conn.D, SMI, Po.Dm y EF son de VOXELES en los dos modos:
                    # son topologicos o de transformada de distancia, y no
                    # tienen version de malla. La superficie interna/externa NO
                    # se copia: el modo malla calcula la suya con su propia BS.
                    mx = morfometria(BW, sp, do_mil=False, extra=extra or ef,
                                     do_poro=poro, do_ef=ef)
                    for k in ("ConnD", "Conn", "euler", "SMI", "dBSdr",
                              "PoDm", "PoDm_mediana", "PoDm_sd", "PoDm_p05",
                              "PoDm_p95", "PoDm_n", "PoDm_frac_ventana",
                              "PoDm_interior", "PoDm_interior_mediana",
                              "PoDm_interior_n", "EF", "EF_mediana", "EF_sd",
                              "EF_p05", "EF_p95", "EF_frac_placa",
                              "EF_frac_barra", "EF_relleno",
                              "EF_n_elipsoides"):
                        if k in mx:
                            m[k] = mx[k]
                return m
            m = morfometria(BW, sp, extra=extra or ef, do_poro=poro, do_ef=ef)
            m["modo"] = "voxel"
            return m

        def tarea():
            ms, infos = {}, {}
            for f in fams:
                if f in mascaras:
                    BW = mascaras[f]
                    info = {"esquema": esquemas.get(f)}
                    sp = sps[f]
                else:
                    BW, _x, info = Visor._generar(ps[f])
                    sp = sp_spin
                ms[f] = medir_una(BW, sp)
                infos[f] = info
            m_v = medir_una(VOI, sp_voi) if medir_voi else None
            return ms, m_v, infos, res, bool(mascaras)

        self.hilo = Trabajador(tarea)
        self.hilo.listo.connect(self._medida_lista)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    def _medida_lista(self, r):
        ms, m_v, infos, res, sobre_mascara = r
        for f, m_f in ms.items():
            self._activar(f)
            self.m_spin = m_f
            self._m_spin_vigente = True
        f0 = self._familias()[0]
        self._activar(f0)
        m_s, info = ms[f0], infos[f0]
        if m_v is not None:
            self.m_voi = m_v
        self._res["morfometria"] = {
            "resolucion": int(res),
            # Si se midio la mascara del ajuste o una regenerada desde los
            # deslizadores: no son la misma estructura.
            "origen": "mascara" if sobre_mascara else "deslizadores",
            "esquema": info.get("esquema"),
            "familias": list(ms),
            "semilla": int(self.spin_semilla.value()),
            "extra": bool(self.chk_extra.isChecked()),
            "ef": bool(self.chk_ef.isChecked()),
            "modo": self._modo_medida(),
        }
        self._rellenar_tabla()
        self._ocupado(False)
        self.statusBar().showMessage(
            _("Medido en {t} s   [superficie: {modo}{extra}]   "
              "[{esq}, semilla {sem}]").format(
                t=f"{time.time()-self._t0:.1f}",
                modo=_(self._modo_medida()),
                extra=(_(", tapas excluidas {tap} mm2, perdida vs superficie "
                         "cruda {per} %").format(
                           tap=f"{m_s['BS_tapas']:.1f}",
                           per=f"{m_s['perdida_vs_mc_pct']:+.1f}")
                       if m_s.get("modo") == "malla"
                       and np.isfinite(m_s.get("perdida_vs_mc_pct", np.nan))
                       else ""),
                esq=info["esquema"], sem=self.spin_semilla.value()))

    def _reconstruir_tabla(self):
        """Columnas segun las familias activas: el VOI y, por cada familia, su
        valor y su diferencia relativa frente al VOI."""
        fams = self._familias()
        n = 2 + 2 * len(fams)
        self.tabla.setColumnCount(n)
        cab = [_("Metrica"), _("VOI")]
        for f in fams:
            cab += [self._etq_fam(f), _("Dif. relativa")]
        self.tabla.setHorizontalHeaderLabels(cab)
        for i in range(len(TABLA)):
            for j in range(1, n):
                if self.tabla.item(i, j) is None:
                    self.tabla.setItem(i, j, QtWidgets.QTableWidgetItem("—"))

    def _rellenar_tabla(self):
        fams = self._familias()
        if self.tabla.columnCount() != 2 + 2 * len(fams):
            self._reconstruir_tabla()
        ms = {f: self._de(f, "m") for f in fams}
        normal = QtGui.QBrush(QtGui.QColor("#000000"))
        for i, (clave, _x, fmt, _y) in enumerate(TABLA):
            vv = self.m_voi.get(clave, np.nan) if self.m_voi else np.nan
            vv = float(vv) if np.isscalar(vv) else np.nan
            it = self.tabla.item(i, 1)
            it.setText(fmt.format(vv) if np.isfinite(vv) else "—")
            it.setForeground(normal)
            it.setToolTip("")

            for k, f in enumerate(fams):
                m = ms[f]
                vs = m.get(clave, np.nan) if m else np.nan
                vs = float(vs) if np.isscalar(vs) else np.nan
                iv = self.tabla.item(i, 2 + 2 * k)
                iv.setText(fmt.format(vs) if np.isfinite(vs) else "—")
                iv.setForeground(normal)
                iv.setToolTip("")

                idf = self.tabla.item(i, 3 + 2 * k)
                if np.isfinite(vv) and np.isfinite(vs) and vv != 0:
                    d = 100.0 * (vs - vv) / abs(vv)
                    idf.setText(f"{d:+.2f} %")
                    # Solo colorea el orden de magnitud del desajuste; no es
                    # un criterio de aceptacion. El error que minimizan los
                    # optimizadores es una suma ponderada, no metrica a metrica.
                    a = abs(d)
                    idf.setForeground(QtGui.QBrush(QtGui.QColor(
                        "#1a7f37" if a < 5 else "#9a6700" if a < 15
                        else "#b62324")))
                else:
                    idf.setText("—")
                    idf.setForeground(QtGui.QBrush(QtGui.QColor("#888")))

                # DA por debajo del suelo de ruido del MIL: no es una
                # diferencia real
                if clave == "DA" and m and m.get("MIL_clamped"):
                    iv.setToolTip(_(
                        "MIL acotado: estructura laminar, el DA es una COTA "
                        "INFERIOR."))

            n = self.tabla.columnCount()
            if clave == "SMI":
                for j in [1] + [2 + 2 * k for k in range(len(fams))]:
                    self.tabla.item(i, j).setToolTip(AVISO_SMI)
            sin = ((clave in ("ConnD", "SMI") and not self.chk_extra.isChecked())
                   or (clave.startswith("EF") and not self.chk_ef.isChecked()))
            if sin:
                for j in range(1, n):
                    self.tabla.item(i, j).setText(_("(sin calcular)"))
                    self.tabla.item(i, j).setForeground(
                        QtGui.QBrush(QtGui.QColor("#aaa")))

    # -- ajuste al VOI -----------------------------------------------------

    def ajustar(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None:
            QtWidgets.QMessageBox.information(
                self, _("Falta el VOI"),
                _("El ajuste necesita un VOI de referencia. Cargalo "
                  "primero."))
            return

        VOI, sp = self.VOI, self.VOI_spacing
        esquema = "rechazo" if self.cmb_esq.currentIndex() == 0 else "equitativo"
        nw = int(self.sl["nw"].valor())
        semilla = int(self.spin_semilla.value())
        m_voi = self.m_voi
        peso_mec = (float(self.spin_peso_mec.value())
                    if self.chk_mec.isChecked() else 0.0)
        res_mec = int(self.spin_res_mec.value())

        fams = self._familias()
        N_EVAL = {"spinodoide": 53, "dual-lattice": 54}
        total = sum(N_EVAL[f] for f in fams)
        self._ocupado(True, _("Ajustando…")
                      + (_(" (con desempate mecanico a {r}³)").format(r=res_mec)
                         if peso_mec else ""),
                      determinada=True, total=total)
        self._t0 = time.time()

        hilo = Trabajador(lambda: None)   # placeholder, se sustituye abajo

        def tarea():
            # En serie y en el mismo hilo: con las dos familias, la barra de
            # progreso avanza de forma monotona sobre el total en vez de
            # saltar entre dos ajustes a la vez.
            out, base = {}, 0
            for f in fams:
                def pr(i, n, etapa, _b=base, _f=f):
                    hilo.informar(_b + int(i), total, f"{_f} · {etapa}")
                if f == "dual-lattice":
                    out[f] = ajustar_dual_lattice(
                        VOI, sp, modo="completo", m_voi=m_voi, seed=semilla,
                        peso_mecanico=peso_mec, res_mec=res_mec,
                        E_s=E_S_PA, nu_s=NU_S, progreso=pr)
                else:
                    out[f] = ajustar_spinodoide(
                        VOI, sp, modo="completo", m_voi=m_voi, esquema=esquema,
                        num_waves=nw, seed=semilla, precision="f32",
                        peso_mecanico=peso_mec, res_mec=res_mec,
                        E_s=E_S_PA, nu_s=NU_S, progreso=pr)
                base += N_EVAL[f]
            return out

        hilo._fn = tarea
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._ajustes_listos)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _ajustes_listos(self, out):
        """Coloca el resultado de cada familia ajustada, en orden."""
        for r in out.values():
            self._ajuste_listo(r)
        self._activar(self._familias()[0])

    def _ajuste_listo(self, r):
        p = r["parametros"]
        fam = r.get("familia", p.get("familia", "spinodoide"))
        self._activar(fam)
        # Se colocan los deslizadores SIN disparar una regeneracion: la mascara
        # buena ya viene calculada dentro del resultado y volver a generarla
        # daria otra realizacion, con otras metricas que las de la tabla.
        self._aplicando = True
        try:
            if fam == "dual-lattice":
                self.sl["d_dens"].fijar(p["densidad"] * 100.0)
                self.sl["d_celdas"].fijar(p["celdas"])
                for k, clave in zip(range(3), ("d_ex", "d_ey", "d_ez")):
                    self.sl[clave].fijar(p["estiramiento"][k])
                self.sl["d_irr"].fijar(p.get("irregularidad", 0.5))
            else:
                self.sl["dens"].sld.setValue(int(round(p["densidad"] * 100)))
                self.sl["wave"].sld.setValue(
                    int(round(p["wave_number_pi"] * 10)))
                self.sl["nw"].sld.setValue(int(p["num_waves"]))
                for k, clave in zip(range(3), ("thx", "thy", "thz")):
                    self.sl[clave].sld.setValue(int(round(p["thetas"][k])))
            self.R_forzada = np.asarray(p["R"], dtype=float)
        finally:
            self._aplicando = False

        rx, ry, rz = p["euler_deg"]
        self.lab_R.setText(_(
            "Orientacion tomada del eje principal del VOI (Euler {e}°). "
            "Manda sobre los deslizadores de rotacion hasta que muevas uno."
        ).format(e=f"{rx:.0f}/{ry:.0f}/{rz:.0f}"))

        self.m_voi = r["metricas_voi"]
        self.m_spin = r["metricas_spin"]
        self._m_spin_vigente = True
        self.BW_vista = r["mascara"]
        # El candidato cambio: su mapa de von Mises ya no le corresponde.
        self._vm_comparado.pop(fam, None)
        self.malla_spin = malla_de_mascara(
            r["mascara"], self._spacing(r["mascara"].shape[0]))
        self._redibujar_spin(reset=True)
        self._rellenar_tabla()

        self._ocupado(False)
        self._res[self._clave_res("ajuste")] = {
            "procedencia": r.get("procedencia"),
            "parametros": r["parametros"], "error": r["error"],
            "diagnostico": r["diagnostico"], "k2": r["k2"],
            "mecanico": r.get("mecanico", {}),
            "objetivo": r.get("objetivo"),
            "incertidumbre": r.get("incertidumbre"),
            "seleccion": r.get("seleccion"),
            "alineacion": r.get("alineacion"),
            "n_evaluaciones": r["n_evaluaciones"], "tiempo_s": r["tiempo_s"],
        }

        k2 = r["k2"]
        nota = ""
        if k2.get("activo") and not k2.get("aplicado"):
            nota = "  ·  " + _("desempate por orientacion activo pero no "
                                "aplicado ({n} candidato/s dentro del 5%)"
                                ).format(n=k2.get("n_elegibles", 0))
        if voi_no_trabecular(r["metricas_voi"])["aviso"]:
            nota += "  ·  " + _("VOI no trabecular (BV/TV {bv})").format(
                bv=f"{r['metricas_voi']['BVTV']:.3f}")
        des = desalineacion(r["metricas_voi"], r["metricas_spin"])
        if des["aviso"]:
            nota += "  ·  " + _("eje a {ang}° del VOI: revisalo antes de "
                                "ensayar").format(ang=f"{des['angulo_deg']:.0f}")
        alin = r.get("alineacion") or {}
        if alin.get("eje_comparado") and not alin.get("alineada"):
            nota += "  ·  " + _("orientacion impuesta sin llegar a la "
                                "tolerancia: eje a {ang}°").format(
                ang=f"{alin.get('angulo_despues_deg', float('nan')):.0f}")
        inc = r.get("incertidumbre") or {}
        if inc.get("K"):
            nota += "  ·  " + _("error con {K} semillas nuevas: {m} ± {s} "
                                "(suelo {f})").format(
                K=inc["K"], m=f"{inc['error']['media']:.5f}",
                s=f"{inc['error']['sd']:.5f}",
                f=f"{inc['suelo_autoconsistente']['media']:.5f}")
        mec = r.get("mecanico") or {}
        if mec.get("traza"):
            nota += ("  ·  " + _("desempate mecanico:") + " "
                     + (_("REORDENO") if mec["reordeno"]
                        else _("confirmo la morfometria")))
        self.statusBar().showMessage(
            _("Ajuste terminado: error={err}, {n} evaluaciones en {t} s"
              ).format(err=f"{r['error']:.5f}", n=r["n_evaluaciones"],
                       t=f"{r['tiempo_s']:.0f}") + nota)

        if mec.get("traza"):
            self._informe_mecanico(mec)

    def ajustar_todos(self):
        """Ajusta con todos los metodos, compara y deja elegir."""
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None:
            QtWidgets.QMessageBox.information(
                self, _("Falta el VOI"),
                _("El ajuste necesita un VOI de referencia. Cargalo "
                  "primero."))
            return
        from dialogo_metodos import DialogoMetodos

        m_voi = self.m_voi
        if m_voi is None:
            self.statusBar().showMessage(_("Midiendo el VOI antes de ajustar…"))
            QtWidgets.QApplication.processEvents()
            m_voi = morfometria(self.VOI, self.VOI_spacing)
            self.m_voi = m_voi

        # Una ventana por familia activa, una detras de otra. Cada familia
        # tiene sus propios metodos, y compararlos en la misma ventana
        # mezclaria errores de dos espacios de busqueda distintos.
        for f in self._familias():
            self._activar(f)
            d = DialogoMetodos(self.VOI, self.VOI_spacing, m_voi, self,
                               semilla=int(self.spin_semilla.value()),
                               num_waves=int(self.sl["nw"].valor()),
                               registro=REGISTROS[f])
            if d.exec_() != QtWidgets.QDialog.Accepted or d.elegido is None:
                self.statusBar().showMessage(_("Ajuste comparado cancelado."))
                continue

            # Se reutiliza la misma ruta que el ajuste sencillo: coloca los
            # deslizadores sin regenerar, fija la R del ajuste y rellena la
            # tabla con las metricas del propio resultado, no con unas nuevas.
            self._ajuste_listo(d.elegido)
            self._res[self._clave_res("ajuste_comparado")] = {
                "elegido": d.elegido["metodo"],
                "comparativa": [
                    {k: v for k, v in fl.items() if np.isscalar(v) or v is None}
                    for fl in comparar_resultados(
                        [r for r in d.resultados.values()], m_voi)],
            }
            if d.replicas:
                self._guardar_replicas(d.replicas, d.elegido)
        self._activar(self._familias()[0])

    def _guardar_replicas(self, rep, elegido):
        """Vuelca las replicas a CSV y resume la dispersion lograda."""
        SAL = DATOS / "resultados"
        SAL.mkdir(parents=True, exist_ok=True)
        sello = time.strftime("%Y%m%d_%H%M%S")
        dual = (elegido.get("parametros") or {}).get("familia") == "dual-lattice"
        prefijo = "replicas_dual" if dual else "replicas"
        destino = SAL / f"{prefijo}_{elegido['metodo']}_{sello}.csv"
        propios = (["replica", "semilla", "densidad", "celdas", "rho_obtenida"]
                   if dual else
                   ["replica", "semilla", "densidad", "wave_number_pi",
                    "wave_number_rad", "rho_obtenida"])
        # Procedencia en cada fila (`spinpy.procedencia`), al final para no
        # mover las columnas que ya leen las hojas hechas.
        pc = procedencia.columnas(
            procedencia.desde_parametros(elegido["parametros"]))
        cols_proc = ["spinpy_version", "formato_version", "familia", "esquema"]
        cols = propios + [c for c in COMPARADAS_MET
                          if c in (rep["filas"][0] if rep["filas"] else {})]
        # La ultima columna es la del vector que no cabe en una celda: los
        # thetas del spinodoide o el estiramiento del dual-lattice.
        vector, fmt_v = (("estiramiento", "%.2f") if dual
                         else ("thetas", "%.1f"))
        with open(destino, "w", encoding="utf-8-sig") as f:
            f.write(";".join(cols) + f";{vector};" + ";".join(cols_proc)
                    + "\n")
            for fila in rep["filas"]:
                f.write(";".join(
                    ("%.10g" % fila[c]).replace(".", ",") if c in fila else ""
                    for c in cols)
                    + ";" + "/".join(fmt_v % t for t in fila[vector])
                    + ";" + ";".join(self._txt(pc[c]) for c in cols_proc)
                    + "\n")

        # Las mascaras NO se guardan: son cientos de MB y se regeneran exacto
        # desde la semilla y los parametros, que si estan en el CSV.
        t = [_("<b>{n} replicas generadas</b> — {como}<br>").format(
                 n=len(rep["filas"]),
                 como=(_("identicas (solo cambia la semilla)")
                       if rep["tecnicas"] else
                       _("variando {mets} con amplitud {amp} %").format(
                           mets=", ".join(rep["variar"]),
                           amp="%.0f" % (100 * rep["amplitud"])))),
             "<table cellpadding=4 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>",
             _("<tr><th>metrica</th><th>media</th><th>sd</th><th>CV</th>"
               "</tr>")]
        for k, d2 in rep["logrado"].items():
            t.append("<tr><td><b>%s</b></td><td>%.5g</td><td>%.3g</td>"
                     "<td>%.2f %%</td></tr>" % (k, d2["media"], d2["sd"],
                                                d2["cv_pct"]))
        t.append("</table>")
        for a in rep["avisos"]:
            t.append("<p style='color:#8a6d00;font-size:10px'>%s</p>" % a)
        t.append(
            "<p style='color:#b62324;font-size:11px'>"
            + _("<b>Estas {n} filas no son {n} especimenes.</b> Son replicas "
                "de UN ajuste: no anaden grados de libertad a ninguna "
                "comparacion biologica. Para eso esta el N efectivo del panel "
                "de lote.").format(n=len(rep["filas"]))
            + "</p>")
        t.append("<p style='color:#666;font-size:10px'>"
                 + _("Las mascaras no se guardan: se regeneran exactas desde "
                     "la semilla y los parametros, que si estan en el CSV."
                     "<br>-> {archivo}").format(archivo=destino.name)
                 + "</p>")

        self._res[self._clave_res("replicas")] = {"n": len(rep["filas"]),
                                 "variar": rep["variar"],
                                 "amplitud": rep["amplitud"],
                                 "logrado": rep["logrado"],
                                 "csv": destino.name}
        self.statusBar().showMessage(
            _("{n} replicas en {archivo} — recuerda que no son {n} "
              "especimenes").format(n=len(rep["filas"]),
                                    archivo=destino.name))
        d3 = QtWidgets.QMessageBox(self)
        d3.setWindowTitle(_("Replicas generadas"))
        d3.setTextFormat(QtCore.Qt.RichText)
        d3.setText("".join(t))
        d3.exec_()

    def _informe_mecanico(self, mec):
        """Que hizo la etapa D, con el orden antes y despues a la vista.

        Lo importante no es el numero final sino SI el desempate decidio algo.
        Si el ganador sigue siendo el mejor por morfometria, la mecanica
        confirmo y no aporto; si cambio, la morfometria por si sola habria
        elegido una estructura con la rigidez equivocada, y eso es un resultado
        que merece verse en vez de quedar enterrado en el JSON.
        """
        voi = mec.get("voi") or {}
        t = [_("<b>Etapa D — desempate mecanico</b><br>"),
             _("<p>VOI: E<sub>z</sub>/E<sub>s</sub> = <b>{ezes}</b> · "
               "E<sub>z</sub>/E<sub>x</sub> = <b>{ezex}</b> "
               "(malla {res}³, peso {peso})</p>").format(
                 ezes=f"{voi.get('Ez_rel', float('nan')):.5f}",
                 ezex=f"{voi.get('Ez_Ex', float('nan')):.4f}",
                 res=mec["res_mec"], peso=f"{mec['peso']:g}"),
             "<table cellpadding=5 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>",
             _("<tr><th>orden morfo.</th><th>error morfo.</th>"
               "<th>error + mecanica</th><th>E<sub>z</sub>/E<sub>s</sub></th>"
               "<th>desvio vs VOI</th></tr>")]
        mejor = min(range(len(mec["traza"])),
                    key=lambda i: mec["traza"][i]["err_con_mec"])
        for i, f in enumerate(mec["traza"]):
            marca = " style='background:#e8f5e9'" if i == mejor else ""
            ez = f.get("Ez_rel")
            dv = ("—" if not ez or not voi.get("Ez_rel")
                  else f"{100*(ez/voi['Ez_rel']-1):+.0f} %")
            t.append(
                f"<tr{marca}><td>{f['orden_morfo']}</td>"
                f"<td>{f['err_morfo']:.6f}</td>"
                f"<td><b>{f['err_con_mec']:.6f}</b></td>"
                f"<td>{('%.5f' % ez) if ez else _('no homogeneizo')}</td>"
                f"<td>{dv}</td></tr>")
        t.append("</table>")
        t.append(
            "<p><b>" + (_("Reordeno") if mec["reordeno"] else _("Confirmo"))
            + ":</b> "
            + (_("el mejor por morfometria NO es el mejor al sumar la "
                 "mecanica. Sin esta etapa el ajuste habria elegido una "
                 "estructura con la rigidez equivocada.")
               if mec["reordeno"] else
               _("el mejor por morfometria lo sigue siendo. La etapa no "
                 "decidio nada, que tambien es informacion: la morfometria "
                 "bastaba."))
            + "</p>")
        t.append(
            "<p style='color:#666;font-size:10px'>"
            + _("E<sub>z</sub>/E<sub>s</sub> depende mucho de la malla de "
                "homogeneizacion y <b>no es citable como rigidez</b>. Vale "
                "aqui porque VOI y candidatos se miden con la misma malla y "
                "el sesgo se cancela al comparar.<br>"
                "Solo se homogeneizaron los finalistas: si el optimo mecanico "
                "estaba en una region que la busqueda morfometrica descarto "
                "pronto, esta etapa no lo encuentra.")
            + "</p>")

        d = QtWidgets.QMessageBox(self)
        d.setWindowTitle(_("Desempate mecanico"))
        d.setTextFormat(QtCore.Qt.RichText)
        d.setText("".join(t))
        self._mostrar(d)

    # -- ensayo de compresion ----------------------------------------------

    def ensayo_fe(self):
        """Ensayo de compresion sobre cada familia activa, una tras otra."""
        self._en_cada_familia(self._ensayo_fe_una)

    def _ensayo_fe_una(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.BW_vista is None and self.VOI is None:
            return
        if not self._confirmar_orientacion():
            return
        n = int(self.spin_res_fe.value())
        apoyo = apoyo_de_combo(self.cmb_apoyo)
        BW_s = self.BW_vista
        sp_s = self._spacing(BW_s.shape[0]) if BW_s is not None else None
        VOI, sp_v = self.VOI, self.VOI_spacing

        # Indice 0 -> Z, 1 -> X, 2 -> Y, 3 -> los tres. El desplegable pone Z
        # primero porque es la direccion axial del VOI y el caso habitual.
        i = self.cmb_eje_fe.currentIndex()
        ejes = (0, 1, 2) if i == 3 else ((2,), (0,), (1,))[i]
        etq_ejes = _("los tres ejes") if i == 3 else EJES[ejes[0]]

        self._ocupado(True, _("Resolviendo la compresion a {n}³ en {eje}…"
                              ).format(n=n, eje=etq_ejes))
        self._t0 = time.time()

        def uno(BW, sp):
            bw, spr = remuestrear_bw(BW, sp, n)
            ejes_ok, primero = {}, None
            for e in ejes:
                r = ensayo_compresion_eje(bw, spr, eje=e, E_s=E_S_PA,
                                          nu_s=NU_S, apoyo=apoyo)
                if not r["ok"]:
                    ejes_ok[EJES[e]] = {"ok": False, "msg": r["msg"]}
                    continue
                p = criterio_pistoia(r)
                p["E_app"] = r["E_app"]
                p["n_elem"] = r["n_elem"]
                p["frac_portante"] = r.get("frac_portante", 1.0)
                p["solver"] = r.get("solver", "")
                p["residuo"] = r.get("residuo_rel", 0.0)
                # `residuo_rel` es la clave que lee `informe.comprobar` para
                # la puerta F10; con solo `residuo` el informe no la veia.
                p["residuo_rel"] = p["residuo"]
                ejes_ok[EJES[e]] = p
                if primero is None:
                    # Los campos que se pintan son los del PRIMER eje que
                    # resolvio. Superponer los tres no tendria sentido: cada
                    # uno corresponde a un caso de carga distinto, y mezclarlos
                    # daria un mapa que no es la respuesta a ninguna carga.
                    primero = (EJES[e], r["campo_eps_eff"], r["campo_vm"], spr)
            if primero is None:
                m = next((v["msg"] for v in ejes_ok.values() if v.get("msg")),
                         "El ensayo no se resolvio.")
                return {"ok": False, "msg": m}
            eje_pintado, campo_eps, campo_vm, spr = primero
            return {"ok": True, "ejes": ejes_ok, "eje_pintado": eje_pintado,
                    "campo": campo_eps, "campo_vm": campo_vm, "spacing": spr}

        def tarea():
            out = {}
            if VOI is not None:
                out["voi"] = uno(VOI, sp_v)
            if BW_s is not None:
                out["spin"] = uno(BW_s, sp_s)
            return out

        self.hilo = Trabajador(tarea)
        self.hilo.listo.connect(self._fe_listo)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    def _fe_listo(self, out):
        transcurrido = time.time() - self._t0
        self._ocupado(False)
        self._fe = {k: (v["campo"], v["spacing"]) for k, v in out.items()
                    if v.get("ok")}
        self._vm = {k: (v["campo_vm"], v["spacing"]) for k, v in out.items()
                    if v.get("ok")}
        # Sin los campos 3D: son cientos de megas y no caben en un CSV.
        self._res[self._clave_res("resistencia")] = {
            "resolucion": int(self.spin_res_fe.value()),
            "apoyo": apoyo_de_combo(self.cmb_apoyo),
            "E_s_Pa": E_S_PA, "nu_s": NU_S,
            "por_estructura": {
                k: {"ejes": {n: {kk: vv for kk, vv in p.items()
                                 if np.isscalar(vv)}
                             for n, p in v["ejes"].items()},
                    "eje_pintado": v["eje_pintado"]}
                for k, v in out.items() if v.get("ok")},
        }

        t = [_("<b>Ensayo de compresion</b> — criterio de Pistoia<br>"),
             "<table cellpadding=5 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>",
             "<tr><th></th><th>" + _("eje") + "</th><th>E<sub>app</sub></th>"
             "<th>E<sub>app</sub>/E<sub>s</sub></th>"
             "<th>&sigma;<sub>fallo</sub></th><th>&epsilon;<sub>eff</sub> p98"
             "</th><th>" + _("&sigma;<sub>vM</sub> p99 superficie al fallo")
             + "</th><th>" + _("&sigma;<sub>vM</sub> max al fallo") + "</th>"
             "<th>" + _("elem.") + "</th></tr>"]
        anis = []
        for cual, etq in (("voi", _("VOI")), ("spin", self._etq_fam())):
            r = out.get(cual)
            if r is None:
                continue
            if not r.get("ok"):
                t.append(f"<tr><td><b>{etq}</b></td>"
                         f"<td colspan=8>{r['msg']}</td></tr>")
                continue
            ejes = r["ejes"]
            for j, (nom, p) in enumerate(ejes.items()):
                cab = f"<b>{etq}</b>" if j == 0 else ""
                if not p.get("ok"):
                    t.append(f"<tr><td>{cab}</td><td><b>{nom}</b></td>"
                             f"<td colspan=7>{p['msg']}</td></tr>")
                    continue
                vmf = p.get("vm_max_fallo")
                # El p99 de la capa superficial va PRIMERO porque es el que se
                # cita: el maximo no converge con la resolucion. Puede faltar
                # (dominio macizo, sin superficie libre) y entonces no se
                # inventa nada.
                vms = p.get("vm_p99_superficie_fallo")
                cel_s = ("—" if vms is None or not np.isfinite(vms)
                         else f"{vms/1e6:.0f} MPa")
                t.append(
                    f"<tr><td>{cab}</td><td><b>{nom}</b></td>"
                    f"<td>{p['E_app']/1e6:.0f} MPa</td>"
                    f"<td>{p['E_app']/E_S_PA:.4f}</td>"
                    f"<td>{p['sigma_fallo']/1e6:.2f} MPa</td>"
                    f"<td>{p['eps_eff_p']:.3e}</td>"
                    f"<td>{cel_s}</td>"
                    f"<td>{vmf/1e6:.0f} MPa</td>"
                    f"<td>{p['n_elem']:,}</td></tr>")
            E = [p["E_app"] for p in ejes.values()
                 if p.get("ok") and np.isfinite(p.get("E_app", np.nan))
                 and p["E_app"] > 0]
            if len(E) >= 2:
                mx, mn = max(E), min(E)
                nom_max = max((k for k, p in ejes.items() if p.get("ok")),
                              key=lambda k: ejes[k]["E_app"])
                anis.append(
                    _("{etq}: E<sub>max</sub>/E<sub>min</sub> = {r} (mas "
                      "rigido en <b>{eje}</b>)").format(
                        etq=etq, r=f"{mx/mn:.2f}", eje=nom_max))
        t.append("</table>")

        if anis:
            t.append("<p><b>" + _("Anisotropia mecanica") + "</b> — "
                     + " &nbsp;·&nbsp; ".join(anis) + "</p>")
            t.append(
                "<p style='color:#666;font-size:10px'>"
                + _("Es una anisotropia MECANICA y no tiene por que coincidir "
                    "con el DA del tensor MIL, que es puramente geometrico. "
                    "Contrastarla con E<sub>x</sub>/E<sub>y</sub>/E<sub>z</sub> "
                    "del boton «Tensor elastico»: alli el contorno es "
                    "PERIODICO y aqui son platos reales sobre una probeta "
                    "finita, asi que solo convergen si la estructura es grande "
                    "frente a la trabecula. Que difieran es efecto de tamano o "
                    "de borde, no un fallo.")
                + "</p>")

        t.append(
            "<p style='color:#666;font-size:10px'>"
            + _("&sigma;<sub>fallo</sub> es tension APARENTE sobre la seccion "
                "bruta; &sigma;<sub>vM</sub> es la del TEJIDO, y compararla "
                "con el limite elastico del hueso mineralizado "
                "(~150–200 MPa) dice si el criterio de Pistoia esta "
                "prediciendo algo fisicamente coherente. Apoyo: {apoyo}.<br>"
                "El 2% del tejido y el 0.7% de deformacion son convenciones "
                "calibradas en radio distal humano, no constantes fisicas.<br>"
                "<b>E<sub>app</sub> no converge con la malla</b> en "
                "estructuras poco densas: a BV/TV 0.28 se midio 700, 797 y "
                "365 MPa a 32³, 48³ y 64³. Haz un estudio de convergencia "
                "antes de citar un valor.<br>Resuelto en {t} s.").format(
                    apoyo=_(apoyo_de_combo(self.cmb_apoyo)),
                    t=f"{transcurrido:.0f}")
            + "</p>")

        d = QtWidgets.QMessageBox(self)
        d.setWindowTitle(_("Resistencia estimada"))
        d.setTextFormat(QtCore.Qt.RichText)
        d.setText("".join(t))
        hay = []
        pintado = ""
        for k, v in out.items():
            if not v.get("ok"):
                continue
            pintado = v["eje_pintado"]
            p = v["ejes"].get(pintado, {})
            if p.get("ok"):
                hay.append(f"{k} ({pintado}): {p['sigma_fallo']/1e6:.2f} MPa")
        self.statusBar().showMessage(
            _("Carga de fallo aparente — {hay}   [{t} s; los campos que se "
              "pintan son los del eje {eje}: cada eje es un caso de carga "
              "distinto]").format(hay="  ·  ".join(hay),
                                  t=f"{transcurrido:.0f}", eje=pintado))
        self._mostrar(d)

    def analisis_comparado(self):
        """Analisis comparado del VOI frente a cada familia activa."""
        self._en_cada_familia(self._analisis_comparado_una)

    def _analisis_comparado_una(self, ajustar=None):
        """Protocolo del articulo, aplicado al VOI y al spinodoide a la vez.

        Un solo boton porque la comparacion solo significa algo si las dos
        estructuras reciben EXACTAMENTE el mismo ensayo: mismo material, misma
        carga, mismo apoyo, misma resolucion y la misma escala de color al
        pintarlas. Separarlo en dos pulsaciones invita a comparar dos ensayos
        que no son el mismo.

        `ajustar=None` obedece a la casilla «Ajustar antes»; el informe
        automatico pasa False porque ya ajusto en su primera etapa, y volver a
        ajustar aqui sustituiria el candidato que ya se midio y homogeneizo.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None or self.BW_vista is None:
            self._avisar(
                _("Faltan estructuras"),
                _("El analisis comparado necesita las DOS: carga un VOI y "
                  "genera el spinodoide."))
            return

        if ajustar is None:
            ajustar = self.chk_ajustar_antes.isChecked()
        # Si se ajusta antes, la orientacion la decide ese ajuste y aun no
        # existe: solo se comprueba sobre el candidato que ya esta en pantalla.
        if not ajustar and not self._confirmar_orientacion():
            return
        n = int(self.spin_res_fe.value())
        VOI, sp_v = self.VOI, self.VOI_spacing
        BW_s, sp_s = self.BW_vista, self._spacing(self.BW_vista.shape[0])
        # El aviso se prepara aqui, en el hilo de la interfaz, y viaja con el
        # resultado: si el spinodoide no viene de un ajuste, la comparacion es
        # contra una estructura que nadie ha acercado al VOI.
        sin_ajuste = (self.R_forzada is None) and not ajustar
        m_voi = self.m_voi
        semilla = int(self.spin_semilla.value())
        nw = int(self.sl["nw"].valor())

        # Los metodos comparables DE ESTA FAMILIA: el dual-lattice no tiene
        # el muestreo equitativo, que es un esquema de ondas del spinodoide.
        reg = REGISTROS[self._fam]
        claves = [c for c in METODOS_COMPARABLES if c in reg]
        if ajustar:
            self._ocupado(True,
                          _("Ajustando con los dos metodos comparables…")
                          if len(claves) == 2 else
                          _("Ajustando con los tres metodos comparables…"),
                          determinada=True,
                          total=sum(reg[c]["n_eval"] for c in claves))
        else:
            self._ocupado(True, f"Analisis comparado a {n}³: 100 N, 18 GPa, "
                                f"apoyo empotrado…")
        self._t0 = time.time()

        def uno(BW, sp):
            bw, spr = remuestrear_bw(BW, sp, n)
            # unidad='mm' porque `spr` esta en milimetros y E_s en pascales:
            # sin el factor, 100 N se convierten en 1e-4 N.
            r = ensayo_compresion(bw, spr, E_s=PAPER_E_S, nu_s=PAPER_NU,
                                  apoyo=PAPER_APOYO, carga_N=PAPER_CARGA_N,
                                  unidad="mm")
            if not r["ok"]:
                return {"ok": False, "msg": r["msg"]}
            # Direccion principal del MIL, MEDIDA SOBRE LA MISMA MALLA que se
            # ensaya. Es lo que permite saber si las dos estructuras tienen el
            # eje rigido en el mismo sitio; sin esto, dos DA iguales pueden
            # corresponder a orientaciones perpendiculares y la comparacion
            # mecanica no significa nada.
            #
            # Convencion del tensor de fabrica: 1/MIL(n)^2 = n^T M n, asi que
            # el autovalor MENOR corresponde al MIL mas largo, que es la
            # direccion de alineacion trabecular y la mas rigida.
            fab = tensor_mil(bw, spr)
            lam = np.asarray(fab["eigenvalues"], float)
            vec = np.asarray(fab["eigenvectors"], float)
            dir_ppal = vec[:, int(np.argmin(lam))]
            vm = np.asarray(r["vm_solido"], float)
            vm = vm[np.isfinite(vm)]
            de = np.asarray(r["desp_solido"], float)
            de = de[np.isfinite(de)]
            # Los campos salen en las unidades de `spr`, que son mm; los
            # desplazamientos tambien. Las tensiones estan en Pa -> MPa, que
            # es como los reporta el articulo.
            # El p99 de la capa superficial es el unico sabor del pico que
            # converge con la resolucion; el maximo y el p99 global se siguen
            # dando porque estan en los resultados de siempre, pero el numero
            # que se cita es este (`resistencia.estadisticos_vm`).
            sup = estadisticos_vm(r)
            return {"ok": True, "spacing": spr,
                    "campo_vm": r["campo_vm"], "campo_desp": r["campo_desp"],
                    "vm_media": float(vm.mean()) / 1e6,
                    "vm_max": float(vm.max()) / 1e6,
                    "vm_p99": float(np.percentile(vm, 99)) / 1e6,
                    "vm_p99_superficie": (
                        float(sup["vm_p99_superficie"]) / 1e6
                        if sup.get("vm_p99_superficie") is not None
                        else float("nan")),
                    "vm_n_superficie": int(sup.get("vm_n_superficie", 0)),
                    # La distribucion de la que sale ese p99, en MPa: la
                    # figura de la cola se dibuja desde aqui y no del campo.
                    "cuantiles_vm_superficie": cuantiles_vm_superficie(
                        r, escala=1e-6),
                    "desp_max": float(r["desp_max"]),
                    "desp_media": float(de.mean()),
                    "E_app": float(r["E_app"]),
                    "eps_app": float(r["eps_app"]),
                    "DA": float(fab["DA"]),
                    "dir_ppal": [float(x) for x in dir_ppal],
                    "ang_z_deg": float(np.degrees(np.arccos(
                        min(1.0, abs(float(dir_ppal[2])))))),
                    "rho": float(bw.mean()),
                    "n_elem": int(r["n_elem"]),
                    "sigma_app": float(r["sigma_app"]) / 1e6,
                    "frac_portante": float(r.get("frac_portante", 1.0)),
                    "residuo": float(r.get("residuo_rel", 0.0)),
                    "solver": r.get("solver", "")}

        hilo = Trabajador(lambda: None)     # se sustituye abajo

        def tarea():
            bw_spin, sp_spin, ajuste = BW_s, sp_s, None
            if ajustar:
                # Los tres, en serie y en el hilo de trabajo. En serie y no en
                # paralelo porque aqui no hay una barra por metodo que mirar:
                # es un solo paso automatico, y en serie el progreso avanza de
                # forma monotona en vez de a saltos.
                ajuste, _x = mejor_ajuste(VOI, sp_v, reg, claves, m_voi=m_voi,
                                          semilla=semilla, num_waves=nw,
                                          informar=hilo.informar)
                # Se ensaya la MASCARA del ganador, tal cual la produjo el
                # ajuste: regenerarla desde los deslizadores la redondearia a
                # los pasos del control y ya no seria la que gano.
                bw_spin = ajuste["mascara"]
                sp_spin = [float(ajuste["parametros"]["spacing_mm"])] * 3
                hilo.informar(TOTAL_AJUSTE, TOTAL_AJUSTE,
                              "ensayo de compresion")
            return {"voi": uno(VOI, sp_v), "spin": uno(bw_spin, sp_spin),
                    "sin_ajuste": sin_ajuste, "ajuste": ajuste}

        TOTAL_AJUSTE = sum(reg[c]["n_eval"] for c in claves)
        hilo._fn = tarea
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._comparado_listo)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _comparado_listo(self, out):
        transcurrido = time.time() - self._t0
        self._ocupado(False)
        ajuste = out.get("ajuste")
        if ajuste is not None:
            # Deja la interfaz en el estado del ganador: deslizadores, R y
            # tabla de metricas. Reutiliza la ruta del ajuste sencillo para no
            # tener dos sitios donde colocar un resultado.
            self._ajuste_listo(ajuste)
        ok = {k: v for k, v in out.items()
              if k in ("voi", "spin") and v.get("ok")}
        # Una estructura que no resuelve NO se descarta en silencio: si el
        # solver no convergio, la tabla saldria con una sola fila y pareceria
        # que solo se pidio una. Se listan sus mensajes con el mismo peso que
        # los resultados.
        fallos = {k: (v.get("msg") or _("El ensayo no se resolvio."))
                  for k, v in out.items()
                  if k in ("voi", "spin") and not v.get("ok")}
        if not ok:
            texto = "<br>".join(f"<b>{k}</b>: {m}" for k, m in fallos.items())
            if self._auto is not None:
                self._auto["error"] = texto
            else:
                QtWidgets.QMessageBox.critical(
                    self, _("El analisis no se resolvio"), texto)
            return

        # Los campos quedan disponibles para colorear, con escala comun.
        self._vm = {k: (v["campo_vm"], v["spacing"]) for k, v in ok.items()}
        self._vm_comparado[self._fam] = dict(self._vm)
        self._desp = {k: (v["campo_desp"], v["spacing"])
                      for k, v in ok.items()}
        clave_res = self._clave_res("analisis_comparado")
        self._res[clave_res] = {
            "protocolo": "Tapia et al. Biology 2026;15:722",
            "E_s_Pa": PAPER_E_S, "nu_s": PAPER_NU,
            "carga_N": PAPER_CARGA_N, "apoyo": PAPER_APOYO,
            "resolucion": int(self.spin_res_fe.value()),
            "por_estructura": {k: {kk: vv for kk, vv in v.items()
                                   if np.isscalar(vv)
                                   or kk == "cuantiles_vm_superficie"}
                               for k, v in ok.items()},
        }

        if ajuste is not None:
            self._res[clave_res]["ajuste_previo"] = {
                "elegido": ajuste["metodo"],
                "error": float(ajuste["error"]),
                "metodos": ajuste.get("_todos", {}),
            }

        if self._fam == "dual-lattice":
            nom_spin = (_("Dual-lattice ajustado ({m})").format(
                            m=ajuste["etiqueta"])
                        if ajuste is not None else _("Dual-lattice"))
        else:
            nom_spin = (_("Spinodoide ajustado ({m})").format(
                            m=ajuste["etiqueta"])
                        if ajuste is not None else _("Spinodoide"))
        NOM = {"voi": _("VOI real"), "spin": nom_spin}
        t = [_("<b>Analisis mecanico comparado</b><br>"),
             f"<p style='color:#555;font-size:11px'>{PAPER_REF}</p>"]
        if ajuste is not None:
            todos = ajuste.get("_todos", {})
            t.append(_("<p><b>Ajuste previo</b> — se ensaya el de menor "
                       "error:</p>")
                     + "<table cellpadding=4 cellspacing=0 border=1 "
                       "style='border-collapse:collapse'>"
                     + _("<tr><th>metodo</th><th>error</th><th>tiempo</th>"
                         "</tr>"))
            for clave, d in sorted(todos.items(), key=lambda x: x[1]["error"]):
                gana = clave == ajuste["metodo"]
                marca = " ←" if gana else ""
                estilo = " style='background:#e8f5e9'" if gana else ""
                t.append(f"<tr{estilo}><td>{d['etiqueta']}{marca}</td>"
                         f"<td align=right>{d['error']:.5f}</td>"
                         f"<td align=right>{d['tiempo_s']:.0f} s</td></tr>")
            t.append("</table>")
            t.append(
                "<p style='color:#666;font-size:10px'>"
                + _("El desempate mecanico no compite aqui: su error incluye "
                    "Ez/Es y Ez/Ex, dos terminos que los otros no tienen, asi "
                    "que su cifra no es comparable con estas.")
                + "</p>")
        t += [
             "<table cellpadding=5 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>",
             _("<tr><th></th><th>BV/TV</th>"
               "<th>von Mises media<br>[MPa]</th>"
               "<th>von Mises p99 superficie<br>[MPa]</th>"
               "<th>von Mises p99 global<br>[MPa]</th>"
               "<th>von Mises max<br>[MPa]</th>"
               "<th>Deform. total max<br>[mm]</th>"
               "<th>&epsilon; aparente<br>[%]</th>"
               "<th>E<sub>app</sub><br>[MPa]</th></tr>")]
        for k in ("voi", "spin"):
            v = ok.get(k)
            if v is None:
                continue
            t.append(
                f"<tr><td><b>{NOM[k]}</b></td>"
                f"<td align=right>{v['rho']:.4f}</td>"
                f"<td align=right>{v['vm_media']:.3f}</td>"
                f"<td align=right>{_fmt_num(v.get('vm_p99_superficie'), '.3f')}"
                "</td>"
                f"<td align=right>{v['vm_p99']:.3f}</td>"
                f"<td align=right>{v['vm_max']:.3f}</td>"
                f"<td align=right>{v['desp_max']:.5f}</td>"
                f"<td align=right>{100.0 * v['eps_app']:.2f}</td>"
                f"<td align=right>{v['E_app'] / 1e6:.1f}</td></tr>")
        if len(ok) == 2:
            a, b = ok["voi"], ok["spin"]
            def dif(x, y):
                return "—" if x == 0 else f"{100.0 * (y / x - 1.0):+.1f} %"
            t.append(
                "<tr><td><i>"
                + (_("dual-lattice vs VOI") if self._fam == "dual-lattice"
                   else _("spin vs VOI"))
                + "</i></td>"
                f"<td align=right><i>{dif(a['rho'], b['rho'])}</i></td>"
                f"<td align=right><i>{dif(a['vm_media'], b['vm_media'])}</i></td>"
                f"<td align=right><i>"
                f"{dif(a.get('vm_p99_superficie', 0.0), b.get('vm_p99_superficie', 0.0))}"
                "</i></td>"
                f"<td align=right><i>{dif(a['vm_p99'], b['vm_p99'])}</i></td>"
                f"<td align=right><i>{dif(a['vm_max'], b['vm_max'])}</i></td>"
                f"<td align=right><i>{dif(a['desp_max'], b['desp_max'])}</i></td>"
                f"<td align=right><i>{dif(a['eps_app'], b['eps_app'])}</i></td>"
                f"<td align=right><i>{dif(a['E_app'], b['E_app'])}</i></td></tr>")
        t.append("</table>")

        # ORIENTACION: el dato que decide si la comparacion mecanica significa
        # algo. Medido sobre el VOI proximal de H4: el ajuste de menor error
        # reproduce el DA (1.456 frente a 1.480) pero deja el eje rigido a 86
        # grados del que tiene el hueso, y por eso sale 57 veces mas blando en
        # la direccion de carga. El error de ajuste usa el DA ESCALAR, que no
        # distingue "rigido en Z" de "rigido en Y".
        if len(ok) == 2:
            a, b = ok["voi"], ok["spin"]
            ca = abs(float(np.dot(a["dir_ppal"], b["dir_ppal"])))
            ang = float(np.degrees(np.arccos(min(1.0, ca))))
            self._res.setdefault(clave_res, {})["orientacion"] = {
                "angulo_deg": ang,
                "dir_voi": a["dir_ppal"], "dir_spin": b["dir_ppal"],
                "DA_voi": a["DA"], "DA_spin": b["DA"],
            }
            if self._fam == "dual-lattice":
                plantilla = _(
                    "<p><b>Orientacion del eje principal</b> (MIL, medido "
                    "sobre la malla que se ensaya): VOI a {av}° del eje Z (DA "
                    "{dav}), dual-lattice a {as_}° (DA {das}). "
                    "<b>Angulo entre ambos: {ang}°.</b></p>")
            else:
                plantilla = _(
                    "<p><b>Orientacion del eje principal</b> (MIL, medido sobre "
                    "la malla que se ensaya): VOI a {av}° del eje Z (DA {dav}), "
                    "spinodoide a {as_}° (DA {das}). "
                    "<b>Angulo entre ambos: {ang}°.</b></p>")
            t.append(plantilla.format(
                    av=f"{a['ang_z_deg']:.0f}", dav=f"{a['DA']:.3f}",
                    as_=f"{b['ang_z_deg']:.0f}", das=f"{b['DA']:.3f}",
                    ang=f"{ang:.0f}"))
            if ang > 30.0:
                t.append(
                    "<p style='color:#b62324'>"
                    + _("<b>La comparacion mecanica no es interpretable tal "
                        "cual:</b> los ejes rigidos estan a {ang}° uno de "
                        "otro. La carga va en Z, asi que cada estructura "
                        "responde por una direccion distinta de su propia "
                        "fabrica y la diferencia de E<sub>app</sub> mide sobre "
                        "todo ese desalineamiento, no la microarquitectura."
                        "<br>El DA coincide porque es un ESCALAR: no distingue "
                        "«rigido en Z» de «rigido en Y», y el error de ajuste "
                        "solo mira ese escalar.").format(ang=f"{ang:.0f}")
                    + "</p>")

        if fallos:
            NOMF = {"voi": _("VOI real"), "spin": self._etq_fam()}
            t.append("<p style='color:#b62324'><b>" + _("No se resolvio") + " "
                     + (" " + _("ni") + " ").join(NOMF[k] for k in fallos)
                     + ":</b></p><ul>")
            for k, m in fallos.items():
                t.append(f"<li style='color:#b62324'>{NOMF[k]}: {m}</li>")
            t.append("</ul>")
            t.append(
                "<p style='color:#b62324'>"
                + _("Con una sola estructura resuelta <b>no hay "
                    "comparacion</b>. Si el solver no convergio, la causa "
                    "habitual es una resolucion demasiado baja para el grosor "
                    "trabecular: a 20³ esta estructura no converge y a 24³ si. "
                    "Sube la resolucion del ensayo y repite.")
                + "</p>")

        v0 = next(iter(ok.values()))
        t.append(_(
            "<p>Tension aparente sobre la seccion bruta: <b>{sig} MPa</b> "
            "({carga} N). Resolucion {res}³, {n} elementos.</p>").format(
                sig=f"{v0['sigma_app']:.4f}", carga=f"{PAPER_CARGA_N:.0f}",
                res=int(self.spin_res_fe.value()),
                n=f"{v0['n_elem']:,}".replace(",", " ")))
        t.append(
            "<p style='color:#555;font-size:10px'>"
            + _("La <b>von Mises</b> localiza donde se concentra la tension "
                "en el tejido; la <b>deformacion total</b> es el modulo del "
                "desplazamiento y mide la rigidez global, por lo que su "
                "maximo esta siempre en la cara cargada. No son el mismo "
                "mapa.")
            + "</p>")
        eps_grande = max((v["eps_app"] for v in ok.values()), default=0.0)
        if eps_grande > 0.01:
            t.append(
                "<p style='color:#b06000'>"
                + _("<b>Aviso:</b> la deformacion aparente llega al {e} %. El "
                    "analisis es LINEAL, asi que los campos escalan con la "
                    "carga y siguen siendo correctos, pero el hueso real "
                    "habria fallado mucho antes: los 100 N son un estimulo "
                    "estandarizado de comparacion, no una condicion "
                    "fisiologica. Para leer la carga de rotura esta el "
                    "criterio de Pistoia.").format(
                        e=f"{100.0 * eps_grande:.1f}")
                + "</p>")
        if out.get("sin_ajuste") and self._fam == "dual-lattice":
            t.append(
                "<p style='color:#b06000'>"
                + _("<b>Aviso:</b> el dual-lattice actual no procede de un "
                    "ajuste al VOI. A densidad igual, el estiramiento de la "
                    "malla decide en que eje es rigido: ajusta al VOI antes de "
                    "comparar, o estaras midiendo la orientacion que dejaron "
                    "los deslizadores.")
                + "</p>")
        elif out.get("sin_ajuste"):
            t.append(
                "<p style='color:#b06000'>"
                + _("<b>Aviso:</b> el spinodoide actual no procede de un "
                    "ajuste al VOI. Aqui eso pesa mas que en la morfometria: a "
                    "densidad IGUAL, mover solo los angulos de cono cambia "
                    "E<sub>app</sub> de 93 a 1390 MPa —un factor de 15— porque "
                    "el cono ancho marca el eje BLANDO. Ajusta al VOI antes de "
                    "comparar, o estaras midiendo la orientacion que dejaron "
                    "los deslizadores.")
                + "</p>")
        t.append(
            "<p style='color:#b06000;font-size:10px'>"
            + _("<b>Malla:</b> aqui es hexaedrica de un voxel; en el articulo "
                "son tetraedros SOLID187 de 0.05 mm generados en ANSYS. Los "
                "patrones espaciales y el orden de magnitud son comparables; "
                "las cifras absolutas no lo son.")
            + "</p>")
        t.append(
            "<p style='color:#666;font-size:10px'>"
            + _("Resuelto en {t} s ({solver}, residuo {res}). Los dos campos "
                "quedan en «Colorear por» con la MISMA escala en los dos "
                "paneles.").format(t=f"{transcurrido:.1f}",
                                   solver=v0["solver"],
                                   res=f"{v0['residuo']:.1e}")
            + "</p>")

        self.statusBar().showMessage(
            _("Analisis comparado listo en {t} s — von Mises y deformacion "
              "total disponibles en «Colorear por»").format(
                t=f"{transcurrido:.1f}"))
        d = QtWidgets.QMessageBox(self)
        d.setWindowTitle(_("Analisis mecanico comparado"))
        d.setTextFormat(QtCore.Qt.RichText)
        d.setText("".join(t))
        self._mostrar(d)

        # Se deja pintada la von Mises: es la figura con la que abre el
        # articulo, y ver el resultado es mas util que un dialogo cerrado.
        self.cmb_color.setCurrentIndex(3)

    def distribuciones_fe(self):
        self._en_cada_familia(self._distribuciones_una)

    def _distribuciones_una(self):
        if not self._fe and not self._vm:
            QtWidgets.QMessageBox.information(
                self, _("Falta el ensayo"),
                _("Pulsa antes «Resolver y estimar el fallo»: las "
                  "distribuciones salen de los campos de ese calculo."))
            return
        # Los valores se sacan de los campos 3D ya cacheados, quedandose con lo
        # finito: el NaN marca el material que no se resolvio —el filtrado por
        # no ser portante— y meterlo como cero diria que esta descargado.
        datos = {}
        for cual in ("voi", "spin"):
            d = {}
            for clave, cache in (("eps", self._fe), ("vm", self._vm)):
                if cual in cache:
                    c = cache[cual][0]
                    d[clave] = c[np.isfinite(c)]
            if d:
                datos[cual] = d
        DialogoDistribucionFE(datos, self, etq_spin=self._etq_fam()).exec_()

    def convergencia_fe(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        cual, BW, sp = self._estructura_para_fe()
        if BW is None:
            return

        n_max = int(self.spin_res_fe.value())
        # La serie termina en la resolucion del ensayo y baja desde ahi. Es al
        # reves de lo que parece natural, pero es lo unico honesto: refinar por
        # encima de lo que el usuario dijo que podia pagar convertiria un boton
        # de diagnostico en una espera de media hora sin avisar.
        resols = [n for n in (int(n_max * f) for f in (0.55, 0.7, 0.85, 1.0))
                  if n >= 12]
        resols = sorted(set(resols))
        if len(resols) < 3:
            self._avisar(
                _("Resolucion insuficiente"),
                _("Con {n} voxeles no salen tres mallas distintas por encima "
                  "de 12. Sube la resolucion del ensayo.").format(n=n_max))
            return

        apoyo = apoyo_de_combo(self.cmb_apoyo)
        i = self.cmb_eje_fe.currentIndex()
        eje = 2 if i in (0, 3) else (0, 1)[i - 1]

        self._ocupado(True, _("Convergencia sobre el {cual}: {n} mallas "
                              "hasta {res}³…").format(
                                  cual=cual, n=len(resols), res=n_max),
                      determinada=True, total=len(resols))
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)

        codigo = "voi" if cual == "VOI" else self._fam

        def tarea():
            r = estudio_convergencia(
                BW, sp, resoluciones=resols, eje=eje, E_s=E_S_PA, nu_s=NU_S,
                apoyo=apoyo,
                progreso=lambda f, m: hilo.informar(
                    int(f * len(resols)), len(resols), m))
            # El informe necesita saber DE QUE estructura es la serie (su
            # Tb.Th da el eje Tb.Th/h de la figura); la etiqueta de pantalla
            # esta traducida y no sirve de clave.
            r["estructura_codigo"] = codigo
            return r

        hilo._fn = tarea
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._convergencia_lista)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _estructura_para_fe(self):
        """VOI si esta cargado, y si no el spinodoide. Devuelve (etq, BW, sp).

        El estudio se hace sobre UNA estructura, no sobre las dos: son varios
        ensayos completos y duplicar el coste para una comprobacion de malla no
        se justifica. Se prefiere el VOI porque es el dato real, y es sobre el
        que se citaria un E_app.
        """
        if self.VOI is not None:
            return "VOI", self.VOI, self.VOI_spacing
        if self.BW_vista is not None:
            return (self._fam, self.BW_vista,
                    self._spacing(self.BW_vista.shape[0]))
        self._avisar(
            _("Sin estructura"),
            _("Carga un VOI o genera un spinodoide."))
        return None, None, None

    def _convergencia_lista(self, res):
        self._ocupado(False)
        self._res["convergencia"] = res
        v = [p for p in res["puntos"] if p["ok"]]
        self.statusBar().showMessage(
            _("Convergencia en {t} s — {n} de {tot} mallas resueltas."
              ).format(t=f"{time.time()-self._t0:.0f}", n=len(v),
                       tot=len(res["puntos"]))
            + " " + res.get("veredicto", "")[:110])
        self._mostrar(lambda: DialogoConvergencia(res, self))

    # -- simulaciones in silico --------------------------------------------

    def _estructura_simulacion(self):
        """(etiqueta, BW, spacing) segun el desplegable de la seccion."""
        if self.cmb_sim_estructura.currentIndex() == 0:
            if self.VOI is None:
                self._avisar(
                    _("Falta el VOI"),
                    _("Carga un VOI de referencia o elige el candidato "
                      "activo."))
                return None, None, None
            return _("VOI"), self.VOI, self.VOI_spacing
        if self.BW_vista is None:
            self._avisar(
                _("Sin estructura"),
                _("Genera o ajusta primero un candidato."))
            return None, None, None
        return (self._etq_fam(), self.BW_vista,
                self._spacing(self.BW_vista.shape[0]))

    def _codigo_simulacion(self):
        """Codigo de estructura de la simulacion, para `informe.comprobar`."""
        if self.cmb_sim_estructura.currentIndex() == 0:
            return "voi"
        return self._fam

    def _eje_ensayo(self):
        """Eje del ensayo de compresion; 'los tres ejes' se reduce a Z."""
        i = self.cmb_eje_fe.currentIndex()
        return 2 if i in (0, 3) else (0, 1)[i - 1]

    def _avance_sim(self, i, n, etapa):
        if self.barra.maximum() != n:
            self.barra.setRange(0, max(1, n))
        self.barra.setValue(i)
        self.statusBar().showMessage(
            _("Simulacion — {etapa} ({i} de {n})").format(etapa=etapa, i=i,
                                                          n=n))

    def _lanzar_simulacion(self, fn, total, msg):
        self._ocupado(True, msg, determinada=True, total=total)
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)
        hilo._fn = lambda: fn(hilo.informar)
        hilo.avance.connect(self._avance_sim)
        hilo.listo.connect(self._simulacion_lista)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def simular_perdida_gui(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        etq, BW, sp = self._estructura_simulacion()
        if BW is None:
            return
        protocolo = PROTOCOLOS_SIM[self.cmb_sim_protocolo.currentIndex()]
        pasos = int(self.spin_sim_pasos.value())
        paso = self.spin_sim_paso.value() / 100.0
        if pasos * paso >= 0.95:
            self._avisar(
                _("Demasiada perdida"),
                _("{p} pasos de {q} % retiran el {t} % del hueso: no quedaria "
                  "estructura que ensayar. Baja los pasos o el hueso por "
                  "paso.").format(p=pasos, q=int(round(100 * paso)),
                                  t=int(round(100 * pasos * paso))))
            return
        kw = dict(protocolo=protocolo, pasos=pasos, perdida_paso=paso,
                  n_mec=int(self.spin_res_fe.value()), eje=self._eje_ensayo(),
                  E_s=E_S_PA, nu_s=NU_S, apoyo=apoyo_de_combo(self.cmb_apoyo),
                  semilla=int(self.spin_semilla.value()))
        total = 1 + pasos * (2 if protocolo == "recuperacion" else 1)

        codigo = self._codigo_simulacion()

        def tarea(informar):
            r = simular_perdida(BW, sp, progreso=informar, **kw)
            r["estructura"] = etq
            r["estructura_codigo"] = codigo
            return r

        self._lanzar_simulacion(
            tarea, total, _("Simulando perdida osea sobre el {cual}: {n} "
                            "pasos…").format(cual=etq, n=total - 1))

    def fallo_progresivo_gui(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        etq, BW, sp = self._estructura_simulacion()
        if BW is None:
            return
        pasos = int(self.spin_sim_pasos.value())
        kw = dict(pasos=pasos, n_mec=int(self.spin_res_fe.value()),
                  eje=self._eje_ensayo(), E_s=E_S_PA, nu_s=NU_S,
                  apoyo=apoyo_de_combo(self.cmb_apoyo))

        codigo = self._codigo_simulacion()

        def tarea(informar):
            r = fallo_progresivo(BW, sp, progreso=informar, **kw)
            r["estructura"] = etq
            r["estructura_codigo"] = codigo
            return r

        self._lanzar_simulacion(
            tarea, pasos + 1, _("Fallo progresivo sobre el {cual}: {n} "
                                "pasos…").format(cual=etq, n=pasos))

    def _simulacion_lista(self, res):
        self._ocupado(False)
        fallo = res.get("tipo") == "fallo_progresivo"
        # Sin las mascaras: son un volumen por paso y no caben en el JSON.
        self._res["fallo_progresivo" if fallo else "simulacion_perdida"] = {
            k: v for k, v in res.items() if k != "mascaras"}
        s = res.get("resumen", {})
        if fallo:
            txt = _("Fallo progresivo en {t} s: carga maxima {F} N en el paso "
                    "{k}.").format(t=f"{res['tiempo_s']:.0f}",
                                   F=f"{s.get('F_max', float('nan')):.1f}",
                                   k=s.get("paso_F_max"))
        else:
            txt = _("Simulacion en {t} s: rigidez final {E} de la inicial."
                    ).format(t=f"{res['tiempo_s']:.0f}",
                             E=f"{s.get('E_rel_final', float('nan')):.2f}")
        self.statusBar().showMessage(txt)
        self._mostrar(lambda: DialogoSimulacion(res, self))

    # -- exportar como solido ----------------------------------------------

    def _modo_medida(self):
        return "malla" if self.cmb_modo_med.currentIndex() == 1 else "voxel"

    def _sync_formatos(self):
        """El .stl solo tiene sentido por la via TET10."""
        tet = self.cmb_solido.currentIndex() == 1
        self.chk_stl.setEnabled(tet)
        self.chk_stl.setToolTip(
            _("Superficie suavizada, reparada y cerrada: es la que se "
              "imprime.")
            if tet else
            _("Solo disponible con la malla TET10. La hexaedrica solo puede "
              "dar la piel escalonada de los voxeles, que no interesa "
              "imprimir."))
        self.chk_feb.setEnabled(not tet)
        self.chk_feb.setToolTip(
            _("Ensayo de compresión en z listo para FEBio 4: la misma malla, "
              "el mismo apoyo y la misma carga (1 MPa sobre la sección bruta) "
              "que el ensayo de la app. Solo el hueso portante, como ese "
              "ensayo. Validado contra FEBio en comparativa_febio/.")
            if not tet else
            _("Solo disponible con la malla hexaédrica: es la que resuelve "
              "el ensayo de la app, y la única con la que se validó."))

    def exportar_solido(self):
        """Exporta cada familia activa, una detras de otra, a su archivo."""
        self._en_cada_familia(self._exportar_solido_una)

    def _exportar_solido_una(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        tet = self.cmb_solido.currentIndex() == 1
        quiere = {"vtu": self.chk_vtu.isChecked(),
                  "inp": self.chk_inp.isChecked(),
                  "apdl": self.chk_apdl.isChecked(),
                  "stl": self.chk_stl.isChecked() and tet,
                  "feb": self.chk_feb.isChecked() and not tet}
        if not any(quiere.values()):
            QtWidgets.QMessageBox.information(
                self, _("Ningun formato marcado"),
                _("Marca al menos un formato que exportar."))
            return
        dual = self._fam == "dual-lattice"
        base, _x = QtWidgets.QFileDialog.getSaveFileName(
            self,
            _("Exportar el dual-lattice como solido") if dual
            else _("Exportar el spinodoide como solido"),
            str(DATOS / ("dual_lattice" if dual else "spinodoide")),
            _("Nombre base (sin extension) (*)"))
        if not base:
            return
        base = Path(base)
        base = base.with_suffix("")

        # Se regenera a la resolucion de MEDIDA con la misma semilla: la vista
        # suele estar a 32^3 y exportar eso seria entregar una malla mucho mas
        # gruesa de la que muestran las metricas de la tabla.
        res = int(self.sl["resm"].valor())

        # Con la resolucion de medida hasta 256, un descuido aqui produce una
        # malla de millones de elementos y un .inp de varios GB. Se avisa con
        # la cuenta estimada -un hexaedro por voxel solido- antes de empezar,
        # no despues de veinte minutos.
        n_est = int(res ** 3 * self._densidad_actual() / 100.0)
        if n_est > 2_000_000:
            r = QtWidgets.QMessageBox.question(
                self, _("La malla va a ser muy grande"),
                _("A {res}³ y densidad {dens} % salen del orden de "
                  "<b>{n} elementos</b>.<br><br>"
                  "El .inp de Abaqus ronda los 200 bytes por elemento, asi "
                  "que serian varios GB, y la via TET10 puede no terminar."
                  "<br><br>¿Exportar de todas formas?").format(
                    res=res, dens=f"{self._densidad_actual():.0f}",
                    n=f"{n_est:,}".replace(",", " ")),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No)
            if r != QtWidgets.QMessageBox.Yes:
                return
        p = self._params(res)
        sp = self._spacing(res)

        self._ocupado(True,
                      _("Generando a {res}³ y mallando…").format(res=res))
        self._t0 = time.time()

        hilo = Trabajador(lambda: None)

        def tarea():
            BW, _x, _y = Visor._generar(p)
            salidas = []
            if tet:
                nod, ele, sup, inf = malla_tet10(
                    BW, sp, progreso=lambda f, m: hilo.informar(int(100 * f), 100, m))
                if not inf["orden_ok"]:
                    raise RuntimeError(inf["orden"]["msg"])
                if quiere["stl"]:
                    salidas.append(escribir_stl(sup, base.with_suffix(".stl")))
            else:
                nod, ele, inf = malla_hex(BW, sp)
                V, neg = volumen_hex(nod, ele)
                inf["elementos_invertidos"] = neg
                if neg:
                    raise RuntimeError(_(
                        "{n} hexaedros con volumen no positivo; la "
                        "conectividad estaria mal orientada.").format(n=neg))
            if quiere["inp"]:
                salidas.append(escribir_abaqus(nod, ele,
                                               base.with_suffix(".inp"),
                                               E_s=E_S_PA, nu_s=NU_S,
                                               nombre=base.name))
            if quiere["apdl"]:
                salidas.append(escribir_apdl(nod, ele,
                                             base.with_suffix(".apdl"),
                                             E_s=E_S_PA, nu_s=NU_S))
            if quiere["vtu"]:
                salidas.append(escribir_vtu(nod, ele,
                                            base.with_suffix(".vtu")))
            if quiere["feb"]:
                # El ensayo de la app, no la malla entera: solo el hueso que
                # une base y techo (las islas dejarian a FEBio con un sistema
                # singular) y la seccion BRUTA del cubo, que es la del VOI
                # aunque el filtro vacie una columna del borde.
                n_f, e_f, _i = malla_hex(_solo_portante(BW), sp)
                salidas.append(escribir_febio(
                    n_f, e_f, base.with_suffix(".feb"), E_s=E_S_PA, nu_s=NU_S,
                    sigma_app=1e6,
                    A_bruta=float(BW.shape[0] * sp[0] * BW.shape[1] * sp[1]),
                    apoyo=apoyo_de_combo(self.cmb_apoyo)))
            return inf, salidas, base

        hilo._fn = tarea
        hilo.avance.connect(self._avance_export)
        hilo.listo.connect(self._export_listo)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _avance_export(self, i, n, etapa):
        if self.barra.maximum() != n:
            self.barra.setRange(0, max(1, n))
        self.barra.setValue(i)
        self.statusBar().showMessage(
            _("Exportando — {etapa}").format(etapa=etapa))

    def _export_listo(self, r):
        inf, salidas, base = r
        transcurrido = time.time() - self._t0
        self._ocupado(False)

        t = [_("<b>Malla exportada</b> — sistema de unidades mm · N · MPa"
               "<br>"),
             _("<p>{ne} elementos de tipo <b>{tipo}</b>, {nn} nodos "
               "({ng} grados de libertad)<br>").format(
                 ne=f"{inf['n_elems']:,}", tipo=inf["tipo"],
                 nn=f"{inf['n_nodos']:,}", ng=f"{inf['n_gdl']:,}"),
             _("Volumen mallado {v} mm³ — <b>{pct}%</b> del volumen en "
               "vóxeles</p>").format(v=f"{inf['volumen']:.5g}",
                                     pct=f"{inf['vol_pct_BV']:.1f}")]
        if inf["tipo"] == "tet10":
            # Estanqueidad, comprobada sobre la malla concreta antes de
            # escribir. Una STL con bordes abiertos no se imprime y no avisa.
            if inf.get("estanca"):
                t.append(
                    "<p style='color:#1a7f37'>"
                    + _("<b>Superficie estanca</b>: 0 bordes abiertos, {nc} "
                        "componente(s). Apta para imprimir y para "
                        "tetraedralizar.").format(
                            nc=inf.get("n_componentes", "?"))
                    + "</p>")
            else:
                t.append(
                    "<p style='color:#b62324'>"
                    + _("<b>Superficie NO estanca</b>: {nb} bordes abiertos, "
                        "{nc} componente(s). El .stl no es imprimible tal "
                        "cual; sube la resolución de medida o usa la malla "
                        "hexaédrica.").format(
                            nb=inf.get("bordes_abiertos", "?"),
                            nc=inf.get("n_componentes", "?"))
                    + "</p>")
            # Perdida REAL (suavizado + decimado + reparacion) frente a la
            # superficie cerrada cruda. La diferencia entre el conteo de
            # voxeles y esa superficie es de definicion, no de perdida.
            vmc = inf.get("vol_pct_MC", float("nan"))
            if np.isfinite(vmc):
                t.append(_(
                    "<p>Frente a la superficie cerrada cruda: <b>{pct}%</b> "
                    "del volumen (la superficie cruda es ya un {dif}% "
                    "respecto al conteo de vóxeles, por definición).</p>"
                ).format(pct=f"{vmc:.1f}",
                         dif=f"{100 * (inf['V_mc'] / inf['BV_voxel'] - 1):+.1f}"))
            if np.isfinite(vmc) and vmc < 95:
                t.append(
                    "<p style='color:#b06000'>"
                    + _("<b>Aviso:</b> el mallado perdió más del 5% del "
                        "volumen de la superficie cruda. Con trabéculas de "
                        "uno o dos vóxeles el suavizado y la reparación se "
                        "comen material: sube la resolución de medida, o usa "
                        "la malla hexaédrica, que es exacta.")
                    + "</p>")
        t += [
             "<table cellpadding=4 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>"]
        for s in salidas:
            nom = Path(s["ruta"]).name
            extra = s.get("tipo", "")
            if "n_base" in s:
                extra += _(" · base {b} / techo {t} nodos").format(
                    b=s["n_base"], t=s["n_techo"])
            if "triangulos" in s:
                extra = _("{n} triángulos").format(n=f"{s['triangulos']:,}")
            t.append(f"<tr><td>{nom}</td><td align=right>{s['MB']:.1f} MB</td>"
                     f"<td>{extra}</td></tr>")
        t.append("</table>")
        pie = [_("Carpeta: {d}").format(d=base.parent)]
        if any(Path(s["ruta"]).suffix == ".inp" for s in salidas):
            pie.append(_("El .inp lleva los conjuntos de nodos BASE y TECHO "
                         "para aplicar el apoyo y la carga con un clic."))
        if any(Path(s["ruta"]).suffix == ".feb" for s in salidas):
            pie.append(_("El .feb es el ensayo completo (apoyo, 1 MPa, "
                         "salidas): se corre tal cual con febio4 -i. FEBio es "
                         "no lineal geométricamente; con voladizos en el "
                         "techo se aparta del ensayo lineal de la app."))
        pie.append(_("Generado en {t} s.").format(t=f"{transcurrido:.1f}"))
        t.append("<p style='color:#666;font-size:10px'>"
                 + "<br>".join(pie) + "</p>")

        d = QtWidgets.QMessageBox(self)
        d.setWindowTitle(_("Exportación completada"))
        d.setTextFormat(QtCore.Qt.RichText)
        d.setText("".join(t))
        self.statusBar().showMessage(
            _("Exportado {tipo}: {n} elementos, {pct}% del volumen, en {t} s"
              ).format(tipo=inf["tipo"], n=f"{inf['n_elems']:,}",
                       pct=f"{inf['vol_pct_BV']:.1f}",
                       t=f"{transcurrido:.1f}"))
        d.exec_()

    # -- tensor elastico ---------------------------------------------------

    def calcular_elastico(self):
        """Tensor elastico de cada familia activa, una tras otra."""
        self._en_cada_familia(self._calcular_elastico_una)

    def _calcular_elastico_una(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.BW_vista is None:
            self._avisar(
                _("Sin estructura"), _("Genera primero un spinodoide."))
            return
        if not self._confirmar_orientacion():
            return

        n = int(self.spin_homog.value())
        BW_s = self.BW_vista
        sp_s = self._spacing(BW_s.shape[0])
        VOI, sp_v = self.VOI, self.VOI_spacing

        self._ocupado(True, _("Homogeneizando a {n}³ (y el VOI si esta "
                              "cargado)…").format(n=n))
        self._t0 = time.time()

        def tarea():
            out = {}
            bw, sp = remuestrear_bw(BW_s, sp_s, n)
            C, info = homogeneizar(bw, E_S_PA, NU_S, vox_size=sp)
            out["spin"] = (C, info, bw.shape)
            if VOI is not None:
                bwv, spv = remuestrear_bw(VOI, sp_v, n)
                Cv, iv = homogeneizar(bwv, E_S_PA, NU_S, vox_size=spv)
                out["voi"] = (Cv, iv, bwv.shape)
            return out

        self.hilo = Trabajador(tarea)
        self.hilo.listo.connect(self._elastico_listo)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    def _elastico_listo(self, out):
        # El tiempo se toma AQUI, antes del dialogo modal: `exec_()` bloquea
        # hasta que el usuario lo cierra, y medir despues contaria ese rato
        # como si fuera calculo.
        transcurrido = time.time() - self._t0
        self._ocupado(False)
        filas = []
        claves = []
        for clave, etiqueta in (("voi", _("VOI")),
                                ("spin", self._etq_fam())):
            if clave not in out:
                continue
            C, info, forma = out[clave]
            claves.append(clave)
            if not info["ok"]:
                filas.append((etiqueta, forma, None, info["msg"]))
                continue
            filas.append((etiqueta, forma, constantes_ingenieria(C),
                          info["solver"]))

        if not filas:
            return

        self._res[self._clave_res("elastico")] = {
            "resolucion": int(self.spin_homog.value()),
            "E_s_Pa": E_S_PA, "nu_s": NU_S,
            "por_estructura": {
                clave: {
                    "rejilla": list(forma),
                    "solver": nota,
                    **({} if ec is None else {k: float(v) for k, v in ec.items()}),
                    # F10: el residuo es lo que dice si C se puede citar. Sin
                    # guardarlo, el informe para publicacion no lo comprueba.
                    "residuo_rel": out[clave][1].get("residuo_rel"),
                    # El tensor completo va al JSON: las constantes de
                    # ingenieria son una lectura de C, no un sustituto. Sin C
                    # no se puede recomputar nada ni comprobar la simetria.
                    "C_Pa": np.asarray(out[clave][0], float).tolist(),
                } for clave, (etq, forma, ec, nota) in zip(claves, filas)},
        }

        txt = [_("<b>Modulos efectivos</b> (material base E<sub>s</sub> = "
                 "20 GPa, &nu;<sub>s</sub> = 0.30)<br>")]
        txt.append("<table cellpadding=4 cellspacing=0 border=1 "
                   "style='border-collapse:collapse'>")
        txt.append(_("<tr><th></th><th>rejilla</th><th>E<sub>x</sub></th>"
                     "<th>E<sub>y</sub></th><th>E<sub>z</sub></th>"
                     "<th>E<sub>z</sub>/E<sub>s</sub></th>"
                     "<th>E<sub>z</sub>/E<sub>x</sub></th></tr>"))
        for etiqueta, forma, ec, nota in filas:
            dims = "×".join(str(s) for s in forma)
            if ec is None:
                txt.append(f"<tr><td><b>{etiqueta}</b></td><td>{dims}</td>"
                           f"<td colspan=5>{nota}</td></tr>")
                continue
            txt.append(
                f"<tr><td><b>{etiqueta}</b></td><td>{dims}</td>"
                f"<td>{ec['Ex']/1e9:.3f} GPa</td>"
                f"<td>{ec['Ey']/1e9:.3f} GPa</td>"
                f"<td>{ec['Ez']/1e9:.3f} GPa</td>"
                f"<td>{ec['Ez']/E_S_PA:.4f}</td>"
                f"<td>{ec['Ez']/ec['Ex']:.3f}</td></tr>")
        txt.append("</table>")
        txt.append("<p style='color:#666;font-size:10px'>"
                   + _("E<sub>z</sub>/E<sub>s</sub> es la rigidez axial "
                       "normalizada y E<sub>z</sub>/E<sub>x</sub> la "
                       "anisotropia elastica. Son los dos terminos mecanicos "
                       "que la app puede sumar al error de ajuste.<br>"
                       "La estructura se remuestreo antes de homogeneizar: el "
                       "coste crece con el cubo del lado.")
                   + "</p>")

        d = QtWidgets.QMessageBox(self)
        d.setWindowTitle(_("Tensor elastico homogeneizado"))
        d.setTextFormat(QtCore.Qt.RichText)
        d.setText("".join(txt))
        self.statusBar().showMessage(
            _("Homogeneizacion terminada en {t} s").format(
                t=f"{transcurrido:.1f}"))
        self._mostrar(d)

    # -- sesion ------------------------------------------------------------

    # Las claves del JSON son las de los controles, no un formato propio: asi
    # el archivo se lee a ojo y se edita a mano sin abrir la app.
    # `cmb_color` queda FUERA a proposito: restaurarlo a "espesor local" o a
    # von Mises lanzaria un calculo pesado nada mas abrir la sesion, o pintaria
    # un panel liso porque los campos no se guardan. Se documenta abajo.
    _COMBOS = {"esquema": "cmb_esq", "eje_recorte": "cmb_eje",
               "apoyo": "cmb_apoyo", "solido": "cmb_solido",
               "eje_ensayo": "cmb_eje_fe", "modo_medida": "cmb_modo_med"}
    _CHECKS = {"auto": "chk_auto", "sync": "chk_sync", "ejes": "chk_ejes",
               "extra": "chk_extra", "mecanico": "chk_mec",
               "fmt_vtu": "chk_vtu", "fmt_inp": "chk_inp",
               "fmt_apdl": "chk_apdl", "fmt_stl": "chk_stl",
               "fmt_feb": "chk_feb",
               "poro": "chk_poro", "ef": "chk_ef"}

    def _procedencia_controles(self):
        """Procedencia de lo que describen los controles, por familia.

        El numero de onda sale del deslizador (multiplo de pi) y el bloque
        guarda las dos lecturas: una sesion no puede volver a perder la pi.
        """
        sem = int(self.spin_semilla.value())
        esq = "rechazo" if self.cmb_esq.currentIndex() == 0 else "equitativo"
        return {
            "spinodoide": procedencia.bloque(
                "spinodoide", esquema=esq, semilla=sem,
                wave_number_pi=float(self.sl["wave"].valor())),
            "dual-lattice": procedencia.bloque("dual-lattice", semilla=sem),
        }

    def _estado(self):
        return {
            "formato": "spinpy.visor/sesion",
            # 2: lleva `procedencia`. Las de version 1 se cargan con aviso.
            "version": 2,
            "procedencia": self._procedencia_controles(),
            "guardado": time.strftime("%Y-%m-%d %H:%M:%S"),
            "voi": self.VOI_ruta,
            "semilla": int(self.spin_semilla.value()),
            "res_homog": int(self.spin_homog.value()),
            "res_fe": int(self.spin_res_fe.value()),
            "peso_mecanico": float(self.spin_peso_mec.value()),
            "res_mec": int(self.spin_res_mec.value()),
            "deslizadores": {k: float(d.valor()) for k, d in self.sl.items()},
            "combos": {k: self._combo(c).currentIndex()
                       for k, c in self._COMBOS.items()},
            "casillas": {k: bool(getattr(self, c).isChecked())
                         for k, c in self._CHECKS.items()},
            # La R del ajuste NO se deriva de los deslizadores (correccion F4):
            # sin guardarla, cargar la sesion daria otra orientacion.
            "R_forzada": (None if self._de("spinodoide", "R") is None
                          else np.asarray(self._de("spinodoide", "R"),
                                          float).tolist()),
            # Familia del selector y R del ajuste del dual-lattice. La R del
            # spinodoide conserva su clave de siempre: una sesion guardada
            # antes de que hubiera familias se sigue cargando igual.
            "familia": self._modo_fam,
            "R_forzada_dual": (None if self._de("dual-lattice", "R") is None
                               else np.asarray(self._de("dual-lattice", "R"),
                                               float).tolist()),
        }

    def _combo(self, nombre):
        return getattr(self, nombre)

    def guardar_sesion(self):
        f, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Guardar la sesion"), str(DATOS / "sesion.json"),
            _("JSON (*.json)"))
        if not f:
            return
        Path(f).write_text(json.dumps(self._estado(), indent=2,
                                      ensure_ascii=False), encoding="utf-8")
        self.statusBar().showMessage(
            _("Sesion guardada en {f} — regenera la MISMA realizacion, porque "
              "la semilla viaja dentro.").format(f=f))

    def cargar_sesion(self):
        f, _x = QtWidgets.QFileDialog.getOpenFileName(
            self, _("Cargar una sesion"), str(DATOS),
            _("JSON (*.json);;Todos (*)"))
        if not f:
            return
        try:
            d = json.loads(Path(f).read_text(encoding="utf-8"))
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self, _("No se pudo leer"), str(e))
            return
        if d.get("formato") != "spinpy.visor/sesion":
            QtWidgets.QMessageBox.critical(
                self, _("No es una sesion"),
                _("El archivo no lleva la marca 'spinpy.visor/sesion'."))
            return

        # `_aplicando` bloquea la regeneracion mientras se colocan los valores:
        # sin el, cada deslizador dispararia una vista y se generarian ocho
        # campos intermedios que nadie va a mirar.
        self._aplicando = True
        try:
            for k, v in (d.get("deslizadores") or {}).items():
                if k in self.sl:
                    self.sl[k].fijar(float(v))
            for k, c in self._COMBOS.items():
                i = (d.get("combos") or {}).get(k)
                if i is not None:
                    w = self._combo(c)
                    w.setCurrentIndex(int(i) if 0 <= int(i) < w.count() else 0)
            for k, c in self._CHECKS.items():
                v = (d.get("casillas") or {}).get(k)
                if v is not None:
                    getattr(self, c).setChecked(bool(v))
            if d.get("semilla") is not None:
                self.spin_semilla.setValue(int(d["semilla"]))
            if d.get("res_homog"):
                self.spin_homog.setValue(int(d["res_homog"]))
            if d.get("res_fe"):
                self.spin_res_fe.setValue(int(d["res_fe"]))
            if d.get("peso_mecanico"):
                self.spin_peso_mec.setValue(float(d["peso_mecanico"]))
            if d.get("res_mec"):
                self.spin_res_mec.setValue(int(d["res_mec"]))
            fam = d.get("familia", "spinodoide")
            if fam in MODOS_FAMILIA:
                self.cmb_familia.setCurrentIndex(MODOS_FAMILIA.index(fam))
            for f_, clave in (("spinodoide", "R_forzada"),
                              ("dual-lattice", "R_forzada_dual")):
                R = d.get(clave)
                self._fijar_R(f_, R)
                self._lab_R[f_].setText(
                    "" if R is None else
                    _("Orientacion FIJADA por una sesion guardada (correccion "
                      "F4). Mueve cualquier rotacion para liberarla."))
        finally:
            self._aplicando = False

        avisos = []
        if int(d.get("version", 1)) < 2:
            avisos.append(_(
                "Sesion en un formato anterior, sin procedencia: no consta "
                "con que version de spinpy se guardo. El numero de onda se "
                "lee como multiplo de pi, que es lo que siempre guardo el "
                "deslizador."))
        else:
            proc = (d.get("procedencia") or {}).get("spinodoide") or {}
            w = (d.get("deslizadores") or {}).get("wave")
            if w is not None and proc.get("wave_number_rad") is not None:
                try:
                    procedencia.numero_onda(float(w), proc["wave_number_rad"])
                except procedencia.ErrorNumeroOnda as e:
                    avisos.append(_("El numero de onda de la sesion no cuadra "
                                    "con su procedencia: {e}").format(e=e))
        ruta = d.get("voi")
        if ruta:
            if Path(ruta).exists():
                self.cargar_voi(ruta)
            else:
                # El VOI puede no estar: el banco vive fuera del repositorio y
                # las rutas son absolutas. Se avisa y se sigue, porque los
                # parametros del spinodoide son utiles igual; fallar entero
                # seria perder tambien lo que si se puede restaurar.
                avisos.append(_("El VOI '{ruta}' no esta donde lo dejaste; "
                                "cargalo a mano si lo necesitas."
                                ).format(ruta=ruta))

        # El coloreado se restaura a material: los campos de espesor, de
        # deformacion y de von Mises no se guardan —son cientos de megas— y
        # dejar el desplegable en uno de ellos mostraria un panel liso sin
        # explicar por que.
        self._esp, self._fe, self._vm = {}, {}, {}
        self._desp, self._clim = {}, None
        for f in FAMILIAS:
            self._cand[f].update(m=None, vigente=False, esp=None, fe=None,
                                 vm=None, desp=None)
        self.cmb_color.setCurrentIndex(0)
        self.m_spin = None
        self._generar_vista()
        self.statusBar().showMessage(
            _("Sesion de {fecha} cargada.").format(
                fecha=d.get("guardado", "?"))
            + ("  " + "  ".join(avisos) if avisos else
               "  " + _("Misma semilla: la estructura es la misma "
                        "realizacion.")))
        if avisos:
            QtWidgets.QMessageBox.warning(
                self, _("Sesion cargada a medias"), "\n".join(avisos))

    # -- dispersion del generador y lote -----------------------------------

    def dispersion(self):
        """Dispersion entre semillas de cada familia activa, una tras otra."""
        self._en_cada_familia(self._dispersion_una)

    def _dispersion_una(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        K = int(self.spin_k.value())
        res = int(self.sl["resm"].valor())
        p = self._params(res)
        sp = self._spacing(res)
        m_voi = self.m_voi
        extra = self.chk_extra.isChecked()
        poro = self.chk_poro.isChecked()

        self._ocupado(True, f"{K} realizaciones a {res}³…",
                      determinada=True, total=K)
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)

        # Semillas seed+1..seed+K, como las replicas del ajuste (U1): la de la
        # busqueda es la que gano, y meterla sesga la media hacia el VOI.
        # `informe.comprobar` marca con reservas una dispersion que la use.
        base = int(p.get("seed", 20260720)) + 1

        def tarea():
            return dispersion_semillas(p, sp, n_semillas=K, semilla_base=base,
                                       m_ref=m_voi, extra=extra,
                                       progreso=hilo.informar)

        hilo._fn = tarea
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._dispersion_lista)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _dispersion_lista(self, disp):
        self._ocupado(False)
        self._res[self._clave_res("dispersion")] = {"resumen": disp["resumen"],
                                   "n_semillas": disp["n_semillas"],
                                   "semilla_base": disp["semilla_base"]}
        r = disp["resumen"]
        fuera = [m for m, d in r.items()
                 if np.isfinite(d.get("z", np.nan)) and abs(d["z"]) >= 2]
        cv = {m: d["cv_pct"] for m, d in r.items()
              if np.isfinite(d.get("cv_pct", np.nan))}
        msg = (_("Dispersion de {n} realizaciones en {t} s — CV entre "
                 "{a}% y {b}%").format(
                     n=disp["n_semillas"], t=f"{time.time()-self._t0:.0f}",
                     a=f"{min(cv.values()):.2f}", b=f"{max(cv.values()):.2f}")
               if cv else _("Dispersion calculada"))
        if self.m_voi is not None:
            msg += ("   ·  " + _("{n} metrica/s fuera de banda (|z|>=2)"
                                 ).format(n=len(fuera))
                    + (f": {', '.join(fuera)}" if fuera else ""))
        self.statusBar().showMessage(msg)
        self._mostrar(lambda: DialogoDispersion(disp, self))

    def correr_lote_gui(self):
        if self.hilo is not None and self.hilo.isRunning():
            return
        ini = str(VOIDIR if VOIDIR.exists() else DATOS)
        rutas, _x = QtWidgets.QFileDialog.getOpenFileNames(
            self, _("VOIs del lote (se pueden elegir varios)"), ini,
            _("VOI (*.vtk *.mat);;Todos (*)"))
        if not rutas:
            return

        n_rep = int(self.spin_repl.value())
        modo = "rapido" if self.cmb_modo_lote.currentIndex() == 0 else "completo"
        # El coste se declara ANTES, no despues: son minutos por VOI y el
        # usuario no puede estimarlo desde la interfaz. Un lote que arranca sin
        # avisar y bloquea la ventana media hora es un boton hostil.
        fams = self._familias()
        est = len(rutas) * len(fams) * (0.7 if modo == "rapido" else 2.5)
        if QtWidgets.QMessageBox.question(
                self, _("Confirmar el lote"),
                _("<p>{nv} VOI(s) × {nr} replicas, modo <b>{modo}</b>.</p>"
                  "<p>Estimacion muy gruesa: <b>{a}–{b} minutos</b>, segun el "
                  "tamano de los VOIs. La ventana queda inutilizable mientras "
                  "tanto.</p><p>¿Seguir?</p>").format(
                    nv=len(rutas), nr=n_rep, modo=modo,
                    a=f"{est:.0f}", b=f"{3*est:.0f}"),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No) != QtWidgets.QMessageBox.Yes:
            return

        self._ocupado(True, f"Lote: {len(rutas)} VOIs…",
                      determinada=True, total=len(rutas) * len(fams))
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)

        def tarea():
            # Un lote completo por familia, en serie: mismos VOIs, mismas
            # replicas y mismo modo, y tablas separadas que llevan su columna
            # `familia` para poder unirlas sin mezclarlas.
            salidas = []
            for k, f in enumerate(fams):
                base = k * len(rutas)
                df, aj, meta = correr_lote(
                    rutas, n_replicas=n_rep, modo=modo, familia=f,
                    progreso=lambda i, n, msg, _b=base: hilo.informar(
                        _b + i, len(rutas) * len(fams), msg))
                inf = informe(df, [m for m in METRICAS if m in df.columns])
                meta["n_filas"] = len(df)
                meta["familia"] = f
                salidas.append((df, aj, inf, meta))
            return salidas

        hilo._fn = tarea
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._lotes_listos)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _lotes_listos(self, salidas):
        for r in salidas:
            self._lote_listo(r)

    def _lote_listo(self, r):
        df, aj, inf, meta = r
        dual = meta.get("familia") == "dual-lattice"
        prefijo = "lote_dual" if dual else "lote"
        self._ocupado(False)
        self._lote = (df, aj, inf)
        self._res["lote_dual" if dual else "lote"] = {
            "meta": meta, "varianza": inf.to_dict(orient="records")}

        # Se guarda SIEMPRE y sin preguntar: el lote acaba de costar minutos u
        # horas de CPU y perderlo por cerrar un dialogo seria inaceptable. El
        # usuario puede reexportarlo donde quiera desde el menu.
        SAL = DATOS / "resultados"
        SAL.mkdir(parents=True, exist_ok=True)
        sello = time.strftime("%Y%m%d_%H%M%S")
        for nombre, tabla in (("especimenes", df), ("ajustes", aj),
                              ("varianza", inf)):
            tabla.to_csv(SAL / f"{prefijo}_{nombre}_{sello}.csv", index=False,
                         sep=";", decimal=",", encoding="utf-8-sig")

        n_ef = _n_efectivo_min(inf)
        self.statusBar().showMessage(
            _("Lote terminado en {t} min — {n} filas, pero N efectivo = "
              "{nef}. CSV en resultados/lote_*_{sello}.csv").format(
                t=f"{(time.time()-self._t0)/60:.1f}", n=len(df),
                nef=f"{n_ef:.1f}", sello=sello))
        DialogoVarianza(inf, meta, self).exec_()

    # -- exportar los resultados numericos ---------------------------------

    def exportar_resultados(self):
        """Vuelca a CSV y JSON todo lo calculado en la sesion.

        DOS ARCHIVOS, PORQUE SIRVEN PARA COSAS DISTINTAS. El CSV es plano y
        largo (una fila por magnitud) para pegarlo en una hoja o leerlo con
        pandas; el JSON conserva lo que no cabe en una tabla —el tensor C
        entero, los parametros, la traza— y es el que permite reproducir.

        Se exporta lo que YA se calculo, sin recalcular nada. Con un generador
        estocastico, recalcular al exportar daria numeros distintos de los que
        el usuario vio en pantalla, y el archivo dejaria de ser el registro de
        lo ocurrido.
        """
        if not self._res and all(self._de(f, "m") is None for f in FAMILIAS):
            QtWidgets.QMessageBox.information(
                self, _("Nada que exportar"),
                _("Mide, homogeneiza o ensaya algo primero."))
            return

        base, _x = QtWidgets.QFileDialog.getSaveFileName(
            self, _("Exportar los resultados"),
            str(DATOS / "resultados_visor"),
            _("Nombre base (sin extension) (*)"))
        if not base:
            return
        base = Path(str(base).rsplit(".", 1)[0])

        doc = self._documento_resultados()
        m_spin = self._de("spinodoide", "m")
        m_dual = self._de("dual-lattice", "m")
        ctrl = doc["procedencia"]["controles"]

        # La columna del dual-lattice solo aparece si hay algo suyo: el CSV de
        # una sesion solo de spinodoide conserva exactamente el formato de
        # antes, que es el que ya leen las hojas de calculo hechas.
        hay_dual = (m_dual is not None
                    or any(str(k).endswith("_dual") for k in self._res))

        def fila(bloque, mag, v, s, dl, uni):
            return ((bloque, mag, v, s, dl, uni) if hay_dual
                    else (bloque, mag, v, s, uni))

        filas = [fila("bloque", "magnitud", "VOI", "spinodoide",
                      "dual-lattice", "unidad")]
        for clave, etq, _x, uni in TABLA:
            v = (self.m_voi or {}).get(clave)
            s = (m_spin or {}).get(clave)
            dl = (m_dual or {}).get(clave)
            if v is None and s is None and dl is None:
                continue
            filas.append(fila("morfometria", etq, self._num(v), self._num(s),
                              self._num(dl), uni))

        ela = self._res.get("elastico", {}).get("por_estructura", {})
        ela_d = self._res.get("elastico_dual", {}).get("por_estructura", {})
        voi_ela = ela.get("voi") or ela_d.get("voi") or {}
        for k in ("Ex", "Ey", "Ez", "Gyz", "Gxz", "Gxy", "anisotropia_E"):
            if any(k in d for d in list(ela.values()) + list(ela_d.values())):
                filas.append(fila("elastico", k,
                                  self._num(voi_ela.get(k)),
                                  self._num(ela.get("spin", {}).get(k)),
                                  self._num(ela_d.get("spin", {}).get(k)),
                                  "Pa" if k[0] in "EG" else ""))

        res = self._res.get("resistencia", {}).get("por_estructura", {})
        res_d = self._res.get("resistencia_dual", {}).get("por_estructura", {})
        voi_r = res.get("voi") or res_d.get("voi") or {}
        ejes = sorted({n for d in list(res.values()) + list(res_d.values())
                       for n in d["ejes"]})
        for eje in ejes:
            for k, uni in (("E_app", "Pa"), ("sigma_fallo", "Pa"),
                           ("factor", ""), ("eps_eff_p", ""),
                           ("vm_p99_superficie_fallo", "Pa"),
                           ("vm_n_superficie", ""),
                           ("vm_max_fallo", "Pa"), ("frac_portante", "")):
                filas.append(fila(
                    "resistencia", f"{k} [{eje}]",
                    self._num(voi_r.get("ejes", {}).get(eje, {}).get(k)),
                    self._num(res.get("spin", {}).get("ejes", {})
                              .get(eje, {}).get(k)),
                    self._num(res_d.get("spin", {}).get("ejes", {})
                              .get(eje, {}).get(k)),
                    uni))

        # Procedencia al FINAL: las hojas ya hechas leen las filas de arriba por
        # posicion y no se mueven.
        for k, uni in (("spinpy", ""), ("version_formato", ""),
                       ("esquema", ""), ("semilla", ""),
                       ("wave_number_pi", "x pi"), ("wave_number_rad", "rad")):
            filas.append(fila("procedencia", k, "",
                              self._txt(ctrl["spinodoide"].get(k)),
                              self._txt(ctrl["dual-lattice"].get(k)), uni))

        # Separador ';' y coma decimal: es lo que Excel en configuracion
        # espanola abre sin pasar por el asistente de importacion. Con ','
        # cada fila cae entera en la primera celda.
        csv = "\n".join(";".join(str(c) for c in fila) for fila in filas)
        (base.with_suffix(".csv")).write_text(csv, encoding="utf-8-sig")
        (base.with_suffix(".json")).write_text(
            json.dumps(doc, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8")

        self.statusBar().showMessage(
            _("Resultados en {csv} y {json} ({n} magnitudes; el JSON lleva "
              "ademas el tensor C completo y la semilla)").format(
                csv=base.name + ".csv", json=base.name + ".json",
                n=len(filas) - 1))

    def _documento_resultados(self):
        """El documento de la sesion: lo escriben «Exportar resultados» y el
        informe para publicacion.

        Uno solo para los dos. Si el informe armara el suyo, podria describir
        una sesion distinta de la que queda exportada al lado.
        """
        estado = self._estado()
        ctrl = estado["procedencia"]
        proc_doc = procedencia.bloque(
            self._modo_fam, esquema=ctrl["spinodoide"]["esquema"],
            semilla=int(self.spin_semilla.value()))
        # La de cada familia segun los controles, y la de cada ajuste tal como
        # lo devolvio el ajuste: si difieren, manda la del ajuste.
        proc_doc["controles"] = ctrl
        proc_doc["ajustes"] = {k: v["procedencia"] for k, v in self._res.items()
                               if isinstance(v, dict) and v.get("procedencia")}
        voi = None
        if self.VOI is not None:
            # La huella del archivo es lo que permite a otro comprobar que
            # parte del mismo VOI. Una carpeta de rebanadas no tiene una sola
            # huella y queda en None, declarado asi en el informe.
            voi = {"nombre": self.VOI_nombre, "ruta": self.VOI_ruta,
                   "forma": [int(s) for s in self.VOI.shape],
                   "spacing_mm": [float(x) for x in
                                  np.ravel(self.VOI_spacing)],
                   "sha256": sha256_archivo(self.VOI_ruta)}
            if self.VOI_recorte:
                # Figura 0 desde la CLI: con esto se relee la pila y se
                # vuelve a situar el cubo.
                voi["recorte"] = self.VOI_recorte
        return {"procedencia": proc_doc,
                "estado": estado, "resultados": dict(self._res),
                "voi": voi,
                "morfometria_voi": self._limpiar(self.m_voi),
                "morfometria_spin": self._limpiar(self._de("spinodoide", "m")),
                "morfometria_dual": self._limpiar(self._de("dual-lattice", "m"))}

    def informe_publicacion(self):
        """Metodos en ES/EN, citabilidad de cada resultado y reproduccion.

        Describe lo que YA se calculo, igual que la exportacion. Lo unico que
        ejecuta es el generador, para anotar la huella de cada estructura
        ajustada, y por eso va en un hilo de trabajo.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        if (not self._res and self.m_voi is None
                and all(self._de(f, "m") is None for f in FAMILIAS)):
            QtWidgets.QMessageBox.information(
                self, _("Nada que informar"),
                _("Mide, ajusta u homogeneiza algo primero: el informe "
                  "describe lo que ya se calculo."))
            return
        carpeta = self._carpeta_informe()
        if carpeta is None:
            return
        self._lanzar_informe(carpeta)

    def _carpeta_informe(self):
        """Pide la carpeta del informe y confirma si pisa uno anterior.

        Devuelve la ruta, o None si el usuario se echa atras.
        """
        carpeta = QtWidgets.QFileDialog.getExistingDirectory(
            self, _("Elige la carpeta del informe para publicacion"),
            str(DATOS))
        if not carpeta:
            return None
        carpeta = Path(carpeta)
        previos = [n for n in informe_pub.ARCHIVOS.values()
                   if (carpeta / n).exists()]
        if previos and QtWidgets.QMessageBox.question(
                self, _("La carpeta ya tiene un informe"),
                _("Se sobrescribiran {n} archivo(s) de un informe anterior "
                  "en esa carpeta. ¿Continuar?").format(n=len(previos))
        ) != QtWidgets.QMessageBox.Yes:
            return None
        return carpeta

    def _lanzar_informe(self, carpeta):
        """Las tres etapas del informe sobre la sesion tal como esta.

        Las figuras 3D salen con las ultimas opciones elegidas en la ventana
        del informe automatico (`opciones_figuras`): el informe manual no
        tiene ventana propia, y asi los dos producen la misma presentacion.
        """
        doc = self._documento_resultados()
        # Figura 8: VOI una vez (el de la primera familia que lo tenga) y el
        # candidato de cada familia, todos del mismo protocolo comparado.
        campos = {}
        for f in FAMILIAS:
            d = self._vm_comparado.get(f) or {}
            if "voi" in d and "voi" not in campos:
                campos["voi"] = d["voi"]
            if "spin" in d:
                campos[f] = d["spin"]
        self._ocupado(True, _("Informe: regenerando las estructuras y "
                              "dibujando las figuras…"))
        self.hilo = Trabajador(self._informe_etapa_datos, doc, carpeta,
                               self.VOI, self.VOI_spacing,
                               opciones_figuras(), campos or None,
                               self.VOI_contexto)
        self.hilo.listo.connect(self._informe_3d)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    @staticmethod
    def _informe_etapa_datos(doc, carpeta, VOI, spacing, opciones_3d=None,
                             campos_vm=None, contexto_voi=None):
        prep = informe_pub.preparar(doc, carpeta, VOI=VOI, spacing=spacing,
                                    opciones_3d=opciones_3d,
                                    campos_vm=campos_vm,
                                    contexto_voi=contexto_voi)
        return informe_pub.figuras_datos(prep)

    def _informe_3d(self, prep):
        """Renders 3D en el hilo de la interfaz y PDF en otro hilo.

        VTK necesita el contexto OpenGL del hilo que ya tiene los paneles 3D:
        renderizar fuera de el revienta el proceso de forma intermitente.
        """
        self.statusBar().showMessage(_("Informe: renderizando las estructuras "
                                       "en 3D…"))
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtWidgets.QApplication.processEvents()
        try:
            informe_pub.figuras_3d(prep)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.statusBar().showMessage(_("Informe: componiendo los PDF…"))
        # Se conserva la referencia al hilo anterior: acaba de emitir `listo` y
        # puede no haber salido aun de `run`. Soltarlo ahora destruiria un
        # QThread en marcha.
        self._hilo_informe_previo = self.hilo
        self.hilo = Trabajador(informe_pub.componer, prep)
        self.hilo.listo.connect(self._informe_listo)
        self.hilo.fallo.connect(self._error)
        self.hilo.start()

    def _informe_listo(self, r):
        auto = self._auto
        if auto is not None:
            t_inf = time.time() - auto.get("t_informe0", time.time())
            guardar_factor_tiempo("informe", t_inf,
                                  self._auto_estimado(auto, "informe", None))
            auto["hechos"].append({
                "clave": "informe", "etapa": _("Informe y figuras"),
                "familia": None, "estado": "hecha", "t_s": t_inf,
                "t_estimado_s": self._auto_estimado(auto, "informe", None)})
            ajustes = QtCore.QSettings("spinpy", "visor")
            try:
                n = int(ajustes.value("tiempos/ejecuciones", 0))
            except (TypeError, ValueError):
                n = 0
            ajustes.setValue("tiempos/ejecuciones", n + 1)
            self._auto_terminar()
        self._ocupado(False)
        idioma = "en" if str(QtCore.QSettings("spinpy", "visor").value(
            "idioma", "es")) == "en" else "es"
        c = r["resumen"]
        t = [_("<b>Informe para publicacion</b>"),
             _("{c} citables, {r} con reservas, {n} no citables.").format(
                 c=c["citable"], r=c["reservas"], n=c["no_citable"])]
        malos = [it for it in r["items"] if it["estado"] == "no_citable"]
        if malos:
            t.append(_("<b>No citables</b> (no deben aparecer como "
                       "resultado):"))
            t.append("<ul>" + "".join(
                f"<li>{html.escape(etiqueta_item(it, idioma))}: "
                f"{html.escape(texto_motivos(it, idioma))}</li>"
                for it in malos) + "</ul>")
        t.append(_("Cada resultado, con su motivo, esta en el informe."))
        n_fig = sum(1 for a in r["archivos"]
                    if Path(a).parent.name == informe_pub.CARPETA_FIGURAS)
        t.append(_("{n} archivos: informe en PDF y Markdown (ES/EN), {f} "
                   "figuras y el paquete de reproduccion.").format(
                       n=len(r["archivos"]), f=n_fig))
        for av in r["avisos"]:
            t.append("⚠ " + html.escape(texto_aviso(av, idioma)))
        t.append(_("Carpeta: {c}").format(c=html.escape(r["carpeta"])))
        if auto is not None:
            t.append(self._auto_resumen_html(auto))
        caja = QtWidgets.QMessageBox(self)
        caja.setWindowTitle(_("Informe para publicacion"))
        caja.setTextFormat(QtCore.Qt.RichText)
        caja.setText("<br>".join(t))
        abrir = caja.addButton(_("Abrir la carpeta"),
                               QtWidgets.QMessageBox.ActionRole)
        caja.addButton(QtWidgets.QMessageBox.Ok)
        caja.exec_()
        if caja.clickedButton() is abrir:
            QtGui.QDesktopServices.openUrl(
                QtCore.QUrl.fromLocalFile(r["carpeta"]))
        self.statusBar().showMessage(_("Informe escrito en {d}").format(
            d=r["carpeta"]))

    @staticmethod
    def _num(v):
        """Numero en formato espanol, o vacio si no lo hay."""
        if v is None or not np.isscalar(v) or not np.isfinite(v):
            return ""
        return repr(float(v)).replace(".", ",")

    @staticmethod
    def _txt(v):
        """Celda de CSV para un valor que puede ser texto o numero."""
        if v is None:
            return ""
        if isinstance(v, (bool, np.bool_, str)):
            return str(v)
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        return Visor._num(v)

    @staticmethod
    def _limpiar(m):
        """Deja un dict de morfometria en algo serializable a JSON.

        Los campos del MIL son arrays (direccion principal, autovectores, la
        propia nube de MIL) y `json` no sabe volcarlos. Se convierten en vez de
        descartarse: la direccion principal es justo lo que hace falta para
        comprobar si el ajuste oriento bien el candidato.
        """
        if not m:
            return None
        out = {}
        for k, v in m.items():
            if isinstance(v, np.ndarray):
                out[k] = v.tolist()
            elif isinstance(v, (np.floating, np.integer, np.bool_)):
                out[k] = v.item()
            else:
                out[k] = v
        return out

    def _error(self, tb):
        a = self._auto
        if a is not None and not a.get("en_informe"):
            # Una etapa que revienta no para la cadena: se anota y se sigue
            # con la siguiente, que suele no depender de ella (el tensor no
            # necesita el ensayo). El informe final lista la etapa con su
            # error, y lo que no llego a calcularse simplemente no aparece.
            a["error"] = tb
            self._ocupado(False, _("Error en la etapa; se sigue con la "
                                   "siguiente"))
            return
        if a is not None:
            self._auto_terminar()
        self._ocupado(False, "Error")
        QtWidgets.QMessageBox.critical(self, _("Error en el calculo"), tb)

    def _mostrar(self, dialogo):
        """Abre un dialogo de resultados, salvo durante el informe automatico.

        En la cadena automatica no hay nadie delante para cerrarlo, y un
        `exec_()` la pararia hasta que volviera alguien. Lo que el dialogo
        ensena ya esta en `_res` y acaba en el informe. `dialogo` puede ser el
        dialogo o una funcion que lo construye; la segunda forma es para los
        que llevan un panel 3D, que no merece la pena crear para tirarlo.
        """
        if self._auto is not None:
            if not callable(dialogo):
                dialogo.deleteLater()
            return
        (dialogo() if callable(dialogo) else dialogo).exec_()

    def _avisar(self, titulo, texto):
        """Aviso de que una accion no se puede lanzar (falta algo, resolucion
        insuficiente...). Durante el informe automatico no se abre: se anota,
        y la etapa sale como «omitida» en el resumen con este motivo."""
        if self._auto is not None:
            self._auto["avisos"].append(f"{titulo}: {texto}")
            return
        QtWidgets.QMessageBox.information(self, titulo, texto)

    # -- informe automatico ------------------------------------------------

    def informe_auto(self):
        """Todo el recorrido sobre el VOI cargado y, al final, el informe.

        Encadena los MISMOS metodos que los botones del panel, una etapa
        detras de otra, en vez de tener un camino de calculo propio: un
        informe automatico que midiera distinto que el manual describiria
        una sesion que nadie puede repetir a mano. Lo unico que cambia es que
        los dialogos de resultados no se abren y los errores se anotan.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None:
            QtWidgets.QMessageBox.information(
                self, _("Falta el VOI"),
                _("El informe automatico parte del VOI de referencia. "
                  "Cargalo primero."))
            return
        d = DialogoInformeAuto(self)
        if d.exec_() != QtWidgets.QDialog.Accepted:
            return
        op = d.opciones()
        carpeta = Path(op["carpeta"])
        previos = [n for n in informe_pub.ARCHIVOS.values()
                   if (carpeta / n).exists()]
        if previos and QtWidgets.QMessageBox.question(
                self, _("La carpeta ya tiene un informe"),
                _("Se sobrescribiran {n} archivo(s) de un informe anterior "
                  "en esa carpeta. ¿Continuar?").format(n=len(previos))
        ) != QtWidgets.QMessageBox.Yes:
            return

        pasos = self._auto_pasos(op)
        self._auto = {"pasos": pasos, "total": len(pasos), "hechos": [],
                      "avisos": [], "error": None, "detener": False,
                      "t0": time.time(), "carpeta": carpeta, "opciones": op,
                      "inicio": time.strftime("%Y-%m-%d %H:%M:%S")}
        self.lab_auto.setVisible(True)
        self.btn_auto_detener.setEnabled(True)
        self.btn_auto_detener.setVisible(True)
        self._ocupado(True, _("Informe automatico: empezando…"))
        QtCore.QTimer.singleShot(0, self._auto_siguiente)

    def _auto_pasos(self, op):
        """Lista de etapas `(clave, rotulo, familia, accion)`, en orden.

        El orden es el del trabajo: sin ajuste no hay candidato que medir, y
        el analisis comparado va despues del ensayo propio para que los campos
        que quedan pintados al final sean los del protocolo del articulo.
        """
        fams = self._familias()
        pasos = []

        def por_familia(clave, rotulo, accion):
            for f in fams:
                pasos.append((clave, rotulo, f, accion))

        if op["ajuste"]:
            por_familia("ajuste", _("Mejor ajuste al VOI"),
                        self._ajuste_mejor_una)
        if op["morfometria"]:
            pasos.append(("morfometria", _("Morfometria completa"), None,
                          lambda: self.medir(sobre_mascara=True)))
        if op["dispersion"]:
            por_familia("dispersion", _("Dispersion entre semillas"),
                        self._dispersion_una)
        if op["elastico"]:
            por_familia("elastico", _("Tensor elastico"),
                        self._calcular_elastico_una)
        if op["ensayo"]:
            por_familia("ensayo", _("Ensayo de compresion y fallo"),
                        self._ensayo_fe_una)
        if op["convergencia"]:
            pasos.append(("convergencia", _("Convergencia de malla (VOI)"),
                          None, self.convergencia_fe))
        if op["comparado"]:
            por_familia("comparado", _("Analisis comparado"),
                        lambda: self._analisis_comparado_una(ajustar=False))
        if op["perdida"]:
            pasos.append(("perdida", _("Simulacion de perdida osea"),
                          fams[0], self.simular_perdida_gui))
        if op["fallo"]:
            pasos.append(("fallo", _("Fallo progresivo"), fams[0],
                          self.fallo_progresivo_gui))
        return pasos

    def _auto_siguiente(self):
        """Lanza la etapa siguiente y se engancha al final de su hilo.

        Mismo mecanismo que `_en_cada_familia`: el `*_listo` de cada etapa lee
        la familia activa, asi que la siguiente no se lanza hasta que el hilo
        de la anterior ha terminado y su resultado ya esta colocado.
        """
        a = self._auto
        if a is None:
            return
        if a["detener"]:
            self._auto_detenido()
            return
        if not a["pasos"]:
            if a.get("tipo") == "fem":
                self._fem_auto_fin()
            else:
                self._auto_informe()
            return
        clave, rotulo, fam, accion = a["pasos"].pop(0)
        if fam is not None:
            self._activar(fam)
            if a.get("tipo") != "fem":
                rotulo = f"{rotulo} · {self._etq_fam(fam)}"
        i = a["total"] - len(a["pasos"])
        from spinpy import tiempos
        if a.get("tipo") == "fem":
            # Una estimacion por etapa, en el orden de las etapas: varias
            # comparten clave y familia y `_auto_estimado` no las distingue.
            t_est = a["estimados"].pop(0) if a["estimados"] else None
            resta = (t_est or 0.0) + sum(x or 0.0 for x in a["estimados"])
            self.lab_auto.setText(
                _("FEM automático {i}/{n}: {etapa}").format(
                    i=i, n=a["total"], etapa=rotulo)
                + "  ·  " + _("quedan {t}").format(t=tiempos.texto(resta)))
        else:
            t_est = self._auto_estimado(a, clave, fam)
            # Lo que queda: esta etapa, las pendientes y el informe final.
            resta = (t_est or 0.0) + (self._auto_estimado(a, "informe", None)
                                      or 0.0) + sum(
                self._auto_estimado(a, p[0], p[2]) or 0.0 for p in a["pasos"])
            self.lab_auto.setText(
                _("Informe automatico {i}/{n}: {etapa}").format(
                    i=i, n=a["total"], etapa=rotulo)
                + "  ·  " + _("quedan {t}").format(t=tiempos.texto(resta)))
        a["error"] = None
        fila = {"clave": clave, "etapa": rotulo, "familia": fam,
                "t_estimado_s": t_est}
        t0 = time.time()
        previo = self.hilo
        try:
            accion()
        except Exception:
            a["error"] = traceback.format_exc()
        h = self.hilo
        lanzado = h is not None and h is not previo

        def cerrar():
            if "estado" in fila or self._auto is not a:
                return
            err = a.get("error")
            fila["estado"] = ("error" if err
                              else "hecha" if lanzado else "omitida")
            if err:
                # La ultima linea de la traza basta para el resumen; la traza
                # entera va al JSON de la sesion.
                fila["error"] = err
            fila["t_s"] = time.time() - t0
            if fila["estado"] == "hecha":
                guardar_factor_tiempo(clave, fila["t_s"], t_est)
            a["hechos"].append(fila)
            QtCore.QTimer.singleShot(0, self._auto_siguiente)

        if lanzado:
            h.finished.connect(cerrar)
            # Si el hilo ya termino antes de conectar, `finished` se perdio.
            if h.isFinished():
                QtCore.QTimer.singleShot(0, cerrar)
        else:
            cerrar()

    def _auto_registro(self, a):
        """Lo que queda en la sesion sobre como se produjo el informe."""
        op = {k: v for k, v in a["opciones"].items() if k != "carpeta"}
        return {"inicio": a["inicio"],
                "duracion_s": time.time() - a["t0"],
                "familias": self._familias(),
                "opciones": op,
                "etapas": [dict(f) for f in a["hechos"]],
                "pendientes": [p[1] for p in a["pasos"]],
                "avisos": list(a["avisos"])}

    @staticmethod
    def _auto_estimado(a, clave, fam):
        """Segundos estimados de una etapa, o None si no hay estimacion. Las
        claves de familia se guardaron como texto ("" = sin familia) para que
        la estimacion quepa tal cual en el JSON de la sesion."""
        d = (a["opciones"].get("estimacion") or {}).get(clave) or {}
        return d.get("" if fam is None else fam)

    def _auto_informe(self):
        a = self._auto
        a["t_informe0"] = time.time()
        self._activar(self._familias()[0])
        self._redibujar()
        # Antes de escribir: asi `resultados_sesion.json` lleva el registro
        # de las etapas junto a los numeros que produjeron.
        self._res["informe_auto"] = self._auto_registro(a)
        a["en_informe"] = True
        self.lab_auto.setText(_("Informe automatico: escribiendo el "
                                "informe…"))
        self.btn_auto_detener.setEnabled(False)
        self._lanzar_informe(a["carpeta"])

    def _auto_detener(self):
        # FEBio se detiene de verdad: el evento lo consulta `febio.correr`,
        # que mata el proceso; lo terminado se guarda.
        if self._febio_evento is not None:
            self._febio_evento.set()
        if self._auto is None:
            self.btn_auto_detener.setEnabled(False)
            return
        self._auto["detener"] = True
        self.btn_auto_detener.setEnabled(False)
        self.lab_auto.setText(self.lab_auto.text() + "  ·  "
                              + _("se detiene al acabar esta etapa"))

    def _auto_detenido(self):
        a = self._auto
        if a.get("tipo") == "fem":
            self._fem_auto_fin(detenido=True)
            return
        self._activar(self._familias()[0])
        self._redibujar()
        self._res["informe_auto"] = self._auto_registro(a)
        self._auto_terminar()
        self._ocupado(False, _("Informe automatico detenido. Lo calculado "
                               "sigue en la sesion."))
        caja = QtWidgets.QMessageBox(self)
        caja.setWindowTitle(_("Informe automatico detenido"))
        caja.setTextFormat(QtCore.Qt.RichText)
        caja.setText(_("No se ha escrito el informe. Lo calculado hasta aqui "
                       "sigue en la sesion: «Informe para publicación…» lo "
                       "describe tal como esta.") + "<br>"
                     + self._auto_resumen_html(a))
        caja.exec_()

    def _auto_terminar(self):
        self._auto = None
        self.lab_auto.setVisible(False)
        self.btn_auto_detener.setVisible(False)

    @staticmethod
    def _duracion(s):
        s = int(round(float(s)))
        if s < 60:
            return f"{s} s"
        if s < 3600:
            return f"{s // 60} min {s % 60:02d} s"
        return f"{s // 3600} h {(s % 3600) // 60:02d} min"

    def _auto_resumen_html(self, a):
        estados = {"hecha": _("hecha"), "omitida": _("omitida"),
                   "error": _("con error")}
        colores = {"hecha": "#1a7f37", "omitida": "#9a6700",
                   "error": "#b62324"}
        from spinpy import tiempos
        est_total = sum(sum(d.values()) for d in
                        (a["opciones"].get("estimacion") or {}).values())
        t = ["<br><b>" + _("Recorrido automatico") + "</b> — "
             + _("{t} en total").format(t=self._duracion(time.time()
                                                         - a["t0"]))
             + (" (" + _("estimado {t}").format(t=tiempos.texto(est_total))
                + ")" if est_total else ""),
             "<table cellpadding=3 cellspacing=0 border=1 "
             "style='border-collapse:collapse'>",
             "<tr><th>" + _("etapa") + "</th><th>" + _("estado")
             + "</th><th>" + _("real") + "</th><th>" + _("estimado")
             + "</th><th></th></tr>"]
        for f in a["hechos"]:
            nota = ""
            if f.get("error"):
                ultima = [x for x in str(f["error"]).strip().splitlines()
                          if x.strip()]
                nota = html.escape(ultima[-1][:160] if ultima else "")
            est = f.get("t_estimado_s")
            t.append(f"<tr><td>{html.escape(f['etapa'])}</td>"
                     f"<td style='color:{colores[f['estado']]}'>"
                     f"{estados[f['estado']]}</td>"
                     f"<td align=right>{self._duracion(f['t_s'])}</td>"
                     f"<td align=right style='color:#666'>"
                     f"{tiempos.texto(est) if est else '—'}</td>"
                     f"<td>{nota}</td></tr>")
        no_lanzada = _("no lanzada")
        for p in a["pasos"]:
            t.append(f"<tr><td>{html.escape(p[1])}</td>"
                     f"<td colspan=4 style='color:#888'>{no_lanzada}</td>"
                     "</tr>")
        t.append("</table>")
        for av in a["avisos"]:
            t.append("<br>⚠ " + html.escape(av))
        return "".join(t)

    # -- FEBio -------------------------------------------------------------
    #
    # Dos puertas —«Analizar con FEBio…» del panel y «FEM automático
    # (FEBio)…» de la barra— y un solo camino: las dos construyen «tareas»
    # (estructura x protocolo x malla x ejes) y las pasan a `_febio_una`, que
    # las resuelve en un hilo con `febio.analizar` / `febio.homogeneizar`. El
    # automatico las recorre con `_auto_siguiente`, una etapa por tarea.

    def _carpeta_febio_def(self):
        nombre = Path(str(self.VOI_nombre or "sesion")).stem
        return DATOS / f"FEM_{nombre}_{time.strftime('%Y-%m-%d')}"

    def _mascara_estructura(self, est):
        if est == "voi":
            return self.VOI, self.VOI_spacing
        BW = self._de(est, "BW")
        return (None, None) if BW is None else (BW, self._spacing(BW.shape[0]))

    def _procedencia_estructura(self, est):
        """La del ajuste si lo hay (es la estructura ensayada); si no, la de
        los controles. El VOI no tiene generador."""
        if est == "voi":
            return procedencia.bloque("voi")
        rec = self._res.get("ajuste" if est == "spinodoide" else "ajuste_dual")
        if isinstance(rec, dict) and rec.get("procedencia"):
            return rec["procedencia"]
        return self._procedencia_controles()[est]

    def _febio_una(self, tareas, agregar=True):
        """Resuelve `tareas` en FEBio, en un hilo, una detras de otra."""
        if self.hilo is not None and self.hilo.isRunning():
            return
        import threading
        trabajos = []
        for t in tareas:
            BW, sp = self._mascara_estructura(t["estructura"])
            if BW is None:
                self._avisar(_("Falta la estructura"),
                             _("{e} no está generado; se omite.").format(
                                 e=t["estructura"]))
                continue
            trabajos.append((dict(t), BW, sp,
                             self._procedencia_estructura(t["estructura"])))
        if not trabajos:
            return
        ev = threading.Event()
        self._febio_evento = ev
        self.btn_auto_detener.setVisible(True)
        self.btn_auto_detener.setEnabled(True)
        self._ocupado(True, _("FEBio: preparando…"))
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)
        hilo._fn = lambda: Visor._febio_tarea(trabajos, ev, hilo.informar)
        hilo.avance.connect(self._febio_avance)
        hilo.listo.connect(lambda out, ag=agregar: self._febio_listo(out, ag))
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _febio_avance(self, i, n, etapa):
        if self.barra.maximum() != n:
            self.barra.setRange(0, max(1, n))
        self.barra.setValue(i)
        self.statusBar().showMessage("FEBio · " + etapa)

    @staticmethod
    def _febio_tarea(trabajos, ev, informar):
        """Sin Qt: corre en el hilo de trabajo."""
        import tempfile

        from spinpy import febio
        out = {"registros": [], "mapas": [], "cancelado": False}
        for t, BW, sp, proc in trabajos:
            p, est, m = t["protocolo"], t["estructura"], t["malla"]
            carpeta = (Path(t["carpeta"]) if t.get("carpeta")
                       else Path(tempfile.mkdtemp(prefix="febio_")))
            for eje in t["ejes"]:
                rot = f"{est} · {p['nombre']} · {m} · {'XYZ'[eje]}"

                def pr(i, n, et, _r=rot):
                    informar(i, n, f"{_r} · {et}")
                informar(0, 1, rot)
                try:
                    if p["tipo"] == "homogeneizacion":
                        r = febio.homogeneizar(
                            BW, sp, p, malla=m, carpeta=carpeta, n=t["n"],
                            exe=t["exe"], hilos=t["hilos"], cancelar=ev,
                            progreso=pr, opciones_malla=t["opciones_malla"],
                            etiqueta=est, conservar=t["conservar"])
                    else:
                        r = febio.analizar(
                            BW, sp, p, malla=m, analisis=t["analisis"],
                            eje=eje, carpeta=carpeta, n=t["n"], exe=t["exe"],
                            hilos=t["hilos"], cancelar=ev,
                            material=t["material"], pasos=t["pasos"],
                            progreso=pr, conservar=t["conservar"],
                            comparar_app=t["comparar_app"],
                            opciones_malla=t["opciones_malla"], etiqueta=est,
                            conv_malla=t["conv_malla"])
                except febio.Cancelado:
                    out["cancelado"] = True
                    return out
                except (febio.ErrorMalla, febio.ErrorFEBio, ValueError,
                        RuntimeError) as e:
                    # Una etapa que falla se registra y no para las demas.
                    r = {"nombre": p["nombre"], "protocolo": p["clave"],
                         "modificado": p["modificado"], "malla": m,
                         "lineal": None, "fallos": [{"corrida": "malla",
                                                      "msg": str(e)}],
                         "carpeta": str(carpeta / f"{est}_fallo")}
                r.update(estructura=est, estructura_codigo=est,
                         tipo=r.get("tipo", p["tipo"]), procedencia=proc,
                         eje_nombre=r.get("eje_nombre", "XYZ"[eje]))
                cl, mm = r.get("_campos_lineal"), r.get("_malla")
                if cl is not None and mm is not None:
                    out["mapas"].append({
                        "titulo": rot, "nodos": mm["nodos"],
                        "elems": mm["elems"],
                        "vm": febio.von_mises(cl["sigma"]) / 1e6})
                out["registros"].append(febio.registro_json(r))
        return out

    def _febio_registrar(self, registros):
        """Registros en la sesion: `_res["febio"]["registros"]`, donde los
        recoge el informe de publicacion. Una repeticion reemplaza a la
        anterior de la misma estructura, protocolo, malla y eje."""
        def clave(r):
            return (r.get("estructura"), r.get("nombre"), r.get("malla"),
                    r.get("eje_nombre"), r.get("tipo"))
        rec = self._res.setdefault("febio", {"registros": []})
        nuevas = {clave(r) for r in registros}
        rec["registros"] = [r for r in rec["registros"]
                            if clave(r) not in nuevas] + list(registros)

    @staticmethod
    def _febio_guardar(carpeta, registros):
        """`resultados_fem.json` y `tabla_fem.csv` en la carpeta, sumando a
        lo que ya hubiera de corridas anteriores."""
        carpeta = Path(carpeta)
        carpeta.mkdir(parents=True, exist_ok=True)
        ruta = carpeta / "resultados_fem.json"
        previos = []
        if ruta.exists():
            try:
                previos = json.loads(ruta.read_text(encoding="utf-8")).get(
                    "registros", [])
            except (ValueError, OSError):
                previos = []

        def clave(r):
            return (r.get("estructura"), r.get("nombre"), r.get("malla"),
                    r.get("eje_nombre"), r.get("tipo"))
        nuevas = {clave(r) for r in registros}
        todos = [r for r in previos if clave(r) not in nuevas] + list(registros)
        doc = {"formato": "spinpy/resultados_fem", "version": 1,
               "spinpy": procedencia.version_spinpy(),
               "escrito": time.strftime("%Y-%m-%d %H:%M:%S"),
               "registros": todos}
        ruta.write_text(json.dumps(doc, ensure_ascii=False, indent=1,
                                   default=float), encoding="utf-8")
        escribir_csv(todos, carpeta / "tabla_fem.csv")
        return todos

    def _febio_listo(self, out, agregar=True):
        self._febio_evento = None
        a = self._auto
        if a is None:
            self.btn_auto_detener.setVisible(False)
        self._ocupado(False)
        regs = out["registros"]
        self._febio_mapas = (self._febio_mapas if a is not None else []) \
            + out["mapas"]
        carpetas = {r.get("carpeta") for r in regs if r.get("carpeta")}
        # La carpeta de la tarea es la raiz; cada registro guarda la suya.
        raiz = None
        for c in carpetas:
            raiz = Path(c).parent
            break
        if raiz is not None:
            self._febio_guardar(raiz, regs)
        if agregar:
            self._febio_registrar(regs)
        if a is not None:
            a.setdefault("registros", []).extend(regs)
            if out["cancelado"]:
                a["detener"] = True
            return
        if out["cancelado"]:
            self.statusBar().showMessage(_("FEBio detenido; lo terminado se "
                                           "conserva."))
        self._mostrar(lambda: DialogoResultadosFEBio(
            self, regs, raiz, pendiente=not agregar, mapas=self._febio_mapas))

    def analizar_febio(self):
        """«Analizar con FEBio…»: la estructura activa y el VOI."""
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None and self.BW_vista is None:
            return
        d = DialogoFEBio(self, str(self._carpeta_febio_def()))
        if d.exec_() != QtWidgets.QDialog.Accepted:
            return
        op = d.opciones()
        p = op["protocolo"]

        def tareas(est):
            return [{"estructura": est, "protocolo": p, "malla": m,
                     "analisis": op["analisis"] or ["lineal"],
                     "ejes": op["ejes"] if p["tipo"] == "compresion" else [2],
                     "n": op["n_hex"] if m == "hex8" else op["n_tet"],
                     "opciones_malla": op["opciones_malla"],
                     "conv_malla": op["conv_malla"] and m == "tet10",
                     "material": op["material"], "pasos": op["pasos"],
                     "hilos": op["hilos"], "exe": op["exe"],
                     "carpeta": op["carpeta"], "conservar": True,
                     "comparar_app": True}
                    for m in op["mallas"]]

        voi_hecho = [False]

        def una():
            ts = []
            if self.VOI is not None and not voi_hecho[0]:
                ts += tareas("voi")
                voi_hecho[0] = True
            if self.BW_vista is not None and self._confirmar_orientacion():
                ts += tareas(self._fam)
            if ts:
                self._febio_una(ts)

        self._en_cada_familia(una)

    def fem_auto(self):
        """«FEM automático (FEBio)…»: todo desatendido, como el informe
        automatico, encadenando las mismas `_febio_una` que el panel."""
        if self.hilo is not None and self.hilo.isRunning():
            return
        if self.VOI is None and all(self._de(f, "BW") is None
                                    for f in FAMILIAS):
            QtWidgets.QMessageBox.information(
                self, _("Nada que analizar"),
                _("Carga el VOI o genera una estructura primero."))
            return
        d = DialogoFEMAuto(self, str(self._carpeta_febio_def()))
        if d.exec_() != QtWidgets.QDialog.Accepted:
            return
        op = d.opciones()
        from dialogo_febio import etiqueta_estructura, etiqueta_malla
        pasos, estimados = [], []
        for est in op["por_ajustar"]:
            pasos.append(("ajuste", _("Mejor ajuste al VOI") + " · "
                          + self._etq_fam(est), est, self._ajuste_mejor_una))
            estimados.append(1100.0)
        for t, s in zip(op["tareas"], op["estimados"]):
            fam = None if t["estructura"] == "voi" else t["estructura"]
            rot = " · ".join((etiqueta_estructura(t["estructura"]),
                              t["protocolo"]["nombre"],
                              etiqueta_malla(t["malla"])))
            pasos.append((f"febio_{t['malla']}", rot, fam,
                          lambda t=t: self._febio_una(
                              [t], agregar=op["sesion"])))
            estimados.append(s)
        self._auto = {"tipo": "fem", "pasos": pasos, "total": len(pasos),
                      "hechos": [], "avisos": [], "error": None,
                      "detener": False, "t0": time.time(),
                      "carpeta": Path(op["carpeta"]),
                      "opciones": {k: v for k, v in op.items()
                                   if k not in ("tareas", "estimados")},
                      "inicio": time.strftime("%Y-%m-%d %H:%M:%S"),
                      "estimados": estimados, "registros": []}
        self._febio_mapas = []
        self.lab_auto.setVisible(True)
        self.btn_auto_detener.setEnabled(True)
        self.btn_auto_detener.setVisible(True)
        self._ocupado(True, _("FEM automático: empezando…"))
        QtCore.QTimer.singleShot(0, self._auto_siguiente)

    def _fem_auto_fin(self, detenido=False):
        a = self._auto
        op = a["opciones"]
        reg = {"inicio": a["inicio"], "duracion_s": time.time() - a["t0"],
               "configuracion": op.get("configuracion"),
               "etapas": [dict(f) for f in a["hechos"]],
               "pendientes": [p[1] for p in a["pasos"]],
               "avisos": list(a["avisos"]), "detenido": bool(detenido),
               "carpeta": str(a["carpeta"])}
        if op.get("sesion"):
            self._res["fem_auto"] = reg
        regs = a.get("registros", [])
        carpeta = a["carpeta"]
        figura = None
        try:
            carpeta.mkdir(parents=True, exist_ok=True)
            (carpeta / "fem_auto.json").write_text(
                json.dumps(reg, ensure_ascii=False, indent=1, default=str),
                encoding="utf-8")
            if op.get("figura") and regs:
                figura = figura_validacion(
                    regs, carpeta / "figura_validacion_fem.png")
        except Exception:
            a["avisos"].append(traceback.format_exc().splitlines()[-1])
        mapas = self._febio_mapas if op.get("mapas") else []
        resumen = self._auto_resumen_html(a)
        self._auto_terminar()
        self._febio_evento = None
        self._ocupado(False, _("FEM automático detenido; lo terminado se "
                               "conserva.") if detenido
                      else _("FEM automático terminado."))
        caja = QtWidgets.QMessageBox(self)
        caja.setWindowTitle(_("FEM automático (FEBio)"))
        caja.setTextFormat(QtCore.Qt.RichText)
        caja.setText(resumen)
        caja.exec_()
        DialogoResultadosFEBio(self, regs, carpeta,
                               pendiente=not op.get("sesion"), mapas=mapas,
                               figura=figura).exec_()

    def _ajuste_mejor_una(self):
        """Etapa de ajuste del informe automatico, para la familia activa.

        El de menor error entre los metodos que comparten definicion de error
        (`mejor_ajuste`), igual que «Ajustar antes» del analisis comparado. El
        desempate mecanico queda fuera por lo mismo que alli: su cifra suma
        dos terminos y no compite en la misma cola.
        """
        if self.hilo is not None and self.hilo.isRunning():
            return
        VOI, sp = self.VOI, self.VOI_spacing
        reg = REGISTROS[self._fam]
        claves = [c for c in METODOS_COMPARABLES if c in reg]
        m_voi = self.m_voi
        semilla = int(self.spin_semilla.value())
        nw = int(self.sl["nw"].valor())
        self._ocupado(True, _("Ajustando con los metodos comparables…"),
                      determinada=True,
                      total=sum(reg[c]["n_eval"] for c in claves))
        self._t0 = time.time()
        hilo = Trabajador(lambda: None)
        hilo._fn = lambda: mejor_ajuste(VOI, sp, reg, claves, m_voi=m_voi,
                                        semilla=semilla, num_waves=nw,
                                        informar=hilo.informar)
        hilo.avance.connect(self._avance)
        hilo.listo.connect(self._mejor_ajuste_listo)
        hilo.fallo.connect(self._error)
        self.hilo = hilo
        hilo.start()

    def _mejor_ajuste_listo(self, out):
        ganador, hechos = out
        self._ajuste_listo(ganador)
        # Misma forma que deja «Ajustar con TODOS los metodos», para que un
        # JSON no tenga que saber por que boton se eligio el ajuste.
        self._res[self._clave_res("ajuste_comparado")] = {
            "elegido": ganador["metodo"],
            "criterio": "menor error entre metodos comparables",
            "comparativa": [
                {k: v for k, v in fl.items() if np.isscalar(v) or v is None}
                for fl in comparar_resultados(list(hechos.values()),
                                              ganador["metricas_voi"])],
        }

    def closeEvent(self, ev):
        for v in (self.vis_voi, *self._vis.values()):
            try:
                v.close()
            except Exception:
                pass
        super().closeEvent(ev)


def _replicar(clave, rapido=False):
    """Ejecuta `Test/replicar_<clave>.py` sin abrir ninguna ventana.

    El guion se carga por RUTA y no por `import`: empaquetado, `Test` no esta
    en el camino de importacion de Python, esta en los datos de PyInstaller.
    Ademas se anade su carpeta a sys.path, porque las replicas se reutilizan
    entre si -la de Zheng usa la superficie elastica ya verificada de la de
    Kumar- y sin eso ese import falla solo dentro del ejecutable.
    """
    import importlib.util
    try:
        from dialogo_validacion import replica_por_clave, ruta_test
    except Exception:
        traceback.print_exc()
        return 2
    ruta = ruta_test()
    rep = replica_por_clave(clave)
    if ruta is None or rep is None:
        print("No se encuentra la carpeta Test o la replica pedida.")
        return 2
    guion = ruta / rep["guion"]
    if not guion.is_file():
        print(f"No se encuentra {guion.name}.")
        return 2
    if str(ruta) not in sys.path:
        sys.path.insert(0, str(ruta))
    spec = importlib.util.spec_from_file_location(guion.stem, guion)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    argv = list(sys.argv)
    sys.argv = [guion.name] + (["--rapido"] if rapido else [])
    try:
        return mod.main()
    finally:
        sys.argv = argv


def _cerrar_obreros():
    """Cierra los procesos obreros que joblib deja vivos para reutilizarlos.

    POR QUE HACE FALTA, Y NO ES COSMETICO
      El backend `loky` mantiene sus obreros arrancados despues de que
      `Parallel` devuelva, para no pagar el arranque la proxima vez. Se quedan
      unos minutos y HEREDAN LA SALIDA ESTANDAR del proceso padre.

      Eso tiene dos consecuencias feas. La primera se vio compilando: el guion
      de publicacion canaliza la salida de `--autocomprobacion`, el padre
      termina, pero los obreros siguen con el extremo de escritura abierto, la
      tuberia nunca ve el final y PowerShell se queda esperando. La
      compilacion se colgo 27 minutos en el paso 3 con un obrero girando al
      100 % de una CPU. La segunda es peor de cara al usuario: la seccion de
      lote de la aplicacion dejaria obreros vivos despues de terminar.

      Se cierran a mano. Si el backend no es loky -o joblib cambia de sitio
      esta funcion- no pasa nada: no habia nada que cerrar.
    """
    try:
        from spinpy.lote import cerrar_obreros
        cerrar_obreros()
    except Exception:
        pass


def autocomprobacion():
    """Ejercita sin interfaz cada camino de calculo y escribe un informe.

    PARA QUE SIRVE
      La aplicacion empaquetada arranca aunque le falte una biblioteca: los
      modulos pesados -tetgen, pymeshfix, pyamg- no se importan hasta que el
      usuario pulsa el boton correspondiente. Un fallo de empaquetado en el
      mallador no se manifiesta al abrir la ventana, sino media hora despues,
      en mitad de una exportacion, en forma de traza que quien la recibe no
      sabe leer. Esto lo destapa en unos segundos y en la maquina donde
      importa: la del usuario, no la del que compilo.

      El informe se guarda ademas en un archivo, para poder enviarlo.

    Uso:  spinpy.exe --autocomprobacion       (o  python visor.py --auto...)
    Devuelve 0 si todo pasa, 1 si algo falla.
    """
    import tempfile

    lineas = []
    fallos = 0

    def prueba(nombre, fn):
        nonlocal fallos
        t0 = time.time()
        try:
            detalle = fn() or ""
            lineas.append(f"  [OK   ] {nombre:38s} {time.time()-t0:6.2f} s  "
                          f"{detalle}")
        except Exception as e:
            fallos += 1
            lineas.append(f"  [FALLA] {nombre:38s} {time.time()-t0:6.2f} s  "
                          f"{type(e).__name__}: {e}")

    lineas.append("=" * 78)
    lineas.append(" spinpy - autocomprobacion")
    lineas.append(f" {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lineas.append(f" empaquetado: {bool(getattr(sys, 'frozen', False))}")
    lineas.append(f" ejecutable : {sys.executable}")
    lineas.append(f" datos      : {DATOS}")
    lineas.append("=" * 78)
    lineas.append("")
    lineas.append(" MODULOS")

    for mod in ("numpy", "scipy", "skimage", "pyvista", "vtkmodules.all",
                "pyamg", "tetgen", "pymeshfix", "pandas", "joblib",
                "matplotlib", "matplotlib.backends.backend_qt5agg", "PyQt5",
                "pyvistaqt"):
        def _imp(m=mod):
            import importlib
            x = importlib.import_module(m)
            return getattr(x, "__version__", "")
        prueba(f"import {mod}", _imp)

    lineas.append("")
    lineas.append(" CALCULO  (mallas pequenas: esto NO valida los numeros,")
    lineas.append("           solo que cada camino se puede recorrer)")

    estado = {}

    def _generar():
        BW, _, info = generar_mascara(resolution=32, wave_number=12 * np.pi,
                                      num_waves=400, thetas=[15., 15., 15.],
                                      rho=0.3, seed=1)
        estado["BW"] = BW
        estado["sp"] = np.full(3, 1.0 / 32)
        return f"BV/TV pedido 0.30, obtenido {info['rho_obtenida']:.3f}"
    prueba("generar spinodoide 32^3", _generar)

    def _dual():
        # Delaunay y el arbol k-d viven en scipy.spatial, que envuelve Qhull:
        # una DLL que el empaquetado puede dejarse sin que nada lo avise hasta
        # que alguien elige la familia dual-lattice.
        from spinpy.dual_lattice import generar_dual_lattice
        _BW, _, info = generar_dual_lattice(32, 4.0, 0.3, seed=1)
        return (f"BV/TV pedido 0.30, obtenido {info['rho_obtenida']:.3f}, "
                f"{info['n_nodos']} nodos")
    prueba("generar dual-lattice 32^3 (Qhull)", _dual)

    def _morfo():
        m = morfometria(estado["BW"], estado["sp"])
        estado["m"] = m
        return (f"BV/TV {m['BVTV']:.4f}  Tb.Th {m['TbTh']:.5f} mm  "
                f"DA {m['DA']:.3f}")
    prueba("morfometria + tensor MIL", _morfo)

    def _malla_sup():
        m = morfometria_malla(estado["BW"], estado["sp"], do_mil=False)
        return f"BS {m['BS']:.4f} mm2, estanca={m['estanca']}"
    prueba("superficie cerrada (pymeshfix)", _malla_sup)

    def _homog():
        BW16, sp16 = remuestrear_bw(estado["BW"], estado["sp"], 16)
        estado["BW16"], estado["sp16"] = BW16, sp16
        C, _inf = homogeneizar(BW16, E_s=E_S_PA, nu_s=NU_S,
                               vox_size=float(sp16[0]))
        e = constantes_ingenieria(C)
        estado["E"] = e
        return (f"Ex {e['Ex']/1e6:.1f}  Ey {e['Ey']/1e6:.1f}  "
                f"Ez {e['Ez']/1e6:.1f} MPa")
    prueba("homogeneizacion 16^3 (pyamg)", _homog)

    def _ensayo():
        BW16, sp16 = remuestrear_bw(estado["BW"], estado["sp"], 16)
        r = ensayo_compresion(BW16, sp16, E_s=E_S_PA,
                              nu_s=NU_S, sigma0=1e6, apoyo="deslizante")
        if not r["ok"]:
            raise RuntimeError(r["msg"])
        return f"E_app {r['E_app']/1e6:.1f} MPa"
    prueba("ensayo de compresion 16^3", _ensayo)

    def _simulacion():
        from spinpy.simulacion import fallo_progresivo, simular_perdida
        r = simular_perdida(estado["BW"], estado["sp"], "trabeculas_finas",
                            pasos=1, perdida_paso=0.05, n_trabajo=32,
                            n_mec=12, guardar_mascaras=False)
        a, b = r["pasos"][0], r["pasos"][-1]
        rf = fallo_progresivo(estado["BW"], estado["sp"], pasos=1, n_mec=12,
                              guardar_mascaras=False)
        return (f"BV/TV {a['BVTV']:.3f} -> {b['BVTV']:.3f}, E/E0 "
                f"{b['E_rel']:.3f}; fallo {len(rf['pasos'])} pasos")
    prueba("simulacion de perdida y fallo 12^3", _simulacion)

    def _hex():
        n, e, _inf = malla_hex(estado["BW"], estado["sp"])
        estado["hex"] = (n, e)
        return f"{len(e)} hexaedros, {len(n)} nodos"
    prueba("malla hexaedrica", _hex)

    def _tet():
        n, e, sup, inf = malla_tet10(estado["BW"], estado["sp"])
        estado["tet"] = (n, e, sup)
        return f"{len(e)} TET10, {inf.get('vol_pct_MC', float('nan')):.1f} % del volumen"
    prueba("malla TET10 (tetgen)", _tet)

    def _escribir():
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            n, e = estado["hex"]
            escribir_vtu(n, e, d / "p.vtu")
            escribir_abaqus(n, e, d / "p.inp", E_s=E_S_PA, nu_s=NU_S)
            escribir_apdl(n, e, d / "p.dat", E_s=E_S_PA, nu_s=NU_S)
            escribir_febio(n, e, d / "p.feb", E_s=E_S_PA, nu_s=NU_S)
            escribir_stl(estado["tet"][2], d / "p.stl")
            tam = {f.suffix: f.stat().st_size for f in d.iterdir()}
        return " ".join(f"{k}{v//1024}kB" for k, v in sorted(tam.items()))
    prueba("escribir VTU / INP / DAT / FEB / STL", _escribir)

    def _febio():
        # En el ejecutable FEBio tiene que ser la copia EMPAQUETADA: si se
        # usara una instalacion de FEBio Studio de la maquina que compilo, la
        # prueba pasaria aqui y fallaria en la de cualquier otro. El TET10
        # pasa por el proceso hijo de `mallar_aislado`, que en un ejecutable
        # congelado depende de `freeze_support`: se prueba lo que se usa.
        from spinpy import febio
        exe = febio.localizar()
        if exe is None:
            if getattr(sys, "frozen", False):
                raise RuntimeError("no se encontro la copia empaquetada de "
                                   "FEBio")
            return "FEBio no instalado (desde el codigo no es obligatorio)"
        org = febio.origen(exe)
        if getattr(sys, "frozen", False) and org != "empaquetado":
            raise RuntimeError(f"FEBio no es el empaquetado: {exe}")
        out = []
        # Cubo con una cavidad esferica: malla de forma fiable a 16^3 (un
        # spinodoide a esa resolucion tiene puntales de 1-2 voxeles, que es
        # justo donde tetgen se cae, y la prueba es del camino, no de eso).
        c = (np.arange(16) - 7.5) / 16
        X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
        cav = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) > 0.2
        with tempfile.TemporaryDirectory() as d:
            for tipo in ("hex8", "tet10"):
                r = febio.analizar(cav, np.full(3, 1.0 / 16),
                                   febio.protocolo("app"), malla=tipo,
                                   carpeta=d, n=16, exe=exe,
                                   comparar_app=False, etiqueta="auto")
                if not r.get("lineal"):
                    raise RuntimeError(f"{tipo}: {r.get('fallos')}")
                out.append(f"{tipo} E_app {r['lineal']['E_app'] / 1e6:.1f} MPa")
        return f"FEBio {febio.version(exe)} ({org}); " + ", ".join(out)
    prueba("FEBio: ensayo lineal hex8 y TET10 16^3", _febio)

    def _lote():
        # Se ejercita con la MISMA politica que usa el lote de verdad
        # (`opciones_paralelo`), no con una llamada de juguete: lo que hay que
        # comprobar es justamente que esa politica funciona aqui. Con
        # `prefer="processes"` dentro del ejecutable esta prueba falla, que es
        # como se descubrio el problema.
        import pandas as pd
        import joblib
        from spinpy.lote import opciones_paralelo
        pd.DataFrame({"a": [1, 2]}).to_csv(os.devnull, index=False)
        op = opciones_paralelo(2)
        r = joblib.Parallel(**op)(joblib.delayed(abs)(x) for x in (-3, -4))
        assert r == [3, 4], r
        _cerrar_obreros()
        return f"pandas y joblib responden (prefer={op['prefer']})"
    prueba("lote (pandas + joblib)", _lote)

    lineas.append("")
    lineas.append("=" * 78)
    lineas.append(f"  FALLOS: {fallos}")
    lineas.append("=" * 78)

    texto = "\n".join(lineas)
    print(texto)
    try:
        destino = DATOS / "autocomprobacion.txt"
        destino.parent.mkdir(parents=True, exist_ok=True)
        destino.write_text(texto, encoding="utf-8")
        print(f"\nInforme guardado en: {destino}")
    except Exception as e:
        print(f"\nNo se pudo guardar el informe: {e}")
    return 1 if fallos else 0


def main():
    # PRIMERO DE TODO, y no es opcional en un ejecutable empaquetado.
    #
    # La seccion de lote usa joblib, que lanza procesos hijos reejecutando el
    # PROPIO ejecutable. Sin `freeze_support()` cada hijo vuelve a entrar por
    # aqui, crea otra ventana y lanza mas hijos: la aplicacion se multiplica
    # sola hasta agotar la maquina. Corriendo desde el codigo fuente no pasa
    # -Python arranca por el interprete, no por el script- asi que es un fallo
    # que solo aparece despues de empaquetar.
    import multiprocessing
    multiprocessing.freeze_support()

    # La autocomprobacion va ANTES de crear la QApplication: tiene que poder
    # correr en una maquina sin sesion grafica -un servidor, una consola
    # remota- que es justo donde uno acaba diagnosticando.
    if "--autocomprobacion" in sys.argv:
        sys.exit(autocomprobacion())

    # Las replicas tambien se pueden pedir por linea de ordenes, y la ventana
    # de Validacion las lanza ASI, como proceso aparte.
    #
    # No es una comodidad: es la unica forma de que no reviente. La replica
    # termina componiendo una figura con `pyvista.Plotter(off_screen=True)`, y
    # VTK no puede crear su contexto OpenGL fuera del hilo principal. Hecho
    # desde un QThread, el proceso MUERE de golpe -0xC0000409, sin excepcion
    # que capturar- justo al llegar a la figura, despues de haber calculado
    # bien todas las comprobaciones. En un proceso nuevo la figura se compone
    # en su propio hilo principal y no hay conflicto; de paso, un fallo en la
    # replica ya no puede llevarse por delante la aplicacion.
    for _clave, _opcion in (("kumar2020", "--replicar-kumar2020"),
                            ("zheng2021", "--replicar-zheng2021"),
                            ("guo2024", "--replicar-guo2024")):
        if _opcion in sys.argv:
            sys.exit(_replicar(_clave, "--rapido" in sys.argv))

    # Identidad propia en la barra de tareas de Windows, ANTES de crear
    # ninguna ventana. Sin esto, al ejecutar desde el codigo Windows agrupa la
    # ventana bajo `python.exe` y le pone el icono generico de Python, por muy
    # bien que se llame a `setWindowIcon`: el icono de la barra lo decide el
    # AppUserModelID del proceso, no el de la ventana.
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                "uv.spinpy.visor")
        except Exception:
            pass

    pv.set_plot_theme("document")
    app = QtWidgets.QApplication(sys.argv)
    _ico = _ruta_icono()
    if _ico is not None:
        app.setWindowIcon(QtGui.QIcon(str(_ico)))

    # La portada va ANTES de instanciar el Visor: construir la ventana entera
    # -con su render de VTK- tarda unos segundos y la pantalla se queda vacia
    # mientras tanto.
    #
    # El baile con el idioma no es un adorno. `DialogoBienvenida` se construye
    # ya traducido, porque se crea una sola vez y no vuelve a pasar por
    # `aplicar`. El `Visor`, en cambio, captura sus textos ORIGINALES en
    # espanol y solo despues se traduce a si mismo, asi que hay que devolver
    # el idioma global a "es" antes de crearlo. Sin esto, arrancar en ingles
    # dejaria grabadas las cadenas inglesas como originales y el diccionario
    # -indexado por el espanol- ya no encontraria ninguna.
    if DialogoBienvenida.mostrar_al_arrancar():
        _cod = QtCore.QSettings("spinpy", "visor").value("idioma", "es")
        if _cod in IDIOMAS:
            fijar_idioma(_cod)
        DialogoBienvenida().exec_()
        fijar_idioma("es")

    v = Visor()
    # Maximizada, no `showFullScreen()`: en pantalla completa de verdad
    # desaparecen la barra de titulo y la de tareas, y no habria forma obvia de
    # cerrar la ventana ni de volver a otra aplicacion.
    v.showMaximized()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
