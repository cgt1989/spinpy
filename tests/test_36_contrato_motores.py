"""
Bloque 36 - El contrato de un motor FEM: lo que debe cumplir cualquier motor
registrado en `spinpy.motores.REGISTRO`, este instalado o no.

POR QUE EXISTE
--------------
Anadir un motor es anadir una ficha al registro y un archivo `m_<clave>.py`
(`docs/AGREGAR_UN_MOTOR_FEM.md`). Esta prueba comprueba la parte ESTATICA del
contrato, la que no necesita la biblioteca instalada: que la ficha y el
archivo existen, que el modulo se importa sin la biblioteca (para que
`disponibles()` pueda sondear sin cargar nada), que ofrece `resolver` y
`version`, y que lo que la ficha declara es coherente. La parte numerica
(bloque macizo, plato, no lineal, convergencia) la hace el bloque 29, que
recorre `motores.MOTORES` y por tanto prueba automaticamente un motor nuevo.

TOLERANCIAS: ninguna; son comprobaciones exactas de estructura.
"""
import importlib
import inspect
from pathlib import Path

import pytest

from spinpy import idioma_textos, motores

RAIZ_MOTORES = Path(motores.__file__).parent
CLAVES = [m.clave for m in motores.REGISTRO]


def test_las_tablas_derivadas_salen_del_registro():
    assert motores.MOTORES == tuple(CLAVES)
    assert len(set(CLAVES)) == len(CLAVES)
    assert motores.CAPACIDADES.keys() == motores.ETIQUETAS.keys() \
        == set(CLAVES)
    assert motores.REGISTRO[0].clave == "app"      # la linea base va primero
    assert motores.RECOMENDADO in CLAVES


@pytest.mark.parametrize("clave", CLAVES)
def test_ficha_y_archivo_coherentes(clave):
    ficha = next(m for m in motores.REGISTRO if m.clave == clave)
    assert (RAIZ_MOTORES / f"m_{clave}.py").is_file()
    assert set(ficha.mallas) <= {"hex8", "tet10"} and ficha.mallas
    assert set(ficha.nl) <= {"svk", "neohookeano"}
    # un motor sin sonda esta siempre; con sonda, se detecta sin importarlo
    assert (ficha.sonda is None) == (clave == "app")
    assert isinstance(motores.instalado(clave), bool)


@pytest.mark.parametrize("clave", CLAVES)
def test_modulo_importable_sin_la_biblioteca(clave):
    """El modulo NO puede importar la biblioteca del motor al nivel del
    modulo: lo hace dentro de `resolver`, en el proceso hijo."""
    mod = importlib.import_module(f"spinpy.motores.m_{clave}")
    assert callable(mod.resolver)
    assert callable(mod.version)
    assert isinstance(mod.NOMBRE, str) and mod.NOMBRE
    assert mod.NOMBRE == motores.ETIQUETAS[clave]
    assert isinstance(mod.LICENCIA, str) and mod.LICENCIA
    # `resolver(problema)` recibe un unico argumento
    assert len(inspect.signature(mod.resolver).parameters) == 1


@pytest.mark.parametrize("clave", CLAVES)
def test_la_etiqueta_esta_en_el_diccionario_de_idioma(clave):
    """`idioma_revisar.py` no ve las etiquetas de motor (llegan a `_()` como
    variable): hay que declararlas a mano en el diccionario."""
    assert motores.ETIQUETAS[clave] in idioma_textos.EN


def test_puede_es_coherente_con_las_capacidades():
    assert motores.puede("app", "hex8", "lineal", control="fuerza")
    assert not motores.puede("app", "tet10")
    assert not motores.puede("app", "hex8", "nl")
    assert not motores.puede("app", "hex8", control="plato")
    assert motores.puede("sfepy", "tet10", "nl", "svk")
    assert not motores.puede("sfepy", "tet10", "nl", "neohookeano")


def test_la_plantilla_sigue_la_forma_del_contrato():
    """`_plantilla.py` es el punto de partida que se copia: si el contrato
    cambia y la plantilla no, la guia miente. No esta en el registro."""
    from spinpy.motores import _plantilla
    assert "plantilla" not in CLAVES
    assert callable(_plantilla.resolver) and callable(_plantilla.version)
    assert len(inspect.signature(_plantilla.resolver).parameters) == 1
    with pytest.raises(motores.NoDisponible):
        _plantilla.resolver({})
    # y cada paso de la guia nombra algo que existe
    guia = (RAIZ_MOTORES.parent.parent / "docs"
            / "AGREGAR_UN_MOTOR_FEM.md").read_text(encoding="utf-8")
    for nombre in ("REGISTRO", "_plantilla.py", "test_29_motores",
                   "idioma_revisar.py", "spinpy.spec"):
        assert nombre in guia, nombre
