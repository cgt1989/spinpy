"""
familias.py — Fuente unica de las familias de microestructura.

Hasta la tercera familia (level-set-foam) las dos que habia estaban repetidas
como tuplas literales `(("spinodoide", ""), ("dual-lattice", "_dual"))` en el
informe, las figuras, los estudios y la interfaz. Aqui vive lo que NO depende
de la familia concreta: el codigo, el sufijo de sus claves de resultados, su
etiqueta y de donde viene.

Anadir una cuarta familia es anadir una entrada en `INFO` (y su generador,
ajuste y registro de metodos); lo que recorre `FAMILIAS` la recoge sola.
"""

from __future__ import annotations

# El orden es el de aparicion en la interfaz, la tabla y el informe.
FAMILIAS = ("spinodoide", "dual-lattice", "level-set-foam")

# Modos del selector: cada familia sola, o todas a la vez.
MODOS_FAMILIA = FAMILIAS + ("todas",)

INFO = {
    "spinodoide": {
        "sufijo": "",                 # el spinodoide conserva las claves de siempre
        "morfometria": "morfometria_spin",
        "etiqueta": {"es": "Spinodoide", "en": "Spinodoid"},
        "ref": "Kumar et al. 2020 (npj Comput. Mater.)",
        "libreria": "spinpy (propia, sobre GIBBON)",
    },
    "dual-lattice": {
        "sufijo": "_dual",
        "morfometria": "morfometria_dual",
        "etiqueta": {"es": "Dual-lattice", "en": "Dual-lattice"},
        "ref": "Vafaeefar et al. 2022",
        "libreria": "spinpy (propia)",
    },
    "level-set-foam": {
        "sufijo": "_ls",
        "morfometria": "morfometria_ls",
        "etiqueta": {"es": "Level-set foam", "en": "Level-set foam"},
        "ref": "Sonon, Francois y Massart 2015; Kilingar et al. 2019",
        "libreria": "morphoGen (GPL-3.0-or-later), opcional",
    },
}


def sufijo(familia):
    """Sufijo de las claves de resultados de la familia ('' para el spinodoide)."""
    return INFO[familia]["sufijo"]


def clave(base, familia):
    """`clave('ajuste', 'dual-lattice')` -> 'ajuste_dual'."""
    return base + INFO[familia]["sufijo"]


def clave_morfometria(familia):
    return INFO[familia]["morfometria"]


def etiqueta(familia, idioma="es"):
    return INFO[familia]["etiqueta"].get(idioma, INFO[familia]["etiqueta"]["es"])


def con_sufijos(base):
    """[(familia, clave)] de una clave de resultados, en orden de FAMILIAS."""
    return [(f, base + INFO[f]["sufijo"]) for f in FAMILIAS]
