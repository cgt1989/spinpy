"""
Bloque 15 — Simulaciones in silico (`spinpy.simulacion`).

QUE SE VERIFICA
  Que cada protocolo retira LO QUE DICE retirar, sobre estructuras construidas
  a mano donde la respuesta se conoce antes de ejecutar. Los numeros de
  rigidez no se contrastan contra nada externo —no hay solucion cerrada para
  una estructura trabecular que pierde hueso—; lo que se comprueba es el
  ORDEN y la CONTABILIDAD, que es donde un error no daria sintomas.

  contabilidad    cada paso de perdida retira exactamente round(f * N0)
                  voxeles, en los cuatro protocolos. Es la regla que hace
                  comparables las curvas.
  adelgazamiento  sobre una losa gruesa, la primera capa que sale es la de
                  superficie: el nucleo sigue intacto.
  trab. finas     con una columna gruesa y otra fina, sale entera la fina
                  antes de tocar la gruesa.
  desuso          con una columna portante en z y una placa colgada que no
                  une los platos, sale antes la placa (no trabaja).
  recuperacion    termina con el mismo numero de voxeles que el inicio.
  fallo           el dano acumulado crece y la rigidez no sube; los rotos del
                  paso 0 son ~frac del tejido (el percentil de Pistoia).
  determinismo    misma entrada y semilla, mismas mascaras.
"""
import numpy as np
import pytest

from spinpy.simulacion import (PROTOCOLOS, fallo_progresivo, simular_perdida)

BLOQUE = "15 Simulaciones in silico"
REF = "construccion a mano"
SP = np.full(3, 0.05)


def _columnas(n=24, gruesa=6, fina=2):
    """Dos columnas en z que unen base y techo: una gruesa y una fina."""
    BW = np.zeros((n, n, n), bool)
    BW[3:3 + gruesa, 3:3 + gruesa, :] = True
    BW[16:16 + fina, 16:16 + fina, :] = True
    return BW


def test_contabilidad(registro):
    BW = _columnas()
    N0 = int(BW.sum())
    ok = True
    notas = []
    for prot in PROTOCOLOS:
        r = simular_perdida(BW, SP, protocolo=prot, pasos=2, perdida_paso=0.05,
                            n_trabajo=24, n_mec=12, guardar_mascaras=False)
        k = int(round(0.05 * N0))
        perdida = [f for f in r["pasos"] if f["fase"] == "perdida"]
        ok &= all(f["cambio_voxeles"] == -k for f in perdida)
        if prot == "recuperacion":
            ok &= r["pasos"][-1]["voxeles"] == N0
        notas.append(f"{prot}:{[f['cambio_voxeles'] for f in r['pasos']]}")
    registro.anotar(BLOQUE, "voxeles retirados por paso", REF, None, None,
                    "exacto", "regla de comparabilidad", ok, nota="; ".join(notas))
    assert ok, notas


def test_adelgazamiento_quita_superficie(registro):
    BW = np.zeros((20, 20, 20), bool)
    BW[4:16, 4:16, :] = True                      # losa gruesa que une platos
    nucleo = np.zeros_like(BW); nucleo[6:14, 6:14, :] = True
    r = simular_perdida(BW, SP, "adelgazamiento", pasos=1, perdida_paso=0.10,
                        n_trabajo=20, n_mec=10)
    m1 = r["mascaras"][1]
    ok = bool(np.all(m1[nucleo]))
    registro.anotar(BLOQUE, "adelgazamiento respeta el nucleo", REF, True, ok,
                    "nucleo intacto", "10 % retirado", ok)
    assert ok


def test_trabeculas_finas_salen_enteras(registro):
    BW = _columnas()
    fina = np.zeros_like(BW); fina[16:18, 16:18, :] = True
    gruesa = np.zeros_like(BW); gruesa[3:9, 3:9, :] = True
    frac_fina = fina.sum() / BW.sum()               # ~0.1
    r = simular_perdida(BW, SP, "trabeculas_finas", pasos=1,
                        perdida_paso=float(frac_fina), n_trabajo=24, n_mec=12)
    m1 = r["mascaras"][1]
    ok = (not m1[fina].any()) and bool(m1[gruesa].all())
    registro.anotar(BLOQUE, "sale la columna fina entera", REF, True, ok,
                    "fina vacia, gruesa intacta", f"f = {frac_fina:.3f}", ok)
    assert ok


def test_desuso_retira_lo_que_no_trabaja(registro):
    n = 24
    BW = np.zeros((n, n, n), bool)
    BW[8:16, 8:16, :] = True                        # columna portante
    placa = np.zeros_like(BW)
    placa[2:6, 2:22, 4:14] = True                   # placa colgada: no toca platos
    BW |= placa
    frac_placa = placa.sum() / BW.sum()
    r = simular_perdida(BW, SP, "desuso", pasos=1,
                        perdida_paso=float(frac_placa) * 0.5, n_trabajo=n,
                        n_mec=12)
    m1 = r["mascaras"][1]
    col = np.zeros_like(BW); col[8:16, 8:16, :] = True
    perdido_placa = placa.sum() - m1[placa].sum()
    perdido_col = col.sum() - m1[col].sum()
    ok = perdido_col == 0 and perdido_placa > 0
    registro.anotar(BLOQUE, "desuso: sale la placa sin carga", REF, 0,
                    int(perdido_col), "columna intacta", "mecanostato", ok,
                    nota=f"placa -{int(perdido_placa)} vox")
    assert ok


def test_fallo_progresivo(registro):
    rng = np.random.default_rng(3)
    BW = _columnas(n=16, gruesa=5, fina=3)
    BW[:, :, 6:9] |= rng.random((16, 16, 3)) < 0.3  # algo de ruido
    r = fallo_progresivo(BW, SP, pasos=3, n_mec=16)
    f = [x for x in r["pasos"] if x["ok"]]
    dano = [x["dano_acumulado"] for x in r["pasos"]]
    E = [x["E_app"] for x in f]
    frac0 = f[0]["rotos"] / f[0]["voxeles"]
    s = r["resumen"]
    # La carga de fallo nunca se lee despues del colapso: ahi la sostiene el
    # tejido ablandado (medido en H4: sube de 165 a 215 N con E en 0.07).
    antes_del_colapso = (s["colapso_en_paso"] is None
                         or s["paso_F_max"] is None
                         or s["paso_F_max"] < s["colapso_en_paso"])
    ok = (all(b >= a for a, b in zip(dano, dano[1:]))
          and all(b <= a * (1 + 1e-9) for a, b in zip(E, E[1:]))
          and 0.0 < frac0 <= 0.05 and antes_del_colapso)
    registro.anotar(BLOQUE, "fallo progresivo: dano crece, E no sube", REF,
                    0.02, float(frac0), "rotos paso 0 en (0, 0.05]",
                    "percentil de Pistoia", ok,
                    nota=f"dano {['%.3f' % d for d in dano]}")
    assert ok


def test_determinismo(registro):
    BW = _columnas()
    a = simular_perdida(BW, SP, "adelgazamiento", pasos=2, n_trabajo=24, n_mec=10)
    b = simular_perdida(BW, SP, "adelgazamiento", pasos=2, n_trabajo=24, n_mec=10)
    ok = all(np.array_equal(x, y) for x, y in zip(a["mascaras"], b["mascaras"]))
    registro.anotar(BLOQUE, "misma semilla, mismas mascaras", REF, True, ok,
                    "identicas", "desempate con semilla", ok)
    assert ok


def test_protocolo_desconocido():
    with pytest.raises(ValueError):
        simular_perdida(_columnas(), SP, protocolo="magia")
