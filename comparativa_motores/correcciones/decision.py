"""
decision.py: Aplica la regla de PRERREGISTRO.md a los resultados y escribe
resultados/decision.json. No mide nada: solo lee referencia, hex8 y tet10.
"""

import json
from pathlib import Path

RES = Path(__file__).resolve().parent / "resultados"
ref = json.loads((RES / "referencia.json").read_text())
hx = json.loads((RES / "hex8.json").read_text())
tt = json.loads((RES / "tet10.json").read_text())

#: Metodo -> (control del ensayo del VOI, region donde se mide, E usado)
METODOS = {"B0": ("fuerza", "voi", "E_app"), "F1": ("plato", "voi", "E_app"),
           "F2": ("fuerza", "nucleo", "E"), "F12": ("plato", "nucleo", "E")}


def errores(d, r, region, clave_E):
    E = d["E_app"] if clave_E == "E_app" else d[region]["E"]
    return {"E": E, "p99": d[region]["p99"],
            "e_E": abs(E / r[region]["E"] - 1),
            "s_E": E / r[region]["E"] - 1,
            "e_p99": abs(d[region]["p99"] / r[region]["p99"] - 1),
            "s_p99": d[region]["p99"] / r[region]["p99"] - 1}


out = {"hex8": {}, "veredicto": {}}
for n in (32, 48, 64):
    r = ref[f"n{n}"]
    out["hex8"][n] = {m: errores(hx[f"n{n}_{c}"], r, reg, kE)
                      for m, (c, reg, kE) in METODOS.items()}
for m in ("F1", "F2", "F12"):
    filas = []
    for n in (32, 48, 64):
        b, x = out["hex8"][n]["B0"], out["hex8"][n][m]
        c1 = x["e_E"] <= 0.5 * b["e_E"] and b["e_E"] - x["e_E"] >= 0.05
        c2 = x["e_p99"] - b["e_p99"] <= 0.02
        filas.append({"n": n, "criterio_E": c1, "criterio_p99": c2})
    ok = all(f["criterio_E"] and f["criterio_p99"] for f in filas)
    out["veredicto"][m] = {"acepta": ok, "detalle": filas}

# A3: brecha TET10 / hex8 bajo cada protocolo (misma region y definicion)
a3 = {}
for n in (32, 48):
    r = ref[f"n{n}"]
    fila = {}
    for m, (c, reg, kE) in METODOS.items():
        k = f"n{n}_omision_{c}"
        if k not in tt:
            continue
        t = errores(tt[k], r, reg, kE)
        h = out["hex8"][n][m]
        fila[m] = {"E_tet": t["E"], "E_hex": h["E"],
                   "brecha_E": t["E"] / h["E"] - 1, "s_E_ref": t["s_E"],
                   "p99_tet": t["p99"], "p99_hex": h["p99"],
                   "brecha_p99": t["p99"] / h["p99"] - 1,
                   "s_p99_ref": t["s_p99"]}
    a3[n] = fila
out["a3"] = a3
cand = {k: v for k, v in tt.items() if k.startswith("n32_")}
out["a3_candidatos_32"] = {
    k: ({"error": v["error"]} if "error" in v else
        {"E_nucleo_plato": v["nucleo"]["E"], "E_app": v["E_app"]})
    for k, v in cand.items()}
F12_48 = a3[48].get("F12", {})
out["veredicto"]["A3_protocolo"] = {
    "brecha_48": F12_48.get("brecha_E"),
    "resuelto": F12_48.get("brecha_E") is not None
    and abs(F12_48["brecha_E"]) <= 0.10}
out["veredicto"]["A3_malla"] = {
    "acepta": False,
    "motivo": "ningun candidato mejora a la malla por omision en 32^3 "
              "(taubin5 no genera malla valida; sin_suavizado es mas blanda)"}
(RES / "decision.json").write_text(json.dumps(out, indent=1))

for n, d in out["hex8"].items():
    print(n, {m: f"E {v['s_E']:+.3f} p99 {v['s_p99']:+.3f}" for m, v in d.items()})
for n, d in a3.items():
    print("A3", n, {m: f"brecha {v['brecha_E']:+.3f} ref {v['s_E_ref']:+.3f} p99 {v['brecha_p99']:+.3f}" for m, v in d.items()})
print(json.dumps(out["veredicto"], indent=1))
