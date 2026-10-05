"""
validar_app.py: Validacion de extremo a extremo del codigo INTEGRADO.

Corre `fem.ensayo` (el camino de la GUI y de la CLI) con los analisis
'lineal' y 'lineal_plato' sobre el espinodoide de referencia y compara lo que
la app publica (`lineal.E_app`, linea base; `corregido`, con las
correcciones) con la referencia embebida de `resultados/referencia.json`.
Escribe resultados/validacion_app.json.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))
sys.path.insert(0, str(AQUI.parent.parent))
import casos                                                    # noqa: E402
from spinpy import fem                                          # noqa: E402

RES = AQUI / "resultados"
ref = json.loads((RES / "referencia.json").read_text())
prot = fem.protocolo("app", E_s=20e9)
out = {}
for tipo, ns in (("hex8", (32, 48, 64)), ("tet10", (32, 48))):
    for n in ns:
        BW, sp = casos.espinodoide(n)
        malla = fem.mallar(BW, sp, tipo)
        t0 = time.perf_counter()
        reg = fem.ensayo(malla, prot, analisis=["lineal", "lineal_plato"],
                         motor="ngsolve", aislado=False)
        t = time.perf_counter() - t0
        r = ref[f"n{n}"]
        lin, cor, lp = reg["lineal"], reg["corregido"], reg["lineal_plato"]
        sref = prot["sigma_app"]
        d = {"t_s": t, "ok": reg["ok"],
             "E_base": lin["E_app"] / 1e6,
             "E_corregido": cor["E_app"] / 1e6, "metodo_E": cor["metodo_E"],
             "E_ref_nucleo": r["nucleo"]["E"], "E_ref_voi": r["voi"]["E"],
             "p99_base": lin["pistoia"]["vm_p99_superficie"] / sref,
             "p99_corregido": cor["vm_p99_superficie"] / sref,
             "metodo_p99": cor["metodo_p99"],
             "p99_ref_voi": r["voi"]["p99"],
             "E_nucleo_traccion": lin["nucleo"]["E_app"] / 1e6,
             "p99_nucleo_traccion": lin["nucleo"]["vm_p99_superficie"] / sref,
             "p99_ref_nucleo": r["nucleo"]["p99"],
             "cociente_plato": lp["cociente_plato_fuerza"],
             "sigma_fallo_base_MPa": lin["pistoia"]["sigma_fallo"] / 1e6,
             "sigma_fallo_plato_MPa": lp["pistoia"]["sigma_fallo"] / 1e6}
        d["err_E_base"] = d["E_base"] / d["E_ref_voi"] - 1
        d["err_E_corregido"] = d["E_corregido"] / d["E_ref_nucleo"] - 1
        d["err_p99_base"] = d["p99_base"] / d["p99_ref_voi"] - 1
        d["err_p99_corregido"] = d["p99_corregido"] / d["p99_ref_voi"] - 1
        out[f"{tipo}_n{n}"] = d
        print(tipo, n, f"{t:.0f}s", f"E base {d['E_base']:.1f} "
              f"({d['err_E_base']:+.3f}) -> {d['E_corregido']:.1f} "
              f"({d['err_E_corregido']:+.3f}) | p99 {d['p99_base']:.1f} "
              f"({d['err_p99_base']:+.3f}) -> {d['p99_corregido']:.1f} "
              f"({d['err_p99_corregido']:+.3f})", flush=True)
(RES / "validacion_app.json").write_text(json.dumps(out, indent=1))
