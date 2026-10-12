#!/bin/bash
# Espera a que acabe la campana de tiempos y luego ajusta y ensaya.
while pgrep -f "correr.py exactos" >/dev/null; do sleep 20; done
date -Is > ../tiempos_sesion/FIN_CAMPANA.txt
cd "$(dirname "$0")"
export PYTHONWARNINGS=ignore
python ajustar.py C2 C4 > ajustar.log 2>&1
python replicar_tapia.py C2 C4 --motores app ngsolve skfem sfepy referencia > replicar.log 2>&1
/tmp/claude-0/mm/fx/bin/python replicar_tapia.py C2 C4 --motores fenicsx > replicar_fx.log 2>&1
date -Is > FIN_TAPIA.txt
