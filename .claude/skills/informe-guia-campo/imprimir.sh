#!/bin/sh
# Imprime el informe HTML a PDF con Chromium sin interfaz.
# Variables: ENTRADA (por omision informe.html) y SALIDA (informe.pdf); las
# rutas relativas lo son a la carpeta de este guion.
# Uso: ENTRADA=informe_x.html SALIDA=../INFORME_X.pdf sh imprimir.sh [ruta/a/chrome]
set -e
cd "$(dirname "$0")"
CHROME="${1:-${CHROME:-chromium}}"
abs() { case "$1" in /*) echo "$1" ;; *) echo "$(pwd)/$1" ;; esac; }
"$CHROME" --headless=new --no-sandbox --disable-gpu --allow-file-access-from-files \
  --no-pdf-header-footer --run-all-compositor-stages-before-draw \
  --virtual-time-budget=20000 \
  --print-to-pdf="$(abs "${SALIDA:-informe.pdf}")" "file://$(abs "${ENTRADA:-informe.html}")"
