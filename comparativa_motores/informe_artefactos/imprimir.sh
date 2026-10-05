#!/bin/sh
# Imprime informe_artefactos.html a PDF con Chromium sin interfaz.
# Uso: sh imprimir.sh [ruta/a/chrome]
set -e
cd "$(dirname "$0")"
CHROME="${1:-${CHROME:-chromium}}"
"$CHROME" --headless=new --no-sandbox --disable-gpu --allow-file-access-from-files \
  --no-pdf-header-footer --run-all-compositor-stages-before-draw \
  --virtual-time-budget=20000 \
  --print-to-pdf="$(pwd)/../INFORME_ARTEFACTOS.pdf" "file://$(pwd)/informe_artefactos.html"
