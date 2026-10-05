#!/bin/sh
# Imprime el informe HTML a PDF con Chrome o Chromium sin interfaz.
# Variables: ENTRADA (por omision informe.html) y SALIDA (informe.pdf); las
# rutas relativas lo son a la carpeta de este guion.
# Uso: ENTRADA=informe_x.html SALIDA=INFORME_X.pdf sh imprimir.sh [ruta/a/chrome]
set -e
cd "$(dirname "$0")"

buscar_chrome() {
  for c in "$CHROME" chromium chromium-browser google-chrome google-chrome-stable \
           "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
           "/Applications/Chromium.app/Contents/MacOS/Chromium" \
           "/c/Program Files/Google/Chrome/Application/chrome.exe" \
           "/c/Program Files (x86)/Google/Chrome/Application/chrome.exe" \
           /opt/pw-browsers/chromium-*/chrome-linux/chrome; do
    [ -n "$c" ] || continue
    if command -v "$c" >/dev/null 2>&1 || [ -x "$c" ]; then echo "$c"; return; fi
  done
  echo "No se encontro Chrome ni Chromium; pasa la ruta como argumento." >&2
  exit 1
}

CHROME="${1:-$(buscar_chrome)}"
abs() { case "$1" in /*) echo "$1" ;; *) echo "$(pwd)/$1" ;; esac; }
"$CHROME" --headless=new --no-sandbox --disable-gpu --allow-file-access-from-files \
  --no-pdf-header-footer --run-all-compositor-stages-before-draw \
  --virtual-time-budget=20000 \
  --print-to-pdf="$(abs "${SALIDA:-informe.pdf}")" "file://$(abs "${ENTRADA:-informe.html}")"
