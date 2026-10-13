#!/bin/bash
# construir_mac.sh - empaqueta la version de PRUEBA de spinpy para macOS
# (Apple Silicon) y la comprime.
#
#   bash instalador/construir_mac.sh [CARPETA_DE_TRABAJO] [VERSION]
#
# Equivale a `construir.ps1 -SoloZip` de Windows: entorno fijado aparte,
# autocomprobacion desde el codigo, PyInstaller, autocomprobacion DEL .app y
# zip. No hay instalador (.dmg) ni firma. Salida:
#   <trabajo>/spinpy-V<version>-macos-arm64-PRUEBA.zip
set -euo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$AQUI")"
TRABAJO="${1:-$HOME/spinpy_build}"
VERSION="${2:-2.2.1}"
PY="${PYTHON:-python3}"
ZIP="$TRABAJO/spinpy-V$VERSION-macos-arm64-PRUEBA.zip"

echo "== spinpy V$VERSION (macOS, prueba) =="
echo "codigo  : $RAIZ"
echo "trabajo : $TRABAJO"
echo "arq.    : $(uname -m)"
[ "$(uname -m)" = "arm64" ] || echo "AVISO: esta receta es para Apple Silicon (arm64)."

mkdir -p "$TRABAJO"

echo "-> entorno con versiones fijadas"
[ -x "$TRABAJO/venv/bin/python" ] || "$PY" -m venv "$TRABAJO/venv"
VPY="$TRABAJO/venv/bin/python"
"$VPY" -m pip install --upgrade pip --quiet
"$VPY" -m pip install -r "$AQUI/requisitos_build_mac.txt" --quiet

echo "-> icono"
"$VPY" - <<PYEOF
from PIL import Image
im = Image.open("$AQUI/spinpy_icono.png").convert("RGBA")
im = im.resize((512, 512), Image.LANCZOS)
im.save("$AQUI/spinpy.icns")
PYEOF

# Sin pantalla: el runner de CI no tiene sesion grafica de usuario.
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

echo "-> autocomprobacion desde el codigo"
cd "$RAIZ"
"$VPY" visor.py --autocomprobacion

echo "-> empaquetando"
SPINPY_VERSION="$VERSION" "$VPY" -m PyInstaller --noconfirm --clean \
    --distpath "$TRABAJO/dist" --workpath "$TRABAJO/work" \
    "$AQUI/spinpy_mac.spec"
APP="$TRABAJO/dist/spinpy.app"
[ -d "$APP" ] || { echo "No se genero $APP"; exit 1; }

echo "-> autocomprobacion del .app"
"$APP/Contents/MacOS/spinpy" --autocomprobacion
codesign --verify --deep --verbose=1 "$APP" || echo "AVISO: firma ad hoc no verificable"

echo "-> zip"
STAGE="$TRABAJO/stage/spinpy-V$VERSION-macos-arm64-PRUEBA"
rm -rf "$TRABAJO/stage" "$ZIP"
mkdir -p "$STAGE"
ditto "$APP" "$STAGE/spinpy.app"
cp "$AQUI/LEEME_MAC.txt" "$AQUI/LICENCIA_BINARIO.txt" "$STAGE/"
ditto -c -k --keepParent "$STAGE" "$ZIP"
echo "   $ZIP ($(du -m "$ZIP" | cut -f1) MB)"
echo "listo"
