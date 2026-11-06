#!/bin/bash
# construir_mac.sh - empaqueta la version de PRUEBA de spinpy para macOS y la
# comprime. Sirve para Apple Silicon (arm64) y para Mac con Intel (x86_64): la
# arquitectura es la del equipo en el que se ejecuta, no se cruza.
#
#   bash instalador/construir_mac.sh [CARPETA_DE_TRABAJO] [VERSION]
#
# Equivale a `construir.ps1 -SoloZip` de Windows: entorno fijado aparte,
# autocomprobacion desde el codigo, PyInstaller, autocomprobacion DEL .app y
# zip. No hay instalador (.dmg) ni firma. Salida:
#   <trabajo>/spinpy-V<version>-macos-<arm64|x86_64>-PRUEBA.zip
set -euo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$AQUI")"
TRABAJO="${1:-$HOME/spinpy_build}"
VERSION="${2:-2.3.0}"
PY="${PYTHON:-python3}"
ARQ="$(uname -m)"
case "$ARQ" in
    arm64|x86_64) ;;
    *) echo "Arquitectura no prevista: $ARQ (se esperaba arm64 o x86_64)"; exit 1 ;;
esac
# Bajo Rosetta, `uname -m` dice x86_64 en un Mac con Apple Silicon: el
# ejecutable resultante seria de Intel. Es lo que se pide, pero se avisa.
if [ "$ARQ" = "x86_64" ] && [ "$(sysctl -n sysctl.proc_translated 2>/dev/null || echo 0)" = "1" ]; then
    echo "AVISO: Python corre bajo Rosetta; se construira la version para Intel."
fi
ZIP="$TRABAJO/spinpy-V$VERSION-macos-$ARQ-PRUEBA.zip"

echo "== spinpy V$VERSION (macOS, prueba) =="
echo "codigo  : $RAIZ"
echo "trabajo : $TRABAJO"
echo "arq.    : $ARQ"

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
STAGE="$TRABAJO/stage/spinpy-V$VERSION-macos-$ARQ-PRUEBA"
rm -rf "$TRABAJO/stage" "$ZIP"
mkdir -p "$STAGE"
ditto "$APP" "$STAGE/spinpy.app"
cp "$AQUI/LEEME_MAC.txt" "$AQUI/LICENCIA_BINARIO.txt" "$STAGE/"
ditto -c -k --keepParent "$STAGE" "$ZIP"
echo "   $ZIP ($(du -m "$ZIP" | cut -f1) MB)"
echo "listo"
