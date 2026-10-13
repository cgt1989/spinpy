#!/bin/bash
# construir_linux.sh - empaqueta spinpy para Linux (x86_64) y lo comprime.
#
#   bash instalador/construir_linux.sh [CARPETA_DE_TRABAJO] [VERSION]
#
# Equivale a `construir_mac.sh`: entorno fijado aparte, autocomprobacion desde
# el codigo, PyInstaller, autocomprobacion DEL ejecutable empaquetado y
# .tar.gz. No hay instalador (.deb/.rpm/AppImage). Salida:
#   <trabajo>/spinpy-V<version>-linux-x86_64.tar.gz
set -euo pipefail

AQUI="$(cd "$(dirname "$0")" && pwd)"
RAIZ="$(dirname "$AQUI")"
TRABAJO="${1:-$HOME/spinpy_build}"
VERSION="${2:-2.2.1}"
PY="${PYTHON:-python3}"
PAQUETE="spinpy-V$VERSION-linux-x86_64"
TGZ="$TRABAJO/$PAQUETE.tar.gz"

echo "== spinpy V$VERSION (Linux) =="
echo "codigo  : $RAIZ"
echo "trabajo : $TRABAJO"
echo "arq.    : $(uname -m)"
[ "$(uname -m)" = "x86_64" ] || echo "AVISO: esta receta es para x86_64."

mkdir -p "$TRABAJO"

echo "-> entorno con versiones fijadas"
[ -x "$TRABAJO/venv/bin/python" ] || "$PY" -m venv "$TRABAJO/venv"
VPY="$TRABAJO/venv/bin/python"
"$VPY" -m pip install --upgrade pip --quiet
"$VPY" -m pip install -r "$AQUI/requisitos_build_linux.txt" --quiet

# Sin pantalla: el runner de CI no tiene sesion grafica. xvfb-run (si esta)
# aporta un servidor X virtual para VTK/OpenGL.
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"
EJEC=()
command -v xvfb-run >/dev/null 2>&1 && EJEC=(xvfb-run -a)

echo "-> autocomprobacion desde el codigo"
cd "$RAIZ"
"${EJEC[@]}" "$VPY" visor.py --autocomprobacion

echo "-> empaquetando"
SPINPY_VERSION="$VERSION" "$VPY" -m PyInstaller --noconfirm --clean \
    --distpath "$TRABAJO/dist" --workpath "$TRABAJO/work" \
    "$AQUI/spinpy_linux.spec"
DIST="$TRABAJO/dist/spinpy"
[ -x "$DIST/spinpy" ] || { echo "No se genero $DIST/spinpy"; exit 1; }

echo "-> autocomprobacion del ejecutable empaquetado"
"${EJEC[@]}" "$DIST/spinpy" --autocomprobacion

echo "-> tar.gz"
STAGE="$TRABAJO/stage/$PAQUETE"
rm -rf "$TRABAJO/stage" "$TGZ"
mkdir -p "$STAGE"
cp -a "$DIST/." "$STAGE/"
cp "$AQUI/LEEME_LINUX.txt" "$AQUI/LICENCIA_BINARIO.txt" "$STAGE/"
tar -C "$TRABAJO/stage" -czf "$TGZ" "$PAQUETE"
echo "   $TGZ ($(du -m "$TGZ" | cut -f1) MB)"
echo "listo"
