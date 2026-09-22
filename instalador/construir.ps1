<#
    construir.ps1 - monta el entorno, empaqueta spinpy y produce el instalador.

    Reconstruir esto a mano son media docena de pasos que hay que dar en el
    orden correcto y en el entorno correcto; el error tipico -empaquetar desde
    la instalacion de Anaconda del dia a dia en vez de desde un entorno
    fijado- produce un ejecutable que arranca bien aqui y falla en otra
    maquina. Por eso esta guionizado.

    Uso:
        .\construir.ps1                 # todo
        .\construir.ps1 -SaltarEntorno  # reutiliza el venv ya montado
        .\construir.ps1 -SoloZip        # no intenta el instalador

    Salidas:
        <Trabajo>\dist\spinpy\               carpeta ejecutable
        <Trabajo>\spinpy-<v>-windows.zip     paquete portable
        instalador\salida\spinpy-<v>-instalador.exe
#>
param(
    [string]$Trabajo = "$env:USERPROFILE\spinpy_build",
    [string]$Version = "0.1.0",
    [switch]$SaltarEntorno,
    [switch]$SoloZip
)

$ErrorActionPreference = "Stop"
$AQUI = Split-Path -Parent $MyInvocation.MyCommand.Path
$RAIZ = Split-Path -Parent $AQUI

Write-Host "== spinpy $Version ==" -ForegroundColor Cyan
Write-Host "codigo   : $RAIZ"
Write-Host "trabajo  : $Trabajo"

# --- 1. entorno -------------------------------------------------------------
# Un venv APARTE y con versiones fijadas, no el entorno de trabajo. Dos
# razones: el de trabajo arrastra media distribucion de Anaconda dentro del
# paquete, y sus versiones cambian sin avisar - el ejecutable tiene que dar
# los mismos numeros que se validaron. Ver requisitos_build.txt.
$py = Join-Path $Trabajo "venv\Scripts\python.exe"
if (-not $SaltarEntorno) {
    if (-not (Test-Path $Trabajo)) { New-Item -ItemType Directory -Force $Trabajo | Out-Null }
    if (-not (Test-Path $py)) {
        Write-Host "-> creando el entorno" -ForegroundColor Yellow
        python -m venv (Join-Path $Trabajo "venv")
    }
    Write-Host "-> instalando dependencias fijadas" -ForegroundColor Yellow
    & $py -m pip install --upgrade pip --quiet
    & $py -m pip install -r (Join-Path $AQUI "requisitos_build.txt") --quiet
}
if (-not (Test-Path $py)) { throw "No hay entorno en $py" }

# --- 2. icono ---------------------------------------------------------------
if (-not (Test-Path (Join-Path $AQUI "spinpy.ico"))) {
    Write-Host "-> generando el icono" -ForegroundColor Yellow
    & $py (Join-Path $AQUI "hacer_icono.py")
}

# --- 3. autocomprobacion ANTES de empaquetar --------------------------------
# Si algo esta roto en el codigo, mejor saberlo ahora que despues de diez
# minutos de compresion LZMA.
Write-Host "-> autocomprobacion desde el codigo" -ForegroundColor Yellow
Push-Location $RAIZ
& $py "visor.py" "--autocomprobacion"
$rc = $LASTEXITCODE
Pop-Location
if ($rc -ne 0) { throw "La autocomprobacion fallo desde el codigo. No se empaqueta." }

# --- 4. empaquetar ----------------------------------------------------------
Write-Host "-> PyInstaller" -ForegroundColor Yellow
& $py -m PyInstaller --noconfirm --clean `
    --distpath (Join-Path $Trabajo "dist") `
    --workpath (Join-Path $Trabajo "build") `
    (Join-Path $AQUI "spinpy.spec")
if ($LASTEXITCODE -ne 0) { throw "PyInstaller fallo" }

$dist = Join-Path $Trabajo "dist\spinpy"
$exe = Join-Path $dist "spinpy.exe"
if (-not (Test-Path $exe)) { throw "No se genero $exe" }
$mb = [math]::Round((Get-ChildItem $dist -Recurse | Measure-Object Length -Sum).Sum / 1MB, 0)
Write-Host "   carpeta: $mb MB"

# --- 5. autocomprobacion DEL EJECUTABLE -------------------------------------
# Esta es la que importa. Que compile no significa que funcione: los modulos
# pesados no se importan hasta que se usan, asi que un fallo de empaquetado en
# el mallador o en el solver no se ve al arrancar la ventana.
Write-Host "-> autocomprobacion del ejecutable" -ForegroundColor Yellow
$informe = Join-Path ([Environment]::GetFolderPath('MyDocuments')) "spinpy\autocomprobacion.txt"
if (Test-Path $informe) { Remove-Item $informe -Force }
$p = Start-Process $exe -ArgumentList "--autocomprobacion" -PassThru -Wait
if (Test-Path $informe) { Get-Content $informe | Select-Object -Last 20 }
if ($p.ExitCode -ne 0) { throw "El EJECUTABLE no pasa su autocomprobacion (codigo $($p.ExitCode))" }

# --- 6. zip portable --------------------------------------------------------
Write-Host "-> zip portable" -ForegroundColor Yellow
$zip = Join-Path $Trabajo "spinpy-$Version-windows.zip"
if (Test-Path $zip) { Remove-Item $zip -Force }
Copy-Item (Join-Path $AQUI "LEEME.txt") $dist -Force
Copy-Item (Join-Path $AQUI "LICENCIA_BINARIO.txt") $dist -Force
Compress-Archive -Path $dist -DestinationPath $zip -CompressionLevel Optimal
Write-Host "   $zip  ($([math]::Round((Get-Item $zip).Length/1MB,0)) MB)"

if ($SoloZip) { Write-Host "listo (solo zip)" -ForegroundColor Green; exit 0 }

# --- 7. instalador ----------------------------------------------------------
$iscc = @(
    "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
    "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
) | Where-Object { Test-Path $_ } | Select-Object -First 1

if (-not $iscc) {
    Write-Host "Inno Setup no esta instalado; no se genera el instalador." -ForegroundColor Yellow
    Write-Host "  winget install JRSoftware.InnoSetup"
    Write-Host "El zip de arriba ya es utilizable: descomprimir y ejecutar spinpy.exe."
    exit 0
}

Write-Host "-> Inno Setup (la compresion tarda varios minutos)" -ForegroundColor Yellow
Push-Location $AQUI
& $iscc "/DFUENTE=$dist" "/DVERSION=$Version" "spinpy.iss"
$rc = $LASTEXITCODE
Pop-Location
if ($rc -ne 0) { throw "Inno Setup fallo" }

$inst = Join-Path $AQUI "salida\spinpy-$Version-instalador.exe"
Write-Host "listo: $inst ($([math]::Round((Get-Item $inst).Length/1MB,0)) MB)" -ForegroundColor Green
