; ---------------------------------------------------------------------------
; spinpy - guion de Inno Setup 6
;
; Produce  spinpy-<version>-instalador.exe : asistente, acceso directo en el
; menu de inicio, icono y desinstalador en "Aplicaciones instaladas".
;
; SE INSTALA PARA EL USUARIO, NO PARA LA MAQUINA
;   PrivilegesRequired=lowest instala en la carpeta del usuario y NO pide
;   contrasena de administrador. Es deliberado: en un ordenador de
;   universidad, quien va a probar esto normalmente no es administrador de su
;   propio equipo, y un instalador que exige elevacion se queda sin instalar.
;   El precio es que la instalacion vale solo para quien la ejecuta.
;
; LICENCIA
;   El ejecutable incluye PyQt5, que es GPL v3. El codigo fuente de spinpy es
;   MIT, pero este binario se distribuye bajo GPL-3.0. La pagina de licencia
;   del asistente lo dice antes de instalar nada.
;
; Compilar:
;   ISCC.exe /DFUENTE="C:\...\dist\spinpy" /DVERSION=0.1.0 spinpy.iss
; ---------------------------------------------------------------------------

#ifndef FUENTE
  #define FUENTE "..\..\..\..\spinpy_build\dist\spinpy"
#endif
#ifndef VERSION
  #define VERSION "0.1.0"
#endif

#define NOMBRE "spinpy"
#define DESCRIPCION "Visor y ajuste de spinodoides a VOIs de hueso trabecular"
#define AUTOR "Carlos Gonzalez-Torres - Universidad de Valparaiso"
#define URL "https://github.com/carlosgonzalezt/spinpy"

[Setup]
AppId={{7B3F2A64-9C1D-4E58-A0F2-5C9D8E1B4A37}
AppName={#NOMBRE}
AppVersion={#VERSION}
AppVerName={#NOMBRE} {#VERSION}
AppPublisher={#AUTOR}
AppPublisherURL={#URL}
AppSupportURL={#URL}
VersionInfoDescription={#DESCRIPCION}
VersionInfoVersion={#VERSION}

; Instalacion por usuario: sin UAC, sin administrador. Ver cabecera.
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
DefaultDirName={autopf}\{#NOMBRE}
DefaultGroupName={#NOMBRE}
DisableProgramGroupPage=yes
AllowNoIcons=yes

LicenseFile=LICENCIA_BINARIO.txt
InfoBeforeFile=ANTES_DE_INSTALAR.txt
SetupIconFile=spinpy.ico
UninstallDisplayIcon={app}\spinpy.exe
UninstallDisplayName={#NOMBRE} {#VERSION}

OutputDir=.\salida
OutputBaseFilename={#NOMBRE}-{#VERSION}-instalador
; LZMA2 al maximo: son ~650 MB de bibliotecas cientificas, muy comprimibles.
; Comprimir tarda varios minutos; descargar 300 MB en vez de 650 los ahorra
; en cada persona que lo instale.
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Languages]
Name: "es"; MessagesFile: "compiler:Languages\Spanish.isl"

[Tasks]
Name: "escritorio"; Description: "Crear un acceso directo en el escritorio"; \
  GroupDescription: "Accesos directos:"

[Files]
Source: "{#FUENTE}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "LICENCIA_BINARIO.txt"; DestDir: "{app}"; Flags: ignoreversion
Source: "LEEME.txt"; DestDir: "{app}"; Flags: ignoreversion isreadme

[Icons]
Name: "{group}\{#NOMBRE}"; Filename: "{app}\spinpy.exe"
; Acceso directo a la autocomprobacion: si algo va mal en otra maquina, esto
; escribe un informe en Documentos\spinpy que se puede enviar por correo, en
; vez de tener que describir el sintoma por telefono.
Name: "{group}\{#NOMBRE} - autocomprobacion"; Filename: "{app}\spinpy.exe"; \
  Parameters: "--autocomprobacion"; \
  Comment: "Comprueba que todos los calculos funcionan en este equipo"
Name: "{group}\Desinstalar {#NOMBRE}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#NOMBRE}"; Filename: "{app}\spinpy.exe"; Tasks: escritorio

[Run]
Filename: "{app}\spinpy.exe"; Description: "Abrir {#NOMBRE} ahora"; \
  Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Los resultados del usuario viven en Documentos\spinpy y NO se tocan al
; desinstalar: son sus datos, no los de la aplicacion.
Type: filesandordirs; Name: "{app}\_internal\__pycache__"
