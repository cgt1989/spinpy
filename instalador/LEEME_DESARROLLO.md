# Empaquetado de spinpy para Windows

Notas para quien tenga que reconstruir o arreglar el ejecutable. No hacen
falta para usarlo.

```powershell
.\construir.ps1
```

Eso monta el entorno, genera el icono, comprueba el codigo, empaqueta,
comprueba **el binario**, hace el zip y, si Inno Setup esta instalado, el
instalador.

---

## Las tres trampas del empaquetado

Ninguna da error al compilar. Las tres producen un ejecutable que se construye
limpio y falla despues, que es la peor forma de fallar.

### 1. `__file__` deja de ser un sitio donde escribir

`visor.py` derivaba de `Path(__file__).parent` la carpeta de resultados, la de
exportacion y la de sesiones. Empaquetado, eso apunta a la carpeta de
instalacion —de solo lectura si esta bajo Archivos de programa— o, en modo
archivo unico, al directorio temporal que Windows **borra al cerrar**. El
usuario exportaria sus resultados, la aplicacion diria que los ha guardado, y
no estarian.

Arreglo: se separan `RAIZ` (donde vive el codigo) y `DATOS` (donde se
escribe). Empaquetado, `DATOS` es `Documentos\spinpy`, y la carpeta se
localiza preguntandole a Windows por el identificador de carpeta conocida, no
componiendo `~/Documents`: con OneDrive activo —el caso de esta maquina— la
ruta real es `~/OneDrive/Documentos` y la ingenua no existe.

### 2. Las DLL del interprete de Anaconda no viajan

Los `.pyd` de la biblioteca estandar de Anaconda enlazan contra DLL propias de
conda que viven en `anaconda3\Library\bin`, una carpeta que PyInstaller no
mira. Resultado:

```
ImportError: DLL load failed while importing _ctypes
```

sin decir cual falta. Medido aqui: `_ctypes.pyd` necesita `ffi.dll` y el propio
`python313.dll` necesita `zlib.dll`; con el cierre transitivo salen 12 DLL, 4,4 MB
en total.

`dlls_conda.py` las resuelve recorriendo la tabla de importaciones de cada
binario recogido, y el `.spec` lo llama. Se hace asi y no con una lista escrita
a mano porque la lista cambia con la version de Anaconda y con lo que haya
instalado: una lista fija caduca en silencio.

Esto desaparece construyendo desde un CPython de python.org en vez de desde
Anaconda. Si el empaquetado se muda a integracion continua, alli lo correcto es
el interprete oficial y `dlls_conda.py` sobra.

### 3. joblib se multiplica sola

La seccion de lote usa `joblib` con `prefer="processes"`. `loky` arranca a sus
obreros reejecutando `sys.executable` con `-c "from loky ..."`; en un
ejecutable congelado `sys.executable` es la propia aplicacion, que no acepta
`-c`. Los obreros mueren nada mas nacer, joblib levanta
`TerminatedWorkerError` y el proceso padre se queda colgado —medido: 24 s hasta
el error, y el padre siguio vivo—.

Y sin `multiprocessing.freeze_support()` es peor: cada hijo vuelve a entrar por
`main()`, abre otra ventana y lanza mas hijos.

Arreglo: `freeze_support()` como primera linea de `main()`, y
`spinpy.lote.opciones_paralelo()` devuelve hilos cuando `sys.frozen`. Se pierde
parte del paralelismo, pero no todo: generar la mascara, la morfometria y la
homogeneizacion pasan casi todo su tiempo dentro de numpy, scipy y pyamg, que
sueltan el GIL. **El lote tarda mas en el ejecutable que desde el codigo.**

---

## La autocomprobacion no es un adorno

`spinpy.exe --autocomprobacion` existe porque un ejecutable **arranca aunque le
falte media aplicacion**: `tetgen`, `pymeshfix` y `pyamg` no se importan hasta
que el usuario pulsa el boton correspondiente. Un fallo de empaquetado en el
mallador no se ve al abrir la ventana, sino media hora despues, en mitad de una
exportacion, en forma de traza que quien la recibe no sabe leer.

Recorre los catorce modulos y los nueve caminos de calculo, y escribe un
informe. `construir.ps1` la ejecuta dos veces —desde el codigo antes de
empaquetar, y sobre el binario despues— y aborta si alguna falla.

Fue esa segunda ejecucion la que destapo el problema de joblib. Sin ella el
instalador habria salido con la seccion de lote rota.

## Por que las versiones van fijadas

`requirements.txt` declara minimos, correcto para quien instala la biblioteca.
`requisitos_build.txt` fija versiones exactas, y es a proposito: todo el valor
del proyecto esta en que los numeros estan verificados contra soluciones
cerradas y contra la app de MATLAB, y esa verificacion se hizo con esas
versiones. Quien recibe el ejecutable no tiene Python y no puede comprobar
nada.

Montando el entorno sin fijar nada salieron numpy 2.5.2, pandas 3.0.5 y
matplotlib 3.11.1 frente a las 2.1.3 / 2.2.3 / 3.10.0 validadas.

Comprobado: la autocomprobacion del binario devuelve exactamente los mismos
numeros que desde el codigo —BV/TV 0.3156, Tb.Th 0.05038 mm, DA 1.168,
Ex 532.9 / Ey 1052.4 / Ez 1381.0 MPa, E_app 428.8 MPa, 66374 TET10 al 97.4 %
del volumen—. El empaquetado no cambia ningun resultado.

## Modo carpeta, no archivo unico

`--onefile` descomprime ~650 MB en el temporal **en cada arranque**: 20-40 s
antes de ver la ventana, todas las veces. En modo carpeta el arranque baja a
unos segundos, y el usuario no ve la carpeta porque el instalador la deja en su
sitio y lo que el toca es un acceso directo.

## Tamanos medidos

| | |
|---|---|
| Carpeta `dist\spinpy` | 657 MB |
| Zip portable | 219 MB |
| `vtk.libs` | 265 MB |
| `PyQt5` | 131 MB |
| `scipy` | 63 MB |

## Si el ejecutable no arranca

Reconstruir con consola para ver la traza, que en modo ventana se pierde tras
un cuadro que solo dice *«Unhandled exception in script»*:

```powershell
$env:SPINPY_CONSOLA = "1"
.\construir.ps1 -SaltarEntorno -SoloZip
```

Y para saber que DLL le faltan a una carpeta ya construida, sin recompilar:

```powershell
python dlls_conda.py C:\ruta\dist\spinpy\_internal
```
