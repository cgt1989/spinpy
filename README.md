<p align="center">
  <img src="docs/media/banner.png" alt="spinpy — microestructuras espinodales ajustadas a hueso trabecular de micro-CT" width="100%">
</p>

<p align="center">
  <img alt="versión 0.2.0" src="https://img.shields.io/badge/versi%C3%B3n-0.2.0-e8a33d">
  <img alt="Python 3.10+" src="https://img.shields.io/badge/python-3.10%2B-3776ab?logo=python&logoColor=white">
  <img alt="Windows" src="https://img.shields.io/badge/ejecutable-Windows%2064%20bits-0078d6?logo=windows&logoColor=white">
  <img alt="Código MIT" src="https://img.shields.io/badge/c%C3%B3digo-MIT-2ea44f">
  <img alt="Ejecutable GPL-3.0" src="https://img.shields.io/badge/ejecutable-GPL--3.0-8a8a8a">
  <img alt="Verificación: 22 bloques" src="https://img.shields.io/badge/verificaci%C3%B3n-22%20bloques-5b6b7f">
  <img alt="Interfaz ES/EN" src="https://img.shields.io/badge/interfaz-ES%20%7C%20EN-5b6b7f">
</p>

<p align="center">
  <b>Del micro-CT al modelo mecánico, diciendo cuánto vale cada número.</b><br>
  <a href="#qué-hace">Qué hace</a> ·
  <a href="#la-aplicación-en-movimiento">Demo</a> ·
  <a href="#instalación">Instalación</a> ·
  <a href="#uso">Uso</a> ·
  <a href="#validación-contra-los-resultados-publicados">Validación</a> ·
  <a href="#citar">Citar</a>
</p>

---

**spinpy** ajusta microestructuras **espinodales** (y *dual-lattice*) a
volúmenes de interés de **hueso trabecular** obtenidos por microtomografía, y
después los mide y los ensaya: morfometría, homogeneización elástica, ensayo
de compresión por elementos finitos e informe listo para publicación. Es el
port a Python de la mitad numérica de `AppFinal_V2.m`, verificado contra
soluciones analíticas y contra la implementación MATLAB original.

Desarrollado sobre huesos sesamoideos equinos escaneados a 51,489 µm, con
soporte también para volúmenes porcinos, de vértebra de ratón y de fémur.

> *In English:* spinpy fits spinodoid and dual-lattice surrogate
> microstructures to trabecular-bone VOIs from micro-CT, then measures them
> (BV/TV, Tb.Th, DA, Conn.D, SMI, Ellipsoid Factor…), homogenizes their
> elastic tensor and runs FE compression tests. The GUI is available in
> Spanish and English.

<table>
<tr>
<td>

```
pila TIFF de micro-CT
 → segmentación, encuadre PCA
 → VOI cúbico
 → morfometría
     BV/TV, Tb.Th, DA, Conn.D,
     SMI, Ellipsoid Factor
 → ajuste: spinodoide
     o dual-lattice
 → homogeneización elástica
 → ensayo de compresión FE
     von Mises, Pistoia
 → informe y exportación
     PDF, STL, VTU, Abaqus, ANSYS, FEBio
```

</td>
<td align="center">
<img src="docs/media/giro.gif" alt="Spinodoide ajustado girando" width="300">
</td>
</tr>
</table>

---

## La aplicación en movimiento

Las animaciones son capturas de la propia interfaz pulsando sus botones, y se
regeneran con `python docs/media/hacer_medios.py`: no son maquetas.

### 1 · Explorar el espacio de diseño

La densidad, el número de onda y los tres ángulos de cono definen la
estructura. Cada cambio se regenera al momento, de la clase isótropa a la
columnar, la lamelar o la ajustada a un hueso real.

<p align="center"><img src="docs/media/01_explorar.gif" alt="Explorando clases de spinodoide en el visor" width="100%"></p>

### 2 · Ajustar al VOI de micro-CT

Se carga el VOI (`.vtk`, `.mat`, carpeta de rebanadas TIFF o TIFF
multipágina), se ajusta la microestructura por búsqueda escalonada y la tabla
compara métrica a métrica el hueso con su candidato. Las dos vistas giran
sincronizadas.

<p align="center"><img src="docs/media/02_ajuste.gif" alt="Ajuste de un spinodoide a un VOI equino" width="100%"></p>

### 3 · Ensayar a compresión

Ensayo de elementos finitos en X, Y, Z o los tres ejes, con campo de
deformación efectiva y de tensión de von Mises, **en la misma escala de color
para el hueso y el candidato**, y carga de fallo por el criterio de Pistoia.

<p align="center"><img src="docs/media/03_mecanica.gif" alt="Campo de von Mises del ensayo de compresión" width="100%"></p>

---

## Qué hace

| Etapa | Qué incluye |
|---|---|
| **Lee** | VOIs `.vtk` y `.mat`, y pilas TIFF de micro-CT (carpeta de rebanadas o multipágina) con la escala leída del *log* SkyScan. Umbraliza, encuadra por PCA y recorta el cubo. |
| **Genera** | Espinodales por campo aleatorio gaussiano (Kumar et al. 2020), con las dos convenciones de muestreo de ondas de la literatura, y *dual-lattice* (Vafaeefar et al. 2022). |
| **Mide** | BV/TV, BS/BV, BS/PV, Tb.Th, Tb.Sp, Tb.N, tensor MIL y DA, Conn.D, SMI, espesor local, tamaño de poro, **Ellipsoid Factor** y curvaturas principales de la interfaz. |
| **Ajusta** | Los parámetros que mejor reproducen un VOI real: búsqueda escalonada, réplicas con semillas nuevas, suelo de ruido autoconsistente y desempate mecánico opcional. |
| **Homogeneiza** | El tensor elástico por celda unidad periódica sobre la rejilla de vóxeles, con multigrid algebraico y comprobación del residuo. |
| **Ensaya** | Compresión en X, Y o Z con campos de von Mises, deformación efectiva y desplazamiento; fallo por Pistoia; estudio de convergencia de malla. |
| **Simula** | Pérdida ósea *in silico* (adelgazamiento, trabéculas finas, desuso, recuperación) y fallo progresivo sobre el gemelo digital del VOI. |
| **Informa** | Un informe para publicación en español e inglés (Markdown y PDF): métodos redactados con los valores usados, figuras a 600 ppp, tabla de citabilidad y huella SHA-256 para reproducir cada máscara. |
| **Exporta** | Sólido hexaédrico o TET10 a Abaqus, ANSYS APDL, VTU y STL, y el ensayo de compresión completo a FEBio 4 (`.feb`). |
| **Procesa lotes** | Descompone la varianza entre y dentro de especímenes, con ICC, N efectivo y equivalencia por TOST. |

Todo con interfaz gráfica (`visor.py`) o desde Python.

---

## Qué lo distingue

La mayoría de las herramientas de morfometría ósea dan un número. Esta
además dice **cuánto vale ese número**:

- **Verificación contra soluciones manufacturadas**, no solo contra otro
  programa: esferas, cilindros y láminas de respuesta analítica conocida, el
  promedio de Backus para el tensor elástico y el tensor MIL contra geometrías
  de anisotropía exacta.
- **Los sesgos están medidos y se reportan pegados al número.** El espesor
  local subestima un ~30 % a 3-4 vóxeles; el SMI arrastra un −15 % de
  discretización y se confunde con la concavidad; el MIL tiene un suelo de
  ruido de DA ≈ 1,07. La interfaz lo muestra donde se lee el valor.
- **Lo que no es citable, se dice.** El máximo de un campo de von Mises no
  converge con la malla, así que el informe cita el p99 de la capa de
  superficie y marca cada valor como *citable*, *con reservas* o *no citable*.
- **El generador es estocástico y eso se mide**, no se supone: la misma
  parametrización con K semillas da el suelo de ruido por debajo del cual
  ninguna diferencia significa nada.
- **Las réplicas no inflan el N.** El módulo estadístico calcula el N efectivo
  y lo pone por delante del número de filas, porque tratar K réplicas como K
  especímenes es pseudorreplicación.
- **Cada número guardado lleva su procedencia**: versión, familia, esquema,
  semilla y el número de onda con sus dos lecturas, para que un resultado
  exportado se pueda regenerar bit a bit.

---

## Instalación

### Sin Python: ejecutable de Windows

Para quien solo quiere **usar** la aplicación hay un ejecutable de 64 bits que
no necesita Python ni ninguna dependencia. Se instala en la carpeta del
usuario y **no pide contraseña de administrador**, que es lo que hace falta en
un equipo de universidad.

| Ejecutable | |
|---|---|
| Instalado | ~650 MB |
| Instalador / zip | ~220-300 MB |
| Guarda en | `Documentos\spinpy` (no se borra al desinstalar) |
| Licencia | **GPL-3.0** (ver abajo) |

Si algo no funciona en ese equipo, el menú de inicio trae un acceso directo
**«spinpy - autocomprobación»**: ejercita uno por uno todos los caminos de
cálculo y deja un informe en `Documentos/spinpy/autocomprobacion.txt`. Es lo
que hay que enviar para diagnosticar. También se lanza desde la consola:

```
spinpy.exe --autocomprobacion
```

Para reconstruirlo, `instalador/construir.ps1` hace el proceso entero, y
`instalador/LEEME_DESARROLLO.md` explica las trampas del empaquetado.

### Con Python: como biblioteca

```
pip install -e .
```

El núcleo solo necesita `numpy`, `scipy`, `scikit-image` y `pyvista`. Los
extras se instalan según haga falta:

```
pip install -e ".[elastico]"   # multigrid: homogeneización y ensayo FE
pip install -e ".[malla]"      # mallado TET10
pip install -e ".[lote]"       # lote y estadística
pip install -e ".[gui]"        # interfaz gráfica
pip install -e ".[informe]"    # PDF y figuras del informe
pip install -e ".[todo]"       # todo
```

`pyamg` figura como extra, pero no es opcional en la práctica: sin él, el
solver cae a gradiente conjugado con Jacobi, que con el contraste de 10⁻⁶
entre hueso y vacío no converge.

---

## Uso

Interfaz gráfica:

```
python visor.py
```

Los paneles siguen el orden del trabajo: VOI de referencia, microestructura,
visualización, morfometría, ajuste, análisis mecánico, simulaciones *in
silico*, exportación y lote. El botón **«Informe para publicación (auto)…»**
encadena todas las etapas y estima de antemano cuánto va a tardar.

### Idioma

La interfaz está en **español e inglés** y se cambia desde el menú *Idioma*
sin reconstruir la ventana: el VOI cargado, las métricas y los campos se
conservan. No se traducen, a propósito, los nombres de métrica (BV/TV, Tb.Th,
DA, SMI…), que son notación internacional (Bouxsein et al. 2010).

Después de tocar cualquier texto hay que pasar el comprobador del diccionario:

```
python idioma_revisar.py
```

### Desde Python

```python
from spinpy import leer_voi, morfometria, ajustar_spinodoide

VOI, spacing = leer_voi("VOI_medio_cubico.vtk")
m = morfometria(VOI, spacing, extra=True)
print(m["BVTV"], m["TbTh"], m["DA"], m["ConnD"])

r = ajustar_spinodoide(VOI, spacing, modo="completo")
print(r["parametros"], r["error"])
```

### Línea de comandos

```
spinpy-lote carpeta_con_VOIs --replicas 20 --salida resultados
spinpy-informe resultados_sesion.json        # rehace el informe desde un JSON
```

---

## Validación contra los resultados publicados

El método que implementa `spinpy` **no es nuestro**: es el de Kumar et al.
(2020). La carpeta `Test/` existe para demostrarlo, poniendo la figura
publicada al lado de la misma figura regenerada con este código, y evaluando
predicciones falsables de cada artículo:

| Guion | Qué comprueba |
|---|---|
| `replicar_kumar2020.py` | El método: unión de conos, conjunto de nivel y las cuatro clases. |
| `replicar_zheng2021.py` | Cotas: cada superficie elástica dentro de Voigt y Hashin-Shtrikman, y la curva de ρ ≥ 0,3. |
| `replicar_guo2024.py` | Curvatura: el perfil (κ₁, κ₂) frente a una superficie nodal de curvaturas exactas. |

Se abren desde la aplicación en el menú **Validación**, o desde la consola:

```
python Test/replicar_kumar2020.py            # completo, ~7 min
python Test/replicar_kumar2020.py --rapido   # ~1,5 min
```

> Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
> spinodoid metamaterials. *npj Computational Materials, 6*, 73.
> https://doi.org/10.1038/s41524-020-0341-6 — acceso abierto bajo
> CC BY 4.0, que es lo que permite reproducir aquí su figura con atribución.
> Procedencia de cada figura de referencia en `Test/referencia/FUENTES.md`.

---

## Verificación

```
python -m pytest tests/ -q
```

22 bloques con tolerancias **declaradas antes de medir**: topología, SMI,
espesor, mecánica, Pistoia, pilas TIFF, curvatura, Ellipsoid Factor,
*dual-lattice*, procedencia, función objetivo, muestreo MIL, informe, capa de
superficie de von Mises… Un fallo aquí es un hallazgo, no un error de la
suite: los bloques 04 y 05 tienen fallos conocidos y documentados.

Contra la implementación MATLAB original (`validar_*.py`, salidas de
referencia en `resultados/`):

| Qué se comprueba | Resultado |
|---|---|
| Morfometría sobre VOIs reales | 5,9 × 10⁻¹⁴ % en lo que no pasa por marching cubes |
| Tensor elástico frente al promedio de Backus | 2,7 × 10⁻¹⁵ de diferencia relativa |
| Función de error, sobre 1173 pares | 3,55 × 10⁻¹⁶ |
| Rejillas de búsqueda del ajuste | exactas |
| von Mises en estado uniaxial | 1,000000000000 |

---

## Comparación con BoneJ

Contrastado con BoneJ 7.2.2 sobre los mismos volúmenes:

- **BV/TV coincide a precisión de máquina**: las dos herramientas ven la misma
  máscara.
- **BS difiere hasta un 46 %, y es una diferencia de convención**: BoneJ
  cierra la superficie contra el borde del volumen y spinpy no, porque las
  caras del cubo son el corte artificial del VOI. Igualando la convención, la
  diferencia baja a ±6 %.
- **DA usa escalas distintas** en las dos herramientas, aunque ambas coinciden
  en la ordenación de los especímenes.
- **SMI no se puede comparar**: BoneJ lo retiró en su versión 2, por la misma
  crítica de confusión con la concavidad que spinpy documenta.

Y algo que cambia cómo leer todo lo anterior: **la incertidumbre de
segmentación es un orden de magnitud mayor que las diferencias entre
implementaciones.** Mover el umbral un ±15 % mueve Tb.Th un 62 %, frente al
±6 % que separa a las dos herramientas.

---

## Comparación con FEBio

El ensayo de compresión de spinpy, exportado a FEBio 4.5 y resuelto allí, sobre
tres VOIs de hueso (BV/TV 0,28 a 0,77, a 32³ y 48³) y un espinodoide ajustado:

- **El cálculo coincide en siete u ocho cifras.** Módulo aparente, campo de
  desplazamientos, tensiones, p99 de von Mises en la capa superficial y carga
  de fallo de Pistoia difieren entre 10⁻⁶ y 10⁻¹⁰, dentro de tolerancias
  declaradas antes de medir. Los dos programas plantean el mismo problema
  discreto y lo resuelven con métodos distintos (multigrid frente a
  factorización directa).
- **FEBio es no lineal geométricamente y spinpy es lineal**, y eso tiene
  consecuencias en el hueso más poroso. En los VOIs de BV/TV 0,55 y 0,77 la
  respuesta no lineal a 1 MPa se aparta un 0,06 % de la lineal. En el de
  BV/TV 0,28 se aparta un 16 % a 32³ y un 48 % a 48³. Los tres VOIs son de
  un solo espécimen (el sesamoideo H4): es un caso medido, no una tasa
  general.
- **La causa es la forma de cargar, no la estructura.** Con fuerza impuesta,
  las trabéculas cortadas por la cara cargada del VOI trabajan como
  voladizos. Con un plato rígido el desvío baja al 0,8 %, y en ese VOI el
  propio módulo lineal es 2,1 veces mayor. En VOIs muy porosos, el módulo
  aparente y la carga de fallo deben citarse declarando la condición de
  carga.

Lo que esta comparación **no** demuestra: coincidir con FEBio prueba que
spinpy resuelve bien el problema que plantea, no que ese problema represente
el hueso real. La malla de vóxeles, el tejido homogéneo de 20 GPa y las
condiciones de contorno son supuestos que los dos programas comparten; eso
solo lo contrasta un ensayo físico.

`tests/test_25_febio.py` reproduce la comparación en pequeño y se salta si
FEBio no está instalado.

---

## Documentación

`docs/MANUAL_spinpy.pdf` documenta cada módulo y cada función, y se genera
del propio código (`python docs/generar_manual.py`). Los docstrings de este
proyecto no son un resumen: llevan las decisiones de diseño, los sesgos
medidos y las trampas que costó encontrar.

---

## Limitaciones que conviene conocer

- **La segmentación es un umbral global.** Es la mayor fuente de
  incertidumbre de todo el proceso, y por eso la aplicación declara cuál se
  aplicó en vez de esconderlo.
- **La familia espinodal tiene un límite propio**: una sola longitud
  característica da un espesor trabecular uniforme, y el hueso no lo es. Un
  spinodoide que iguala BV/TV, Tb.Th y DA puede ser bastante más blando que el
  hueso al que se ajustó; spinpy lo mide y lo dice.
- Los optimizadores de Pareto, bayesiano y MOBO siguen solo en MATLAB.
- Por debajo de ρ ≈ 0,25 (clase isótropa) la homogeneización periódica no
  alcanza la tolerancia declarada: es una propiedad del régimen cercano al
  umbral de rigidez, y la aplicación lo reporta.

---

## Citar

Si usas este software, cítalo con los metadatos de [`CITATION.cff`](CITATION.cff).

**Autores:** Carlos González-Torres
([ORCID](https://orcid.org/0000-0002-7765-6388)), David Ortiz-Puerta
([ORCID](https://orcid.org/0000-0001-6285-3066)) y Mauricio A.
Sarabia-Vallejos ([ORCID](https://orcid.org/0000-0001-5128-796X)) —
Universidad de Valparaíso.

Métodos implementados: Kumar et al. 2020 (generador espinodal), Vafaeefar et
al. 2022 (*dual-lattice*), Andreassen & Andreasen 2014 (homogeneización),
Harrigan & Mann 1984 (tensor MIL), Pistoia et al. 2002 (criterio de fallo),
Odgaard & Gundersen 1993 (Conn.D), Hildebrand & Rüegsegger 1997 (espesor local
y SMI), Doube 2015 (Ellipsoid Factor), Parfitt (modelo de placas).

## Licencia

**El código fuente es MIT.** Ver [`LICENSE`](LICENSE). Quien importe el
paquete `spinpy` desde un script o un cuaderno no arrastra ninguna obligación
de copyleft: el núcleo no depende de Qt.

**El ejecutable de Windows es GPL-3.0**, porque incluye PyQt5, que es GPL v3.
Es una obra combinada y su redistribución queda sujeta a esa licencia; el uso,
en cambio, no tiene restricción alguna, y lo que produzcas con él es tuyo. El
detalle, con la lista de componentes de terceros, está en
`instalador/LICENCIA_BINARIO.txt`.

Migrar la interfaz a PySide6 (LGPL) no bastaría por sí solo para tener un
binario sin copyleft: `pymeshfix` es GPL v3 y `tetgen` deriva de código AGPL.
