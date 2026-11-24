# Contribuir a spinpy

*English summary: bug reports, FEM engines and morphometry fixes are welcome.
Run `python -m pytest tests/ -q -m "not lento"` before opening a pull request.
Adding a FEM engine means one file plus one registry entry, see
`docs/AGREGAR_UN_MOTOR_FEM.md`. The code and its comments are in Spanish; issues
and pull requests in English are fine.*

Gracias por querer mejorar `spinpy`. Esta guía dice cómo empezar según lo que
sepas hacer. No hace falta entender todo el programa: casi cualquier aporte
cabe en un módulo.

## Probar el código sin instalar nada

- **Navegador**: abre `ejemplos/spinpy_inicio.ipynb` en Colab (botón en el
  README). Instala el paquete, carga un VOI porcino abierto, mide su
  morfometría y ensaya a compresión.
- **Contenedor de desarrollo**: en VS Code o en GitHub Codespaces, «Reabrir en
  contenedor». El entorno queda listo, con las pruebas verdes.
- **Ejecutable de Windows, macOS o Linux**: sin Python, en la
  [página de releases](https://github.com/cgt1989/spinpy/releases).

## Instalar para desarrollar

```
git clone https://github.com/cgt1989/spinpy.git
cd spinpy
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[elastico,fem,malla,lote,test,informe]" tifffile
python -m pytest tests/ -q -m "not lento"
```

La interfaz (`python visor.py`) necesita además `pip install -e ".[gui]"`. El
paquete en sí no necesita Qt.

## Qué puedes aportar y por dónde empezar

| Si sabes de… | Aporte típico | Guía |
|---|---|---|
| **Elementos finitos** | un motor nuevo (otra biblioteca o un resolvedor propio); un elemento; un material no lineal | `docs/AGREGAR_UN_MOTOR_FEM.md` |
| **Morfometría ósea** | una métrica, una corrección de sesgo, una validación frente a BoneJ o CTAn | `morphometry.py` y su cabecera |
| **Procesamiento de imagen** | segmentación, lectura de otros formatos de micro-CT | `voi.py`, `io.py` |
| **Python científico** | rendimiento, pruebas, empaquetado, tipos | `docs/ARQUITECTURA.md` |
| **Un lenguaje que no es español** | traducir la interfaz | `spinpy/idioma_textos.py` |
| **Nada de lo anterior, pero encontraste un error** | un *issue* con el formulario de errores | abajo |

Para orientarte: `docs/ARQUITECTURA.md` es el mapa del código, y la cabecera de
cada módulo explica por qué existe y qué decisiones lo moldearon.

## Reportar un error

Abre un *issue* con el formulario «Error». Lo que más ayuda:

1. El resultado de `python visor.py --autocomprobacion` o, en el ejecutable, el
   archivo `autocomprobacion.txt` (el acceso directo del menú de inicio lo
   genera).
2. Un VOI mínimo que lo reproduzca. Los diez VOIs porcinos de
   `datos/porcino_koria2020/` valen si el error aparece con ellos.
3. Lo que esperabas, con la fuente: una definición, una solución analítica, el
   resultado de otro programa.

## Reglas del código

1. **Cada cambio numérico lleva su prueba**, en `tests/test_NN_tema.py` o en un
   bloque nuevo. La tolerancia se **declara antes de medir**, con su
   justificación, en la cabecera del módulo de pruebas. No se relaja una
   tolerancia para que pase una prueba: un fallo es un hallazgo.
2. **Se verifica contra una respuesta conocida** (solución analítica,
   definición publicada, invariante), no contra el número que dio ayer el mismo
   código.
3. **Una sola medición**: las métricas pasan por `morphometry.morfometria`
   (o `morfometria_dominio`), para que VOI y candidato recorran el mismo camino.
4. **El paquete no importa Qt.** La lógica va en `spinpy/`; `visor.py` solo
   orquesta y dibuja.
5. **Texto de la interfaz**: escríbelo en español dentro de `_("…")` y añade su
   versión inglesa a `spinpy/idioma_textos.py`. `python idioma_revisar.py` debe
   dar 0 sin traducir.
6. **Si cambia el número de bloques de `tests/`**, actualiza la insignia y el
   texto del README (la prueba 19 lo comprueba).
7. **Documenta el porqué**, no el qué. Un comentario que explica una decisión
   (con el valor medido que la justifica) vale más que diez que repiten el
   código.

## Antes de abrir el pull request

```
python -m pytest tests/ -q -m "not lento"
python idioma_revisar.py        # solo si tocaste textos de la interfaz
ruff check .
```

El flujo de GitHub Actions (`.github/workflows/pruebas.yml`) corre la suite
completa en Linux, Windows y macOS, y la construcción del instalador exige que
pase.

## Licencias

El código es MIT. El ejecutable que incluye PyQt5 se distribuye bajo GPL-3.0
(ver `LICENSE` y el README). Al aportar, aceptas que tu contribución se
publique bajo la licencia del repositorio.
