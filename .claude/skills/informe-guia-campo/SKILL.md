---
name: informe-guia-campo
description: Redacta un informe técnico en formato de guía de campo, en PDF maquetado desde HTML con Chromium. Lleva portada, resumen ejecutivo de una página con cifras destacadas y tabla semáforo, una ficha por hallazgo (qué es, por qué aparece, qué medimos, cómo resolverlo), hoja de ruta priorizada, glosario, referencias verificadas y lista de chequeo de la guía de reporte. Úsala cuando el usuario pida un informe «como el de artefactos», «en formato de guía de campo», «con fichas», o un informe técnico explicativo en PDF que no sea el Markdown convertido de los informes anteriores.
---

# Informe técnico en formato de guía de campo

Informe de diagnóstico pensado para que lo entienda alguien que no hizo el
trabajo. Cada hallazgo se presenta siempre igual: qué es, por qué aparece,
qué se midió y cómo resolverlo. El ejemplo completo es
`comparativa_motores/informe_artefactos/` (PDF en
`comparativa_motores/INFORME_ARTEFACTOS.pdf`). Léelo antes de empezar.

## Archivos de esta skill

| Archivo | Para qué |
|---|---|
| `plantilla.html` | Esqueleto con cada componente y huecos en MAYÚSCULAS |
| `estilo.css` | Maquetación A4, tipografía, cajas, fichas, insignias de gravedad |
| `fuentes/` | Source Serif 4, IBM Plex Sans y Mono en local (OFL 1.1) |
| `estilo_figuras.py` | `aplicar()` y `guardar()` para figuras matplotlib coherentes |
| `imprimir.sh` | HTML a PDF con Chromium sin interfaz |

## Procedimiento

1. **Datos primero.** El informe no genera cifras: todas salen de scripts
   versionados (p. ej. `estudio.py` → `resultados/*.json`). Incluye un
   control que reproduzca un valor ya publicado, para que el banco de
   pruebas quede validado.
2. **Carpeta del informe.** Crea `<carpeta>/informe_<tema>/` y copia ahí
   `plantilla.html` (renómbralo), `estilo.css`, `imprimir.sh` y `fuentes/`.
   Cambia en `estilo.css` el texto del pie (`TÍTULO CORTO DEL INFORME`).
3. **Figuras.** Escribe un `figuras_<tema>.py` que use `estilo_figuras.py` y
   guarde PNG y SVG en `informe_<tema>/figs/`. Reglas: paleta en orden fijo
   (azul, naranja, aqua; más series se resuelven con facetas), un solo eje y
   por gráfico, leyenda siempre que haya dos series o más, texto en tinta y
   nunca en el color de la serie, coma decimal. Abre cada PNG y revísalo
   (solapamientos, etiquetas cortadas) antes de seguir.
4. **Referencias.** Búscalas en PubMed (y Consensus si hace falta) y
   comprueba con el resumen que cada una respalda la frase donde se cita.
   Formato Vancouver con DOI. No inventes autores, DOI ni resultados.
5. **Redacción.** Rellena la plantilla en este orden:
   - **Portada**: tipo de informe, título, subtítulo de una frase, contexto,
     recuadro «Preguntas que responde» y pie con versión, commit y archivos.
   - **El informe en una página**: entradilla, cuatro cifras destacadas,
     tabla semáforo (código, hallazgo, a qué afecta, magnitud medida,
     gravedad) y caja «Lo esencial» con las tres acciones prioritarias.
   - **Concepto**: definición operativa y un diagrama SVG propio que sitúe
     cada hallazgo en la cadena del proceso.
   - **Método**: datos, ensayos (tabla con cifras) y una caja «Cómo se lee»
     por métrica.
   - **Panorama**: figuras de conjunto con caja «Cómo leer la figura».
   - **Una ficha por hallazgo**, cada una en página nueva: cabecera con
     código, insignia de gravedad y etiquetas de lo que afecta; «Qué es»,
     «Por qué aparece», «Qué medimos»; columna lateral con un esquema SVG
     propio y un bloque «En cifras»; figura con pie; «Cómo resolverlo» con
     opciones numeradas, cada una con su esfuerzo y lo que corrige, y la
     ubicación en el código cuando exista.
   - **Síntesis** (qué es robusto y qué no), **hoja de ruta** (acción,
     corrige, dónde en el código, esfuerzo), **experimentos que faltan**,
     **limitaciones**, **glosario**, **reproducir** (órdenes exactas),
     **referencias** y, como anexo, la **lista de chequeo** de la guía de
     reporte que corresponda al diseño (identifica el diseño y verifica la
     versión vigente de la guía).
6. **Imprimir.** `ENTRADA=informe_<tema>.html SALIDA=../INFORME_<TEMA>.pdf sh imprimir.sh /opt/pw-browsers/chromium-1194/chrome-linux/chrome`
   (en otro equipo, la ruta de Chrome o Chromium). Comprueba con `pdffonts`
   que las fuentes están incrustadas.
7. **Revisión página a página.** Renderiza con `pdftoppm -r 70 -png` y mira
   todas las páginas: sin solapamientos, sin cajas cortadas, sin páginas
   casi vacías (si una figura salta de página, reduce `max-height` o mueve
   texto). Después contrasta **cada cifra del texto** con los JSON de
   resultados y corrige las que no coincidan exactamente.
8. **Versionar.** Se versionan scripts, JSON de resultados, HTML, figuras,
   fuentes y PDF. Los campos grandes regenerables van a `.gitignore`.

## Redacción

- Español académico natural, sin el carácter em dash (U+2014). Gravedad siempre con su
  etiqueta además del color.
- Distingue asociación de causalidad y declara lo que el diseño no permite
  separar. No presentes como cota lo que no lo es en sentido estricto.
- Las cajas explican; el cuerpo argumenta. Cada figura tiene un pie que se
  entiende sin leer el texto.
- Prefiere cifras concretas («el 21 % del hueso, a menos de 0,3 mm de una
  cara») a calificativos («mucho», «significativo»).
