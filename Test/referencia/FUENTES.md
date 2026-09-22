# Procedencia y licencia de las figuras de referencia

Las imágenes de esta carpeta **no son nuestras**. Son figuras publicadas,
reproducidas aquí para poder poner la réplica al lado del original. Cada una
lleva abajo su cita completa, su licencia y de dónde se obtuvo exactamente.

Si alguna vez se añade una figura cuya licencia **no** permita
redistribuirla, no se pega aquí: se deja la referencia y se pide al lector que
la consulte en el artículo. Reproducir una figura con copyright cerrado en un
repositorio público es una infracción, por mucho que sea "para validar".

---

## `Kumar2020_Fig1.png` y `Kumar2020_Fig2.png`

**Cita (APA 7)**

> Kumar, S., Tan, S., Zheng, L., & Kochmann, D. M. (2020). Inverse-designed
> spinodoid metamaterials. *npj Computational Materials, 6*(1), 1–10.
> Artículo 73. https://doi.org/10.1038/s41524-020-0341-6

**Licencia**
Creative Commons Attribution 4.0 International (**CC BY 4.0**),
https://creativecommons.org/licenses/by/4.0/
© The Author(s) 2020. Publicado por Springer Nature en asociación con el
Shanghai Institute of Ceramics of the Chinese Academy of Sciences.

La declaración de licencia aparece literalmente en la última página del
artículo:

> *"Open Access This article is licensed under a Creative Commons Attribution
> 4.0 International License, which permits use, sharing, adaptation,
> distribution and reproduction in any medium or format, as long as you give
> appropriate credit to the original author(s) and the source, provide a link
> to the Creative Commons license, and indicate if changes were made."*

**Qué es cada figura**

| Archivo | Figura | Contenido |
|---|---|---|
| `Kumar2020_Fig1.png` | Fig. 1 | Diseño de metamateriales espinodoides: (a) vía clásica por campo de fase, (b) los ángulos de cono θ₁, θ₂, θ₃ y los N vectores de onda muestreados dentro de los conos, (c) la estrategia GRF + homogeneización |
| `Kumar2020_Fig2.png` | Fig. 2 | Seis topologías anisótropas (a–f) con su superficie elástica E(d), para ρ = 0.5 |

**Cómo se obtuvieron**
Extraídas del PDF de acceso abierto depositado por los propios autores en el
repositorio institucional de TU Delft
(https://repository.tudelft.nl/record/uuid:d78e6f8a-50ec-4f45-89bd-a32a60145801),
con PyMuPDF, sin recortar, escalar ni retocar. **No se ha modificado nada**,
que es lo que la licencia pide declarar.

**Nota sobre las etiquetas**
Los rótulos de los paneles (a–f) y las ternas (θ₁, θ₂, θ₃) de la Fig. 2 están
en la capa vectorial de la página, no en el mapa de bits, así que no aparecen
en el PNG extraído. Por eso `replicar_kumar2020.py` **no** intenta reproducir
la Fig. 2 panel por panel: reproduce las cuatro clases que el texto nombra,
con las ternas que el artículo sí publica en la leyenda de la Fig. 3.

---

## `Zheng2021_Fig1.jpeg`

**Cita (APA 7)**

> Zheng, L., Kumar, S., & Kochmann, D. M. (2021). Data-driven topology
> optimization of spinodoid metamaterials with seamlessly tunable anisotropy.
> *Computer Methods in Applied Mechanics and Engineering, 383*, 113894.
> https://doi.org/10.1016/j.cma.2021.113894

**Licencia**
Creative Commons Attribution (**CC BY**),
http://creativecommons.org/licenses/by/4.0/
© 2021 The Author(s). Publicado por Elsevier B.V.

La declaración aparece literalmente en la primera página del artículo:

> *"© 2021 The Author(s). Published by Elsevier B.V. This is an open access
> article under the CC BY license (http://creativecommons.org/licenses/by/4.0/)."*

**Qué es**

| Archivo | Figura | Contenido |
|---|---|---|
| `Zheng2021_Fig1.png` (y el `.jpeg` original) | Fig. 1 | (a) los ángulos de cono θ₁, θ₂, θ₃ y los vectores de onda muestreados; (b–e) las cuatro clases —lamelar, columnar, cúbica, isótropa— con su superficie elástica dentro de la esfera de Voigt y, en la isótropa, la de Hashin-Shtrikman |

**Por qué esta figura y no otra**
Es la única de los tres artículos que publica **las cuatro ternas (θ₁, θ₂, θ₃)
completas** en los pies de sus paneles, incluida la cúbica que la Fig. 2 de
Kumar et al. (2020) no da y que `replicar_kumar2020.py` tuvo que deducir. Y es
la única que dibuja las **cotas**, que es lo que permite comprobar magnitudes
absolutas y no solo proporciones.

**Cómo se obtuvo, y la única modificación que hay, declarada**
Extraída con PyMuPDF del PDF de acceso abierto (página 6, objeto de imagen
xref 107), sin recortar, escalar ni retocar. Se conservan **dos archivos**:

- `Zheng2021_Fig1.jpeg` — los bytes tal como vienen incrustados en el PDF, sin
  tocar. Es la evidencia de procedencia.
- `Zheng2021_Fig1.png` — los **mismos píxeles** re-codificados en PNG, que es
  el archivo que muestra la ventana de Validación.

La conversión es necesaria y no cosmética: el PyQt5 de este entorno **no trae
el complemento de imagen JPEG** (`QImageReader.supportedImageFormats()` no
incluye `jpeg`), de modo que la figura no se vería. Se comprobó que la
conversión es sin pérdida comparando los dos mapas de bits píxel a píxel:
son idénticos. No se ha recortado, escalado ni retocado nada; lo único que
cambia es el contenedor.

---

## `Guo2024_Fig7.png`

**Cita (APA 7)**

> Guo, Y., Sharma, S., & Kumar, S. (2024). Inverse designing surface curvatures
> by deep learning. *Advanced Intelligent Systems, 6*(6), 2300789.
> https://doi.org/10.1002/aisy.202300789

**Licencia**
Creative Commons Attribution (**CC BY 4.0**),
https://creativecommons.org/licenses/by/4.0/
© 2024 The Authors. Advanced Intelligent Systems publicado por Wiley-VCH GmbH.

La declaración aparece literalmente en la primera página del artículo:

> *"This is an open access article under the terms of the Creative Commons
> Attribution License, which permits use, distribution and reproduction in any
> medium, provided the original work is properly cited."*

**Qué es**

| Archivo | Figura | Contenido |
|---|---|---|
| `Guo2024_Fig7.png` | Fig. 7 | Tres casos de diseño inverso por curvatura: (a) hueso trabecular de Tozzi et al., (b) superficie espinodal con β = 15π, Q = 1000, ρ = 0.3, θ = (60°, 30°, 10°), (c) superficie nodal periódica sin x sin 1.8y + sin y sin 1.8z + sin z sin 1.8x = 0.5. De cada uno, el perfil de curvaturas «Target» y el «Reconstruction» de su red |

**Qué se replica de ella**
Solo las columnas **«Target»** de (b) y (c): son geometrías completamente
especificadas, con sus parámetros en la propia figura, y medibles con
`spinpy.curvatura`. El diseño inverso por redes neuronales y la muestra de
hueso de (a) no se replican, y `replicar_guo2024.py` dice por qué.

**Cómo se obtuvo**
Extraída con PyMuPDF del PDF de acceso abierto (página 9, objeto de imagen
xref 103), sin recortar, escalar ni retocar. **No se ha modificado nada.**
