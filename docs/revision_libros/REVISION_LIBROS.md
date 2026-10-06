# Revisión de la biblioteca de elementos finitos y análisis numérico para spinpy

**Versión evaluada:** spinpy V2.1.1 (rama `claude/fem-numerical-analysis-books-objtnc`)
**Fecha:** 6 de octubre de 2026
**Fuentes:** 23 libros y apuntes de la biblioteca privada del grupo (repositorio `cgt1989/LibrosFEM`, que no se redistribuye). Se citan por autor, sección y página impresa.

---

## Resumen

Se revisaron los 23 textos de la biblioteca, buscando dos cosas: soluciones para los problemas que la aplicación declara abiertos y mejoras para métodos que ya implementa. Cada propuesta se contrastó con el código de la V2.1.1 antes de incluirla, y se indica el archivo y la línea afectados. Los tres textos en formato DjVu (Bathe 1996, Bathe y Wilson 1976, Reddy 2002) se revisaron en una segunda fase (sección 6). Siete propuestas (R1 a R7) se implementaron y probaron; los guiones y los resultados están en `comparativa_motores/libros/`, y el informe, con extractos de los libros, se conserva solo en local.

Los libros con más aplicación directa son Zienkiewicz, Taylor y Zhu (2005), Nocedal y Wright (2006), Brenner y Scott (2008) y Šolín, Segeth y Doležel (2004). Tartar (2009), Atkinson y Han (2009), Modersitzki (2004) y los dos textos de cálculo de variaciones aportan justificación teórica para decisiones ya tomadas o para cambios concretos. El resto tiene poca relación con lo que hace spinpy (sección 5).

Cinco hallazgos se pueden aplicar con poco esfuerzo y merecen ir primero:

1. **El multigrid de la homogeneización no recibe los modos rígidos.** `elastic.py:273` construye el precondicionador sin el espacio casi nulo, a diferencia del ensayo de compresión (`resistencia.py:659`) y de los motores (`motores/_comun.py:113`). Es la explicación más barata de comprobar para la falta de convergencia por debajo de ρ ≈ 0,25.
2. **El Newton no lineal no tiene búsqueda lineal ni control de paso.** `motores/m_ngsolve.py:283-309` aplica el paso completo de Newton y, por omisión, toda la carga en un solo incremento (`fem.py:1045`, `pasos=1`). Nocedal y Wright describen la globalización que falta.
3. **El material no lineal por omisión es St. Venant-Kirchhoff**, que no es estable en compresión (Kamensky 2022, §4.3.4). El neo-Hookeano ya está implementado.
4. **El hexaedro de un vóxel es demasiado rígido en flexión.** Zienkiewicz et al. (§9.8-9.9) muestran el origen (deformaciones de corte espurias) y la corrección clásica con modos incompatibles, que no añade grados de libertad globales. Es la mejora con más impacto potencial sobre trabéculas de 1 a 3 vóxeles de espesor, el régimen en que trabajan los VOIs reales a 40³.
5. **Los criterios de parada de los iterativos no acotan el error.** Quarteroni et al. (§4.6.2) dan la cota que el informe de motores ya usó para explicar el error oculto de scikit-fem; la aplicación podría estimar κ durante el propio CG y declararla.

---

## 1. Alcance y método

**Diseño.** Revisión narrativa dirigida de textos técnicos, sin pretensión de exhaustividad bibliográfica. No es una revisión sistemática ni un estudio empírico, por lo que no corresponde una guía de reporte de la red EQUATOR (PRISMA, STROBE u otras). Para que el trabajo sea reproducible se declaran las fuentes, el procedimiento y las limitaciones.

**Procedimiento.**

1. Inventario de lo que implementa spinpy: los 45 módulos de `spinpy/` y `spinpy/motores/`, sus docstrings y los informes `comparativa_*`. Los problemas abiertos se tomaron de la sección «Limitaciones» del README, de `comparativa_febio_tet/ESTADO.md` y de los informes de artefactos, correcciones y VOIs reales.
2. Extracción del texto de los 20 PDF (`pdftotext`) y de sus índices (marcadores del PDF o página de contenidos).
3. Búsqueda dirigida en el texto completo, por términos asociados a cada problema (por ejemplo *arc length*, *superconvergent patch recovery*, *condition number*, *stopping criteria*, *Signorini*, *polyconvex*, *boundary layer*), y lectura de las páginas encontradas.
4. Contraste de cada propuesta con el código. Solo se incluye una propuesta si el código actual no la implementa ya; las que ya están cubiertas se señalan como tales (sección 4).

**Limitaciones.** La búsqueda por términos puede pasar por alto secciones pertinentes que usen otra terminología. No se leyeron los textos completos, sino sus índices y los pasajes localizados. El texto de Reddy (2002) no tenía capa de texto y se leyó por reconocimiento óptico de caracteres, con errores tipográficos que pueden ocultar algún término. Las propuestas R8 a R17 no se han implementado ni medido: su efecto es una hipótesis.

---

## 2. Tabla de prioridades

Esfuerzo: **bajo**, un cambio local con prueba en la suite; **medio**, un módulo nuevo o un cambio de formulación; **alto**, una vía de cálculo nueva.

| # | Propuesta | Dónde | Problema que ataca | Fuente principal | Esfuerzo | Prioridad |
|---|---|---|---|---|---|---|
| R1 | Modos rígidos en el AMG de la homogeneización | `elastic.py:273` | ρ < 0,25 sin alcanzar la tolerancia | Dolean y Tabeart §6.3; Brenner y Scott cap. 6-7 | bajo | alta |
| R2 | Búsqueda lineal sobre la energía e incremento de carga adaptativo | `motores/m_ngsolve.py:283-309`; `fem.py:1045` | no lineal que diverge | Nocedal y Wright cap. 3 y §11.2 | bajo a medio | alta |
| R3 | Neo-Hookeano como material no lineal por omisión | `fem.py:907, 1144, 1536` | no lineal en compresión grande | Kamensky §4.3.4; Dacorogna §3.5 | bajo | alta |
| R4 | Cota de error en la parada de los iterativos (κ estimado en el CG) | `motores/_comun.py:96`; `elastic.py:284`; `resistencia.py:659` | error oculto con residuo pequeño | Quarteroni et al. §4.6.2; Brenner y Scott §9.8 | bajo | alta |
| R5 | Aviso de ejes PCA degenerados | `voi.py:477` | orientación arbitraria del recorte | Modersitzki cap. 5 | bajo | media |
| R6 | Hexaedro con modos incompatibles para la malla de vóxeles | `elastic.hex8_ke` (l. 79); `resistencia.py` | rigidez excesiva con Tb.Th/h bajo | Zienkiewicz et al. §9.8-9.9 | medio | alta |
| R7 | Hexaedros cuadráticos sobre la rejilla de vóxeles, con condensación estática | `motores/m_ngsolve.py`; `fem.mallar` | TET10: memoria a 48³, astillas, tetgen | Šolín et al. §2.2.4, §3.5.9 | medio | alta |
| R8 | Indicador de error por recuperación de tensiones (ZZ/SPR) | `fem.estadisticos` | pico y cola de tensión sin convergencia (A4) | Zienkiewicz et al. cap. 13 | medio | media |
| R9 | Margen del núcleo a partir de la capa límite medida | `fem.py:654` | margen fijo de 0,625 mm | Tartar cap. 18 | medio | media |
| R10 | Plato con contacto unilateral sin fricción | análisis `lineal_plato` | plato pegado al techo | Atkinson y Han §11.5.2; Nocedal y Wright cap. 16-17 | medio | media |
| R11 | Calidad de la malla TET10 y escalado diagonal | `solido.malla_tet10`; motores | astillas y mal condicionamiento | Brenner y Scott §9.6-9.7; Zienkiewicz et al. §8.2.4 | medio | media |
| R12 | Optimización sin derivadas para funciones con ruido en el ajuste | `fit.py`; `fit_dual.py` | optimizadores solo en MATLAB | Nocedal y Wright cap. 9 | medio | media |
| R13 | Seguimiento de trayectoria (longitud de arco) | motores no lineales | carga de fallo con fuerza impuesta | Nocedal y Wright §11.3; Bathe §8.4.3 | alto | baja |
| R14 | Estimación de error orientada al E aparente | `resistencia.estudio_convergencia` | convergencia verificada solo en el cubo | Zienkiewicz et al. §13.9; Šolín et al. §6.3 | alto | baja |
| R15 | Detección de mecanismos por autovalores bajos de K | `fem.mallar` | bisagras (A5) y fragmentos casi sueltos | Quarteroni et al. §5.11; Sun y Zhou cap. 9 | medio | baja |
| R16 | Formulación solo sólido en la homogeneización | `elastic.homogeneizar` | contraste 10⁻⁶ entre hueso y vacío | Tartar cap. 16; Brenner y Scott §5.2 | medio | baja |
| R17 | Métodos inmersos (celda finita) con el conjunto de nivel | vía nueva | escalera de vóxeles y tetgen | Kamensky §8.2; Cottrell et al. | alto | baja |

---

## 3. Fichas

### R1. Modos rígidos en el multigrid de la homogeneización

**Qué hace hoy spinpy.** `elastic.homogeneizar` resuelve los seis problemas de celda con `pyamg.smoothed_aggregation_solver(Kff_csr, max_coarse=500)` (`elastic.py:273`). No pasa el argumento `B`, de modo que pyamg toma como espacio casi nulo un vector constante, el que corresponde a un problema escalar. El ensayo de compresión (`resistencia.py:659`) y el adaptador común de los motores (`motores/_comun.py:113`) sí pasan los modos rígidos. La aplicación declara que por debajo de ρ ≈ 0,25 la homogeneización no alcanza la tolerancia.

**Qué dicen los libros.** El multigrid funciona porque el suavizador elimina las componentes oscilantes del error y la corrección en la malla gruesa elimina las suaves, las de baja energía (Dolean y Tabeart 2026, §6.3; Brenner y Scott 2008, cap. 6). En elasticidad, los modos de energía más baja de una pieza poco conectada son movimientos casi rígidos de cada trabécula. Si el espacio grueso no puede representarlos, el ciclo no los corrige y el CG se estanca, que es lo que se espera precisamente cerca del umbral de conectividad.

**Propuesta.** Construir `B` con los tres modos de traslación y los tres de rotación (la función `motores/_comun.modos_rigidos` ya existe) y pasarlo al constructor. En el problema periódico las rotaciones globales no son compatibles con la periodicidad, pero el espacio casi nulo del AMG es local por agregados, y ahí las rotaciones sí son modos de baja energía.

**Cómo verificarlo.** Repetir la curva de ρ de `Test/replicar_zheng2021.py` con y sin `B`, registrando iteraciones, residuo final y diferencia con el LU directo en los tamaños donde el directo cabe. El criterio de aceptación sería alcanzar la tolerancia declarada en más puntos sin cambiar el tensor más allá de 10⁻⁸ relativo frente al directo.

**Riesgo.** Si la falta de convergencia se debe a la degeneración del propio problema (R16 y Tartar, cap. 16), el cambio mejorará el número de iteraciones sin eliminar el límite. Esa diferencia es en sí misma informativa.

### R2. Globalización del Newton no lineal

**Qué hace hoy spinpy.** El motor NGSolve resuelve el problema no lineal con Newton completo (`motores/m_ngsolve.py:290-307`): calcula el incremento con la tangente y lo aplica entero, sin búsqueda lineal; si no converge en 30 iteraciones, lanza `ErrorMotor`. Las cargas se reparten en incrementos iguales (`fem.py:1045`) y por omisión hay un solo incremento (`pasos=1`). Con FEBio, el ensayo no lineal con fuerza sobre la malla TET10 de un VOI real divergió tras 14 recortes de paso, con incrementos de desplazamiento de 223 mm en un VOI de 1,6 mm (`comparativa_febio_tet/INFORME.md`).

**Qué dicen los libros.** Nocedal y Wright (2006) tratan el problema como minimización de la energía potencial, que es la forma natural en hiperelasticidad. Un paso de Newton solo garantiza descenso si la tangente es definida positiva y la longitud del paso se controla (cap. 3, pp. 30-65): búsqueda lineal con condiciones de Armijo o Wolfe (§3.1) y, si la tangente deja de ser definida positiva, modificación de la matriz (§3.4, p. 48). Para sistemas no lineales generales proponen una función de mérito y búsqueda lineal o región de confianza sobre ella (§11.2, pp. 285-296). La pérdida de definición de la tangente tiene además lectura física. Una segunda variación definida positiva es la condición suficiente de mínimo (Gelfand y Fomin 2000, cap. 5) y, en elasticidad, el mínimo de la energía potencial corresponde al equilibrio estable. Un pivote negativo en la factorización de la tangente apunta, por tanto, a una inestabilidad de la estructura (por ejemplo, el pandeo de una trabécula) más que a un fallo numérico.

**Propuesta.** En este orden: (a) búsqueda lineal por retroceso sobre la energía, que NGSolve evalúa con `BilinearForm.Energy`; (b) incremento de carga adaptativo, que reduzca el paso a la mitad cuando Newton no converge y lo aumente cuando converge en pocas iteraciones; (c) registro de la inercia de la tangente (pivotes negativos) en cada incremento, para declarar en el informe si la estructura perdió estabilidad.

**Cómo verificarlo.** Los dos casos con solución cerrada de la comparativa (compresión de un bloque con SVK y con el neo-Hookeano de FEBio) deben dar las mismas fuerzas, porque la búsqueda lineal no modifica la solución. El caso que hoy falla (TET10 de VOI real con fuerza) es la prueba de que el cambio sirve.

### R3. Neo-Hookeano como material no lineal por omisión

**Qué hace hoy spinpy.** El material por omisión es `svk` en la API, la interfaz y la línea de comandos (`fem.py:907, 1144, 1536`).

**Qué dicen los libros.** Kamensky (2022, §4.3.4, p. 57) señala que el modelo de St. Venant-Kirchhoff «no es estable en compresión»: su energía no diverge cuando J → 0 y solo es adecuado para deformaciones pequeñas con rotaciones grandes. Propone como alternativa el neo-Hookeano compresible. Dacorogna (2004, §3.5, p. 100) da el marco teórico: la existencia de minimizadores en elasticidad no lineal se apoya en la policonvexidad (Ball), propiedad que tiene el neo-Hookeano y no el SVK.

**Propuesta.** Cambiar el valor por omisión a `neohookeano` y conservar `svk` como opción, porque la comparativa con FEBio se hizo con ambos. Con las deformaciones del ensayo lineal los dos coinciden; la diferencia aparece cerca del fallo, que es donde la aplicación usa el no lineal.

**Riesgo.** Cambia los números no lineales publicados con la V2.1.1. Hay que declararlo en las notas de versión y en el informe de publicación.

### R4. Cota de error en la parada de los iterativos

**Qué hace hoy spinpy.** Todos los iterativos paran por residuo relativo (10⁻⁸ en la homogeneización y en el ensayo de la app; 10⁻¹⁰ en `motores/_comun.resolver_scipy`) y comprueban después ese residuo. El informe de motores mostró el límite de ese criterio: el CG + pyamg de scikit-fem se detuvo en un residuo de 10⁻¹⁰ con un error de desplazamiento de 0,9 % en la malla TET10.

**Qué dicen los libros.** Con el residuo normalizado, el error relativo queda acotado por ε·K(A), no por ε (Quarteroni et al. 2007, §4.6.2, p. 174). En la norma de energía, el CG reduce el error en un factor que depende de √κ (Brenner y Scott 2008, §9.8, p. 266). Los coeficientes del CG definen la matriz tridiagonal de Lanczos, cuyos autovalores extremos aproximan los de la matriz precondicionada (Quarteroni et al., §4.4 y §5.11; Dolean y Tabeart, cap. 4).

**Propuesta.** Guardar los coeficientes α y β del CG, estimar κ con los autovalores extremos de la tridiagonal de Lanczos y registrar en cada resolución la cota ε·κ junto al residuo. Cuando la cota supere un umbral, avisar en el informe, como ya se hace con el residuo. pyamg admite una retrollamada por iteración, así que el cambio queda dentro de `resolver_scipy` y de sus equivalentes.

**Cómo verificarlo.** En el caso de scikit-fem a 32³ la cota debe ser del orden del error medido (8,9·10⁻³). En los casos en que el iterativo coincide con el directo, la cota debe quedar por debajo de la tolerancia.

### R5. Aviso de ejes PCA degenerados

**Qué hace hoy spinpy.** `voi.marco_pca` orienta el recorte con la descomposición en valores singulares de la nube de vóxeles óseos (`voi.py:477`) y fija el signo de cada eje. No comprueba si dos valores singulares son casi iguales.

**Qué dicen los libros.** La descomposición de la covarianza es «esencialmente única, salvo el signo de las columnas», solo si sus autovalores son simples (Modersitzki 2004, cap. 5, p. 46). Si dos coinciden, cualquier par de ejes de ese plano es igual de válido y la orientación del recorte queda arbitraria. El mismo capítulo (§5.3, p. 48) describe una variante robusta basada en la distribución t (Kent y Tyler), menos sensible a masas atípicas, como podría ser la cortical de la pieza.

**Propuesta.** Calcular el cociente entre valores singulares consecutivos y avisar cuando sea cercano a 1, con un umbral que se fije midiendo la estabilidad de los ejes al perturbar el umbral de segmentación. La variante robusta es opcional y de menor prioridad.

### R6. Hexaedro con modos incompatibles para la malla de vóxeles

**Qué hace hoy spinpy.** El ensayo de la app y la homogeneización usan el hexaedro trilineal de 8 nodos con integración de Gauss 2×2×2 (`elastic.hex8_ke`, l. 79). En los VOIs reales a 40³ hay entre 1,3 y 2,8 vóxeles por espesor trabecular (`INFORME_VOIS_REALES.pdf`), y el estudio de convergencia del porcino mostró que el E aparente de los candidatos aún cambia entre 34³ y 40³.

**Qué dicen los libros.** El elemento bilineal o trilineal introduce deformaciones de corte espurias en flexión pura, porque su campo de desplazamientos no puede curvarse dentro del elemento. El resultado es una rigidez excesiva que solo desaparece al refinar (Zienkiewicz et al. 2005, §9.8, p. 343, y Fig. 9.10). Los modos incompatibles de Wilson y Taylor añaden por elemento las funciones 1 − ξ² y 1 − η² (en tres dimensiones, también 1 − ζ²), que se condensan estáticamente. Con ellos el elemento reproduce la flexión pura en un elemento rectangular y pasa la prueba de la parcela si se corrigen las derivadas (§9.8, ec. 9.17). El libro lo desarrolla en deformación plana; la extensión al hexaedro es la estándar, pero su comportamiento en la malla de vóxeles hay que medirlo. La prueba de la parcela de orden superior (§9.9, p. 347) muestra la mejora en flexión.

**Por qué encaja con spinpy.** En la malla de vóxeles todos los elementos son cubos idénticos con jacobiano constante, así que la matriz elemental condensada se calcula una sola vez y el número de incógnitas globales no cambia. Una trabécula de dos vóxeles de espesor en flexión es exactamente el caso en que el trilineal falla.

**Cómo verificarlo.** (a) Prueba de la parcela con deformación constante; (b) viga en flexión con solución cerrada, a 2, 3 y 4 elementos por espesor; (c) estudio de convergencia del E aparente en el porcino y en los VOIs reales, frente a la referencia embebida de `comparativa_motores/vois_reales/`. La hipótesis que se pone a prueba es que, a igual resolución, el error con Tb.Th/h bajo disminuye. Si no disminuye, el cambio no se integra.

**Riesgo.** Cambia el problema discreto validado frente a FEBio. Debe entrar como elemento opcional, con su propia validación, y no sustituir al trilineal hasta tener esas pruebas.

### R7. Hexaedros cuadráticos sobre la rejilla de vóxeles

**Qué hace hoy spinpy.** Para tener elementos cuadráticos, la aplicación genera una superficie suavizada y la tetraedraliza con tetgen (TET10). Esa vía tiene tres problemas medidos: la malla a 48³ no cabe en 15 GB con VOIs reales, aparecen astillas y tetraedros diminutos que empeoran el condicionamiento, y en un VOI equino tetgen no pudo tetraedralizar la superficie.

**Qué dicen los libros.** Šolín et al. (2004) construyen elementos jerárquicos de orden arbitrario sobre el hexaedro de referencia (§2.2.4, p. 62) y muestran que las funciones internas de cada elemento se pueden condensar antes del montaje global (§3.5.9, p. 191). La comparativa de motores ya midió que, a igual número de grados de libertad, el TET10 dio un error L2 19 veces menor que el TET4 en el cubo con solución exacta.

**Propuesta.** Resolver la malla de vóxeles con hexaedros de orden 2 en NGSolve, que admite orden alto sobre hexaedros y condensación estática. Esta vía no necesita superficie, tetgen ni corrección de volumen, y no produce astillas. Conserva la escalera de vóxeles (artefacto A4), que la malla suave sí atenúa.

**Cómo verificarlo.** Medir con `fem.tamano_previsto` la memoria a 40³ y 48³ frente a TET10, y comparar el E aparente frente a la referencia embebida. Los datos actuales no permiten anticipar cuál de las dos vías dará menor error; es una pregunta para medir.

### R8. Indicador de error por recuperación de tensiones

**Qué hace hoy spinpy.** La tensión de cada hexaedro se evalúa en su centroide (`resistencia.py:129`), que es el punto óptimo de muestreo del gradiente para elementos lineales (Zienkiewicz et al. 2005, §13.2, p. 459). Esa parte ya es correcta. El informe cita el p99 de la capa superficial porque el pico no converge; el informe de artefactos propuso añadir un indicador de error por recuperación.

**Qué dicen los libros.** El capítulo 13 de Zienkiewicz et al. describe la recuperación por parcelas superconvergente (SPR, §13.4, p. 467): un ajuste local por mínimos cuadrados de los valores óptimos alrededor de cada nodo. La diferencia entre la tensión recuperada y la del elemento da el estimador de Zienkiewicz y Zhu (§13.6, p. 476), con índice de efectividad que tiende a 1 cuando la recuperación es superconvergente (ec. 13.48). Los autores advierten que el error de los valores recuperados sigue sin resolverse (§13.10, p. 494).

**Propuesta.** Calcular el estimador ZZ por elemento y usarlo para dos cosas: un mapa de fiabilidad junto al de von Mises y una fracción de la capa superficial con error estimado alto, declarada junto al p99. No se propone sustituir el p99 por un valor recuperado. En la superficie las parcelas son incompletas y la recuperación es menos precisa, justo donde se mide la tensión.

### R9. Margen del núcleo a partir de la capa límite medida

**Qué hace hoy spinpy.** El E corregido se mide en el núcleo del VOI, a 0,625 mm fijos de sus caras (`fem.py:654`). El informe de VOIs reales midió perfiles de tensión por distancia a la cara lateral.

**Qué dicen los libros.** En homogeneización periódica, la solución cerca de una frontera plana converge exponencialmente al comportamiento interior (Tartar 2009, cap. 18, pp. 195-196, lema de Lions con la demostración variacional de Tartar). El autor añade que la mayor parte del error de sustituir el material oscilante por el homogeneizado parece venir de lo que ocurre cerca de la frontera.

**Propuesta.** Ajustar una exponencial al perfil de tensión media por distancia a la cara y fijar el margen en un número dado de longitudes de decaimiento, con el valor actual como mínimo. El margen pasaría a depender de la arquitectura del VOI y no de un número fijo validado en un espinodoide.

**Riesgo.** El resultado de Tartar es para una ecuación escalar con coeficientes periódicos. Su extrapolación a elasticidad en una arquitectura trabecular aleatoria es razonable, pero no está demostrada. Por eso la propuesta mide el decaimiento en cada VOI en lugar de suponerlo.

### R10. Plato con contacto unilateral sin fricción

**Qué hace hoy spinpy.** El análisis `lineal_plato` impone el mismo desplazamiento vertical a todos los nodos del techo. Un nodo que en un ensayo físico tendería a separarse del plato queda así sujeto a él y transmite tracción.

**Qué dicen los libros.** Atkinson y Han (2009, §11.5.2, p. 465) formulan el problema de Signorini sin fricción como inecuación variacional: el cuerpo no puede penetrar el soporte rígido y solo puede recibir compresión de él. Nocedal y Wright (2006) dan los algoritmos para el problema discreto: programación cuadrática con conjunto activo (§16.5, p. 467) o lagrangiano aumentado (§17.3, p. 514).

**Propuesta.** Primero diagnosticar: calcular la reacción de cada nodo del techo en el análisis actual y medir la fracción de fuerza que es de tracción. Solo si esa fracción es apreciable en los VOIs reales merece la pena implementar el contacto unilateral con conjunto activo, que en elasticidad lineal converge en pocas iteraciones.

### R11. Calidad de la malla TET10 y escalado diagonal

**Qué hace hoy spinpy.** La malla suave se mejora con la razón radio-arista de tetgen (`minratio`, `solido.py:246`). La comparativa de motores encontró 577 tetraedros diminutos a 32³ y astillas a 48³, y relacionó ambos con el error oculto del iterativo.

**Qué dicen los libros.** Brenner y Scott (2008, §9.6-9.7, pp. 261-266) demuestran que, en tres dimensiones, el número de condición se mantiene en O(N^(2/3)) aunque la malla tenga elementos de tamaños muy distintos, siempre que se cumplan dos condiciones: que la malla sea no degenerada (regularidad de forma) y que la base se escale de forma natural, lo que equivale a un escalado diagonal. De ahí se deduce una distinción útil para spinpy: un tetraedro diminuto pero bien formado no daña el condicionamiento si se escala la diagonal, mientras que una astilla viola la hipótesis de regularidad y sí lo daña. Zienkiewicz et al. (§8.2.4, p. 280) describen el suavizado laplaciano con restricciones como paso final de mejora de calidad.

**Propuesta.** (a) Aplicar escalado diagonal simétrico antes de los iterativos que no lo hagan ya, y medir su efecto en el caso de scikit-fem; (b) informar la calidad de forma (diedro mínimo, razón radio-arista) junto al tamaño, porque son problemas distintos; (c) probar un suavizado restringido de los nodos interiores después de tetgen. Si se adopta R7, esta línea pierde prioridad.

### R12. Optimización sin derivadas para el ajuste con ruido

**Qué hace hoy spinpy.** El ajuste es una búsqueda escalonada sobre rejillas, con réplicas por semilla y un suelo de ruido autoconsistente (`fit.py`, `incertidumbre.py`). Los optimizadores de Pareto, bayesiano y MOBO siguen solo en MATLAB.

**Qué dicen los libros.** Nocedal y Wright (2006, cap. 9) tratan justamente funciones objetivo que dependen de una simulación estocástica y tienen un error aleatorio en cada evaluación (§9.1, p. 221). Las diferencias finitas pierden toda precisión cuando el ruido domina el intervalo de diferencia (§9.1). Proponen métodos basados en modelos de interpolación con región de confianza (§9.2, p. 223), Nelder-Mead (§9.5, p. 238) y el filtrado implícito (§9.6, p. 240), diseñado para funciones suaves con ruido superpuesto.

**Propuesta.** Usar el nivel de ruido que spinpy ya mide (el suelo autoconsistente) como el η de Nocedal y Wright, y probar el filtrado implícito o un método de modelo como sustituto en Python de los optimizadores que faltan. La búsqueda escalonada actual serviría de referencia.

**Cómo verificarlo.** Comparar error final, número de evaluaciones y dispersión entre semillas con la búsqueda escalonada, en el banco de VOIs del estudio de familias.

### R13 a R17. Líneas de largo plazo

- **R13. Seguimiento de trayectoria.** Nocedal y Wright (§11.3, pp. 296-300) describen la continuación con longitud de arco como parámetro: carga y desplazamiento avanzan juntos, y la trayectoria puede seguirse más allá de un punto límite. Sería la forma correcta de obtener la carga máxima con fuerza impuesta. El tratamiento estructural clásico es el de Bathe (1996, §8.4.3, pp. 761-764): multiplicador de carga, restricción de longitud de arco esférica o de trabajo externo constante, y un algoritmo que «debe detener la iteración cuando la divergencia es inminente y reiniciarse con nuevos parámetros». Conviene hacerlo después de R2 y R3.
- **R14. Error orientado al E aparente.** El E aparente es una funcional de la solución (la reacción en el techo). Zienkiewicz et al. (§13.9, p. 490) y Šolín et al. (§6.3, p. 324) muestran cómo estimar su error con un problema adjunto resuelto en la misma malla. Eso permitiría declarar el error del E aparente en el espécimen sin repetir el ensayo en cuatro mallas, que es lo que hoy pide la lista de chequeo de Erdemir et al. (2012).
- **R15. Mecanismos por autovalores.** Los autovalores más bajos de la matriz de rigidez, calculados con Lanczos (Quarteroni et al. §5.11; Sun y Zhou 2017, cap. 9), localizan las partes que casi pueden moverse como sólido rígido: bisagras por arista o vértice (artefacto A5) y trabéculas casi sueltas. Hoy el filtro de conectividad por caras (`resistencia._solo_portante`) solo elimina las piezas que no unen base y techo.
- **R16. Homogeneización solo sólido.** Tartar (2009, cap. 16, pp. 177-178) formula la homogeneización con agujeros y condición de Neumann, y su análisis exige una constante de extensión acotada y un dominio conexo. Cerca del umbral de conectividad esa constante no puede mantenerse acotada. Eso sugiere que la falta de convergencia por debajo de ρ ≈ 0,25 es en parte intrínseca, como ya declara la aplicación. Eliminar el vacío (en lugar de darle rigidez 10⁻⁶) quitaría el contraste artificial, pero exige tratar las islas como hace `_solo_portante` y fijar las traslaciones (Brenner y Scott §5.2, problema de Neumann puro).
- **R17. Métodos inmersos.** Kamensky (2022, §8.2, p. 101) resume la idea del método de la celda finita: para conservar el orden de convergencia basta con que la cuadratura se ajuste al dominio, aunque el espacio de funciones no lo haga. Solo hace falta saber si un punto está dentro o fuera. El espinodoide se define por un conjunto de nivel analítico (`grf.py`), así que es un candidato natural. La vía eliminaría la escalera de vóxeles y la dependencia de tetgen, a cambio de cuadratura adaptativa y de una menor robustez, que el propio autor señala. Cottrell, Hughes y Bazilevs (2009) aportan las bases spline para esa línea.

---

## 4. Lo que ya está bien resuelto según los libros

La revisión confirma varias decisiones de diseño que no conviene tocar:

- **Muestreo de la tensión en el centroide del hexaedro**, el punto óptimo para elementos lineales (Zienkiewicz et al. §13.2).
- **BDDC en NGSolve para la malla suave.** Para un problema modelo, Brenner y Scott (§7.8, teorema 7.8.19, p. 210) acotan el número de condición del sistema precondicionado por C[1 + ln(H/h)]², con C independiente del tamaño de malla, del tamaño de subdominio y del número de subdominios. La comparativa confirmó en la práctica que fue fiable donde el AMG genérico no lo fue.
- **Verificación por soluciones manufacturadas y orden de convergencia observado**, que es la práctica que recomiendan Langtangen (2016, cap. 2 y §5.3) y Langtangen y Linge (2017). El estudio de convergencia de la V2.0.2 reproduce las tasas a priori que demuestran Brenner y Scott (cap. 4).
- **Extrapolación de Richardson en el estudio de convergencia del espécimen**, ya implementada (`resistencia.estudio_convergencia`).
- **Cotas de Voigt y Hashin-Shtrikman en la validación** (`Test/replicar_zheng2021.py`), cuyo fundamento está en Tartar (cap. 21 y 25).
- **Laminado de Backus como prueba del tensor**, que coincide con la homogenización de laminados de Tartar (cap. 12 y 27).

---

## 5. Relevancia de cada texto

| Texto | Relevancia | Para qué sirve en spinpy |
|---|---|---|
| Zienkiewicz, Taylor y Zhu 2005 | alta | R6 (modos incompatibles), R8 (SPR y ZZ), R11 (calidad de malla), R14 |
| Nocedal y Wright 2006 | alta | R2 (globalización), R10 (QP y lagrangiano aumentado), R12 (sin derivadas con ruido), R13 |
| Brenner y Scott 2008 | alta | R1 (multigrid), R4 y R11 (condicionamiento), R16; confirma BDDC |
| Šolín, Segeth y Doležel 2004 | alta | R7 (hexaedros jerárquicos, condensación), R14 (adaptatividad orientada) |
| Quarteroni, Sacco y Saleri 2007 | media | R4 (criterios de parada, Lanczos), R15 |
| Tartar 2009 | media | R9 (capa límite), R16 (agujeros con Neumann); respaldo de las cotas ya usadas |
| Kamensky 2022 | media | R3 (SVK frente a neo-Hookeano), R17 (métodos inmersos) |
| Atkinson y Han 2009 | media | R10 (Signorini); teoría de Newton y de formulaciones débiles |
| Modersitzki 2004 | media | R5 (ejes principales y su variante robusta) |
| Dolean y Tabeart 2026 | media | R1 (corrección gruesa), R4 (Lanczos y CG) |
| Dacorogna 2004 | baja | R3 (policonvexidad) |
| Gelfand y Fomin 2000 | baja | R2 (segunda variación y estabilidad) |
| Cottrell, Hughes y Bazilevs 2009 | baja | R17 (bases spline) |
| Sun y Zhou 2017 | baja | R15 (autovalores de matrices grandes) |
| Gatica 2014 | baja | métodos mixtos para obtener tensiones equilibradas; no hay un problema actual que lo pida, porque ν = 0,3 no produce bloqueo volumétrico |
| Allaire 2007 | baja | introducción general; no aporta nada que no cubran los anteriores |
| Langtangen 2016; Langtangen y Linge 2017 | baja | diferencias finitas; sus prácticas de verificación ya están en la suite |
| Langtangen y Logg 2016 (tutorial de FEniCS) | baja | escrito para la versión antigua de FEniCS; FEniCSx es un motor opcional |
| Introducción al Método de Elementos Finitos (apuntes, cap. 1 de A. Brewer) | baja | texto introductorio; útil para formación, no para cambios en el código |
| Bathe 1996 | alta | R2 (§8.4.1-8.4.2: control de paso y búsqueda lineal), R13 (§8.4.3), criterio de energía (§8.4.4), R6 (§4.4.1, modos incompatibles en cubos), R3 (§6.6.1) |
| Reddy 2002 | baja | bloqueo por cortante en vigas y placas (§9.4), apoyo conceptual a R6 |
| Bathe y Wilson 1976 | baja | antecedente de Bathe 1996; no añade nada que este no cubra |

---

## 6. Textos en DjVu (segunda fase)

Se instalaron `djvulibre` y `tesseract` en el entorno. Bathe (1996) y Bathe y Wilson (1976) tienen capa de texto; Reddy (2002) se leyó por reconocimiento óptico.

**Bathe (1996)** es el texto con más aplicación directa de los tres:

- **§8.4.1 y §8.4.2 (pp. 755-761).** Con la tangente exacta, «el procedimiento principal para alcanzar la convergencia es reducir el incremento de carga» (p. 758), y la búsqueda lineal «puede evitar la divergencia», que es la razón principal de su eficacia (p. 761). Es el fundamento de R2.
- **§8.4.3 (pp. 761-764).** Métodos de restricción carga-desplazamiento (longitud de arco) para atravesar puntos límite: fundamento de R13.
- **§8.4.4 (pp. 764-765).** Criterios de convergencia por desplazamiento, por fuerza y por energía. El de energía combina los dos primeros; el motor NGSolve usa hoy solo el de fuerza.
- **§4.4.1 y ejemplo 4.28 (pp. 262-268).** Modos incompatibles: en un elemento cuadrado la integral de la matriz B de los modos es nula y la prueba de la parcela se cumple sin corrección. En la malla de vóxeles todos los elementos son cubos, así que R6 no necesita la corrección de los elementos distorsionados. Al ser no conforme, la energía deja de ser una cota y la convergencia puede no ser monótona.
- **§6.6.1 (pp. 583-589).** El material de St. Venant-Kirchhoff es natural con grandes desplazamientos y rotaciones pero deformaciones pequeñas; con deformaciones grandes la respuesta cambia por completo. Apoya R3.
- **§10.2 (p. 849).** Propiedad de la secuencia de Sturm: el número de elementos negativos de D en la factorización LDLᵀ es el número de autovalores por debajo del desplazamiento. Permite detectar la pérdida de estabilidad en el Newton (R2) y apoya R15.

**Reddy (2002)** trata principios energéticos y métodos variacionales. Su §9.4 explica el bloqueo por cortante de los elementos de viga de Timoshenko y lo corrige con integración reducida, la misma raíz que la rigidez excesiva del hexaedro trilineal en flexión (R6).

**Bathe y Wilson (1976)** es el antecedente de Bathe (1996) y no aporta nada que este no cubra.

**Lo que ningún texto resuelve.** La segmentación por umbral global sigue siendo la mayor fuente de incertidumbre, y ninguno de los libros trata la segmentación de micro-TC (Modersitzki trata registro, no segmentación). La carga de fallo de Pistoia con plato sigue sin referencia para validarla, y la validación del modelo frente a un ensayo físico no puede sustituirse con métodos numéricos.

---

## Referencias

1. Zienkiewicz OC, Taylor RL, Zhu JZ. The finite element method: its basis and fundamentals. 6th ed. Oxford: Elsevier Butterworth-Heinemann; 2005.
2. Nocedal J, Wright SJ. Numerical optimization. 2nd ed. New York: Springer; 2006.
3. Brenner SC, Scott LR. The mathematical theory of finite element methods. 3rd ed. New York: Springer; 2008.
4. Šolín P, Segeth K, Doležel I. Higher-order finite element methods. Boca Raton: Chapman & Hall/CRC; 2004.
5. Quarteroni A, Sacco R, Saleri F. Numerical mathematics. 2nd ed. Berlin: Springer; 2007. (Texts in Applied Mathematics 37).
6. Tartar L. The general theory of homogenization: a personalized introduction. Berlin: Springer; 2009. (Lecture Notes of the Unione Matematica Italiana 7).
7. Kamensky D. Finite element analysis for coupled problems. Lecture notes for MAE 207. University of California San Diego; 2022.
8. Atkinson K, Han W. Theoretical numerical analysis: a functional analysis framework. 3rd ed. New York: Springer; 2009.
9. Modersitzki J. Numerical methods for image registration. Oxford: Oxford University Press; 2004.
10. Dolean V, Tabeart J. Advanced linear algebra with applications, Part I. Lecture notes. arXiv:2608.21234; 2026.
11. Dacorogna B. Introduction to the calculus of variations. London: Imperial College Press; 2004.
12. Gelfand IM, Fomin SV. Calculus of variations. Mineola: Dover; 2000.
13. Cottrell JA, Hughes TJR, Bazilevs Y. Isogeometric analysis: toward integration of CAD and FEA. Chichester: Wiley; 2009.
14. Sun J, Zhou A. Finite element methods for eigenvalue problems. Boca Raton: CRC Press; 2017.
15. Gatica GN. A simple introduction to the mixed finite element method: theory and applications. Cham: Springer; 2014.
16. Allaire G. Numerical analysis and optimization: an introduction to mathematical modelling and numerical simulation. Oxford: Oxford University Press; 2007.
17. Langtangen HP. Finite difference computing with exponential decay models. Cham: Springer; 2016.
18. Langtangen HP, Linge S. Finite difference computing with PDEs: a modern software approach. Cham: Springer; 2017.
19. Langtangen HP, Logg A. Solving PDEs in Python: the FEniCS tutorial I. Cham: Springer; 2016.
20. Erdemir A, Guess TM, Halloran J, Tadepalli SC, Morrison TM. Considerations for reporting finite element analysis studies in biomechanics. J Biomech. 2012;45(4):625-33. doi:10.1016/j.jbiomech.2011.11.038
21. Bathe KJ. Finite element procedures. Englewood Cliffs: Prentice Hall; 1996.
22. Reddy JN. Energy principles and variational methods in applied mechanics. 2nd ed. New York: Wiley; 2002.
23. Bathe KJ, Wilson EL. Numerical methods in finite element analysis. Englewood Cliffs: Prentice-Hall; 1976.
