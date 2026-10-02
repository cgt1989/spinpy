# Qué hizo y qué asumió cada autor

Estado: **V** = leído o confirmado en una base científica en esta sesión (resumen o fragmento de texto); **R** = tomado de las notas del repositorio spinpy, que leyó el artículo completo, pendiente de reconfirmar. Donde solo se leyó el resumen se indica. Antes de citar en un manuscrito, reconfirma la afirmación concreta en el artículo.

## Kumar, Tan, Zheng, Kochmann (2020), npj Comput Mater 6:73
- Define el spinodoide por GRF de ondas planas, conos en tres ejes y conjunto de nivel con phi0 = sqrt(2) erfinv(2 rho − 1). (V, resumen; ecuaciones R)
- Ajusta mediante aprendizaje automático la rigidez anisótropa y la densidad; el resumen afirma que reproduce propiedades del hueso trabecular. (V, resumen) Objetivo: rigidez, no morfometría.
- theta_i en {0} ∪ [theta_min, 90°] con theta_min = 15°; la unión de conos se muestrea de forma uniforme (rechazo). (R)
- beta = 15π en la leyenda de la Fig. 5; beta de la Fig. 2 no publicado. (R)

## Zheng, Kumar, Kochmann (2021), CMAME 383:113894
- Optimización topológica multiescala con sustituto neuronal; homogeneización con condiciones afines (cota superior). (V, resumen del preprint; BC R)
- Fig. 1 publica las cuatro ternas; rho >= 0.3 para evitar dominios disjuntos; theta_min = pi/6 en texto y 15° en pies de figura (inconsistencia interna). (R)

## Guo, Sharma, Kumar (2024), Adv Intell Syst 6(6):2300789
- Perfil de curvaturas (k1, k2) ponderado por área como descriptor; diseño inverso por aprendizaje profundo; PNS de referencia y espinodal con beta = 15π, Q = 1000, rho = 0.3, theta = (60, 30, 10). (R)

## Deng, Kumar, Vallone, Kochmann, Greer (2024), Adv Mater 36(34)
- Mismo marco GRF y conos ("three mutually orthogonal pairs of cones", theta en [0, pi/2]); redacción "uniformly and randomly sampled from within these cones", ambigua entre muestreos. (V, texto)
- Modelo inverso sobre el tensor de rigidez completo, aplicado a 127 muestras trabeculares de cabeza femoral de 32 pacientes; manufactura aditiva. (V, texto)

## Vafaeefar, Moerman, Kavousi, Vaughan (2022 en línea; J Mech Behav Biomed Mater 138:105584)
- Compara gyroide, spinodoide y dual-lattice calibrados a BV/TV; gyroide y spinodoide difieren en otros parámetros morfométricos y topológicos y muestran propiedades mecánicas efectivas menores que el hueso; el dual-lattice se acerca más. Software en GIBBON. (V, PubMed, resumen)
- Consistente con la brecha de rigidez observada en spinpy (spinodoide ajustado a morfometría más blando que el hueso). No es contradicción de Kumar: el objetivo de ajuste es distinto.

## Golnary et al. (2024), Int J Mech Mater Des
- Anisotropía por un índice universal; conos grandes dan baja anisotropía; densidades bajas tienden a mayor anisotropía; diseño inverso con múltiples candidatos y distancia de Mahalanobis a la distribución de tensores. (V, resumen)

## Raßloff et al. (2024), Comput Mech
- Diseño inverso por optimización bayesiana en régimen de pocos datos (rigidez). (V, resumen)

## Otto, Rosenkranz, Kalina, Kästner (2025), GAMM-Mitteilungen 48(4)
- Homogeneización con solver FFT, sustituto neuronal, optimización en el espacio de descriptores, pérdida logarítmica, método para determinar la clase de anisotropía de un tensor; validado con fémur. (V, resumen y texto)

## Röding et al. (2022), Sci Rep
- Spinodoides con anisotropía ajustable por GRF; red convolucional para difusividad direccional y diseño inverso bayesiano aproximado. (V, resumen)

## Park et al. (2026), Mater Des
- Los descriptores de cono idénticos pueden dar morfologías distintas y dispersión de propiedades; regresión de proceso gaussiano heterocedástica; el óptimo determinista es sensible a la incertidumbre. (V, resumen) Respalda tratar la realización como fuente de varianza.

## Mandolesi et al. (2025), Eur J Mech A/Solids (dos artículos)
- Resuelven la ecuación de Cahn-Hilliard adimensional y señalan que los GRF son válidos solo en las etapas iniciales de la descomposición espinodal; anisotropía por movilidad direccional. (V, resumen)

## Risthaus y Schneider (2024), PAMM 24(4)
- Condiciones de contorno de Dirichlet en homogeneización FFT sobre microestructuras bicontinuas por conjunto de nivel de GRF; el umbral parte de erfinv y se refina con Newton para lograr la fracción de volumen en la rejilla; comparan con valores de Soyarslan et al. (V, texto)
- Implica distinguir densidad nominal (erfinv) y realizada (en la máscara a esa resolución).

## Yıldız et al. (2026), J Phys Mater
- Optimización multifísica de spinodoides bifásicos con CNN; la variación inherente del GRF limita optimizadores basados en gradiente. (V, resumen)

## Citados en la bibliografía de otros y no verificados directamente
Thakolkaran y col. (spinodoides con curvas tensión-deformación; mencionado en Otto 2025), Rosenkranz y col. (sustituto equivariante), Soyarslan y col. (microestructuras espinodales isótropas), Vidyasagar y col. (origen anisótropo por energía de superficie o movilidad, citado en Deng 2024). Buscarlos y verificarlos antes de usarlos.

## Puntos donde spinpy difiere o podría diferir, y cómo declararlos

| Punto | Autores | spinpy | Tratamiento |
|---|---|---|---|
| Objetivo de ajuste | rigidez (Kumar, Deng, Otto, Raßloff) | morfometría (como Vafaeefar) | desviación declarada |
| Condiciones de contorno | afines (Zheng) | periódicas | desigualdades, no cifras |
| Muestreo | uniforme en la unión (ec. 2) | `rechazo` por defecto, `equitativo` opcional | declarar siempre |
| theta_min | 15° o 30° según la parte de Zheng | barrido de presets | declarar y no mezclar |
| Densidad | nominal | impuesta con erfinv exacta; la realizada debe medirse | reportar ambas |
| Superficie | no especificada | marching cubes sobre binario, sesgo +8.5 a +8.8 % del área | declarar y usar campo continuo si se compara con curvaturas |
| Régimen de validez | GRF, etapas iniciales | igual | no extrapolar a morfologías tardías |
