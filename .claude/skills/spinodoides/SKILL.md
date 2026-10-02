---
name: spinodoides
description: Trabajo con microestructuras espinodales (spinodoides, Kumar et al. 2020) y su uso como sustitutos sintéticos del hueso trabecular en el proyecto spinpy: definición del generador (GRF, conos, conjunto de nivel), clases de anisotropía, convenciones de muestreo de ondas, ajuste a VOIs de micro-CT, morfometría, homogeneización, cotas de validación (Voigt, Hashin-Shtrikman), percolación, curvaturas, comparación con dual-lattice, pseudorreplicación y redacción de métodos y discusión con las referencias ya verificadas del proyecto. Úsala siempre que el usuario mencione spinodoide, espinodal, spinodoid, Kumar, Zheng, Guo, GRF, conos theta, TPMS, dual-lattice, ajuste a VOI, spinpy o hueso trabecular sintético, aunque no pida explícitamente la skill.
---

# Spinodoides

Esta skill reúne lo que el proyecto spinpy ya estableció sobre los spinodoides: el método publicado, las decisiones de implementación, los hallazgos medidos y las referencias. Su objetivo es que cualquier trabajo nuevo parta de ese conocimiento y no lo rederive ni lo contradiga.

## Flujo de trabajo

1. Identifica qué se pide: generar, ajustar, medir, validar, interpretar o redactar.
2. Lee el archivo de `references/` que corresponda (tabla de abajo). No cargues todos.
3. Antes de afirmar algo sobre el método, distingue si es **citado** (está en el artículo), **deducido** (se infiere de él) o **medido** (resultado del proyecto). Decláralo así en el texto.
4. Si la tarea es de redacción científica, aplica además la skill `redaccion-articulos-cientificos` (estilo, verificación en PubMed y Consensus, guía de reporte).
5. Las ejecuciones largas (réplica de Zheng ~90 min) se lanzan en segundo plano y no se repiten sin motivo.

| Necesidad | Leer |
|---|---|
| Ecuaciones, parámetros, clases, muestreo, ternas de ángulos | `references/metodo-kumar2020.md` |
| Ajuste a VOI, métricas, invariantes, función de error, N efectivo | `references/ajuste-y-morfometria.md` |
| Cotas, percolación, curvaturas, hallazgos y errores ya corregidos | `references/validacion-y-hallazgos.md` |
| Citas completas con DOI y qué respalda cada una | `references/bibliografia.md` |

## Núcleo del método (resumen operativo)

Un spinodoide es el conjunto de nivel de un campo aleatorio gaussiano construido como suma de ondas planas de un único número de onda:

```
GRF(x) = sqrt(2/N) * SUM_i cos( beta * <n_i, x> + gamma_i ),  gamma_i ~ U[0, 2pi)
sólido = { x : GRF(x) <= phi0 },   phi0 = sqrt(2) * erfinv(2*rho - 1)
```

- La **densidad se impone** con `phi0`, no se busca: BV/TV objetivo queda exacto para un campo N(0,1).
- La **anisotropía** viene de restringir las `n_i` a la unión de conos de semiángulos (theta1, theta2, theta3) alrededor de los ejes. Un cono estrecho en un eje produce láminas apiladas a lo largo de ese eje, de modo que ese eje es el **blando**.
- `beta` fija la escala (Tb.Th y Tb.Sp a la vez, no por separado). `N` fija la suavidad estadística. La semilla fija la realización.
- Clases: lamelar (un cono), columnar (dos conos), cúbica (tres conos iguales), isótropa (90, 90, 90).

Módulos del repositorio: `spinpy/grf.py` (generador y umbral), `morphometry.py`, `elastic.py` (homogeneización periódica), `resistencia.py` (Pistoia), `curvatura.py`, `dual_lattice.py`, `fit.py`, `estadistica.py`. Réplicas en `Test/replicar_*.py`.

## Reglas que no se negocian

- **Declarar siempre el muestreo de ondas.** `rechazo` (GIBBON, ecuación 2 del artículo) y `equitativo` (TPMS-Scaffolds-generator) dan anisotropías distintas con conos desiguales (27 % de diferencia en DA con (15, 45, 90) grados). Solo coinciden con conos iguales. Ninguno es "el correcto", pero solo `rechazo` corresponde al artículo.
- **Declarar semilla, número de ondas y resolución** en cualquier cifra publicada. El generador es estocástico y entre N = 400 y 1000 la desconexión a rho = 0.30 cambia de 2 % a 48 %.
- **Una realización no es el ensemble.** Las afirmaciones de simetría del artículo valen para la distribución de direcciones. Las pruebas deben formularse sobre promedios o en relativo, no sobre una semilla.
- **K réplicas de un VOI no son N = K.** Son N = 1 con K réplicas técnicas. Promediar divide la varianza del generador por K y deja intacta la biológica. Reportar N efectivo e ICC y no usar los grados de libertad de las réplicas.
- **Equivalencia, no ausencia de diferencia.** Para sostener que el sintético sustituye al real, usar TOST con margen declarado antes de ver los datos (razonable: el suelo de ruido del generador).
- **Tb.Th = 2·BV/BS** (placas, Parfitt), no 4·BV/BS. **DA del tensor MIL**, con suelo de ruido de ~1.07. **Superficie sin las seis tapas del cubo.** Mismo motor de morfometría para VOI y candidato.
- **Condiciones de contorno:** Kumar y Zheng usan afines (cota superior); spinpy usa periódicas (valor más bajo). Comparar con desigualdades, no cifra a cifra.
- **No extrapolar.** Los resultados proceden de sesamoideos equinos a 51.489 µm (97³ vóxeles) y de réplicas de artículos concretos. Pistoia (2 % del tejido, 0.7 % de deformación) se calibró en radio distal humano y no son constantes físicas.
- **Inconsistencia interna del artículo de Zheng:** theta_min es 30° en el texto y 15° en los pies de figura. No son intercambiables para la regla rho >= 0.3 (ver hallazgos).
- Cuando haya conflicto entre lo publicado y lo medido, se informa ambos y se explica la diferencia; no se ajusta un parámetro no publicado hasta que la figura se parezca.

## Escritura

Si se redacta para un manuscrito: estilo académico natural, sin em dash (U+2014), sin muletillas ni tríadas artificiales, distinguiendo asociación, predicción y efecto. Las referencias de `references/bibliografia.md` provienen de las citas del repositorio; **antes de incluirlas en un manuscrito, verifica cada una en PubMed o Consensus** y comprueba que respalda la frase concreta. Los spinodoides son un diseño computacional y de validación, no observacional ni experimental clásico: no asumas STROBE ni CONSORT. Evalúa si corresponde TRIPOD (modelos predictivos, p. ej. sustitutos estructura a propiedad) o una guía de modelado y simulación, y verifica la versión vigente en EQUATOR antes de usarla.
