---
name: spinodoides
description: Guía para cualquier análisis con microestructuras espinodales (spinodoides, Kumar et al. 2020): generación (GRF, conos, conjunto de nivel), ajuste a VOIs de hueso trabecular, morfometría, homogeneización, validación, curvaturas, percolación, diseño inverso y redacción, asegurando coherencia con lo publicado por Kumar, Zheng, Guo, Vafaeefar, Deng, Otto y otros autores y con los hallazgos previos del proyecto spinpy. Impone un protocolo de convenciones (ficha antes de analizar, comparación con cada autor, clasificación de desviaciones) para que ningún resultado nuevo contradiga o se compare de forma indebida con trabajos previos. Úsala siempre que el usuario mencione spinodoide, espinodal, spinodoid, Kumar, Zheng, Guo, GRF, conos theta, dual-lattice, TPMS, ajuste a VOI, spinpy o hueso trabecular sintético, aunque no pida la skill explícitamente.
---

# Spinodoides: guía de análisis coherente con la literatura

Esta skill sirve para **guiar cualquier análisis posterior** sobre spinodoides, no para repetir uno anterior. Su principio rector: **ningún resultado nuevo se produce ni se compara con un autor previo sin haber declarado antes la convención con la que se generó y se midió**. Casi todas las aparentes contradicciones con la literatura de spinodoides vienen de convenciones distintas (muestreo de ondas, condiciones de contorno, definición de métricas, objetivo del ajuste), no de errores.

## Flujo de trabajo obligatorio

1. **Clasifica la tarea**: generar, ajustar a un VOI, medir, homogeneizar, validar, diseño inverso o sustituto, interpretar o redactar.
2. **Rellena la ficha de convenciones** (plantilla en `references/protocolo-coherencia.md`) con los valores que usará el análisis. Si un valor no está definido, pregunta o declara el supuesto; no lo dejes implícito.
3. **Compara la ficha con las convenciones de cada autor relevante** (`references/convenciones-autores.md`). Clasifica cada diferencia como: equivalente, desviación declarada o conflicto.
4. **Aplica las reglas de diseño** (abajo). Cuando el análisis toque un tema específico, lee solo el archivo que corresponda.
5. **Antes de afirmar algo contra o a favor de un autor**, comprueba si compararon lo mismo (mismo objetivo de ajuste, misma convención, mismo régimen). Si no, formula la afirmación como complementaria, no como contradicción.
6. **Reporta** con la sección de convenciones y las desviaciones declaradas. Si se redacta un manuscrito, usa además la skill `redaccion-articulos-cientificos` y verifica cada cita en PubMed, Consensus o Scholar Gateway.

| Tema | Leer |
|---|---|
| Ficha, clasificación de desviaciones, redacción de contrastes | `references/protocolo-coherencia.md` |
| Qué hizo y qué asumió cada autor (con estado de verificación) | `references/convenciones-autores.md` |
| Ecuaciones, ternas de ángulos, muestreo | `references/metodo-kumar2020.md` |
| Ajuste a VOI, métricas, N efectivo (implementación spinpy) | `references/ajuste-y-morfometria.md` |
| Cotas, percolación, curvaturas, errores ya corregidos | `references/validacion-y-hallazgos.md` |
| Citas con DOI y su estado de verificación | `references/bibliografia.md` |

## Núcleo del método

```
GRF(x) = sqrt(2/N) * SUM_i cos( beta * <n_i, x> + gamma_i ),   gamma_i ~ U[0, 2pi)
sólido = { x : GRF(x) <= phi0 },   phi0 = sqrt(2) * erfinv(2*rho - 1)
```

Las `n_i` se restringen a la unión de tres conos de semiángulos (theta1, theta2, theta3) alrededor de los ejes; un cono estrecho en un eje produce láminas apiladas a lo largo de él y ese eje es el blando. `beta` fija la escala, `N` la suavidad estadística, la semilla la realización. Clases: lamelar, columnar, cúbica, isótropa. Este es el marco de Kumar et al. (2020) y lo reutilizan Zheng, Deng, Golnary, Raßloff, Otto y Röding.

## Reglas de coherencia que no se negocian

- **Convención declarada, siempre**: muestreo de ondas (`rechazo` o `equitativo`), theta_min, beta, N, semilla, resolución, tamaño y periodicidad del dominio, y densidad nominal frente a densidad realizada. La redacción de Deng et al. ("uniformly sampled from within these cones") es ambigua entre ambos muestreos, de modo que citar un autor no sustituye declarar la convención.
- **Compara solo lo comparable.** Kumar y Deng ajustan el tensor de rigidez (objetivo mecánico); Vafaeefar y spinpy calibran morfometría (BV/TV, Tb.Th, DA). Un spinodoide ajustado a morfometría puede ser más blando que el hueso sin contradecir a quienes ajustan la rigidez directamente.
- **Una realización no es el ensemble.** Las simetrías y las afirmaciones de ortotropía valen para la distribución de direcciones. Park et al. (2026) muestran que parámetros de cono idénticos dan dispersión de propiedades por la naturaleza estocástica del GRF. Probar sobre promedios o con varias semillas.
- **K réplicas no son N = K.** Son N = 1 con K réplicas técnicas (pseudorreplicación). Reportar N efectivo e ICC; usar TOST con margen previo para afirmar sustitución.
- **Condiciones de contorno**: afines dan cota superior (Kumar, Zheng), periódicas dan valores menores, Dirichlet es otra opción (Risthaus y Schneider). Comparar con desigualdades y declarar cuál se usó.
- **Definiciones de métricas**: Tb.Th = 2·BV/BS (placas), DA del tensor MIL (suelo de ruido ~1.07), superficie sin las seis tapas, mismo motor para VOI y candidato. Cambiar una definición invalida la comparación con cifras previas.
- **Régimen de validez**: el GRF aproxima solo las etapas iniciales de la descomposición espinodal (Mandolesi et al., 2025); no extrapolar a morfologías tardías ni a otros tejidos sin justificarlo. Los resultados de spinpy proceden de sesamoideos equinos a 51.489 µm.
- **Citado, deducido o medido**: etiquetar cada parámetro y afirmación así. No ajustar un parámetro no publicado hasta que una figura se parezca.
- **Discrepancias internas de los autores** (p. ej. theta_min 30° en el texto de Zheng y 15° en pies de figura) se reportan, no se resuelven en silencio.
- **Hallazgos que fallan se listan primero**; las pruebas se reformulan solo si el error era de la prueba, y se documenta por qué.

## Escritura

Estilo académico natural, sin em dash (U+2014), sin muletillas ni tríadas artificiales; distinguir asociación, predicción y efecto. Las referencias de `references/bibliografia.md` indican su estado: **verificada** en una base científica o **tomada del repositorio, pendiente**. Solo las verificadas pueden entrar sin comprobación adicional, y aun así debe confirmarse que respaldan la frase concreta. El trabajo con spinodoides es computacional y de validación: no asumas STROBE ni CONSORT; evalúa TRIPOD para sustitutos predictivos u otra guía apropiada y verifica la versión vigente en EQUATOR.
