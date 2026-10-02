# Protocolo de coherencia con autores previos

## 1. Ficha de convenciones (completar antes de analizar)

Copiar, rellenar y conservar con los resultados. Cada campo lleva su origen: C (citado de un artículo, indicar cuál), D (deducido), M (medido o elegido aquí).

```
GENERADOR
- Familia: spinodoide GRF / dual-lattice / otra
- Muestreo de ondas: rechazo | equitativo            origen:
- Ternas (theta1, theta2, theta3) en grados y eje de cada una   origen:
- theta_min usado (15° o 30°) y por qué             origen:
- beta (múltiplos de pi) y relación con el tamaño de dominio    origen:
- N de ondas                                          origen:
- Semilla(s) y número de realizaciones
- Umbral: erfinv exacta o aproximada; densidad NOMINAL y densidad REALIZADA (medida en la máscara)
DOMINIO
- Resolución (vóxeles), tamaño físico, periodicidad, orientación (rotaciones)
MEDICIÓN
- Definición de cada métrica (Tb.Th, DA, BS con o sin tapas, Conn.D, SMI)
- Método de superficie (marching cubes sobre binario o sobre el campo continuo)
MECÁNICA
- Condiciones de contorno (afines, periódicas, Dirichlet, ensayo empotrado o deslizante)
- E del sólido, nu; solver
OBJETIVO DEL ANÁLISIS
- ¿Se ajusta a rigidez (como Kumar, Deng, Otto) o a morfometría (como Vafaeefar, spinpy)?
- Criterio de éxito y tolerancias declaradas ANTES de medir
ESTADÍSTICA
- Unidad de análisis (animal, sitio, réplica), N efectivo, ICC, equivalencia (margen TOST)
```

## 2. Clasificar cada diferencia con un autor

| Clase | Definición | Qué hacer |
|---|---|---|
| Equivalente | misma definición con distinta notación o implementación, diferencia numérica acotada | citar al autor y dar la correspondencia |
| Desviación declarada | convención distinta pero explícita y justificada (p. ej. BC periódicas, métricas de morfometría) | declararla en métodos y escribir las comparaciones como desigualdades o como complementarias |
| Conflicto | resultado propio que contradice una afirmación del autor en condiciones comparables | revisar primero la ficha (convención, régimen, N, semilla); si persiste, reportar ambos con la discrepancia y no suavizarla |

Un conflicto solo existe si coinciden objetivo, convención y régimen. Si falta uno, es una desviación declarada.

## 3. Comprobaciones antes de comparar con un artículo

1. ¿El artículo publica los parámetros de la figura o resultado que se quiere replicar? Si no (p. ej. beta de la Fig. 2 de Kumar), no inventarlos.
2. ¿Usó el mismo objetivo de ajuste?
3. ¿Misma condición de contorno y misma definición de rigidez (tensor completo, módulos direccionales)?
4. ¿Comparan promedios sobre ensemble o una realización?
5. ¿Misma resolución y N? La percolación a rho = 0.30 cambió de 2 % a 48 % entre N = 400 y 1000.
6. ¿El texto y las figuras del artículo coinciden entre sí? Anotar las inconsistencias.

## 4. Cómo redactar un contraste con la literatura

- Fórmulas modelo: "En condiciones comparables a las de X, nuestros resultados son consistentes con..."; "A diferencia de X, que ajustó el tensor de rigidez, calibramos morfometría, por lo que no es esperable que..."; "La diferencia puede deberse a la convención de muestreo (rechazo frente a reparto equitativo), que no está especificada en X con el detalle necesario".
- Distinguir asociación, predicción y efecto; no extrapolar de hueso equino a otras especies ni de un régimen de densidad a otro sin justificarlo.
- Citar solo lo verificado. Marcar `[VERIFICAR: autor, año]` cuando una cita no se haya comprobado.

## 5. Registro de decisiones

Al terminar un análisis, dejar junto a los resultados la ficha, la lista de desviaciones declaradas y las comprobaciones que fallaron. Eso permite que un análisis posterior parta de lo establecido sin reinterpretarlo.
