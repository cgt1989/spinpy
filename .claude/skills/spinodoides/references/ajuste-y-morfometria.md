# Ajuste a VOIs de micro-CT, morfometría y diseño estadístico

## Contexto de datos del proyecto

VOIs cúbicos de 97³ vóxeles a 51.489 µm (4.994 mm de lado) de sesamoideos equinos (H1 a H8), tres sitios por hueso: proximal, medio y distal. En H4:

| Métrica | proximal | medio | distal |
|---|---|---|---|
| BV/TV | 0.2804 | 0.5430 | 0.7652 |
| Tb.Th [mm] | 0.1694 | 0.3654 | 0.6415 |
| Tb.Sp [mm] | 0.4348 | 0.3076 | 0.1969 |
| DA (MIL) | 1.509 | 1.367 | 1.208 |

BV/TV y DA varían acopladas entre sitios. Con tres puntos correlacionados no se separa el efecto de cada variable sobre la rigidez; eso justifica el generador sintético. A BV/TV alto (distal) la descripción trabecular pierde sentido para hueso y spinodoide.

## Métricas y definiciones (motor único `morfometria`)

BV/TV; Po.tot = (1 − BV/TV)·100; BS por marching cubes sin las seis tapas; BS/BV; Tb.Th = 2·BV/BS (Parfitt, placas); Tb.Sp = Tb.Th·(1/BV/TV − 1); Tb.N = BV/TV / Tb.Th; DA del tensor MIL (Harrigan y Mann 1984), DA2 como diagnóstico columnar frente a laminar; fracción portante (material en caminos que cruzan la probeta); Conn.D (Odgaard y Gundersen), SMI (Hildebrand y Rüegsegger) opcionales. El SMI no es fiable a BV/TV alto y no mide barras y placas en hueso (Salmon et al. 2015). Elipsoide Factor: Doube 2015.

Invariantes protegidos (versiones anteriores los rompieron): Tb.Th = 2·BV/BS; DA del MIL y no de la covarianza de la nube de puntos; superficie sin tapas; el elipsoide MIL degenera en estructuras muy laminares (se acota el autovalor y se marca como cota inferior).

## Función de error

Promedio ponderado de diferencias relativas al cuadrado: peso 3 BV/TV, 2 DA, 1 BS/BV, Tb.Th, Tb.Sp, Tb.N, 0.5 Po.tot. Se divide por el peso efectivamente usado y los términos no calculables no se omiten en silencio, porque eso premiaba geometrías degeneradas. El error **no es comparable entre métodos con distinto número de términos**.

## Búsqueda escalonada

A) rejilla densidad por número de onda; B) presets de ángulos con permutaciones (el ángulo grande marca el eje blando); C) refinado. Número de onda barrido en [8, 25] y ángulos barridos (antes el candidato salía casi isótropo, DA 1.1 a 1.2, y Tb.Th se sobreestimaba ~43 % en trabéculas finas). Semilla fija 20260720 en los optimizadores; con semilla fija los errores se reprodujeron idénticos en dos ejecuciones.

Resultado de referencia (proximal H4, 700 ondas): método completo error 0.00283, BV/TV −2.6 %, Tb.Th +5.7 %, Tb.Sp +9.6 %, DA 1.446 frente a 1.509 (diferencia del orden del suelo de ruido del MIL, ~0.07). Tb.Th y Tb.Sp se mueven en sentidos opuestos porque comparten la escala del campo.

## Qué compran y qué no las realizaciones sintéticas

Compran: experimentos numéricos controlados (variar BV/TV con DA fijo), cubrir huecos del muestreo, reducir la varianza del generador por promedio de K realizaciones, alimentar un modelo sustituto estructura a propiedad.

No compran tamaño muestral: `sigma2_total = sigma2_entre + sigma2_dentro`; promediar K divide solo la segunda. Reportar N efectivo e ICC; los grados de libertad los ponen los animales. Pseudorreplicación según Hurlbert (1984). Guardar columnas animal, sitio y réplica separadas para poder descomponer la varianza.

## Brecha de rigidez (problema abierto)

El spinodoide ajustado a BV/TV, Tb.Th y DA resultó ~8 veces más blando que el hueso equino (ver `Estudio_Discriminadores` en el repositorio). Hipótesis a contrastar: es propio del spinodoide o de cualquier estructura calibrada solo con esas tres métricas; el dual-lattice es el contraste previsto. Además el ajuste cae en el "rincón frágil" de la lamelar (rho 0.339, 13.8 % de hueso en islas, cerca del umbral de percolación). Tratar ambas cosas como abiertas, no como resueltas.
