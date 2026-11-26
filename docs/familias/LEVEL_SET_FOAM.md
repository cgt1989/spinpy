# Level-set foam: por qué está en spinpy y qué se midió

Esta nota recoge la evidencia con la que se incorporó la tercera familia de
microestructura (`spinpy.level_set`, `spinpy.fit_levelset`). Es lo que el
README resume en «Familias de microestructura». Los scripts para repetir cada
paso están en [`experimento/`](experimento/LEEME.md).

## Qué es y de dónde viene

Espuma de celda abierta definida como el **conjunto de nivel** de campos de
distancia de elipsoides empaquetados. Los puntales siguen las aristas de la
teselación de las semillas y su sección se controla con tres parámetros: la
curvatura de la sección (K_c) y el engrosamiento hacia los nodos (c₁, c₂). La
densidad no es un parámetro libre: el generador resuelve el espesor del puntal
para la fracción volumétrica pedida.

- **Implementación:** [morphoGen](https://github.com/jdirr/morphoGen)
  (Dirrenberger, GPL-3.0-or-later; DOI de todas las versiones
  [10.5281/zenodo.24340315](https://doi.org/10.5281/zenodo.24340315)). spinpy la
  **llama**, no copia nada de ella.
- **Método:** Sonon, François y Massart (2015), *Computational Mechanics*,
  [10.1007/s00466-015-1168-8](https://doi.org/10.1007/s00466-015-1168-8);
  Kilingar, Kamel, Sonon, Massart y Noels (2019), *European Journal of
  Mechanics A/Solids*,
  [10.1016/j.euromechsol.2019.103847](https://doi.org/10.1016/j.euromechsol.2019.103847).
  Ambas se comprobaron (registro completo con DOI) en SciSpace, en noviembre de
  2026.
- **No comprobadas:** la documentación de morphoGen cita además a Sonon et al.
  (2012) y a Degdeg y Dirrenberger (2026). No se localizaron en ninguna base
  bibliográfica consultada (SciSpace; Consensus no estaba disponible por
  agotar su cupo), así que spinpy **no las cita**.

## Qué se midió

Cuatro VOI porcinos (C1, C5, V3, V2; Koria et al. 2020, CC BY 4.0,
`datos/porcino_koria2020/`), cinco configuraciones y el mismo protocolo para
todas, con el código de spinpy:

- candidato a `min(lado, 96)` vóxeles sobre el cubo físico del VOI;
- misma morfometría y misma función de error (`spinpy.error`);
- búsqueda escalonada de 53 a 62 evaluaciones (spinodoide y dual-lattice con su
  ajuste propio; las dos de morphoGen con la misma estructura de etapas);
- el ganador se **reevalúa con 10 semillas nuevas** (1001 a 1010), distintas de
  la de la búsqueda, para no premiar la mejor tirada;
- las propiedades que **no** entran en el objetivo (SMI, número de poros,
  rigidez elástica periódica) se miden aparte: son la comprobación
  independiente.

### 1. Con el objetivo de siempre, el error no discrimina

Con BV/TV, Tb.Th, Tb.Sp, Tb.N, BS/BV, DA y DA₂ el error de las cuatro familias
queda en torno a 1e-3: casi saturado. El dual-lattice gana ese error en 3 de 4
VOI, pero no dice cuál reproduce mejor lo demás.

### 2. Ampliar el objetivo

Se añadieron los tres términos que spinpy ya admitía en `pesos` (tamaño de
poro PoDm, variación del espesor local TbTh_CV y desviación de la curvatura
H_desv; la propia sugerencia `PESOS_DISCRIMINANTES`, con TbSp a 0) y **Conn.D**,
que se registró como término opcional (`error.TERMINOS_OPCIONALES["ConnD"]`).
Los términos de forma del VOI se miden al **mismo tamaño de vóxel** que los del
candidato: curvatura, espesor, poro y conectividad dependen del vóxel.

Con la densidad libre, el ajuste la sacrificó hasta un 15 % para mejorar los
términos de forma y la rigidez empeoró en tres familias. **Fijar la densidad al
BV/TV medido** lo corrigió, y es lo que hace `ajustar_level_set`.

### 3. Resultado con densidad fija (media de los 4 VOI; menor = más cerca del VOI)

| Familia | Error ampliado | SMI | Nº de poros (\|ln razón\|) | Rigidez 16³ | Rigidez 32³ |
|---|---|---|---|---|---|
| Spinodoide | 0,0216 | 1,97 | 1,86 | 0,37 | 1,77 |
| Dual-lattice | 0,0274 | 1,88 | 1,02 | 0,42 | 1,07 |
| `open_foam` (morphoGen) | 0,0183 | 1,59 | 0,31 | 0,39 | 0,54 |
| `level_set_foam` | 0,0151 | 0,25 | 0,37 | 0,29 | 0,38 |
| `level_set_foam`, plantilla de espuma de cobre | 0,0308 | 1,00 | 1,97 | 0,26 | 0,32 |

La rigidez es la distancia log-euclídea entre el tensor elástico periódico
(ejes materiales) del candidato y el del VOI, al mismo `res_mec`.

### 4. Ajustar la forma (K_c, c₁, c₂) con la densidad fija

| Familia | C1 | C5 | V3 | V2 | Media |
|---|---|---|---|---|---|
| `level_set_foam`, forma por defecto | 0,0113 | 0,0162 | 0,0119 | 0,0212 | 0,0151 |
| `level_set_foam`, **forma ajustada** | **0,0053** | **0,0081** | **0,0117** | **0,0144** | **0,0099** |

El error baja un tercio y las propiedades independientes no empeoran
(SMI 0,32; nº de poros 0,30; rigidez 16³ 0,26 y 32³ 0,36). Es el ajuste que
usa la aplicación.

## Qué no dice esta evidencia

- **Una sola población.** Los cuatro VOI son del talus porcino. No dice nada
  sobre el sesamoideo equino, que es el objetivo principal de spinpy.
- **La mecánica es poco fiable.** Una sola realización a 32³, con trabéculas
  apenas resueltas, y el VOI se submuestrea desde 188³ mientras los candidatos
  lo hacen desde 96³. A 16³ ninguna familia se distingue de otra por rigidez.
  Spinpy ya advierte que esos valores sirven para comparar, no para citar.
- **Sin referencia mecánica externa.** El único criterio es el propio VOI.
- **Sin prueba estadística.** Son cuatro VOI; los recuentos «mejor en k de 4» son
  descriptivos.
- **Presupuestos distintos.** El ajuste del spinodoide y del dual-lattice está
  calibrado con VOI equinos; el de morphoGen es una primera pasada y usa los
  presets del dual-lattice. Entre el 16 % y el 31 % de las evaluaciones de
  `level_set_foam` fueron combinaciones que el modelo no admite (el
  empaquetado denso no cabe con pocas semillas muy alargadas).
- **SMI.** Está confundido por la concavidad (Salmon et al. 2015, ya citado por
  spinpy). Que sea la única familia con SMI negativo, como el VOI, indica
  superficies cóncavas parecidas, no necesariamente la misma proporción de
  placas y barras.
- **El sobreajuste es moderado.** El error de la semilla de búsqueda frente al
  de las semillas nuevas pasa de un cociente de 1,06 a 1,14 al ajustar la forma;
  en V3 la mejora de la forma sobre la de por defecto es nula.
- Los parámetros de forma ganadores **no son los mismos** en los cuatro VOI
  (K_c de -0,85 a +0,30; engrosamiento de nodos moderado en 3 de 4), así que no
  hay una forma universal que fijar.

## Decisiones de diseño derivadas

| Decisión | Razón |
|---|---|
| Dependencia **opcional** (`spinpy[level-set]`), importada de forma perezosa | morphoGen es GPL-3.0-or-later y exige Python ≥ 3.12; el núcleo de spinpy es MIT y admite 3.10 |
| La densidad **no se busca**: se fija al BV/TV medido | con la densidad libre la rigidez empeoró (apartado 2) |
| Objetivo ampliado por omisión (`error.PESOS_LEVEL_SET`) | el objetivo base no discrimina (apartado 1) |
| Sin rotación `R`; la orientación la da el estiramiento por eje | la teselación es periódica y sus ejes son los de la caja |
| Orientación por **permutación de ejes** del estiramiento tras la búsqueda | el objetivo no ve la orientación; es exacta (no interpola) |
| Combinaciones inviables = no medibles, no error | el empaquetado denso exige tres diámetros de semilla por eje |
| Cuarto color de las figuras: ámbar `#f2a600` | pasa el validador con pares adyacentes; no con todos los pares (ver `figuras.py`) |
