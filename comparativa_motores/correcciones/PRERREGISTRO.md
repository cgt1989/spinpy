# Correcciones de A1, A2 y A3: criterio de aceptacion (escrito antes de medir)

Fecha: 2026-10-05. Rama local `experimento/a1-a3` (worktree temporal; no se
sube). Solo lo que cumpla este criterio pasa a la rama de trabajo.

## Referencia (la "verdad" sin artefactos de borde)

Configuracion embebida: el mismo espinodoide (`casos.ESPINODOIDE`, misma
semilla) evaluado en un bloque de 7,5 mm que contiene el VOI de 5 mm rodeado
por 1,25 mm de hueso por cada lado (n/4 voxeles, identicos voxel a voxel en
el interior). Bloque con plato rigido en el techo y base deslizante; todo se
mide SOLO en la region del VOI (o de la ROI que reporte cada metodo).
Sensibilidad de la referencia: traccion frente a plato en el bloque a 48^3.

Definiciones comunes a referencia y metodos (region R = caja de lados
paralelos a los del VOI):
  sigma_R = -sum_{e en R} sigma_zz,e V_e / V_caja(R)
  eps_R   = extensometro de franjas: (u_z medio - u_z medio) / (z medio -
            z medio) entre la franja superior y la inferior de R, de 0,25 mm
            de espesor cada una, con u_z y z en el centroide de cada elemento
            de hueso y medias ponderadas por volumen. Se calcula igual con
            hex8 y TET10 (enmienda anterior a cualquier medida: el
            extensometro por nodos en un plano no existe en TET10).
  E_R     = sigma_R / eps_R
  p99_R   = p99 de von Mises (capa superficial, ponderado) en R / sigma_R

## Metodos evaluados (hex8 a 32^3, 48^3 y 64^3)

  B0  linea base de la app: traccion uniforme en el techo; R = VOI; E es
      el que publica hoy la app (definicion por area del techo).
  F1  (A1) plato rigido lineal en el techo; R = VOI; E = F/A / (u_plato/H),
      lo que publicaria la app con plato.
  F2  (A2) traccion; R = nucleo a >= 0,6 mm de las seis caras del VOI
      (medido dentro del ensayo del VOI; extensometro interno).
  F12 F1 + F2.

## Metricas

  e_E   = |E_R(metodo) / E_R(referencia) - 1|
  e_p99 = |p99_R(metodo) / p99_R(referencia) - 1|
  campo: mediana de |vm - vm_ref| / vm_ref en la capa superficial de R
         (voxel a voxel) y fraccion de la cola de la referencia (> p99)
         recuperada en la cola del metodo.

## Regla de decision

Una correccion se ACEPTA si, en las TRES resoluciones:
  1. reduce e_E al menos a la mitad respecto de B0 Y en al menos 5 puntos
     porcentuales, y
  2. no empeora e_p99 en mas de 2 puntos porcentuales.
Si una correccion mejora pero no cumple 1 en las tres resoluciones, la mejora
se considera LEVE y se descarta. Si empeora, se descarta.

## A3 (tipo de malla)

TET10 a 32^3 (seleccion) y 48^3 (validacion), bajo el mejor protocolo
aceptado para A1/A2 (o B0 si ninguno se acepta). Referencia: bloque embebido
hex8 a la misma resolucion. Candidatos de malla (opciones de
`fem.opciones_malla`): por omision (Taubin 20, decimado 0,5); Taubin 5 sin
decimado; sin suavizado ni decimado. Se elige en 32^3 el de menor e_E y se
ACEPTA solo si en 48^3 reduce e_E respecto de la malla por omision al menos
a la mitad y en al menos 5 puntos, sin empeorar e_p99 en mas de 2 puntos.
Ademas: A3 se da por resuelto con el protocolo si la brecha TET10/hex8 en E
bajo ese protocolo es <= 10 % en 48^3.

## Integracion en la app (solo si se acepta)

Codigo minimo en `spinpy/fem.py` (y lo que dependa), con pruebas nuevas en
`tests/`, y la suite de motores (bloque 29) sin fallos.
