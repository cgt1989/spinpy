# Verificación de spinpy contra literatura y soluciones cerradas

**Invariantes topológicos, definiciones publicadas, soluciones exactas de la elasticidad y cotas universales — con las tolerancias declaradas antes de medir**

## Qué se verifica y con qué criterio

Cada comprobación compara una función de `spinpy` con una respuesta que se conoce **de antemano** y que no depende del propio código: un número de Euler, el valor que una definición publicada asigna a una figura ideal, la solución cerrada de un problema elástico, o una desigualdad que toda microestructura cumple. La tolerancia de cada prueba está escrita en el módulo de pruebas, con su justificación, **antes** de ejecutarla; no se ajusta después para que pase.

Una comprobación que no pasa no es un fallo de la suite: es un hallazgo sobre el código o sobre su documentación, y se lista primero.

> **En palabras sencillas.** No comparamos el programa con otro programa, sino con cosas que se saben seguras: una rosquilla tiene exactamente un agujero, un bloque macizo se comprime exactamente como dice la ley de Hooke, y ningún material poroso puede ser más rígido que el sólido del que está hecho. Si el programa acierta esas cosas, sus números merecen confianza.

Resumen: **61 de 66 comprobaciones dentro de tolerancia**. Ejecutado el 2026-09-01 18:15:14 con Python 3.13.5, numpy 2.1.3, scipy 1.15.3, scikit-image 0.25.0.

## Lo que NO pasó

| Bloque | Prueba | Esperado | Obtenido | Error | Tolerancia |
|---|---|---|---|---|---|
| 04 Espesor local (esferas inscritas) | losa t=3: mediana Th | 3 | 4 | +33,33 % | ¦Th - t¦ <= 0.5 vox |
| 04 Espesor local (esferas inscritas) | losa t=5: mediana Th | 5 | 6 | +20,00 % | ¦Th - t¦ <= 0.5 vox |
| 04 Espesor local (esferas inscritas) | losa t=9: mediana Th | 9 | 10 | +11,11 % | ¦Th - t¦ <= 0.5 vox |
| 04 Espesor local (esferas inscritas) | losa t=15: mediana Th | 15 | 16 | +6,67 % | ¦Th - t¦ <= 0.5 vox |
| 05 Elasticidad (cerradas y cotas) | Gibson-Ashby: exponente n de E/E_s ~ rho^n | — | 3,87899 | — | 1 <= n <= 3 |

Cada uno de estos casos se discute en su bloque.

## 02 Topologia (Euler, Conn.D)

Referencias: Odgaard & Gundersen 1993, Bone 14:173; Euler-Poincare chi = b0-b1+b2.

11 de 11 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| bola maciza: chi | 1 | 1 | +0,00 % | exacta | sí |
| bola maciza: Conn = 1 - chi | 0 | 0 | 0 | exacta | sí |
| toro macizo: chi | 0 | 0 | 0 | exacta | sí |
| toro macizo: Conn = 1 - chi | 1 | 1 | +0,00 % | exacta | sí |
| reticulo 27 nodos / 54 barras: chi | -27 | -27 | +0,00 % | exacta | sí |
| reticulo 27 nodos / 54 barras: Conn = 1 - chi | 28 | 28 | +0,00 % | exacta | sí |
| dos bolas separadas: chi | 2 | 2 | +0,00 % | exacta | sí |
| dos bolas separadas: Conn = 1 - chi | -1 | -1 | +0,00 % | exacta | sí |
| cascara esferica hueca: chi | 2 | 2 | +0,00 % | exacta | sí |
| cascara esferica hueca: Conn = 1 - chi | -1 | -1 | +0,00 % | exacta | sí |
| Conn.D con spacing/2 = 8 x Conn.D | 8 | 8 | +0,00 % | 1e-12 | sí |

- *reticulo 27 nodos / 54 barras: chi*: b1 = E - V + 1 = 28 conexiones redundantes
- *reticulo 27 nodos / 54 barras: Conn = 1 - chi*: b1 = E - V + 1 = 28 conexiones redundantes
- *dos bolas separadas: chi*: Conn negativo: la formula presupone b0 = 1
- *dos bolas separadas: Conn = 1 - chi*: Conn negativo: la formula presupone b0 = 1
- *cascara esferica hueca: chi*: Conn negativo: la formula presupone b2 = 0
- *cascara esferica hueca: Conn = 1 - chi*: Conn negativo: la formula presupone b2 = 0

## 03 SMI (Hildebrand & Ruegsegger)

Referencias: Hildebrand & Ruegsegger 1997; Salmon et al. 2015.

5 de 5 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| losa: SMI | 0 | 0 | 0 | ¦SMI¦ <= 0.05 | sí |
| cilindro d=16: SMI | 3 | 2,71595 | -9,47 % | 2.40 <= SMI <= 3.00 | sí |
| esfera d=28: SMI | 4 | 3,44981 | -13,75 % | 3.20 <= SMI <= 4.00 | sí |
| esfera: dBS/dr > 0 (dilata hacia fuera) | — | 356,553 | — | > 0 | sí |
| cavidad esferica: SMI < 0 (concavidad) | — | -48,4919 | — | < 0 | sí |

- *cavidad esferica: SMI < 0 (concavidad)*: es la razon por la que el SMI no es fiable a BV/TV alto

## 04 Espesor local (esferas inscritas)

Referencias: Hildebrand & Ruegsegger 1997, J Microsc 185:67.

9 de 13 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| losa t=3: mediana Th | 3 | 4 | +33,33 % | ¦Th - t¦ <= 0.5 vox | **NO** |
| losa t=5: mediana Th | 5 | 6 | +20,00 % | ¦Th - t¦ <= 0.5 vox | **NO** |
| losa t=9: mediana Th | 9 | 10 | +11,11 % | ¦Th - t¦ <= 0.5 vox | **NO** |
| losa t=15: mediana Th | 15 | 16 | +6,67 % | ¦Th - t¦ <= 0.5 vox | **NO** |
| cilindro d=4: mediana Th | 4 | 2,82843 | -29,29 % | ¦err¦ <= 40% | sí |
| cilindro d=6: mediana Th | 6 | 5,65685 | -5,72 % | ¦err¦ <= 15% | sí |
| cilindro d=10: mediana Th | 10 | 8,94427 | -10,56 % | ¦err¦ <= 15% | sí |
| cilindro d=20: mediana Th | 20 | 18,868 | -5,66 % | ¦err¦ <= 15% | sí |
| bola d=20: mediana Th | 20 | 18,5472 | -7,26 % | ¦err¦ <= 15 % | sí |
| Th escala linealmente con el spacing | 0,05 | 0,05 | -0,00 % | 1e-06 | sí |
| cilindro d=10, eje ENTRE voxeles: mediana Th | 10 | 8,94427 | -10,56 % | ¦err¦ <= 40 % | sí |
| cilindro d=10, eje SOBRE un voxel: mediana Th | 10 | 10,198 | +1,98 % | ¦err¦ <= 40 % | sí |
| signo del error: entre voxeles < 0 < sobre voxel | — | 0,125377 | — | > 0 | sí |

- *losa t=3: mediana Th*: max = 4.00
- *losa t=5: mediana Th*: max = 6.00
- *losa t=9: mediana Th*: max = 10.00
- *losa t=15: mediana Th*: max = 16.00
- *cilindro d=4: mediana Th*: max = 2.83
- *cilindro d=6: mediana Th*: max = 5.66
- *cilindro d=10: mediana Th*: max = 8.94
- *cilindro d=20: mediana Th*: max = 18.87
- *bola d=20: mediana Th*: max = 18.55
- *cilindro d=10, eje ENTRE voxeles: mediana Th*: error -10.6%
- *cilindro d=10, eje SOBRE un voxel: mediana Th*: error +2.0%
- *signo del error: entre voxeles < 0 < sobre voxel*: entre -10.6%, sobre +2.0%

## 05 Elasticidad (cerradas y cotas)

Referencias: Andreassen & Andreasen 2014; Gibson & Ashby 1997; Kumar 2020; Hashin & Shtrikman 1963; Hill 1952; Hooke; Pistoia 2002, def. energetica; def. von Mises; modulo edometrico E(1-nu)/((1+nu)(1-2nu)).

16 de 17 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| homogeneizacion cubo macizo: max¦Ch - C_iso¦ | 0 | 1.65e-16 | 1.65e-16 | 1e-09 | sí |
| constantes de ingenieria: Ez | 1 | 1 | -0,00 % | 1e-09 | sí |
| constantes de ingenieria: Gxy | 0,384615 | 0,384615 | +0,00 % | 1e-09 | sí |
| constantes de ingenieria: nu_xy | 0,3 | 0,3 | -0,00 % | 1e-09 | sí |
| ensayo deslizante cubo macizo: E_app | 1 | 1 | -0,00 % | 1e-09 | sí |
| eps_eff uniforme = sigma0/E_s (Pistoia uniaxial) | 0,001 | 1.00e-03 | -0,00 % | 1e-07 | sí |
| von Mises uniforme = sigma0 | 0,001 | 0,001 | +0,00 % | 1e-07 | sí |
| ensayo empotrado: E_s <= E_app <= M_oed | 1,34615 | 1,03239 | -23,31 % | desigualdad | sí |
| spinodoide rho=0.32: Ch simetrico | 0 | 0 | 0 | 1.346153846153846e-09 | sí |
| spinodoide rho=0.32: Ch >= 0 (Reuss) | — | 0,00621827 | — | >= -1.3e-09 | sí |
| spinodoide rho=0.32: rho*Cs - Ch >= 0 (Voigt) | — | 0,0998245 | — | >= -1.3e-09 | sí |
| spinodoide rho=0.51: Ch simetrico | 0 | 0 | 0 | 1.346153846153846e-09 | sí |
| spinodoide rho=0.51: Ch >= 0 (Reuss) | — | 0,0485472 | — | >= -1.3e-09 | sí |
| spinodoide rho=0.51: rho*Cs - Ch >= 0 (Voigt) | — | 0,0901405 | — | >= -1.3e-09 | sí |
| HS+: K_hill <= 1.05 K_HS+ (rho=0.39) | 0,165012 | 0,0578986 | -64,91 % | 5 % | sí |
| HS+: G_hill <= 1.05 G_HS+ (rho=0.39) | 0,0974806 | 0,0390565 | -59,93 % | 5 % | sí |
| Gibson-Ashby: exponente n de E/E_s ~ rho^n | — | 3,87899 | — | 1 <= n <= 3 | **NO** |

- *eps_eff uniforme = sigma0/E_s (Pistoia uniaxial)*: en estado uniaxial eps_eff = |eps_zz|
- *ensayo empotrado: E_s <= E_app <= M_oed*: la coaccion lateral solo puede rigidizar
- *Gibson-Ashby: exponente n de E/E_s ~ rho^n*: rho=[0.262, 0.35, 0.509]  E/E_s=[0.0167, 0.0661, 0.2236]

## 06 Criterio de Pistoia

Referencias: Pistoia et al. 2002, Bone 30:842.

8 de 8 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| cubo macizo: sigma_fallo = eps_crit * E_s | 0,007 | 0,007 | -0,00 % | 1e-09 | sí |
| cubo macizo: F_fallo = sigma_fallo * A_bruta | 0,7 | 0,7 | -0,00 % | 1e-09 | sí |
| cubo macizo: factor = eps_crit * E_s / sigma0 | 7 | 7 | -0,00 % | 1e-09 | sí |
| sigma_fallo independiente de sigma0 (linealidad) | 0,007 | 0,007 | +0,00 % | 1e-09 | sí |
| factor se divide por 2 al duplicar sigma0 | 2 | 2 | +0,00 % | 1e-09 | sí |
| frac = 2 % selecciona el grupo mas deformado | 0,002 | 0,002 | +0,00 % | 1e-12 | sí |
| frac = 10 % cae en el grupo base | 0,001 | 0,001 | +0,00 % | 1e-12 | sí |
| factor(2 %) = factor(10 %) / 2 | 3,5 | 3,5 | +0,00 % | 1e-12 | sí |

## 07 Medida sobre malla y estanqueidad

Referencias: Taubin 1995; Attene 2010 (PyMeshFix); Parfitt 1987.

12 de 12 dentro de tolerancia.

| Prueba | Esperado | Obtenido | Error | Tolerancia | OK |
|---|---|---|---|---|---|
| esfera r=14: BS malla | 2463,01 | 2478,14 | +0,61 % | 2 % | sí |
| esfera r=14: BV malla | 11494 | 11505,1 | +0,10 % | 1 % | sí |
| esfera: la malla mejora el area de voxeles | — | 0,0714919 | — | < 0.5 | sí |
| losa t=10: BS = 2 L^2 (tapas excluidas) | 1922 | 1938,58 | +0,86 % | 1 % | sí |
| losa t=10: BV = L^2 t | 9610 | 9617,97 | +0,08 % | 1 % | sí |
| TV malla = ((n-1)/n)^3 TV voxeles | 1,52088 | 1,52088 | -0,00 % | 1e-12 | sí |
| spinodoide 40^3: estanca | — | 0 | — | 0 bordes | sí |
| spinodoide 40^3: perdida vs superficie cruda (misma caja) | 0 | -0,02348 | -0,02348 | ¦.¦ <= 3 % | sí |
| spinodoide 40^3: BS malla / BS voxel | 0,915 | 0,910778 | -0,46 % | [0.85, 0.97] | sí |
| spinodoide 40^3: BV/TV malla / BV/TV voxel - 1 | 0 | -0,0452979 | -0,0452979 | ¦.¦ <= 8 % (sanidad) | sí |
| TET10 40^3: superficie estanca antes de tetgen | — | 0 | — | 0 bordes | sí |
| TET10 40^3: volumen frente a superficie cruda | 100 | 96,4861 | -3,51 % | >= 95 % | sí |

- *esfera r=14: BS malla*: voxeles: +8.6 %
- *losa t=10: BS = 2 L^2 (tapas excluidas)*: tapas excluidas: 1193.8
- *spinodoide 40^3: BV/TV malla / BV/TV voxel - 1*: definicion MC vs conteo -4.6 %; perdida vs MC -2.3 %
- *TET10 40^3: superficie estanca antes de tetgen*: 1 componente(s)
- *TET10 40^3: volumen frente a superficie cruda*: frente al conteo de voxeles: 92.0 %
