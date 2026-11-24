# Agregar un motor de elementos finitos

Esta guía es para quien sabe de elementos finitos y quiere que su biblioteca
(o su propio resolvedor) pueda ensayar los VOIs de `spinpy`, sin tener que
entender el resto del programa. Un motor es **un archivo de Python** y **una
fila de una tabla**. Todo lo demás (interfaz, informe, comparativa entre
motores, pruebas) lo recoge solo.

## Qué hace un motor

Recibe un problema de compresión ya mallado y devuelve el desplazamiento de
cada nodo. Nada más. La malla, las condiciones de apoyo, el postproceso
(módulo aparente, von Mises, criterio de Pistoia) y el informe son comunes, de
modo que, cuando se comparan dos motores, la diferencia es del motor.

```
VOI binario ──► fem.mallar ──► problema_de_malla ──► motor.resolver ──► u, F
                 (hex8/TET10)    (un dict, mm-N-MPa)     (tu archivo)      │
                                                                           ▼
                                                  E aparente, von Mises, Pistoia
```

## Pasos

1. **Copia la plantilla.** `spinpy/motores/_plantilla.py` a
   `spinpy/motores/m_<clave>.py`. La clave es un nombre corto en minúsculas.
2. **Escribe `resolver(p)`.** La plantilla explica el contenido de `p` y lo que
   debe devolver; la descripción completa está en la cabecera de
   `spinpy/motores/_comun.py`. `m_skfem.py` es un ejemplo corto en Python puro
   (ensamblaje propio y Newton propio); `m_ngsolve.py`, uno con solver directo
   multihilo.
3. **Registra el motor.** Añade una ficha a `REGISTRO` en
   `spinpy/motores/__init__.py`: clave, etiqueta, módulo que delata la
   instalación, nombre en PyPI, mallas, materiales no lineales y si resuelve el
   plato rígido. Las demás tablas (`MOTORES`, `ETIQUETAS`, `CAPACIDADES`) se
   derivan de ahí.
4. **Traduce la etiqueta.** Añade el nombre a `spinpy/idioma_textos.py` y a
   `DINAMICAS` en `idioma_revisar.py` (las etiquetas llegan a la interfaz como
   variable y el comprobador no las ve).
5. **Declara el motor en el empaquetado**, si va a distribuirse en el
   ejecutable: `ocultos` en `instalador/spinpy.spec` (y en las variantes de
   Linux y macOS). Sin esto el ejecutable funciona en tu máquina y no en la de
   los demás.
6. **Ejecuta las pruebas.**

   ```
   python -m pytest tests/test_29_motores.py tests/test_36_contrato_motores.py -q
   ```

   El bloque 29 recorre `motores.MOTORES`: pone a tu motor frente a soluciones
   cerradas sin que escribas una prueba. Si el motor no está instalado, esas
   pruebas se saltan y se anota.

## Lo que se le exige a un motor

| Prueba | Qué verifica | Tolerancia |
|---|---|---|
| Bloque macizo, lineal | E aparente = E del material; fuerza = carga aplicada | 1e-9 relativa |
| Plato rígido, lineal | reacción del plato = E·A·ε | 1e-9 relativa |
| No lineal, compresión uniaxial hasta 20 % | fuerza por paso frente a la solución cerrada (St. Venant-Kirchhoff y neo-Hookeano) | 1e-8 relativa |
| Entre motores | mismo espinodoide 16³ da el mismo campo que la app | 1e-6 del máx. de \|u\| |
| Combinación no disponible | lanza `NoDisponible`, nunca otro modelo | exacto |
| Proceso hijo | da lo mismo y se puede cancelar | exacto |
| Convergencia | pendiente L2 = k + 1 con solución exacta | 0,15 |

Las tolerancias están **declaradas antes de medir**: no se relajan para que un
motor pase. Si un motor no las cumple, es un hallazgo y se documenta (ver
`comparativa_motores/INFORME.md`), no se ajusta la tolerancia.

## Reglas que ahorran tiempo

- **No importes la biblioteca al nivel del módulo.** La interfaz sondea los
  motores sin cargarlos (sus DLL chocan con las de VTK y Qt); se importa dentro
  de `resolver`, que corre en un proceso hijo cancelable.
- **Unidades**: mm, N, MPa.
- **Orden de nodos**: el de `spinpy`. Si tu motor renumera, vuelve al orden con
  `_comun.emparejar` (por coordenadas).
- **Lo que no resuelves, se declara**: `NoDisponible`. Un resultado de otro
  modelo, aunque se parezca, falsea la comparación.
- **Si no converge**, `ErrorMotor` con la causa (residuo, iteraciones).

## Cómo comparar tu motor con los demás

```
python -m spinpy.fem VOI.vtk --protocolo app --motores ngsolve skfem <tu_clave>
```

La comparativa original (tiempos, memoria, exactitud con TET10) está en
`comparativa_motores/`; sus adaptadores son independientes de los del paquete.
