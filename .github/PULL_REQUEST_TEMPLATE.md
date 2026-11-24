## Qué cambia y por qué

<!-- Una o dos frases. Si corrige un error, di cuál era el resultado malo y cuál el correcto, con la fuente. -->

## Cómo se verificó

- [ ] Hay una prueba nueva o ampliada en `tests/`, con la tolerancia **declarada antes de medir**
- [ ] La prueba compara contra una respuesta conocida (analítica, definición publicada), no contra la salida anterior del mismo código
- [ ] `python -m pytest tests/ -q -m "not lento"` pasa
- [ ] `ruff check .` da 0 avisos

## Si toca la interfaz

- [ ] Los textos nuevos están en `idioma_textos.py` y `python idioma_revisar.py` da 0 sin traducir
- [ ] La lógica de cálculo está en `spinpy/`, no en `visor.py`

## Si añade o cambia un motor FEM

- [ ] Ficha en `REGISTRO` y archivo `m_<clave>.py` (ver `docs/AGREGAR_UN_MOTOR_FEM.md`)
- [ ] `tests/test_29_motores.py` y `tests/test_36_contrato_motores.py` pasan con el motor instalado

## Si cambia el número de bloques de pruebas

- [ ] Insignia y texto de `README.md` y `README.en.md` actualizados
