"""
interfaz: ventanas secundarias de la aplicacion grafica.

La ventana principal sigue en `visor.py`, en la raiz del proyecto, porque es
el punto de entrada: la receta de PyInstaller, `instalador/construir.ps1` y la
orden documentada (`python visor.py`) apuntan a el. Aqui viven los dialogos
que el visor abre:

  dialogo_febio        FEM automatico y «Analizar con FEBio…»
  dialogo_metodos      ajustar con todos los metodos, comparar y elegir
  dialogo_validacion   las replicas de los resultados publicados (`Test/`)

Como `visor.py`, dependen de PyQt5 (GPL v3). El paquete `spinpy` no importa
nada de esta carpeta, y por eso sigue libre de Qt y de sus obligaciones de
licencia.
"""
