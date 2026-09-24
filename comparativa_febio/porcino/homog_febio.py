"""homog_febio.py — La homogeneizacion periodica de spinpy, escrita para FEBio 4.

`elastic.homogeneizar` resuelve la celda unidad sobre la rejilla de voxeles con
conectividad periodica (el nodo n+1 ES el nodo 1) y rellena los poros con un
material de rigidez 1e-6*E_s. FEBio no tiene esa conectividad, pero el MISMO
problema discreto se escribe con restricciones lineales:

    u(nodo imagen en la cara maxima) = u(nodo en la cara minima) + E . (x' - x)

Cada nodo con algun indice igual a n es DEPENDIENTE (FEBio lo elimina del
sistema) de su imagen modulo n, que nunca es dependiente: no hay cadenas. El
desplazamiento resultante es E.x + fluctuacion periodica, exactamente el campo
de spinpy (u0 - chi). La traslacion rigida se quita fijando el nodo (0,0,0),
igual que spinpy fija el nodo 0.

La rigidez homogeneizada sale del promedio de volumen de la tension,
C[:, j] = <sigma> para la deformacion unitaria j. spinpy la calcula por
energia, (u0-chi)' K (u0-chi) / V; en la solucion exacta del problema discreto
las dos son la misma cantidad. El promedio de los 8 puntos de Gauss de un
hexaedro rectangular es su integral / V, asi que la media de las tensiones de
elemento del logfile ES <sigma>.

FEBio es no lineal geometricamente: se aplica la deformacion a dos amplitudes y
se extrapola a amplitud nula, igual que el ensayo de compresion.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")


def tensor_de_voigt(e):
    """Voigt [xx yy zz yz xz xy] con distorsiones de INGENIERIA -> tensor."""
    return np.array([[e[0], e[5] / 2, e[4] / 2],
                     [e[5] / 2, e[1], e[3] / 2],
                     [e[4] / 2, e[3] / 2, e[2]]], float)


def malla_completa(n_xyz, spacing):
    """Nodos (n+1)^3 y hexaedros n^3 en el orden de spinpy (Fortran, i antes).

    Devuelve nodos (N,3) en mm, elems (M,8) base 0, y los indices (i,j,k) de
    cada nodo. El orden de los 8 nodos es el de `_edof_periodico`.
    """
    nx, ny, nz = n_xyz
    dx, dy, dz = spacing
    I, J, K = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1),
                          np.arange(nz + 1), indexing="ij")
    I, J, K = (a.ravel(order="F") for a in (I, J, K))
    nodos = np.stack([I * dx, J * dy, K * dz], axis=1).astype(float)

    def nid(a, b, c):
        return a + (nx + 1) * b + (nx + 1) * (ny + 1) * c

    ex, ey, ez = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                             indexing="ij")
    ex, ey, ez = (a.ravel(order="F") for a in (ex, ey, ez))
    elems = np.stack([nid(ex, ey, ez), nid(ex + 1, ey, ez),
                      nid(ex + 1, ey + 1, ez), nid(ex, ey + 1, ez),
                      nid(ex, ey, ez + 1), nid(ex + 1, ey, ez + 1),
                      nid(ex + 1, ey + 1, ez + 1), nid(ex, ey + 1, ez + 1)],
                     axis=1)
    return nodos, elems, np.stack([I, J, K], axis=1)


def escribir_periodico(BW, spacing, ruta, E_voigt, E_s=20e9, nu_s=0.30,
                       escala_vacio=1e-6):
    """Escribe la celda periodica bajo la deformacion macroscopica E_voigt.

    `BW` en orden (x, y, z). Dos materiales: hueso (E_s) y vacio
    (escala_vacio*E_s), el mismo relleno que `homogeneizar`. Unidades mm-N-MPa.
    """
    ruta = Path(ruta)
    BW = np.asarray(BW, bool)
    nx, ny, nz = BW.shape
    sp = np.asarray(spacing, float)
    nodos, elems, ijk = malla_completa(BW.shape, sp)
    solido = BW.ravel(order="F")
    Emac = tensor_de_voigt(np.asarray(E_voigt, float))
    E_MPa = float(E_s) / 1e6
    stem = ruta.stem
    n_lim = np.array([nx, ny, nz])

    dep = np.where((ijk == n_lim).any(axis=1))[0]
    img = ijk[dep] % n_lim
    ind = img[:, 0] + (nx + 1) * img[:, 1] + (nx + 1) * (ny + 1) * img[:, 2]
    salto = nodos[dep] - nodos[ind]                      # (m, 3) mm
    offset = salto @ Emac.T                              # u_dep - u_ind

    with open(ruta, "w", encoding="utf-8") as fh:
        w = fh.write
        w('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        w("<!-- Celda periodica generada por comparativa_febio/porcino/"
          "homog_febio.py.\n     E_macro (Voigt, ingenieria) = %s -->\n"
          % ", ".join(f"{v:.6g}" for v in E_voigt))
        w('<febio_spec version="4.0">\n')
        w('\t<Module type="solid">\n\t\t<units>mm-N-s</units>\n\t</Module>\n')
        w("\t<Control>\n\t\t<analysis>STATIC</analysis>\n")
        w("\t\t<time_steps>1</time_steps>\n\t\t<step_size>1</step_size>\n")
        w("\t\t<solver>\n\t\t\t<max_refs>50</max_refs>\n")
        w('\t\t\t<qn_method type="BFGS">\n\t\t\t\t<max_ups>0</max_ups>\n'
          "\t\t\t</qn_method>\n\t\t\t<dtol>1e-9</dtol>\n")
        w("\t\t\t<etol>1e-12</etol>\n\t\t\t<rtol>1e-12</rtol>\n")
        w("\t\t\t<lstol>0.9</lstol>\n")
        w('\t\t\t<linear_solver type="pardiso"/>\n\t\t</solver>\n')
        w("\t</Control>\n")
        w("\t<Material>\n")
        for i, (nom, E) in enumerate((("hueso", E_MPa),
                                      ("vacio", E_MPa * escala_vacio)), 1):
            w(f'\t\t<material id="{i}" name="{nom}" type="isotropic elastic">'
              f"\n\t\t\t<E>{E:.17g}</E>\n\t\t\t<v>{nu_s:.9g}</v>\n"
              "\t\t</material>\n")
        w("\t</Material>\n")
        w('\t<Mesh>\n\t\t<Nodes name="todos">\n')
        w("".join(f'\t\t\t<node id="{i}">{p[0]:.12g},{p[1]:.12g},{p[2]:.12g}'
                  "</node>\n" for i, p in enumerate(nodos, start=1)))
        w("\t\t</Nodes>\n")
        ids = np.arange(1, elems.shape[0] + 1)
        for nom, m in (("hueso", solido), ("vacio", ~solido)):
            if not m.any():
                continue
            w(f'\t\t<Elements type="hex8" name="{nom}">\n')
            w("".join(f'\t\t\t<elem id="{i}">' + ",".join(map(str, c))
                      + "</elem>\n"
                      for i, c in zip(ids[m].tolist(),
                                      (elems[m] + 1).tolist())))
            w("\t\t</Elements>\n")
        w('\t\t<NodeSet name="origen">1</NodeSet>\n\t</Mesh>\n')
        w("\t<MeshDomains>\n")
        for nom, m in (("hueso", solido), ("vacio", ~solido)):
            if m.any():
                w(f'\t\t<SolidDomain name="{nom}" mat="{nom}"/>\n')
        w("\t</MeshDomains>\n")
        w("\t<Boundary>\n")
        w('\t\t<bc name="origen" type="zero displacement" node_set="origen">'
          "\n\t\t\t<x_dof>1</x_dof>\n\t\t\t<y_dof>1</y_dof>\n"
          "\t\t\t<z_dof>1</z_dof>\n\t\t</bc>\n")
        partes = []
        for d_, i_, off in zip((dep + 1).tolist(), (ind + 1).tolist(),
                               offset.tolist()):
            for c, dn in enumerate("xyz"):
                partes.append(
                    f'\t\t<bc type="linear constraint"><node>{d_}</node>'
                    f"<dof>{dn}</dof><offset>{off[c]:.17g}</offset>"
                    f"<child_dof><node>{i_}</node><dof>{dn}</dof>"
                    "<value>1</value></child_dof></bc>\n")
        w("".join(partes))
        w("\t</Boundary>\n")
        w('\t<LoadData>\n\t\t<load_controller id="1" name="rampa" '
          'type="loadcurve">\n')
        w("\t\t\t<interpolate>LINEAR</interpolate>\n\t\t\t<points>\n")
        w("\t\t\t\t<point>0,0</point>\n\t\t\t\t<point>1,1</point>\n")
        w("\t\t\t</points>\n\t\t</load_controller>\n\t</LoadData>\n")
        # El plotfile NO es opcional: FEBio 4.5.0 termina con una violacion de
        # acceso (codigo 0xC0000005) si la seccion Output solo tiene logfile.
        # Comprobado en un bloque de 27 elementos, con y sin restricciones.
        w('\t<Output>\n\t\t<plotfile type="febio">\n'
          '\t\t\t<var type="displacement"/>\n\t\t</plotfile>\n\t\t<logfile>\n')
        w(f'\t\t\t<element_data data="sx;sy;sz;syz;sxz;sxy" delim=" " '
          f'file="{stem}_s.txt"/>\n')
        w(f'\t\t\t<node_data data="ux;uy;uz" delim=" " file="{stem}_u.txt"/>\n')
        w("\t\t</logfile>\n\t</Output>\n</febio_spec>\n")

    return {"ruta": str(ruta), "n_nodos": int(nodos.shape[0]),
            "n_elems": int(elems.shape[0]), "n_dependientes": int(dep.size),
            "n_restricciones": int(3 * dep.size),
            "MB": ruta.stat().st_size / 1e6}
