"""
fuentes_libros.py: Trazabilidad de cada cambio de la revision de la
biblioteca hasta el parrafo del libro del que sale.

Cada fuente da el libro, la seccion, la pagina impresa, la pagina del archivo
de la biblioteca (PDF o DjVu del repositorio privado cgt1989/LibrosFEM), el
parrafo (por sus primeras palabras) y la cita textual en el idioma original.
Las citas se transcribieron del texto extraido con pdftotext, djvutxt o, en
Reddy (2002), por reconocimiento optico; se corrigieron solo artefactos de la
extraccion (cortes de palabra con guion, ligaduras como «fi», numeracion de
ecuaciones mal leida). Las omisiones se marcan con [...]. «Respalda» dice que
afirmacion del informe sostiene la cita y, cuando corresponde, que parte es
inferencia de spinpy y no del texto.

`componer_informe.py` usa este modulo para escribir los bloques «De donde
sale» de cada ficha y la matriz de trazabilidad del anexo.
"""

LIBROS = {
    "bathe1996": ("Bathe KJ. Finite element procedures. Englewood Cliffs: "
                  "Prentice Hall; 1996.",
                  "Klaus-Jurgen Bathe - Finite element procedures. Part 1-2-"
                  "Prentice Hall (1996).djvu"),
    "dolean2026": ("Dolean V, Tabeart J. Advanced linear algebra with "
                   "applications, Part I. Lecture notes. arXiv:2608.21234; "
                   "2026.",
                   "Advanced Lineal Algebra with Applications Part 1.pdf"),
    "nocedal2006": ("Nocedal J, Wright SJ. Numerical optimization. 2nd ed. "
                    "New York: Springer; 2006.",
                    "Numerical Optimization by Jorge Nocedal, Stephen Wright"
                    ".pdf"),
    "kamensky2022": ("Kamensky D. Finite element analysis for coupled "
                     "problems. Lecture notes for MAE 207. University of "
                     "California San Diego; 2022.",
                     "Finite_element_analysis_for_coupled_problems.pdf"),
    "dacorogna2004": ("Dacorogna B. Introduction to the calculus of "
                      "variations. London: Imperial College Press; 2004.",
                      "INTRODUCTION TO THE CALCULUS OF VARIATIONS - Bernard "
                      "Dacorogna.pdf"),
    "quarteroni2007": ("Quarteroni A, Sacco R, Saleri F. Numerical "
                       "mathematics. 2nd ed. Berlin: Springer; 2007. "
                       "doi:10.1007/b98885",
                       "Quarteroni, Sacco, Saleri - Numerical Mathematics "
                       "(2007).pdf"),
    "brenner2008": ("Brenner SC, Scott LR. The mathematical theory of finite "
                    "element methods. 3rd ed. New York: Springer; 2008. "
                    "doi:10.1007/978-0-387-75934-0",
                    "Susanne Brenner. The mathematical theory of finite "
                    "element methods.pdf"),
    "modersitzki2004": ("Modersitzki J. Numerical methods for image "
                        "registration. Oxford: Oxford University Press; 2004.",
                        "Numerical_Methods_for_Image_Registration_"
                        "Modersitzki.pdf"),
    "zienkiewicz2005": ("Zienkiewicz OC, Taylor RL, Zhu JZ. The finite element "
                        "method: its basis and fundamentals. 6th ed. Oxford: "
                        "Elsevier Butterworth-Heinemann; 2005.",
                        "Zienkiewicz 2005, The finite element method its "
                        "basis and fundamentals.pdf"),
    "reddy2002": ("Reddy JN. Energy principles and variational methods in "
                  "applied mechanics. 2nd ed. New York: Wiley; 2002.",
                  "J. N. Reddy - Energy Principles and Variational Methods in "
                  "Applied Mechanics -Wiley (2002).djvu"),
    "solin2004": ("Šolín P, Segeth K, Doležel I. Higher-order finite element "
                  "methods. Boca Raton: Chapman & Hall/CRC; 2004.",
                  "Higher-Order Finite Element Methods.pdf"),
    "bevill2009": ("Bevill G, Keaveny TM. Trabecular bone strength "
                   "predictions using finite element analysis of micro-scale "
                   "images at limited spatial resolution. Bone. "
                   "2009;44(4):579-84. doi:10.1016/j.bone.2008.11.020",
                   "PubMed 19135184 (resumen)"),
}

FUENTES = [
    # ------------------------------------------------------------- C1
    dict(id="F1", cambio="C1", libro="dolean2026",
         ubicacion="§6.3.2 «Coarse-grid correction», p. 87",
         archivo_pag="93",
         parrafo="párrafo 1 de §6.3.2 y recuadro siguiente («Coarse-grid "
                 "correction cannot be an iterative method on its own»), con "
                 "el párrafo que lo cierra",
         cita="Section 6.3.1 left us with an iterate whose error is smooth: "
              "the smoother has removed the oscillatory components and stalls "
              "on the rest. [...] Coarse-grid correction cannot be an "
              "iterative method on its own. [...] Coarse-grid correction "
              "removes the part of the error that the coarse space can see, "
              "and by construction can do nothing about the rest. This is the "
              "precise sense in which the two ingredients are complementary, "
              "and it is why neither is optional. The smoother cannot touch "
              "smooth error; the coarse correction cannot touch anything "
              "outside its range.",
         respalda="El multigrid solo elimina el error de baja energía que su "
                  "espacio grueso puede representar. Que en elasticidad esos "
                  "modos sean los movimientos rígidos de cada trabécula es la "
                  "aplicación de spinpy (ya documentada en resistencia.py), "
                  "no una afirmación del texto."),
    dict(id="F2", cambio="C1", libro="dolean2026",
         ubicacion="§6.3.1, recuadro «Takeaways», p. 86", archivo_pag="92",
         parrafo="recuadro «Takeaways», tras la proposición 6.3.3",
         cita="Weighted Jacobi with ω = 2/3 is not a scalable solver, but it "
              "is an effective smoother: it damps oscillatory components "
              "uniformly while leaving smooth components for the coarse grid "
              "to remove.",
         respalda="El suavizador no elimina el error suave; esa tarea "
                  "corresponde al espacio grueso."),
    # ------------------------------------------------------------- C2
    dict(id="F3", cambio="C2", libro="bathe1996",
         ubicacion="§8.4 «Solution of nonlinear equations», p. 755",
         archivo_pag="768",
         parrafo="párrafo 1 de la p. 755, antes de §8.4.1",
         cita="The methods we now present are basic techniques that in "
              "practice would be combined in a self-adaptive procedure that "
              "chooses load steps, iterative method, and convergence criteria "
              "automatically depending on the problem considered and solution "
              "accuracy sought.",
         respalda="Combinar Newton, búsqueda lineal y elección automática del "
                  "incremento de carga en un mismo procedimiento."),
    dict(id="F4", cambio="C2", libro="bathe1996",
         ubicacion="§8.4.1 «Newton-Raphson schemes», pp. 757-758",
         archivo_pag="770-771",
         parrafo="último párrafo de la p. 757 («The practical consequence of "
                 "these properties») y párrafo 2 de la p. 758 («In an "
                 "effective finite element program»)",
         cita="The practical consequence of these properties is that if the "
              "current solution iterate is sufficiently close to the solution "
              "U* and if the tangent stiffness matrix does not change "
              "abruptly, we can expect rapid (i.e., quadratic) convergence. "
              "The assumption is of course that the exact tangent stiffness "
              "matrix is used in the iteration [...]. On the other hand, if "
              "the current solution iterate is not sufficiently close to U* "
              "and/or the stiffness matrix used is not the exact tangent "
              "matrix and/or changes abruptly, then the iteration may diverge. "
              "In an effective finite element program, the exact tangent "
              "stiffness matrix will be used, if possible, and hence the "
              "primary procedure for reaching convergence (if convergence "
              "difficulties are encountered) is to decrease the magnitude of "
              "the load step.",
         respalda="Newton completo diverge lejos de la solución; con la "
                  "tangente exacta (la de NGSolve lo es), el remedio principal "
                  "es reducir el incremento de carga."),
    dict(id="F5", cambio="C2", libro="bathe1996",
         ubicacion="§8.4.2 «The BFGS method», p. 761", archivo_pag="774",
         parrafo="párrafo 3 de la p. 761 («As pointed out above, the line "
                 "search is an integral part»)",
         cita="As pointed out above, the line search is an integral part of "
              "the solution method. Of course, such line searches as "
              "performed in (8.98) and (8.99) can also be used in the "
              "Newton-Raphson methods presented in Section 8.4.1. With the "
              "line search performed within an iteration (i), the expense of "
              "the iteration increases, but fewer iterations may be needed for "
              "convergence. Also, the line search may prevent divergence of "
              "the iterations, and in practice, this increased robustness is "
              "the major reason why a line search can in general be "
              "effective.",
         respalda="Añadir búsqueda lineal al Newton completo para evitar la "
                  "divergencia."),
    dict(id="F6", cambio="C2", libro="nocedal2006",
         ubicacion="§3.1 «Step length», «The Wolfe conditions», p. 33",
         archivo_pag="55",
         parrafo="párrafos 1 y 2 de «The Wolfe conditions»",
         cita="A popular inexact line search condition stipulates that αk "
              "should first of all give sufficient decrease in the objective "
              "function f, as measured by the following inequality: "
              "f(xk + α pk) ≤ f(xk) + c1 α ∇fkᵀ pk, (3.4) for some constant "
              "c1 ∈ (0, 1). In other words, the reduction in f should be "
              "proportional to both the step length αk and the directional "
              "derivative ∇fkᵀ pk. Inequality (3.4) is sometimes called the "
              "Armijo condition. [...] In practice, c1 is chosen to be quite "
              "small, say c1 = 10⁻⁴.",
         respalda="La condición de aceptación del paso sobre la energía "
                  "potencial (f) y el valor c1 = 10⁻⁴ del código."),
    dict(id="F7", cambio="C2", libro="nocedal2006",
         ubicacion="§3.3 «Rate of convergence», p. 41", archivo_pag="63",
         parrafo="párrafo 2 de §3.3 («Algorithmic strategies that achieve "
                 "rapid convergence»)",
         cita="Algorithmic strategies that achieve rapid convergence can "
              "sometimes conflict with the requirements of global "
              "convergence, and vice versa. For example, the steepest descent "
              "method is the quintessential globally convergent algorithm, "
              "but it is quite slow in practice, as we shall see below. On the "
              "other hand, the pure Newton iteration converges rapidly when "
              "started close enough to a solution, but its steps may not even "
              "be descent directions away from the solution. The challenge is "
              "to design algorithms that incorporate both properties: good "
              "global convergence guarantees and a rapid rate of convergence.",
         respalda="El orden elegido: Newton completo primero (rápido cerca de "
                  "la solución) y globalización solo si falla. La versión con "
                  "búsqueda lineal siempre activa ralentizaba casos que Newton "
                  "resolvía (medido)."),
    dict(id="F8", cambio="C2", libro="nocedal2006",
         ubicacion="§3.4 «Newton's method with Hessian modification», p. 48",
         archivo_pag="70", parrafo="párrafo 1 de §3.4",
         cita="Away from the solution, the Hessian matrix ∇²f(x) may not be "
              "positive definite, so the Newton direction pkN defined by "
              "∇²f(xk) pkN = −∇f(xk) (3.38) (see (3.30)) may not be a "
              "descent direction.",
         respalda="El recuento de direcciones de Newton que no son de "
                  "descenso como señal de tangente no definida positiva."),
    dict(id="F9", cambio="C2", libro="bathe1996",
         ubicacion="§8.4.3 «Load-displacement-constraint methods», p. 762",
         archivo_pag="775",
         parrafo="párrafo 2 de la p. 762 («In order to calculate the response "
                 "in Fig. 8.13»)",
         cita="In order to calculate the response in Fig. 8.13 initially "
              "relatively large load increments can be employed, but as the "
              "collapse of the structural model is approached, the load "
              "increments must become smaller and there is also the "
              "difficulty of traversing the collapse point. At that point, "
              "the stiffness matrix is singular (the slope of the "
              "load-displacement response curve is zero), and beyond that "
              "point a special solution procedure that allows for a decrease "
              "in load and an increase in displacement must be used to "
              "calculate the ensuing response.",
         respalda="Que la tangente deje de ser definida positiva cerca del "
                  "colapso es estructural, no numérico; y que pasado el punto "
                  "límite con fuerza impuesta haga falta otro método (R13)."),
    dict(id="F10", cambio="C2", libro="bathe1996",
         ubicacion="§10.2 «Fundamental facts used in the solution of "
                   "eigensystems», p. 849", archivo_pag="862",
         parrafo="párrafo 2 de la p. 849 («An important fact that follows»)",
         cita="The important fact is that in the decomposition of K − μM, "
              "the number of negative elements in D is equal to the number of "
              "eigenvalues smaller than μ.",
         respalda="Mejora pendiente: contar los pivotes negativos de la "
                  "factorización LDLᵀ de la tangente (μ = 0) daría el número "
                  "exacto de modos inestables. No está implementado."),
    # ------------------------------------------------------------- C3
    dict(id="F11", cambio="C3", libro="kamensky2022",
         ubicacion="§4.3.4 (hiperelasticidad), p. 57", archivo_pag="64",
         parrafo="párrafo que empieza «The small-deformation correspondence "
                 "between E and the small strain»",
         cita="The small-deformation correspondence between E and the small "
              "strain from linear elasticity motivates the formally-simplest "
              "hyperelastic model, the St. Venant–Kirchhoff model, "
              "ψ = ½ E : C : E, (4.51) where C is the rank-4 elasticity "
              "tensor. However, this model is not stable under compression, "
              "and it can be seen that the energy density does not diverge as "
              "J → 0, as would be needed to prevent non-physical degenerated "
              "deformations. The St. Venant–Kirchhoff model is only "
              "well-suited to problems with small deformations (but "
              "potentially large rotations, unlike linear elasticity). A "
              "slightly more complex model that retains stability under "
              "compression is the compressible neo-Hookean model, [...]",
         respalda="Cambiar el material por omisión de SVK a neo-Hookeano."),
    dict(id="F12", cambio="C3", libro="bathe1996",
         ubicacion="§6.6.1 «Elastic material behavior. Generalization of "
                   "Hooke's law», p. 584", archivo_pag="597",
         parrafo="párrafo 2 de la p. 584 («Considering this material "
                 "description») y párrafo 4 («The preceding observations»)",
         cita="[...] However, an important observation is that in large "
              "displacement and large rotation but small strain analysis, the "
              "relation in (6.184) provides a natural material description "
              "because the components of the second Piola-Kirchhoff stress "
              "and Green-Lagrange strain tensors do not change under rigid "
              "body rotations [...]. The preceding observations are of "
              "special importance because, in practice, Hooke's law is "
              "applicable only to small strains and because there are many "
              "engineering problems in which large displacements, large "
              "rotations, but only small strain conditions are encountered.",
         respalda="El SVK (ecuación 6.184) es adecuado con rotaciones grandes "
                  "y deformaciones pequeñas."),
    dict(id="F13", cambio="C3", libro="bathe1996",
         ubicacion="§6.6.1, p. 589", archivo_pag="602",
         parrafo="último párrafo de la p. 589 («However, when large strains "
                 "are modeled»)",
         cita="However, when large strains are modeled using (6.184) and "
              "(6.197) with the same elastic material constants, completely "
              "different response predictions must be expected.",
         respalda="Con deformaciones grandes, la respuesta depende de la "
                  "elección de la ley; no da igual SVK que otra."),
    dict(id="F14", cambio="C3", libro="dacorogna2004",
         ubicacion="§3.5 «The vectorial case», observación 3.20 (tras el "
                   "teorema 3.19), p. 100", archivo_pag="113",
         parrafo="apartado (v) de la observación 3.20",
         cita="(v) A function f that can be written in terms of a convex "
              "function F as in the theorem is called polyconvex. The theorem "
              "is due to Morrey (see also Ball [7] for important applications "
              "of such results to non linear elasticity).",
         respalda="El marco de la policonvexidad. Que el neo-Hookeano "
                  "compresible sea policonvexo y el SVK no es un resultado de "
                  "Ball (referencia [7] del texto), no una afirmación de "
                  "Dacorogna."),
    # ------------------------------------------------------------- C4
    dict(id="F15", cambio="C4", libro="quarteroni2007",
         ubicacion="§4.6.2 «A stopping test based on the residual», p. 174",
         archivo_pag="193", parrafo="la subsección completa",
         cita="A different stopping criterion consists of continuing the "
              "iteration until ‖r(k)‖ ≤ ε, ε being a fixed tolerance. Note "
              "that ‖x − x(k)‖ = ‖A⁻¹b − x(k)‖ = ‖A⁻¹r(k)‖ ≤ ‖A⁻¹‖ ε. "
              "Considering instead a normalized residual, i.e. stopping the "
              "iteration as soon as ‖r(k)‖/‖b‖ ≤ ε, we obtain the following "
              "control on the relative error ‖x − x(k)‖/‖x‖ ≤ "
              "‖A⁻¹‖ ‖r(k)‖/‖x‖ ≤ K(A) ‖r(k)‖/‖b‖ ≤ εK(A). In the case of "
              "preconditioned methods, the residual is replaced by the "
              "preconditioned residual, so that the previous criterion "
              "becomes ‖P⁻¹r(k)‖/‖P⁻¹r(0)‖ ≤ ε, where P is the "
              "preconditioning matrix.",
         respalda="La cota ε·κ y su versión con el residuo precondicionado."),
    dict(id="F16", cambio="C4", libro="quarteroni2007",
         ubicacion="§4.4.3, observación 4.5 «The conjugate gradient method», "
                   "p. 168", archivo_pag="187",
         parrafo="observación 4.5, párrafo 1",
         cita="If A is symmetric and positive definite, starting from the "
              "Lanczos method for linear systems it is possible to derive the "
              "conjugate gradient method already introduced in Section 4.3.4 "
              "(see [Saa96]). The conjugate gradient method is a variant of "
              "the Lanczos method where the orthonormalization process "
              "remains incomplete.",
         respalda="Que los coeficientes del CG definen la tridiagonal de "
                  "Lanczos con la que se estimó κ."),
    dict(id="F17", cambio="C4", libro="brenner2008",
         ubicacion="§9.8 «Applications to the conjugate-gradient method», "
                   "pp. 266-267", archivo_pag="277-278",
         parrafo="párrafos 1 y 2 de §9.8",
         cita="The conjugate-gradient method for solving a linear system of "
              "the form AU = F is an iterative method whose convergence "
              "properties can be estimated in terms of the condition number "
              "of A (cf. Luenberger 1973). [...] ‖U − U(k)‖_A ≤ "
              "C exp(−2k/√κ₂(A)) ‖U‖_A [...] This estimate says that to "
              "reduce the relative error ‖u − u(k)‖_a/‖u‖_a to O(ε) requires "
              "at most k = O(√κ₂(A) |log ε|) iterations.",
         respalda="La relación entre κ, iteraciones y error del CG."),
    # ------------------------------------------------------------- C5
    dict(id="F18", cambio="C5", libro="modersitzki2004",
         ubicacion="cap. 5 «Principal axes-based registration», p. 46",
         archivo_pag="56", parrafo="párrafo 1 de la p. 46",
         cita="For normalization purposes, we arrange the columns of D_B "
              "such that for the standard deviations we have σ_B,1 ≥ ... ≥ "
              "σ_B,d ≥ 0. If the eigenvalues of Cov_B are simple, the "
              "decomposition is essentially unique, up to a sign in the "
              "columns of D_B.",
         respalda="Los ejes PCA solo están bien definidos (salvo el signo) "
                  "si las varianzas son distintas."),
    dict(id="F19", cambio="C5", libro="quarteroni2007",
         ubicacion="§5.2.1, propiedad 5.5, p. 190", archivo_pag="208",
         parrafo="último párrafo de §5.2.1 («We conclude this section with a "
                 "stability result») y propiedad 5.5",
         cita="We conclude this section with a stability result for the "
              "approximation of the eigenvector associated with a simple "
              "eigenvalue. [...] Property 5.5 The eigenvectors xk and xk(ε) "
              "of the matrices A and A(ε) = A + εE, with ‖xk(ε)‖₂ = ‖xk‖₂ = "
              "1 for k = 1, ..., n, satisfy ‖xk(ε) − xk‖₂ ≤ ε‖E‖₂ / "
              "min_{j≠k} |λk − λj| + O(ε²), ∀k = 1, ..., n. Analogous to "
              "(5.11), the quantity κ(xk) = 1/min_{j≠k} |λk − λj| can be "
              "regarded as being the condition number of the eigenvector xk. "
              "Computing xk might be an ill-conditioned operation if some "
              "eigenvalues λj are \"very close\" to the eigenvalue λk "
              "associated with xk.",
         respalda="La sensibilidad en grados por 1 % que calcula "
                  "voi.estabilidad_pca, con ‖E‖ = 0,01 λk."),
    # ------------------------------------------------------------- C6
    dict(id="F20", cambio="C6", libro="zienkiewicz2005",
         ubicacion="§9.8 «Application of the patch test to an incompatible "
                   "element», p. 343", archivo_pag="354",
         parrafo="párrafos 1 y 2 de §9.8",
         cita="In order to demonstrate the use of the patch test for a finite "
              "element formulation which violates the usually stated "
              "requirements for shape function continuity, we consider the "
              "plane strain incompatible modes first introduced by Wilson et "
              "al. and discussed by Taylor et al. The specific incompatible "
              "formulation considered uses the element displacement "
              "approximations: û = Na ũa + N1n α1 + N2n α2 (9.12) where Na "
              "(a = 1, ..., 4) are the usual conforming bilinear shape "
              "functions and the last two terms are incompatible modes of "
              "deformation defined by the hierarchical functions "
              "N1n = 1 − ξ² and N2n = 1 − η² (9.13) defined independently for "
              "each element. The shape functions used are illustrated in Fig. "
              "9.10. The first, a set of standard bilinear type, gives a "
              "displacement pattern which, as shown in Fig. 9.10(b), "
              "introduces spurious shear strains in pure bending. The second, "
              "in which the parameters α1 and α2 are strictly associated with "
              "a specific element, therefore introduces incompatibility but "
              "assures correct bending behaviour in an individual rectangular "
              "element.",
         respalda="Las funciones 1 − ξ², 1 − η² y el defecto que corrigen "
                  "(corte espurio en flexión). La extensión a tres "
                  "dimensiones con 1 − ζ² es la de spinpy."),
    dict(id="F21", cambio="C6", libro="bathe1996",
         ubicacion="§4.4.1 «Incompatible displacement-based models», "
                   "pp. 262-263", archivo_pag="275-276",
         parrafo="párrafo 2 de §4.4.1 («Since in finite element analysis "
                 "using incompatible») y párrafo de la p. 263 («When "
                 "considering displacement-based elements with "
                 "incompatibilities»)",
         cita="Since in finite element analysis using incompatible "
              "(nonconforming) elements the requirements presented in Section "
              "4.3.2 are not satisfied, the calculated total potential energy "
              "is not necessarily an upper bound to the exact total potential "
              "energy of the system, and consequently, monotonic convergence "
              "is not ensured. [...] When considering displacement-based "
              "elements with incompatibilities, if the patch test is passed, "
              "convergence is ensured (although convergence may not be "
              "monotonic and convergence may be slow).",
         respalda="hex8i pasa la prueba de la parcela, luego converge, pero "
                  "la energía deja de ser una cota y la convergencia puede no "
                  "ser monótona."),
    dict(id="F22", cambio="C6", libro="bathe1996",
         ubicacion="ejemplo 4.28, pp. 267-268", archivo_pag="280-281",
         parrafo="tras la ecuación (a) de la p. 267 y primer párrafo de la "
                 "p. 268",
         cita="In practice, the incompatible displacement parameters α would "
              "now be statically condensed out to obtain the element "
              "stiffness matrix corresponding to only the u degrees of "
              "freedom. [...] We can now easily check that the condition in "
              "(c) is satisfied for the square element: [...] However, we can "
              "also check that the condition is not satisfied for the general "
              "quadrilateral element.",
         respalda="La condensación estática de los modos y que en elementos "
                  "cuadrados (cubos, en la malla de vóxeles) la prueba de la "
                  "parcela se cumple sin corrección."),
    dict(id="F23", cambio="C6", libro="reddy2002",
         ubicacion="§9.4 «Finite element models of the Timoshenko beam "
                   "theory», p. 467 (texto por reconocimiento óptico)",
         archivo_pag="240 (página doble, mitad derecha)",
         parrafo="párrafo 1 de la p. 467",
         cita="This is equivalent to dφ/dx = 0, which is an incorrect "
              "condition to be satisfied as it forces the curvature and hence "
              "the bending energy to zero. Thus, the finite element equations "
              "[...], in an effort to satisfy the constraints [...], will "
              "yield the trivial solution [...]. This is known in the finite "
              "element literature as shear locking.",
         respalda="El bloqueo por cortante de elementos de bajo orden en "
                  "flexión, la misma raíz que la rigidez excesiva del "
                  "trilineal. El texto trata vigas, no hexaedros."),
    dict(id="F24", cambio="C6", libro="bevill2009",
         ubicacion="resumen (Bone 2009;44(4):579-84)",
         archivo_pag="PubMed 19135184",
         parrafo="frase 4 del resumen",
         cita="However, at a voxel size of 120 microm, the predictive ability "
              "of yield stress was slightly less than that of stiffness, "
              "likely due to the large convergence-related errors that could "
              "develop with larger element sizes.",
         respalda="Contexto: con elementos grandes en hueso trabecular pueden "
                  "aparecer errores de convergencia importantes (Sección 3)."),
    # ------------------------------------------------------------- C7
    dict(id="F25", cambio="C7", libro="solin2004",
         ubicacion="§2.2.4 «Brick master element», p. 62", archivo_pag="80",
         parrafo="párrafo 1 de §2.2.4",
         cita="The first three-dimensional master element of arbitrary "
              "order, K1B, will be associated with the reference brick "
              "domain [...]",
         respalda="Elementos de orden arbitrario sobre el hexaedro de "
                  "referencia, la base del Q2 de NGSolve."),
    dict(id="F26", cambio="C7", libro="solin2004",
         ubicacion="§3.5.9 «Static condensation of internal DOF», "
                   "pp. 191-192", archivo_pag="208-209",
         parrafo="párrafo 1 de §3.5.9 y párrafo 2 de la p. 192",
         cita="This procedure is a great tool for the parallelization of p "
              "and hp finite element codes for linear problems. It consists "
              "of three steps: 1. reduction of the size of the discrete "
              "problem by leaving out all bubble functions, 2. solution of "
              "the reduced linear system and 3. calculation of the remaining "
              "coefficients for the bubble functions by solving elementwise "
              "local linear problems. [...] The global degrees of freedom "
              "associated with basis functions that are nonzero within single "
              "mesh elements (bubble functions) are viewed as internal, while "
              "the remaining DOF (associated with basis functions that are "
              "nonzero within more than one mesh element) are by definition "
              "external.",
         respalda="La condensación estática como vía para reducir el sistema "
                  "global de Q2. En las pruebas de este informe NO se usó: se "
                  "resolvió el sistema completo."),
]

#: Afirmaciones del informe que no salen de los libros sino de las medidas.
MEDIDAS = {
    "C1": "Umbral del respaldo LU (45 000 GDL) y su coste: p1_homogeneizacion.py.",
    "C2": "Que la búsqueda lineal siempre activa ralentizaba casos que Newton "
          "completo resolvía (TET10 a 20³: 119 frente a 19 iteraciones): "
          "p2_newton.py y la prueba documentada en m_ngsolve.py.",
    "C3": "Las cifras de tensión de los dos materiales: soluciones cerradas "
          "de comparativa_motores/comun.py (p3_material.py).",
    "C4": "Que la cota falla en TET10: p4_cota_error.py.",
    "C5": "El umbral de 5 °/% y su relación con el giro medido: p5_pca.py.",
    "C6": "Las mejoras en la viga y en el espinodoide: p6_flexion.py y "
          "p6_espinodoide.py.",
    "C7": "GDL, memoria y error: p7_hex_orden2.py.",
}
