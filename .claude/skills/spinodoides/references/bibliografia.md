# Bibliografía del proyecto

Tomada de las citas del repositorio (`Test/README.md`, `Test/referencia/FUENTES.md`, `spinpy/informe.py`, `CITATION.cff`). **No sustituye la verificación**: antes de usarlas en un manuscrito, confirma cada una en PubMed o Consensus y comprueba que respalda la frase concreta. Las marcadas (v) tienen DOI en el repositorio; las demás faltan o deben comprobarse.

## Núcleo espinodal

- Kumar S, Tan S, Zheng L, Kochmann DM. Inverse-designed spinodoid metamaterials. npj Comput Mater. 2020;6:73. doi:10.1038/s41524-020-0341-6 (CC BY 4.0). Respalda: generador GRF, conos, conjunto de nivel, clases, diseño inverso por red neuronal.
- Zheng L, Kumar S, Kochmann DM. Data-driven topology optimization of spinodoid metamaterials with seamlessly tunable anisotropy. Comput Methods Appl Mech Eng. 2021;383:113894. doi:10.1016/j.cma.2021.113894 (CC BY). Respalda: ternas de las cuatro clases, rho >= 0.3, theta_min, cotas de Voigt y Hashin-Shtrikman, ortotropía, condiciones afines.
- Guo Y, Sharma S, Kumar S. Inverse designing surface curvatures by deep learning. Adv Intell Syst. 2024;6(6):2300789. doi:10.1002/aisy.202300789 (CC BY 4.0). Respalda: perfil de curvaturas (k1, k2) ponderado por área, PNS de referencia, Fig. 7.
- Cahn JW. Phase separation by spinodal decomposition in isotropic systems. J Chem Phys. 1965;42(1):93-99. doi:10.1063/1.1695731. Respalda: origen físico (descomposición espinodal).
- Moerman KM. GIBBON: The geometry and image-based bioengineering add-on. J Open Source Softw. 2018;3(22):506. doi:10.21105/joss.00506. Respalda: implementación `spinodoid.m` del muestreo por rechazo.
- Vafaeefar M, Moerman KM, Kavousi M, Vaughan TJ. A morphological, topological and mechanical investigation of gyroid, spinodoid and dual-lattice algorithms as structural models of trabecular bone. J Mech Behav Biomed Mater. 2023;138:105584. doi:10.1016/j.jmbbm.2022.105584. Respalda: dual-lattice, comparación de familias. El código la cita como 2022 (año de aceptación) y el CFF como 2022; unificar año al citar.

## Morfometría ósea

- Parfitt AM, Drezner MK, Glorieux FH, et al. Bone histomorphometry: standardization of nomenclature, symbols, and units. J Bone Miner Res. 1987;2(6):595-610. doi:10.1002/jbmr.5650020617. Respalda: Tb.Th = 2·BV/BS (modelo de placas).
- Harrigan TP, Mann RW. Characterization of microstructural anisotropy in orthotropic materials using a second rank tensor. J Mater Sci. 1984;19(3):761-767. doi:10.1007/BF00540446. Respalda: tensor MIL, DA.
- Odgaard A, Gundersen HJG. Quantification of connectivity in cancellous bone, with special emphasis on 3-D reconstructions. Bone. 1993;14(2):173-182. doi:10.1016/8756-3282(93)90245-6. Respalda: Conn.D, Euler-Poincaré. (DOC_SPINODOIDES.md cita Odgaard 1997, Bone 20:315; conciliar.)
- Hildebrand T, Rüegsegger P. Quantification of bone microarchitecture with the structure model index. Comput Methods Biomech Biomed Engin. 1997;1(1):15-23. doi:10.1080/01495739708936692. Respalda: SMI.
- Salmon PL, Ohlsson C, Shefelbine SJ, Doube M. Structure model index does not measure rods and plates in trabecular bone. Front Endocrinol. 2015;6:162. doi:10.3389/fendo.2015.00162. Respalda: límites del SMI.
- Doube M. The Ellipsoid Factor for quantification of rods, plates, and intermediate forms in 3D geometries. Front Endocrinol. 2015 (completar volumen y DOI en la verificación).
- Doube M, Kłosowski MM, Arganda-Carreras I, et al. BoneJ: Free and extensible bone image analysis in ImageJ. Bone. 2010;47(6):1076-1079. doi:10.1016/j.bone.2010.08.023.
- Bouxsein ML, Boyd SK, Christiansen BA, et al. Guidelines for assessment of bone microstructure in rodents using micro-computed tomography. J Bone Miner Res. 2010;25(7):1468-1486. doi:10.1002/jbmr.141.
- Callens SJP, Tourolle né Betts DC, Müller R, Zadpoor AA. The local and global geometry of trabecular bone. Acta Biomater. 2021;130:343-361. doi:10.1016/j.actbio.2021.06.013. Respalda: curvatura del hueso trabecular.
- Rusinkiewicz S. Estimating curvatures and their derivatives on triangle meshes. 3DPVT 2004:486-493. doi:10.1109/TDPVT.2004.1335277. Respalda: estimador discreto.

## Mecánica y homogeneización

- Andreassen E, Andreasen CS. How to determine composite material properties using numerical homogenization. Comput Mater Sci. 2014;83:488-495. doi:10.1016/j.commatsci.2013.09.006.
- Hill R. The elastic behaviour of a crystalline aggregate. Proc Phys Soc A. 1952;65(5):349-354. doi:10.1088/0370-1298/65/5/307.
- Hashin Z, Shtrikman S. A variational approach to the theory of the elastic behaviour of multiphase materials. J Mech Phys Solids. 1963;11(2):127-140. doi:10.1016/0022-5096(63)90060-7.
- Gibson LJ, Ashby MF. Cellular solids: structure and properties. 2nd ed. Cambridge University Press; 1997. doi:10.1017/CBO9781139878326.
- Pistoia W, van Rietbergen B, Lochmüller EM, et al. Estimation of distal radius failure load with micro-finite element analysis models based on three-dimensional peripheral quantitative computed tomography images. Bone. 2002;30(6):842-848. doi:10.1016/S8756-3282(02)00736-6. Parámetros calibrados en radio distal humano.

## Diseño estadístico

- Hurlbert SH. Pseudoreplication and the design of ecological field experiments. Ecol Monogr. 1984;54(2):187-211. Respalda: pseudorreplicación.

## Otras

- Tozzi et al. (muestra de hueso trabecular usada por Guo 2024; el repositorio la cita como 2017 sin referencia completa): no replicada; citar solo a través de Guo salvo verificación directa.
- Tapia D, González A, Vidal F, Salinas P. Biology. 2026;15:722 (protocolo micro-CT/FEA de referencia del proyecto). Verificar antes de citar.

## Herramientas

Attene 2010 (PyMeshFix), Taubin 1995 (suavizado), Si 2015 (TetGen), Sullivan y Kaszynski 2019 (PyVista), Virtanen 2020 (SciPy), Harris 2020 (NumPy), van der Walt 2014 (scikit-image): citas completas con DOI en `Test/README.md`, sección "Herramientas".
