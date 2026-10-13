// Construye Tapia_C1_C5.pptx a partir de figs/ y resultados/.
//   node construir_presentacion.js
// Las cifras salen de los JSON, no se escriben a mano.
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");
const SKILL = "/root/.claude/skills/synced/0f07208f-0f93-477c-a320-b385defad481_0e654f15-27ad-4f71-9a8c-1cb528d44549/pptx";
const { applyTheme } = require(path.join(SKILL, "scripts/apply_theme.js"));

const AQUI = __dirname;
const R = JSON.parse(fs.readFileSync(path.join(AQUI, "figs/resumen.json"), "utf8"));
const B = JSON.parse(fs.readFileSync(path.join(AQUI, "figs/bloque.json"), "utf8"));
const S = JSON.parse(fs.readFileSync(path.join(AQUI, "resultados/sanidad_cubo.json"), "utf8"));

const THEME = {
  name: "Hueso trabecular",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2A33", lt1: "FFFFFF", dk2: "1F3A44", lt2: "EEF2F3",
    accent1: "C8553D", accent2: "2E7D7A", accent3: "E0A030",
    accent4: "5B6573", accent5: "8FB8B5", accent6: "6B5B95",
    hlink: "2E7D7A", folHlink: "5B6573",
  },
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";            // 13.33 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Protocolo de Tapia con C1 a C5";
pres.author = "spinpy";
const C = pres.SchemeColor;

// ---------------------------------------------------------------- cifras
const ESTR = ["voi", "spinodoide", "dual-lattice"];
const NOM = { voi: "VOI", spinodoide: "Spinodoide", "dual-lattice": "Dual-lattice" };
const VOIS = ["C1", "C2", "C3", "C4", "C5"];
const T = {};
R.tabla.forEach((t) => (T[t.voi + "|" + t.estructura] = t));
const g = (v, e, k) => T[v + "|" + e][k];
const f1 = (x) => x.toFixed(1).replace(".", ",");
const f2 = (x) => x.toFixed(2).replace(".", ",");
const f3 = (x) => x.toFixed(3).replace(".", ",");
const f4 = (x) => x.toFixed(4).replace(".", ",");
const pct = (x) => (100 * x).toFixed(0);
const razon = (e, k) => VOIS.map((v) => g(v, e, k) / g(v, "voi", k));
const rng = (a) => [Math.min(...a), Math.max(...a)];
const rE_s = rng(razon("spinodoide", "E_app_MPa")), rE_d = rng(razon("dual-lattice", "E_app_MPa"));
const rV_s = rng(razon("spinodoide", "vm_p99_MPa")), rV_d = rng(razon("dual-lattice", "vm_p99_MPa"));
const e_ok = S.filas.filter((f) => f.ok === true).length;
const e_info = S.filas.filter((f) => f.ok === null).length;
const e_mal = S.filas.filter((f) => f.ok === false).length;
const sum = (g) => {
  const fs_ = S.filas.filter((f) => g(f.grupo));
  const ok = fs_.filter((f) => f.ok === true);
  return { n: ok.length, info: fs_.filter((f) => f.ok === null).length,
           max: ok.length ? Math.max(...ok.map((f) => f.error)) : 0 };
};
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const exp = (x) => {
  if (x === 0) return "0";
  const [m, e] = x.toExponential(1).split("e");
  return m.replace(".", ",") + " × 10" + String(+e).split("").map((c) => SUP[c]).join("");
};
const dl = B.deslizante, em = B.empotrado;
const leer = (f) => fs.readFileSync(path.join(AQUI, f), "utf8").trim().split("\n").map((l) => JSON.parse(l));
const CONV = leer("resultados/convergencia.jsonl").filter((r) => r.E_app_MPa !== undefined);
const BASE = leer("resultados/tapia_c1_c5.jsonl").filter((r) => r.clave === "referencia" && r.E_app_MPa !== undefined);
const SENS = leer("resultados/sensibilidad.jsonl").filter((r) => r.E_app_plato_MPa !== undefined);
// serie por (voi, estructura, n)
const MAL = {};
BASE.forEach((r) => (MAL[r.voi + "|" + r.estructura + "|40"] = r));
CONV.forEach((r) => (MAL[r.voi + "|" + r.estructura + "|" + r.n] = r));
const NS = [40, ...[...new Set(CONV.map((r) => r.n))].sort((a, b) => a - b)]
  .filter((n) => VOIS.every((v) => ESTR.every((e) => MAL[v + "|" + e + "|" + n])));
const mal = (v, e, n, k) => MAL[v + "|" + e + "|" + n][k];
const SEN = {};
SENS.forEach((r) => (SEN[r.voi + "|" + r.estructura] = r));
const cocientes = (num, den) => VOIS.map((v) => num(v) / den(v));
const rg = (a) => f0(Math.min(...a)) + " a " + f0(Math.max(...a));
const f0 = (x) => (100 * x).toFixed(0);
const SIGMA = dl.sigma_MPa, LIM_PISTOIA = 0.007 * 18000;

// ---------------------------------------------------------------- layouts
const PIE = (s) => s;
pres.defineSlideMaster({
  title: "PORTADA", background: { color: "1F3A44" },
  objects: [],
});
pres.defineSlideMaster({
  title: "CONTENIDO", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.35, w: 12.1, h: 0.95,
        fontFace: "Cambria", fontSize: 36, bold: true, color: "1F3A44", align: "left", valign: "middle", margin: 0 }, text: "" } },
    { text: { text: "spinpy V2.2.0 · protocolo de Tapia con C1 a C5", options: {
        x: 0.6, y: 7.02, w: 8, h: 0.3, fontFace: "Calibri", fontSize: 10, color: "5B6573", margin: 0 } } },
  ],
  slideNumber: { x: 12.2, y: 7.02, w: 0.55, h: 0.3, fontFace: "Calibri", fontSize: 10, color: "5B6573", align: "right" },
});
pres.defineSlideMaster({
  title: "CIERRE", background: { color: "1F3A44" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.5, w: 12.1, h: 0.95,
        fontFace: "Cambria", fontSize: 40, bold: true, color: "FFFFFF", align: "left", valign: "middle", margin: 0 }, text: "" } },
  ],
});

const txt = (s, text, o) => s.addText(text, Object.assign({ isTextBox: true, fontFace: "Calibri", color: "1B2A33", margin: 0, valign: "top" }, o));
const notas = (s, t) => s.addNotes(t);

// ---------------------------------------------------------------- 1 portada
pres.addSection({ title: "Apertura" });
{
  const s = pres.addSlide({ masterName: "PORTADA", sectionTitle: "Apertura" });
  txt(s, "Protocolo de Tapia con C1 a C5", { x: 0.8, y: 2.0, w: 11.7, h: 1.2, fontFace: "Cambria", fontSize: 48, bold: true, color: "FFFFFF", valign: "middle" });
  txt(s, "Deformación total y tensión de von Mises en cinco VOIs porcinos, sus spinodoides y sus dual-lattice", { x: 0.8, y: 3.35, w: 10.5, h: 1.1, fontSize: 24, color: "CFE3E1" });
  txt(s, "spinpy V2.2.0  ·  VOIs: Koria, Mengoni y Brockett (2020), CC BY 4.0, doi:10.5518/787", { x: 0.8, y: 6.5, w: 11.7, h: 0.4, fontSize: 14, color: "8FB8B5" });
  notas(s, "Presentación breve del protocolo de ensayo de Tapia et al. (2026) aplicado a cinco VOIs porcinos y a sus candidatos espinodal y dual-lattice. Se muestran la deformación total y la tensión de von Mises. Antes de los resultados se verifica la formulación con un cubo macizo.");
}

// ---------------------------------------------------------------- 2 experimento
pres.addSection({ title: "Experimento y método" });
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Experimento y método" });
  s.addText("El experimento", { placeholder: "title" });
  const stats = [["5", "VOIs porcinos\nC1 a C5"], ["15", "estructuras:\nVOI, spinodoide, dual-lattice"], ["100 N", "compresión axial\nsobre 3,0 mm de lado"]];
  stats.forEach((st, i) => {
    const x = 0.6 + i * 4.1;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "tarjeta " + i, x, y: 1.55, w: 3.8, h: 1.9, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.1 });
    txt(s, st[0], { x: x + 0.3, y: 1.7, w: 3.2, h: 0.95, fontFace: "Cambria", fontSize: 48, bold: true, color: i === 2 ? "C8553D" : "2E7D7A", valign: "middle" });
    txt(s, st[1], { x: x + 0.3, y: 2.65, w: 3.3, h: 0.7, fontSize: 16, color: "1B2A33" });
  });
  const pasos = ["VOI de micro-CT\n16 µm, cubo de 3,0 mm", "Ajuste del spinodoide\ny del dual-lattice", "Malla de ladrillos\n40³ elementos", "Compresión axial\nbase adherida, 100 N", "Campos |u|\ny von Mises"];
  pasos.forEach((p, i) => {
    const x = 0.6 + i * 2.5;
    s.addShape(pres.shapes.OVAL, { objectName: "num " + i, x: x, y: 4.1, w: 0.5, h: 0.5, fill: { color: "1F3A44" }, line: { color: "1F3A44" } });
    txt(s, String(i + 1), { x: x, y: 4.1, w: 0.5, h: 0.5, fontSize: 16, bold: true, color: "FFFFFF", align: "center", valign: "middle" });
    txt(s, p, { x: x, y: 4.75, w: 2.25, h: 0.9, fontSize: 15 });
    if (i < 4) s.addShape(pres.shapes.LINE, { objectName: "flecha " + i, x: x + 0.6, y: 4.35, w: 1.5, h: 0, line: { color: "5B6573", width: 1.5, endArrowType: "triangle" } });
  });
  txt(s, "Se replica el protocolo de ensayo de Tapia et al. (2026), no el estudio: los valores absolutos no son comparables con los publicados. E = 18 GPa, ν = 0,30, elasticidad lineal.", { x: 0.6, y: 6.1, w: 12.1, h: 0.7, fontSize: 14, color: "5B6573" });
  notas(s, "Los cinco VOIs son del conjunto público de Koria, Mengoni y Brockett (2020). Cada candidato se ajusta a su VOI con la función objetivo morfométrica de spinpy. La malla es de hexaedros a 40 elementos por lado, la resolución por omisión de la app.");
}

// ---------------------------------------------------------------- 3 objetivo
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Experimento y método" });
  s.addText("El objetivo: la ecuación", { placeholder: "title" });
  txt(s, "Resolver el equilibrio elástico del sólido", { x: 0.6, y: 1.45, w: 6.4, h: 0.5, fontSize: 20, bold: true, color: "2E7D7A" });
  const R_ = (t, o) => ({ text: t, options: o || {} });
  const eq = [[[R_("∇ · σ = 0")], "equilibrio en el sólido, sin fuerzas de volumen"],
              [[R_("σ = C : ε,   ε = ½ (∇u + ∇u"), R_("T", { superscript: true }), R_(")")], "ley de Hooke isótropa, pequeñas deformaciones"],
              [[R_("E = 18 GPa,   ν = 0,30")], "tejido homogéneo e isótropo"],
              [[R_("K u = f")], "problema discreto, hexaedros de 8 nodos, solver directo"]];
  eq.forEach((e, i) => {
    const y = 2.1 + i * 1.12;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "ecuacion " + i, x: 0.6, y, w: 6.5, h: 0.98, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, e[0], { x: 0.85, y: y + 0.08, w: 6.1, h: 0.45, fontFace: "Cambria", fontSize: 22 });
    txt(s, e[1], { x: 0.85, y: y + 0.55, w: 6.1, h: 0.35, fontSize: 14, color: "5B6573" });
  });
  txt(s, "Lo que se obtiene", { x: 7.7, y: 1.45, w: 5, h: 0.5, fontSize: 20, bold: true, color: "C8553D" });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "salida u", x: 7.7, y: 2.1, w: 5.0, h: 2.1, fill: { color: "FFFFFF" }, line: { color: "8FB8B5", width: 1.5 }, rectRadius: 0.08 });
  txt(s, "Deformación total", { x: 7.95, y: 2.2, w: 4.5, h: 0.4, fontSize: 18, bold: true });
  txt(s, "|u| = √(ux² + uy² + uz²)", { x: 7.95, y: 2.65, w: 4.5, h: 0.5, fontFace: "Cambria", fontSize: 20 });
  txt(s, "Es un desplazamiento, no una deformación unitaria. Su máximo es local; la rigidez global sale del desplazamiento medio de la cima.", { x: 7.95, y: 3.25, w: 4.5, h: 0.85, fontSize: 14, color: "5B6573" });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "salida vm", x: 7.7, y: 4.4, w: 5.0, h: 2.1, fill: { color: "FFFFFF" }, line: { color: "E8A18F", width: 1.5 }, rectRadius: 0.08 });
  txt(s, "Tensión de von Mises", { x: 7.95, y: 4.5, w: 4.5, h: 0.4, fontSize: 18, bold: true });
  txt(s, [{ text: "σ" }, { text: "vm", options: { subscript: true } }, { text: " = √(3/2 · s : s),   s = σ − ⅓ tr(σ) I" }], { x: 7.95, y: 4.95, w: 4.7, h: 0.5, fontFace: "Cambria", fontSize: 18 });
  txt(s, "Solo la parte desviadora. Localiza dónde se concentra la tensión.", { x: 7.95, y: 5.55, w: 4.5, h: 0.85, fontSize: 14, color: "5B6573" });
  notas(s, "El campo de von Mises y el de deformación total no se reescalan entre sí: la deformación total acumula lo que se movió la probeta, y su máximo se ubica en la capa superior en la mayoría de los casos; von Mises localiza los nudos trabeculares. La rigidez aparente no usa el máximo de |u|, sino el desplazamiento medio de la cima ponderado por área. En un estado uniaxial von Mises vale exactamente la tensión axial.");
}

// ---------------------------------------------------------------- 4 ejes
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Experimento y método" });
  s.addText("Método: ejes fijos y libres", { placeholder: "title" });
  // diagrama del cubo
  const x0 = 1.1, y0 = 2.9, w = 2.8, h = 2.6;
  s.addShape(pres.shapes.RECTANGLE, { objectName: "cubo", x: x0, y: y0, w, h, fill: { color: "EEF2F3" }, line: { color: "1F3A44", width: 2 } });
  txt(s, "VOI\n3,0 mm", { x: x0, y: y0 + 0.85, w, h: 0.9, fontSize: 18, bold: true, align: "center", color: "1F3A44" });
  s.addShape(pres.shapes.RECTANGLE, { objectName: "apoyo", x: x0 - 0.3, y: y0 + h, w: w + 0.6, h: 0.3, fill: { color: "5B6573" }, line: { color: "5B6573" } });
  for (let i = 0; i < 9; i++) s.addShape(pres.shapes.LINE, { objectName: "rayado " + i, x: x0 - 0.3 + i * 0.4, y: y0 + h + 0.3, w: -0.22, h: 0.26, line: { color: "5B6573", width: 1.5 } });
  txt(s, "Base: ux = uy = uz = 0", { x: x0 - 0.6, y: y0 + h + 0.65, w: w + 1.2, h: 0.4, fontSize: 14, bold: true, align: "center" });
  for (let i = 0; i < 5; i++) s.addShape(pres.shapes.LINE, { objectName: "carga " + i, x: x0 + 0.3 + i * 0.55, y: y0 - 0.85, w: 0, h: 0.75, line: { color: "C8553D", width: 2.5, endArrowType: "triangle" } });
  txt(s, "100 N en −z sobre las caras óseas", { x: x0 - 0.6, y: y0 - 1.35, w: w + 1.2, h: 0.4, fontSize: 14, bold: true, color: "C8553D", align: "center" });
  s.addShape(pres.shapes.LINE, { objectName: "libre izq", x: x0 - 0.85, y: y0 + 0.1, w: 0.7, h: 0, line: { color: "2E7D7A", width: 2, beginArrowType: "triangle", endArrowType: "triangle" } });
  s.addShape(pres.shapes.LINE, { objectName: "libre der", x: x0 + w + 0.15, y: y0 + 0.1, w: 0.7, h: 0, line: { color: "2E7D7A", width: 2, beginArrowType: "triangle", endArrowType: "triangle" } });
  txt(s, "x, y libres", { x: x0 + w + 0.15, y: y0 + 0.2, w: 1.5, h: 0.3, fontSize: 12, color: "2E7D7A", bold: true });
  // tabla
  const hd = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 14, align: "center" } });
  const fx = { text: "fijo", options: { color: "C8553D", bold: true, align: "center", fontSize: 14 } };
  const lb = { text: "libre", options: { color: "2E7D7A", bold: true, align: "center", fontSize: 14 } };
  const cell = (t, o) => ({ text: t, options: Object.assign({ fontSize: 14 }, o || {}) });
  s.addTable([
    [hd("Superficie"), hd("ux"), hd("uy"), hd("uz"), hd("Carga")],
    [cell("Base (z = 0)", { bold: true }), fx, fx, fx, cell("ninguna")],
    [cell("Cara superior", { bold: true }), lb, lb, lb, cell("100 N como presión uniforme sobre las caras óseas, sin plato")],
    [cell("Caras laterales", { bold: true }), lb, lb, lb, cell("ninguna")],
  ], { x: 5.4, y: 1.6, w: 7.3, colW: [1.9, 0.8, 0.8, 0.8, 3.0], rowH: 0.55, border: { type: "solid", color: "CFD6D9", pt: 1 }, fontFace: "Calibri", valign: "middle" });
  txt(s, [
    { text: "Base empotrada: adherida, con fricción infinita contra el plato. Es la condición del protocolo de Tapia aquí replicado.", options: { bullet: true, breakLine: true } },
    { text: "Carga: presión uniforme σ·A_bruta/A_ósea sobre las caras óseas de la cima, no un plato rígido. Con un plato la rigidez sube (ver sensibilidad).", options: { bullet: true, breakLine: true } },
    { text: "Contraste, apoyo deslizante: uz = 0 en la base y solo anclajes contra el movimiento rígido.", options: { bullet: true, breakLine: true } },
    { text: "Malla de 40³ hexaedros (h = 75 µm); NGSolve con solver directo. La app solo actúa como control del solver. Tapia usa otro elemento y 50 µm.", options: { bullet: true } },
  ], { x: 5.4, y: 4.0, w: 7.3, h: 2.9, fontSize: 14, paraSpaceAfter: 5 });
  notas(s, "Revisión de la formulación en el código: en el protocolo de Tapia la base queda con los tres ejes fijos. La cara superior no es un plato rígido: recibe una presión uniforme sobre las caras óseas (fuerza igual por cara de elemento) y sus ejes x e y quedan libres, igual que las caras laterales. El VOI de 188 elementos por lado se lleva a 40 por vecino más próximo y los candidatos sintéticos de 96 a 40. El apoyo deslizante fija solo uz en la base más los anclajes mínimos contra el movimiento de sólido rígido.");
}

// ---------------------------------------------------------------- 5 cubo
pres.addSection({ title: "Chequeo con el cubo macizo" });
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Chequeo con el cubo macizo" });
  s.addText("Chequeo con un cubo macizo", { placeholder: "title" });
  s.addImage({ path: path.join(AQUI, "figs/bloque.png"), objectName: "figura cubo", x: 0.6, y: 1.4, w: 5.5, h: 5.5 * 1277 / 1263 });
  const cajas = [
    ["Apoyo deslizante: exacto", `E aparente / E del tejido = 1,000000000002. Von Mises = σ = ${f2(dl.sigma_MPa)} MPa en todo el cubo. Desplazamientos exactos a 10⁻⁹.`, "2E7D7A"],
    ["Base empotrada: no es uniaxial", `E_app es ${f1(100 * (em.E_app_sobre_E_s - 1))} % mayor. Von Mises máximo = ${f2(em.vm_max_sobre_sigma)} σ en las esquinas de la base (a 40³; con 16³ era 1,27 σ: crece al refinar). El borde superior se expande ${pct(em.ux_borde_techo_mm / em.ux_borde_exacto_libre_mm)} % de lo libre.`, "C8553D"],
    ["Lectura", "La base fija impide que el material se dilate lateralmente: rigidiza y concentra tensión. La cara superior sí tiene x e y libres.", "1F3A44"],
  ];
  cajas.forEach((c, i) => {
    const y = 1.45 + i * 1.8;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "caja " + i, x: 6.5, y, w: 6.2, h: 1.65, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, c[0], { x: 6.75, y: y + 0.12, w: 5.7, h: 0.4, fontSize: 18, bold: true, color: c[2] });
    txt(s, c[1], { x: 6.75, y: y + 0.58, w: 5.7, h: 1.0, fontSize: 14 });
  });
  notas(s, "Cubo macizo de 3,008 mm con 100 N, E de 18 GPa, malla de 40 elementos por lado y NGSolve con solver directo. Con apoyo deslizante la solución es uniaxial y exacta. Con base empotrada no hay solución cerrada simple: se reportan la rigidización, el máximo de von Mises y la reducción de la expansión lateral. El máximo de von Mises de la base crece al refinar la malla (1,27 σ a 16 elementos por lado, 1,76 σ a 40), por eso no se cita.");
}

// ---------------------------------------------------------------- 6 sanidad
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Chequeo con el cubo macizo" });
  s.addText("Sanidad de todos los análisis", { placeholder: "title" });
  const dlz = sum((g) => g.indexOf("deslizante") >= 0 && g.indexOf("empotrado") < 0 && /hex8|tet10/.test(g));
  const emp = sum((g) => g.indexOf("empotrado") >= 0 && /hex8 empotrado|tet10 empotrado/.test(g) && g.indexOf("acuerdo") < 0);
  const acu = sum((g) => g.indexOf("acuerdo") >= 0);
  const hom = sum((g) => g.indexOf("homogeneizacion") >= 0);
  const tri = sum((g) => g.indexOf("triaxial") >= 0);
  const pri = sum((g) => g.indexOf("prisma") >= 0);
  const col = sum((g) => g.indexOf("columnas") >= 0);
  const H = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 13, align: "center" } });
  const c = (t, a) => ({ text: t, options: { fontSize: 13, align: a || "left" } });
  s.addTable([
    [H("Análisis (solución exacta o equilibrio)"), H("Pruebas"), H("Error máx. frente a lo exacto")],
    [c("Compresión lineal con apoyo deslizante: E_app, von Mises, |u|, equilibrio, plato rígido, E corregido, no lineal (SVK), Pistoia; app, NGSolve, scikit-fem; hex8 y TET10"), c(String(dlz.n), "center"), c(exp(dlz.max), "center")],
    [c("Compresión lineal con base empotrada: equilibrio, ejes fijos"), c(String(emp.n), "center"), c(exp(emp.max), "center")],
    [c("Acuerdo entre motores con base empotrada (E_app, von Mises, |u|)"), c(String(acu.n), "center"), c(exp(acu.max), "center")],
    [c("Homogeneización periódica: C₁₁, C₁₂, E, G, ν"), c(String(hom.n), "center"), c(exp(hom.max), "center")],
    [c("Compresión en X, Y, Z y convergencia con la malla (app)"), c(String(tri.n), "center"), c(exp(tri.max), "center")],
    [c("Prisma de lados distintos 12 × 16 × 24: área, altura y ejes correctos"), c(String(pri.n), "center"), c(exp(pri.max), "center")],
    [c("Columnas aisladas (E_app = φ·E del tejido) y malla 24 × 16 × 20 cargada en X, Y y Z; fragmentos flotantes"), c(String(col.n), "center"), c(exp(col.max), "center")],
  ], { x: 0.6, y: 1.4, w: 12.1, colW: [8.7, 1.4, 2.0], border: { type: "solid", color: "CFD6D9", pt: 1 }, fontFace: "Calibri", valign: "middle", rowH: [0.42, 0.75, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4] });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "resumen sanidad", x: 0.6, y: 5.3, w: 5.6, h: 1.55, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
  txt(s, `${e_ok} de ${e_ok + e_mal}`, { x: 0.85, y: 5.38, w: 5.1, h: 0.75, fontFace: "Cambria", fontSize: 40, bold: true, color: "2E7D7A", valign: "middle" });
  txt(s, `pruebas con respuesta exacta pasan. Otras ${e_info} son informativas (sin solución cerrada) y se reportan.`, { x: 0.85, y: 6.13, w: 5.1, h: 0.65, fontSize: 14 });
  txt(s, [
    { text: "Tolerancias declaradas antes de medir: 10⁻⁸ en lineal, 10⁻⁶ en no lineal.", options: { bullet: true, breakLine: true } },
    { text: "La primera versión tuvo 2 fallos: eran un defecto de mi prueba (tetraedros agrupados por centroide). Se acotó a hexaedros y queda documentado.", options: { bullet: true, breakLine: true } },
    { text: "Un cubo es simétrico y no detecta ejes, áreas o alturas confundidos: por eso se añadió un prisma.", options: { bullet: true, breakLine: true } },
    { text: "No aplican a un cubo: morfometría de superficie y fallo progresivo.", options: { bullet: true } },
  ], { x: 6.5, y: 5.3, w: 6.2, h: 1.65, fontSize: 12, paraSpaceAfter: 3 });
  notas(s, "El chequeo recorre cada análisis mecánico de la app con un cubo macizo, que tiene solución exacta: compresión lineal, plato rígido, módulo corregido del núcleo, no lineal de Saint Venant-Kirchhoff frente a su solución cerrada, criterio de Pistoia, homogeneización periódica, compresión en los tres ejes y convergencia con la malla. Los motores usados fueron la app, NGSolve y scikit-fem, en hexaedros y TET10. Guion: sanidad_cubo.py.");
}

// ---------------------------------------------------------------- 7 candidatos
pres.addSection({ title: "Resultados" });
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Resultados" });
  s.addText("Candidatos ajustados a cada VOI", { placeholder: "title" });
  const A = {};
  R.ajustes.forEach((a) => (A[a.voi + "|" + a.familia] = a));
  const H = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 14, align: "center" } });
  const rows = [[H("VOI"), H("BV/TV nativo del VOI"), H("BV/TV spinodoide"), H("BV/TV dual-lattice"), H("Error spinodoide"), H("Error dual-lattice")]];
  let gana = 0;
  VOIS.forEach((v) => {
    const a = A[v + "|spinodoide"], d = A[v + "|dual-lattice"];
    if (d.error < a.error) gana++;
    const best = (x, y) => ({ text: exp(x), options: { fontSize: 14, align: "center", bold: y } });
    rows.push([{ text: v, options: { bold: true, fontSize: 14, align: "center" } },
      { text: f3(a.BVTV_voi), options: { fontSize: 14, align: "center" } },
      { text: f3(a.BVTV_generada), options: { fontSize: 14, align: "center" } },
      { text: f3(d.BVTV_generada), options: { fontSize: 14, align: "center" } },
      best(a.error, a.error < d.error), best(d.error, d.error < a.error)]);
  });
  s.addTable(rows, { x: 0.6, y: 1.55, w: 12.1, colW: [1.3, 2.0, 2.2, 2.2, 2.2, 2.2], rowH: 0.52, border: { type: "solid", color: "CFD6D9", pt: 1 }, fontFace: "Calibri", valign: "middle" });
  const cards = [[`${gana} de 5`, "VOIs donde el dual-lattice ajusta con menor error morfométrico que el spinodoide"], [f2(100 * Math.max(...R.ajustes.map((a) => Math.abs(a.BVTV_generada - a.BVTV_voi)))) + " pp", "diferencia máxima de BV/TV entre un candidato y su VOI (puntos porcentuales)"], ["1", "semilla y un ajuste por VOI y familia: sin réplicas"]];
  cards.forEach((c, i) => {
    const x = 0.6 + i * 4.1;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "dato " + i, x, y: 5.0, w: 3.8, h: 1.75, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, c[0], { x: x + 0.25, y: 5.1, w: 3.3, h: 0.7, fontFace: "Cambria", fontSize: 34, bold: true, color: "2E7D7A", valign: "middle" });
    txt(s, c[1], { x: x + 0.25, y: 5.8, w: 3.3, h: 0.9, fontSize: 14 });
  });
  notas(s, "Los ajustes usan el modo completo, la resolución de comparación por omisión (96 elementos por lado, luego llevada a 40 para el ensayo), semilla 20260720, sin término mecánico ni réplicas. Los de C2 y C4 reproducen los de la campaña anterior. El error es la función objetivo morfométrica de spinpy.");
}

// ---------------------------------------------------------------- 8 y 9 mapas
const mapa = (titulo, img, lectura, nota) => {
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Resultados" });
  s.addText(titulo, { placeholder: "title" });
  s.addImage({ path: path.join(AQUI, img), objectName: "mapa", x: 0.5, y: 1.3, w: 9.0, h: 9.0 * 1474 / 2502 });
  lectura.forEach((l, i) => {
    const y = 1.45 + i * 1.75;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "lectura " + i, x: 9.75, y, w: 3.1, h: 1.6, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, l[0], { x: 9.95, y: y + 0.1, w: 2.8, h: 0.55, fontFace: "Cambria", fontSize: 22, bold: true, color: l[2] });
    txt(s, l[1], { x: 9.95, y: y + 0.68, w: 2.8, h: 0.9, fontSize: 13 });
  });
  txt(s, "Escala común por columna: de 0 al p99 de la estructura más exigida de cada VOI. Amarillo = zona más deformada o tensionada.", { x: 0.5, y: 6.65, w: 9.0, h: 0.3, fontSize: 11, color: "5B6573" });
  notas(s, nota);
};
mapa("Deformación total |u|", "figs/mapa_u.png", [
  [`${pct(rE_s[0])} a ${pct(rE_s[1])} %`, "de la rigidez del VOI conserva el spinodoide a 40³: se deforma más (ver convergencia)", "C8553D"],
  [`${pct(rE_d[0])} a ${pct(rE_d[1])} %`, "de la rigidez del VOI tiene el dual-lattice", "2E7D7A"],
  ["En la cima", `el máximo de |u| cae en la capa superior en ${R.tabla.filter((t) => t.capa_del_max_u === t.n_capas - 1).length} de las ${R.tabla.length} estructuras`, "1F3A44"],
], "Deformación total |u| en milímetros, mapa 3D de cada una de las 15 estructuras con la misma escala por VOI. La rigidez aparente se toma del desplazamiento medio del techo ponderado por área. Las cifras salen de NGSolve con solver directo; la app coincide a " + exp(R.dif_max_app_vs_ngsolve) + " en E_app.");
mapa("Tensión de von Mises", "figs/mapa_vm.png", [
  [`${f1(rV_s[0])} a ${f1(rV_s[1])}×`, "el p99 del VOI alcanza el spinodoide a 40³; sube al refinar", "C8553D"],
  [`${f1(rV_d[0])} a ${f1(rV_d[1])}×`, "el p99 del VOI alcanza el dual-lattice", "2E7D7A"],
  ["p99, no máximo", "el máximo no converge con la malla; el p99 es el valor citable", "1F3A44"],
], "Tensión de von Mises en MPa. El máximo de un campo de von Mises no converge con la malla (con el cubo macizo crece al refinar), por eso las tablas citan el percentil 99. Las zonas amarillas son nudos y puntales finos del spinodoide.");

// ---------------------------------------------------------------- 10 cifras
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Resultados" });
  s.addText("Rigidez y tensión por estructura", { placeholder: "title" });
  const ser = (k) => ESTR.map((e) => ({ name: NOM[e], labels: VOIS, values: VOIS.map((v) => Math.round(g(v, e, k) * 10) / 10) }));
  const cfg = (titulo) => ({
    barDir: "col", barGrouping: "clustered", chartColors: ["2B2F36", "C8553D", "2E7D7A"],
    showTitle: true, title: titulo, titleFontSize: 16, titleColor: "1B2A33", titleFontFace: "+mn-lt",
    showLegend: true, legendPos: "b", legendFontSize: 12, legendFontFace: "+mn-lt", legendColor: "1B2A33",
    showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 10, dataLabelFontFace: "+mn-lt", dataLabelColor: "1B2A33", dataLabelFormatCode: "0",
    catAxisLabelFontSize: 13, catAxisLabelFontFace: "+mn-lt", catAxisLabelColor: "1B2A33",
    valAxisLabelFontSize: 12, valAxisLabelFontFace: "+mn-lt", valAxisLabelColor: "5B6573",
    valGridLine: { color: "E1E6E8", size: 0.5 }, catGridLine: { style: "none" }, barGapWidthPct: 60,
  });
  s.addChart(pres.charts.BAR, ser("E_app_MPa"), Object.assign({ x: 0.5, y: 1.35, w: 6.2, h: 4.5 }, cfg("Módulo aparente E_app [MPa]")));
  s.addChart(pres.charts.BAR, ser("vm_p99_MPa"), Object.assign({ x: 6.75, y: 1.35, w: 6.2, h: 4.5 }, cfg("Von Mises, percentil 99 [MPa]")));
  txt(s, "Valores de NGSolve con solver directo, 40³, base adherida, 100 N. La app coincide en E_app hasta " + exp(R.dif_max_app_vs_ngsolve) + " (relativo). No están convergidos con la malla (diapositiva siguiente).", { x: 0.6, y: 6.05, w: 12.1, h: 0.75, fontSize: 14, color: "5B6573" });
  notas(s, "Con base empotrada, 100 N y E de 18 GPa. E_app se calcula con el desplazamiento medio del techo ponderado por área. El p99 de von Mises es sobre los elementos del sólido. C2 y C4 reproducen los valores de la campaña anterior de la comparativa de motores.");
}


// ---------------------------------------------------------------- convergencia
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Resultados" });
  s.addText("Convergencia con la malla", { placeholder: "title" });
  const H = (t, cs) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 13, align: "center", colspan: cs || 1 } });
  const r1 = [H("VOI"), H("E_app VOI [MPa]", NS.length)];
  const r2 = [H("")];
  NS.forEach(() => {});
  const filas = [[H("VOI"), H("Spinodoide / VOI, rigidez [%]", NS.length), H("Dual-lattice / VOI, rigidez [%]", NS.length), H("Spinodoide / VOI, p99 [×]", NS.length)],
                 [H("n por lado"), ...NS.map((n) => H(String(n))), ...NS.map((n) => H(String(n))), ...NS.map((n) => H(String(n)))]];
  VOIS.forEach((v) => {
    const c = (t, b) => ({ text: t, options: { fontSize: 13, align: "center", bold: b } });
    filas.push([c(v, true),
      ...NS.map((n) => c(f0(mal(v, "spinodoide", n, "E_app_MPa") / mal(v, "voi", n, "E_app_MPa")))),
      ...NS.map((n) => c(f0(mal(v, "dual-lattice", n, "E_app_MPa") / mal(v, "voi", n, "E_app_MPa")))),
      ...NS.map((n) => c(f1(mal(v, "spinodoide", n, "vm_p99_MPa") / mal(v, "voi", n, "vm_p99_MPa"))))]);
  });
  const nc = 1 + 3 * NS.length, w1 = 1.4, wc = (12.1 - w1) / (3 * NS.length);
  s.addTable(filas, { x: 0.6, y: 1.45, w: 12.1, colW: [w1, ...Array(3 * NS.length).fill(wc)], rowH: 0.42, border: { type: "solid", color: "CFD6D9", pt: 1 }, fontFace: "Calibri", valign: "middle" });
  const maxcam = (e) => Math.max(...VOIS.map((v) => Math.abs(mal(v, e, NS[NS.length - 1], "E_app_MPa") / mal(v, e, NS[0], "E_app_MPa") - 1)));
  const cards = [
    [f0(maxcam("voi")) + " %", `cambio máximo de E_app del VOI entre ${NS[0]}³ y ${NS[NS.length - 1]}³: converge razonablemente`],
    [f0(maxcam("spinodoide")) + " %", "cambio máximo en el spinodoide: sus puntales finos siguen resolviéndose, la malla no ha convergido"],
    ["Orden", "VOI > dual-lattice > spinodoide se conserva en cada malla (salvo C5, donde dual y spinodoide quedan cercanos)"],
  ];
  cards.forEach((c, i) => {
    const x = 0.6 + i * 4.1;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "conv " + i, x, y: 4.65, w: 3.8, h: 2.15, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, c[0], { x: x + 0.25, y: 4.75, w: 3.3, h: 0.7, fontFace: "Cambria", fontSize: 32, bold: true, color: i === 1 ? "C8553D" : "2E7D7A", valign: "middle" });
    txt(s, c[1], { x: x + 0.25, y: 5.5, w: 3.35, h: 1.25, fontSize: 13 });
  });
  notas(s, "NGSolve con solver directo, misma carga y misma condición de apoyo, remuestreando cada estructura a distintas mallas. Las razones entre estructuras dependen de la resolución; la que se cita a 40 elementos por lado es una primera aproximación, no un valor convergido. El p99 de von Mises del spinodoide sube con la malla porque los puntales finos quedan mejor resueltos.");
}

// ---------------------------------------------------------------- sensibilidad a la carga
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Resultados" });
  s.addText("Sensibilidad a la carga superior", { placeholder: "title" });
  const H = (t, cs) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 13, align: "center", colspan: cs || 1 } });
  const rel = (e, k) => cocientes((v) => SEN[v + "|" + e][k], (v) => SEN[v + "|voi"][k]);
  const aum = (e) => cocientes((v) => SEN[v + "|" + e].E_app_plato_MPa, (v) => SEN[v + "|" + e].E_app_traccion_MPa);
  const filas = [[H("Protocolo de carga en la cima"), H("Spinodoide / VOI [%]"), H("Dual-lattice / VOI [%]"), H("E_app VOI, cambio frente a presión [%]")]];
  const c = (t, b, a) => ({ text: t, options: { fontSize: 14, align: a || "center", bold: b } });
  const g3 = (k) => [c(rg(rel("spinodoide", k)) , false), c(rg(rel("dual-lattice", k)), false)];
  filas.push([c("Presión sobre las caras óseas (el aquí usado)", true, "left"), ...g3("E_app_traccion_MPa"), c("0")]);
  filas.push([c("Plato rígido (uz uniforme en la cima)", true, "left"), ...g3("E_app_plato_MPa"), c("+" + rg(aum("voi").map((x) => x - 1)))]);
  filas.push([c("Plato rígido, módulo del núcleo corregido", true, "left"), ...g3("E_corregido_MPa"), c("+" + rg(cocientes((v) => SEN[v + "|voi"].E_corregido_MPa, (v) => SEN[v + "|voi"].E_app_traccion_MPa).map((x) => x - 1)))]);
  s.addTable(filas, { x: 0.6, y: 1.5, w: 12.1, colW: [5.0, 2.2, 2.2, 2.7], rowH: [0.55, 0.5, 0.5, 0.5], border: { type: "solid", color: "CFD6D9", pt: 1 }, fontFace: "Calibri", valign: "middle" });
  const fr = R.tabla.map((t) => t.frac_cima), pl = fr.map((x) => SIGMA / x);
  txt(s, [
    { text: `Solo ${pct(Math.min(...fr))} a ${pct(Math.max(...fr))} % de la sección de la cima es hueso: la presión local sobre el hueso vale ${f0(Math.min(...pl) / 100)} a ${f0(Math.max(...pl) / 100)} veces ${f1(SIGMA)} MPa (${f1(Math.min(...pl))} a ${f1(Math.max(...pl))} MPa).`, options: { bullet: true, breakLine: true } },
    { text: "El resultado depende de cómo se aplique la carga: cambia el nivel absoluto y también la razón entre estructuras.", options: { bullet: true, breakLine: true } },
    { text: "El orden VOI > dual-lattice > spinodoide se mantiene en los tres esquemas, salvo C5 con el módulo corregido, donde dual-lattice y spinodoide empatan (79 %).", options: { bullet: true, breakLine: true } },
    { text: "El protocolo original de Tapia et al. (2026) no está replicado al detalle (elemento, malla y contacto), por lo que no se comparan cifras absolutas.", options: { bullet: true } },
  ], { x: 0.6, y: 3.85, w: 12.1, h: 3.0, fontSize: 16, paraSpaceAfter: 8 });
  notas(s, "Sensibilidad a 40 elementos por lado. Tres esquemas de carga en la cima: presión uniforme sobre las caras óseas (usada en todo el documento), plato rígido (uz uniforme en la cima) y plato rígido con el módulo del núcleo corregido. Las razones respecto del VOI varían, lo que muestra que la magnitud de las diferencias entre estructuras no es independiente del protocolo.");
}

// ---------------------------------------------------------------- 11 lectura y límites
pres.addSection({ title: "Cierre" });
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Cierre" });
  s.addText("Lectura y límites", { placeholder: "title" });
    const col = (x, titulo, color, items) => {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "columna " + titulo, x, y: 1.5, w: 5.95, h: 4.6, fill: { color: "EEF2F3" }, line: { color: "EEF2F3" }, rectRadius: 0.08 });
    txt(s, titulo, { x: x + 0.3, y: 1.65, w: 5.4, h: 0.5, fontSize: 22, bold: true, color });
    txt(s, items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1 } })), { x: x + 0.3, y: 2.3, w: 5.4, h: 3.7, fontSize: 16, paraSpaceAfter: 10 });
  };
  col(0.6, "Lo que se puede decir", "2E7D7A", [
    `A 40³ y con esta carga, el spinodoide conserva ${pct(rE_s[0])} a ${pct(rE_s[1])} % de la rigidez del VOI y su p99 de von Mises es ${f1(rV_s[0])} a ${f1(rV_s[1])} veces mayor.`,
    `El dual-lattice conserva ${pct(rE_d[0])} a ${pct(rE_d[1])} % y ${f1(rV_d[0])} a ${f1(rV_d[1])} veces.`,
    `La formulación se verificó con cubo, prisma y columnas macizos: ${e_ok} pruebas con respuesta exacta pasan.`,
  ]);
  col(6.75, "Límites", "C8553D", [
    "Los resultados a 40³ no están convergidos con la malla, sobre todo los spinodoides; las razones dependen de la resolución y de la carga superior.",
    `El p99 supera ${f0(LIM_PISTOIA / 100)} MPa (0,7 % de deformación crítica de Pistoia con E = 18 GPa) en ${R.tabla.filter((t) => t.vm_p99_MPa > LIM_PISTOIA).length} de las ${R.tabla.length} estructuras: la hipótesis lineal es discutible.`,
    "Remuestreo desigual (VOI 188→40, candidatos 96→40), una semilla por ajuste y cinco VOIs de un mismo estudio.",
    "Es simulación: no hay contraste con un ensayo físico de estos huesos.",
  ]);
  notas(s, "El resultado son comparaciones relativas entre estructuras bajo una misma condición de carga, no predicciones de resistencia de un espécimen. El caso C3 spinodoide es el más blando y es el que más se aparta de la linealidad.");
}

// ---------------------------------------------------------------- 12 cierre
{
  const s = pres.addSlide({ masterName: "CIERRE", sectionTitle: "Cierre" });
  s.addText("Conclusiones", { placeholder: "title" });
  const it = [
    ["1", "La formulación es la esperada: base adherida con los tres ejes fijos, cima libre en x e y y presión sobre las caras óseas. El cubo macizo es exacto con apoyo deslizante y rigidiza " + f1(100 * (em.E_app_sobre_E_s - 1)) + " % con la base adherida."],
    ["2", "En estos 5 VOIs, a 40³ y con una semilla, el dual-lattice queda más cerca del VOI que el spinodoide en rigidez y en tensión (" + pct(rE_d[0]) + " a " + pct(rE_d[1]) + " % frente a " + pct(rE_s[0]) + " a " + pct(rE_s[1]) + " %). Las magnitudes dependen de la malla y de la carga."],
    ["3", "Falta convergencia con la malla y contraste con ensayos mecánicos del mismo hueso. Hasta entonces, son comparaciones relativas."],
  ];
  it.forEach((t, i) => {
    const y = 1.9 + i * 1.5;
    s.addShape(pres.shapes.OVAL, { objectName: "numero " + i, x: 0.8, y, w: 0.8, h: 0.8, fill: { color: "2E7D7A" }, line: { color: "2E7D7A" } });
    txt(s, t[0], { x: 0.8, y, w: 0.8, h: 0.8, fontFace: "Cambria", fontSize: 26, bold: true, color: "FFFFFF", align: "center", valign: "middle" });
    txt(s, t[1], { x: 1.9, y: y - 0.05, w: 10.6, h: 1.2, fontSize: 20, color: "FFFFFF" });
  });
  notas(s, "Próximos pasos: obtener los datos mecánicos del conjunto de Koria et al. si están publicados, o ensayar piezas propias, y repetir el contraste con spinpy-validar-exp con criterios declarados antes de mirar los datos.");
}

// ---------------------------------------------------------------- 13 anexo
{
  const s = pres.addSlide({ masterName: "CONTENIDO", sectionTitle: "Cierre" });
  s.addText("Anexo: cifras por estructura", { placeholder: "title" });
  const H = (t, cs) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: "1F3A44" }, fontSize: 12, align: "center", colspan: cs || 1 } });
  const rows = [[H("VOI"), H("Estructura"), H("BV/TV"), H("E_app [MPa]"), H("von Mises p99 [MPa]"), H("von Mises máx. (no citable) [MPa]"), H("|u| máx. [mm]"), H("|u| p99 [mm]")]];
  VOIS.forEach((v) => ESTR.forEach((e, i) => {
    const o = { fontSize: 12, align: "center", fill: { color: i === 0 ? "EEF2F3" : "FFFFFF" } };
    rows.push([{ text: i === 0 ? v : "", options: Object.assign({ bold: true }, o) }, { text: NOM[e], options: Object.assign({}, o, { align: "left" }) },
      { text: f3(g(v, e, "BVTV")), options: o }, { text: f1(g(v, e, "E_app_MPa")), options: o }, { text: f1(g(v, e, "vm_p99_MPa")), options: o },
      { text: f1(g(v, e, "vm_max_MPa")), options: o }, { text: f4(g(v, e, "u_tot_max_mm")), options: o }, { text: f4(g(v, e, "u_p99_mm")), options: o }]);
  }));
  s.addTable(rows, { x: 0.6, y: 1.3, w: 12.1, colW: [0.8, 1.8, 1.1, 1.6, 1.9, 2.3, 1.4, 1.2], rowH: 0.3, border: { type: "solid", color: "CFD6D9", pt: 0.5 }, fontFace: "Calibri", valign: "middle" });
  notas(s, "Cifras de NGSolve con solver directo a 40³. El von Mises máximo se da a título informativo: no converge con la malla y no debe citarse.");
}

(async () => {
  const out = path.join(AQUI, "Tapia_C1_C5.pptx");
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("escrito", out);
})();
