#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
md2pdf.py — Compositor del anexo de validacion (y del informe para publicacion).

Vive en el paquete desde que el visor compone PDF: `docs/` no viaja ni en el
ejecutable ni en la copia de GitHub. `docs/md2pdf.py` queda como puente, asi
que el manual, el anexo y `paper2pdf.py` siguen llamandolo igual. Requiere
reportlab (extra `informe` de pyproject.toml); por eso `spinpy` no lo importa
al cargar, solo `informe.componer` cuando hace falta.

Convierte el informe en Markdown a un PDF con formato de anexo metodologico:
portada, indice, encabezados numerados, tablas con cabecera sombreada, figuras
con pie, bloques de codigo monoespaciados y —elemento distintivo de este
documento— las explicaciones en lenguaje sencillo destacadas en un recuadro
tintado, para que el lector no especialista las localice de un vistazo.

Subconjunto de Markdown soportado (deliberadamente restringido, para que la
composicion sea predecible):
    #, ##, ###, ####      encabezados
    parrafos              texto corrido
    -, *                  listas con vinetas
    1.                    listas numeradas
    | a | b |             tablas GFM (con fila separadora ---)
    ```                   bloques de codigo
    ```math               ecuacion numerada (LaTeX del subconjunto mathtext de
                          matplotlib; una linea `%eq N` fija el numero)
    >                     cita -> recuadro "En palabras sencillas"
    ![pie](ruta)          figura
    ---                   linea horizontal
    **negrita**, *cursiva*, `codigo`
    X_{sub}, X^{sup}      subindice y superindice en el texto corrido

Uso:
    python md2pdf.py entrada.md salida.pdf
"""

import os
import re
import sys
from xml.sax.saxutils import escape as xml_escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, HRFlowable, Image,
                                KeepTogether, ListFlowable, ListItem,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents

# --------------------------------------------------------------------------
# Paleta. Sobria, con un unico color de acento; el recuadro de lenguaje
# sencillo usa un tinte calido para distinguirse del cuerpo tecnico sin
# competir con el.
# --------------------------------------------------------------------------
AZUL      = colors.HexColor('#1F3864')
AZUL_SUAVE= colors.HexColor('#2E5C9A')
GRIS_LIN  = colors.HexColor('#B8BFC9')
GRIS_CAB  = colors.HexColor('#E4E8EE')
GRIS_ALT  = colors.HexColor('#F5F7FA')
CAJA_BG   = colors.HexColor('#FBF4E4')
CAJA_BORDE= colors.HexColor('#C9A227')
COD_BG    = colors.HexColor('#F2F4F7')

ANCHO_UTIL = A4[0] - 40 * mm

# Rotulo del pie de figura. El compositor nacio para documentos en espanol;
# `componer(..., etiqueta_figura="Figure")` lo cambia para un informe en ingles.
ETIQUETA_FIGURA = 'Figura'
# Rotulo del recuadro de lenguaje sencillo, por la misma razon.
TITULO_LLANO = 'EN PALABRAS SENCILLAS'


# ==========================================================================
# Estilos
# ==========================================================================
def construir_estilos():
    ss = getSampleStyleSheet()
    S = {}

    S['cuerpo'] = ParagraphStyle(
        'cuerpo', parent=ss['BodyText'], fontName='Times-Roman', fontSize=10,
        leading=14.2, alignment=TA_JUSTIFY, spaceBefore=0, spaceAfter=6)

    S['h1'] = ParagraphStyle(
        'h1', parent=ss['Heading1'], fontName='Helvetica-Bold', fontSize=17,
        leading=21, textColor=AZUL, spaceBefore=18, spaceAfter=9,
        keepWithNext=1)
    S['h2'] = ParagraphStyle(
        'h2', parent=ss['Heading2'], fontName='Helvetica-Bold', fontSize=13,
        leading=16.5, textColor=AZUL_SUAVE, spaceBefore=14, spaceAfter=6,
        keepWithNext=1)
    S['h3'] = ParagraphStyle(
        'h3', parent=ss['Heading3'], fontName='Helvetica-Bold', fontSize=11,
        leading=14, textColor=colors.HexColor('#333B48'), spaceBefore=11,
        spaceAfter=4, keepWithNext=1)
    S['h4'] = ParagraphStyle(
        'h4', parent=ss['Heading4'], fontName='Helvetica-BoldOblique',
        fontSize=10, leading=13, textColor=colors.HexColor('#4A5262'),
        spaceBefore=9, spaceAfter=3, keepWithNext=1)

    S['codigo'] = ParagraphStyle(
        'codigo', fontName='Courier', fontSize=8.2, leading=10.6,
        textColor=colors.HexColor('#1A1A1A'), spaceBefore=0, spaceAfter=0)

    S['tabla'] = ParagraphStyle(
        'tabla', fontName='Helvetica', fontSize=7.9, leading=9.8,
        alignment=TA_JUSTIFY)
    S['tabla_cab'] = ParagraphStyle(
        'tabla_cab', fontName='Helvetica-Bold', fontSize=7.9, leading=9.8,
        textColor=AZUL)

    S['llano'] = ParagraphStyle(
        'llano', parent=S['cuerpo'], fontName='Times-Roman', fontSize=9.8,
        leading=14, spaceAfter=5, textColor=colors.HexColor('#3B3222'))
    S['llano_tit'] = ParagraphStyle(
        'llano_tit', fontName='Helvetica-Bold', fontSize=9.2, leading=12,
        textColor=colors.HexColor('#8A6D1F'), spaceAfter=4)

    S['pie_fig'] = ParagraphStyle(
        'pie_fig', fontName='Helvetica-Oblique', fontSize=8.3, leading=10.5,
        alignment=TA_CENTER, textColor=colors.HexColor('#4A5262'),
        spaceBefore=4, spaceAfter=10)

    S['portada_tit'] = ParagraphStyle(
        'portada_tit', fontName='Helvetica-Bold', fontSize=25, leading=30,
        alignment=TA_CENTER, textColor=AZUL, spaceAfter=10)
    S['portada_sub'] = ParagraphStyle(
        'portada_sub', fontName='Helvetica', fontSize=13.5, leading=18,
        alignment=TA_CENTER, textColor=AZUL_SUAVE, spaceAfter=8)
    S['portada_txt'] = ParagraphStyle(
        'portada_txt', fontName='Times-Roman', fontSize=10.5, leading=15,
        alignment=TA_CENTER, textColor=colors.HexColor('#333B48'))

    S['toc1'] = ParagraphStyle('toc1', fontName='Helvetica-Bold', fontSize=10.5,
                               leading=16, textColor=AZUL, spaceBefore=5)
    S['toc2'] = ParagraphStyle('toc2', fontName='Helvetica', fontSize=9.5,
                               leading=13.5, leftIndent=14)
    S['toc3'] = ParagraphStyle('toc3', fontName='Helvetica', fontSize=8.8,
                               leading=12, leftIndent=30,
                               textColor=colors.HexColor('#4A5262'))
    return S


# ==========================================================================
# Saneado de caracteres no representables
# ==========================================================================
# Las fuentes base-14 (Times, Helvetica, Courier) usan WinAnsiEncoding.
# ReportLab resuelve las griegas y los operadores matematicos sustituyendolos
# por sus glifos de la fuente Symbol, de modo que 'theta', 'raiz', 'flecha' o
# 'menos' se componen bien. Lo que NO existe en ninguna base-14 son los
# SUPERINDICES y SUBINDICES Unicode: se componian como un recuadro negro.
# Se convierten aqui a marcado <super>/<sub>, que ReportLab compone con la
# fuente en curso y con el tamano correcto.
SUPERINDICES = {
    '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4',
    '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
    '⁺': '+', '⁻': '-', '⁼': '=', '⁽': '(', '⁾': ')',
    'ⁿ': 'n', 'ᵀ': 'T', 'ᵃ': 'a', 'ᵇ': 'b', 'ᶜ': 'c',
}
SUBINDICES = {
    '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
    '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    '₊': '+', '₋': '-', 'ᵢ': 'i', 'ⱼ': 'j', 'ₓ': 'x',
    'ₖ': 'k', 'ₙ': 'n',
}
# Restos sin glifo ni en WinAnsi ni en Symbol.
OTROS = {
    '∥': '||',   # paralelo
    'ł': 'l',    # l polaca (Klosowski)
}

_RE_SUP = re.compile('[' + ''.join(SUPERINDICES) + ']+')
_RE_SUB = re.compile('[' + ''.join(SUBINDICES) + ']+')


def saneado(txt):
    """Sustituye los caracteres que las fuentes base-14 no pueden componer.

    Se aplica DESPUES de xml_escape, de modo que las etiquetas insertadas
    lleguen intactas a ReportLab. Las secuencias contiguas de superindices se
    agrupan en una sola etiqueta: '10⁻¹⁴' -> '10<super>-14</super>',
    y no en una etiqueta por caracter, para que la linea base sea uniforme.
    """
    txt = _RE_SUP.sub(
        lambda m: '<super>%s</super>' % ''.join(SUPERINDICES[c] for c in m.group()),
        txt)
    txt = _RE_SUB.sub(
        lambda m: '<sub>%s</sub>' % ''.join(SUBINDICES[c] for c in m.group()),
        txt)
    for viejo, nuevo in OTROS.items():
        txt = txt.replace(viejo, nuevo)
    return txt


# ==========================================================================
# Texto en linea
# ==========================================================================
def inline(txt):
    """Convierte marcado en linea de Markdown a etiquetas de ReportLab."""
    txt = saneado(xml_escape(txt))
    # codigo en linea primero: su contenido no debe reinterpretarse
    piezas = []
    for i, trozo in enumerate(re.split(r'`([^`]+)`', txt)):
        if i % 2 == 1:
            piezas.append('<font face="Courier" size="8.6" '
                          'backColor="#F2F4F7">%s</font>' % trozo)
        else:
            # Subindices y superindices explicitos con llaves (E_{s}, 10^{-6}):
            # solo con llaves, para no tocar los guiones bajos de nombres de
            # archivo o de claves que aparecen en el texto.
            trozo = re.sub(r'_\{([^{}]+)\}', r'<sub>\1</sub>', trozo)
            trozo = re.sub(r'\^\{([^{}]+)\}', r'<super>\1</super>', trozo)
            trozo = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', trozo)
            trozo = re.sub(r'(?<![\w*])\*([^*\n]+?)\*(?![\w*])', r'<i>\1</i>',
                           trozo)
            trozo = re.sub(r'\[([^\]]+)\]\(([^)]+)\)',
                           r'<link href="\2" color="#2E5C9A">\1</link>', trozo)
            piezas.append(trozo)
    return ''.join(piezas)


# ==========================================================================
# Bloques
# ==========================================================================
def tabla_flowable(filas, S):
    """Construye una tabla GFM con anchos proporcionales al contenido."""
    if not filas:
        return None
    ncol = max(len(f) for f in filas)
    filas = [f + [''] * (ncol - len(f)) for f in filas]

    datos = []
    for r, fila in enumerate(filas):
        est = S['tabla_cab'] if r == 0 else S['tabla']
        datos.append([Paragraph(inline(c), est) for c in fila])

    # Ancho por columna proporcional a la longitud tipica de su contenido,
    # acotado para que ninguna columna colapse ni acapare la pagina.
    pesos = []
    for c in range(ncol):
        largo = max(len(filas[r][c]) for r in range(len(filas)))
        pesos.append(min(max(largo, 6), 46))
    total = float(sum(pesos))
    anchos = [ANCHO_UTIL * p / total for p in pesos]

    t = Table(datos, colWidths=anchos, repeatRows=1, hAlign='LEFT')
    estilo = [
        ('BACKGROUND', (0, 0), (-1, 0), GRIS_CAB),
        ('LINEBELOW', (0, 0), (-1, 0), 0.9, AZUL),
        ('LINEABOVE', (0, 0), (-1, 0), 0.9, AZUL),
        ('LINEBELOW', (0, -1), (-1, -1), 0.9, AZUL),
        ('INNERGRID', (0, 1), (-1, -1), 0.25, GRIS_LIN),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 3.2),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3.2),
    ]
    for r in range(2, len(datos), 2):
        estilo.append(('BACKGROUND', (0, r), (-1, r), GRIS_ALT))
    t.setStyle(TableStyle(estilo))
    return t


def caja_llano(parrafos, S):
    """Recuadro tintado para las explicaciones en lenguaje sencillo."""
    interior = [Paragraph(TITULO_LLANO, S['llano_tit'])]
    primero = True
    for tipo, dato in parrafos:
        if tipo == 'p':
            if primero:
                # El recuadro ya lleva su propio rotulo: si el texto repite la
                # formula de apertura, se suprime para no leerla dos veces.
                dato = re.sub(r'^\*\*(En palabras sencillas|In plain words)'
                              r'\.?\*\*\s*', '', dato)
                primero = False
            if not dato.strip():
                continue
            interior.append(Paragraph(inline(dato), S['llano']))
        elif tipo == 'ul':
            interior.append(ListFlowable(
                [ListItem(Paragraph(inline(x), S['llano']), leftIndent=12)
                 for x in dato],
                bulletType='bullet', start='•', leftIndent=13,
                bulletFontSize=7, spaceBefore=1, spaceAfter=3))
        elif tipo == 'ol':
            interior.append(ListFlowable(
                [ListItem(Paragraph(inline(x), S['llano']), leftIndent=12)
                 for x in dato],
                bulletType='1', leftIndent=15, spaceBefore=1, spaceAfter=3))
    t = Table([[interior]], colWidths=[ANCHO_UTIL], hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), CAJA_BG),
        ('LINEBEFORE', (0, 0), (0, -1), 2.6, CAJA_BORDE),
        ('BOX', (0, 0), (-1, -1), 0.4, colors.HexColor('#E0D2A8')),
        ('LEFTPADDING', (0, 0), (-1, -1), 10),
        ('RIGHTPADDING', (0, 0), (-1, -1), 9),
        ('TOPPADDING', (0, 0), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))
    return t


def caja_codigo(lineas, S):
    interior = [Paragraph(saneado(xml_escape(l).replace(' ', '&nbsp;')) or '&nbsp;',
                          S['codigo']) for l in lineas]
    t = Table([[interior]], colWidths=[ANCHO_UTIL], hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), COD_BG),
        ('BOX', (0, 0), (-1, -1), 0.4, GRIS_LIN),
        ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ('RIGHTPADDING', (0, 0), (-1, -1), 6),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
    ]))
    return t


# Por encima de este peso una figura se reduce ANTES de incrustarla en el PDF.
# Los renders 3D y sus paneles salen a 2400-6500 px y varios MB en PNG, que es
# lo que pide una revista; en un PDF de lectura, a 168 mm de ancho, 2400 px
# son ~360 ppp y no se distinguen. El archivo de `figuras/` no se toca.
PESO_REDUCIR = 1500000
LADO_PDF = 2400


def _para_pdf(ruta):
    """La ruta tal cual, o una version reducida en memoria si pesa mucho."""
    try:
        if os.path.getsize(ruta) <= PESO_REDUCIR:
            return ruta
        import io
        from PIL import Image as PILImage
        im = PILImage.open(ruta)
        if im.mode in ('RGBA', 'LA', 'P'):
            fondo = PILImage.new('RGB', im.size, (255, 255, 255))
            im = im.convert('RGBA')
            fondo.paste(im, mask=im.split()[-1])
            im = fondo
        else:
            im = im.convert('RGB')
        if max(im.size) > LADO_PDF:
            f = LADO_PDF / float(max(im.size))
            im = im.resize((max(1, int(im.size[0] * f)),
                            max(1, int(im.size[1] * f))), PILImage.LANCZOS)
        buf = io.BytesIO()
        im.save(buf, format='JPEG', quality=90, optimize=True)
        buf.seek(0)
        return buf
    except Exception:
        return ruta


# Ecuaciones: cuerpo del PDF en Times, asi que STIX, que es la Times de las
# matematicas. A 300 ppp una ecuacion de 10 pt se lee nitida al imprimir.
PUNTOS_ECUACION = 11
DPI_ECUACION = 300
_RE_NUM_EC = re.compile(r'^\s*%\s*eq\s+(\S+)\s*$')


def texto_ecuacion(lineas):
    """(tex, numero o None) de las lineas de un bloque ```math."""
    numero, cuerpo = None, []
    for l in lineas:
        m = _RE_NUM_EC.match(l)
        if m:
            numero = m.group(1)
        elif l.strip() and not l.strip().startswith('%'):
            cuerpo.append(l.strip())
    return ' '.join(cuerpo), numero


def imagen_ecuacion(tex):
    """PNG en memoria de una ecuacion y su tamano en puntos, o None.

    Sin matplotlib, o si el LaTeX cae fuera del subconjunto de mathtext, se
    devuelve None y el bloque sale como codigo: mejor el fuente legible que un
    hueco en el documento.
    """
    try:
        import io
        import matplotlib
        from matplotlib import mathtext
        from matplotlib.font_manager import FontProperties
        buf = io.BytesIO()
        with matplotlib.rc_context({'mathtext.fontset': 'stix'}):
            mathtext.math_to_image('$' + tex + '$', buf,
                                   prop=FontProperties(size=PUNTOS_ECUACION),
                                   dpi=DPI_ECUACION, format='png')
        buf.seek(0)
        from reportlab.lib.utils import ImageReader
        w, h = ImageReader(buf).getSize()
        buf.seek(0)
        return buf, w * 72.0 / DPI_ECUACION, h * 72.0 / DPI_ECUACION
    except Exception:
        return None


def ecuacion(lineas, S, numero_por_orden):
    tex, numero = texto_ecuacion(lineas)
    img = imagen_ecuacion(tex)
    if img is None:
        return caja_codigo(lineas, S)
    buf, w, h = img
    lado = 22 * mm
    # 4 mm de holgura a cada lado: una ecuacion ancha no toca su numero.
    ancho = min(w, ANCHO_UTIL - 2 * lado - 8 * mm)
    alto = h * ancho / w
    etiqueta = '(%s)' % (numero if numero is not None else numero_por_orden)
    t = Table([[' ', Image(buf, width=ancho, height=alto),
                Paragraph(etiqueta, S['cuerpo'])]],
              colWidths=[lado, ANCHO_UTIL - 2 * lado, lado])
    t.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('ALIGN', (1, 0), (1, 0), 'CENTER'),
        ('ALIGN', (2, 0), (2, 0), 'RIGHT'),
        ('LEFTPADDING', (0, 0), (-1, -1), 0),
        ('RIGHTPADDING', (0, 0), (-1, -1), 0),
        ('TOPPADDING', (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
    ]))
    return t


def figura(ruta, pie, base, S, contador):
    if not os.path.isabs(ruta):
        ruta = os.path.normpath(os.path.join(base, ruta))
    if not os.path.exists(ruta):
        return [Paragraph('<i>[Figura no encontrada: %s]</i>'
                          % xml_escape(ruta), S['pie_fig'])]
    from reportlab.lib.utils import ImageReader
    iw, ih = ImageReader(ruta).getSize()
    ancho = min(ANCHO_UTIL, 168 * mm)
    alto = ancho * ih / float(iw)
    tope = 205 * mm
    if alto > tope:
        alto = tope
        ancho = alto * iw / float(ih)
    partes = [Image(_para_pdf(ruta), width=ancho, height=alto,
                    hAlign='CENTER')]
    if pie:
        partes.append(Paragraph('<b>%s %d.</b> %s' % (ETIQUETA_FIGURA, contador,
                                                     inline(pie)),
                                S['pie_fig']))
    return [KeepTogether(partes)]


# ==========================================================================
# Analizador
# ==========================================================================
def parsear(md, base, S):
    lineas = md.split('\n')
    flow = []
    i = 0
    n_fig = 0
    n_h1 = 0
    n_ec = 0

    def leer_cita(j):
        """Recoge un bloque '>' completo y lo estructura."""
        bloques, buf, lista, tipo_lista = [], [], [], None
        while j < len(lineas) and (lineas[j].startswith('>') or
                                   (lineas[j].strip() == '' and
                                    j + 1 < len(lineas) and
                                    lineas[j + 1].startswith('>'))):
            cont = re.sub(r'^>\s?', '', lineas[j])
            s = cont.strip()
            if not s:
                if buf:
                    bloques.append(('p', ' '.join(buf))); buf = []
                if lista:
                    bloques.append((tipo_lista, lista)); lista = []
                j += 1
                continue
            m_ul = re.match(r'^[-*]\s+(.*)', s)
            m_ol = re.match(r'^\d+\.\s+(.*)', s)
            if m_ul or m_ol:
                if buf:
                    bloques.append(('p', ' '.join(buf))); buf = []
                nuevo = 'ul' if m_ul else 'ol'
                if tipo_lista and tipo_lista != nuevo and lista:
                    bloques.append((tipo_lista, lista)); lista = []
                tipo_lista = nuevo
                lista.append((m_ul or m_ol).group(1))
            elif lista and re.match(r'^\s{2,}\S', cont):
                lista[-1] += ' ' + s          # continuacion de item
            else:
                if lista:
                    bloques.append((tipo_lista, lista)); lista = []
                buf.append(s)
            j += 1
        if buf:
            bloques.append(('p', ' '.join(buf)))
        if lista:
            bloques.append((tipo_lista, lista))
        return bloques, j

    while i < len(lineas):
        ln = lineas[i]
        s = ln.strip()

        if not s:
            i += 1
            continue

        if s.startswith('<!--'):
            while i < len(lineas) and '-->' not in lineas[i]:
                i += 1
            i += 1
            continue

        # --- ecuacion ---
        if s.startswith('```math'):
            i += 1
            buf = []
            while i < len(lineas) and not lineas[i].strip().startswith('```'):
                buf.append(lineas[i])
                i += 1
            i += 1
            n_ec += 1
            flow.append(ecuacion(buf, S, n_ec))
            continue

        # --- bloque de codigo ---
        if s.startswith('```'):
            i += 1
            buf = []
            while i < len(lineas) and not lineas[i].strip().startswith('```'):
                buf.append(lineas[i])
                i += 1
            i += 1
            flow.append(Spacer(1, 3))
            flow.append(caja_codigo(buf, S))
            flow.append(Spacer(1, 7))
            continue

        # --- cita -> recuadro de lenguaje sencillo ---
        if s.startswith('>'):
            bloques, i = leer_cita(i)
            flow.append(Spacer(1, 4))
            flow.append(caja_llano(bloques, S))
            flow.append(Spacer(1, 9))
            continue

        # --- encabezados ---
        m = re.match(r'^(#{1,4})\s+(.*)', s)
        if m:
            niv = len(m.group(1))
            txt = m.group(2).strip()
            if niv == 1:
                n_h1 += 1
                if n_h1 > 1:
                    flow.append(PageBreak())
            p = Paragraph(inline(txt), S['h%d' % niv])
            # Anclaje para el indice (niveles 1-3)
            if niv <= 3:
                clave = 'sec%d' % len(flow)
                p._bookmark = clave
                flow.append(EntradaIndice(p, niv - 1, txt, clave, S))
            else:
                flow.append(p)
            i += 1
            continue

        # --- linea horizontal ---
        if re.match(r'^(-{3,}|\*{3,}|_{3,})$', s):
            flow.append(Spacer(1, 5))
            flow.append(HRFlowable(width='100%', thickness=0.7, color=GRIS_LIN))
            flow.append(Spacer(1, 7))
            i += 1
            continue

        # --- figura ---
        m = re.match(r'^!\[([^\]]*)\]\(([^)]+)\)\s*$', s)
        if m:
            n_fig += 1
            flow.append(Spacer(1, 5))
            flow.extend(figura(m.group(2), m.group(1), base, S, n_fig))
            i += 1
            continue

        # --- tabla ---
        if s.startswith('|') and i + 1 < len(lineas) and \
                re.match(r'^\s*\|[\s:|-]+\|\s*$', lineas[i + 1]):
            filas = []
            while i < len(lineas) and lineas[i].strip().startswith('|'):
                cruda = lineas[i].strip()
                if not re.match(r'^\|[\s:|-]+\|$', cruda):
                    celdas = [c.strip() for c in cruda.strip('|').split('|')]
                    filas.append(celdas)
                i += 1
            flow.append(Spacer(1, 3))
            t = tabla_flowable(filas, S)
            if t:
                flow.append(t)
            flow.append(Spacer(1, 9))
            continue

        # --- listas ---
        m_ul = re.match(r'^[-*]\s+(.*)', s)
        m_ol = re.match(r'^\d+\.\s+(.*)', s)
        if m_ul or m_ol:
            es_ol = m_ol is not None
            items = []
            while i < len(lineas):
                cur = lineas[i]
                cs = cur.strip()
                mm_ = (re.match(r'^\d+\.\s+(.*)', cs) if es_ol
                       else re.match(r'^[-*]\s+(.*)', cs))
                if mm_:
                    items.append(mm_.group(1))
                    i += 1
                elif items and cs and re.match(r'^\s{2,}\S', cur):
                    items[-1] += ' ' + cs
                    i += 1
                else:
                    break
            flow.append(ListFlowable(
                [ListItem(Paragraph(inline(x), S['cuerpo']), leftIndent=13)
                 for x in items],
                bulletType='1' if es_ol else 'bullet',
                start=None if es_ol else '•',
                leftIndent=15 if es_ol else 14,
                bulletFontSize=8 if not es_ol else 9,
                spaceBefore=2, spaceAfter=6))
            continue

        # --- parrafo ---
        buf = [s]
        i += 1
        while i < len(lineas):
            nxt = lineas[i]
            ns = nxt.strip()
            if (not ns or ns.startswith(('#', '>', '|', '```', '![')) or
                    re.match(r'^[-*]\s+', ns) or re.match(r'^\d+\.\s+', ns) or
                    re.match(r'^(-{3,}|\*{3,}|_{3,})$', ns)):
                break
            buf.append(ns)
            i += 1
        flow.append(Paragraph(inline(' '.join(buf)), S['cuerpo']))

    return flow


class EntradaIndice(Paragraph):
    """Parrafo de encabezado que ademas debe aparecer en el indice.

    La notificacion al indice la emite la plantilla en `afterFlowable`, que es
    el punto donde ReportLab conoce ya el numero de pagina definitivo.
    """

    def __init__(self, parrafo, nivel, texto, clave, S):
        Paragraph.__init__(self, parrafo.text, parrafo.style)
        self._nivel = nivel
        self._texto = texto
        self._clave = clave

    def draw(self):
        Paragraph.draw(self)
        self.canv.bookmarkPage(self._clave)


# ==========================================================================
# Plantilla de pagina
# ==========================================================================
class Documento(BaseDocTemplate):
    def __init__(self, ruta, titulo, **kw):
        BaseDocTemplate.__init__(self, ruta, pagesize=A4,
                                 leftMargin=20 * mm, rightMargin=20 * mm,
                                 topMargin=20 * mm, bottomMargin=18 * mm,
                                 title=titulo, author='spinpy',
                                 subject=titulo,
                                 **kw)
        self.titulo_encabezado = titulo
        marco = Frame(self.leftMargin, self.bottomMargin,
                      self.width, self.height, id='normal')
        self.addPageTemplates([
            PageTemplate(id='portada', frames=[marco]),
            PageTemplate(id='cuerpo', frames=[marco], onPage=self.decorar),
        ])
        self.titulo = titulo
        self._ultimo_nivel = -1

    def beforeDocument(self):
        # multiBuild recorre el documento varias veces para resolver el indice,
        # y en cada pasada ReportLab reinicia su arbol de marcadores. El
        # contador de niveles debe reiniciarse con el, o la segunda pasada
        # arrancaria heredando el nivel con el que termino la primera.
        self._ultimo_nivel = -1

    def afterFlowable(self, flowable):
        """Alimenta el indice y el arbol de marcadores del PDF."""
        if not isinstance(flowable, EntradaIndice):
            return
        # El arbol de marcadores no admite saltos de nivel: un documento que
        # empiece por un encabezado de segundo nivel obliga a degradarlo.
        nivel = min(flowable._nivel, self._ultimo_nivel + 1, 3)
        self._ultimo_nivel = nivel
        self.canv.addOutlineEntry(flowable._texto[:110], flowable._clave,
                                  level=nivel)
        self.notify('TOCEntry',
                    (flowable._nivel, flowable._texto, self.page,
                     flowable._clave))

    def decorar(self, canv, doc):
        canv.saveState()
        canv.setFont('Helvetica', 7.4)
        canv.setFillColor(colors.HexColor('#7A828F'))
        # El encabezado sale del TITULO del documento, no de una constante.
        # En el anexo original estaba fijo, lo que hacia que cualquier otro
        # documento compuesto con esta herramienta se rotulara como si fuera
        # el anexo. (Unica diferencia con la copia de Validacion_Anexo/code.)
        canv.drawString(20 * mm, A4[1] - 12.5 * mm, getattr(
            self, 'titulo_encabezado', 'Anexo de validacion metodologica'))
        canv.setStrokeColor(GRIS_LIN)
        canv.setLineWidth(0.4)
        canv.line(20 * mm, A4[1] - 14.5 * mm, A4[0] - 20 * mm, A4[1] - 14.5 * mm)
        canv.line(20 * mm, 13.5 * mm, A4[0] - 20 * mm, 13.5 * mm)
        canv.drawCentredString(A4[0] / 2.0, 9.5 * mm, '%d' % canv.getPageNumber())
        canv.restoreState()


def portada(S, meta):
    e = []
    e.append(Spacer(1, 46 * mm))
    e.append(Paragraph(meta['titulo'], S['portada_tit']))
    e.append(Spacer(1, 4 * mm))
    e.append(HRFlowable(width='58%', thickness=1.4, color=AZUL_SUAVE,
                        hAlign='CENTER'))
    e.append(Spacer(1, 7 * mm))
    e.append(Paragraph(meta['subtitulo'], S['portada_sub']))
    e.append(Spacer(1, 30 * mm))
    for l in meta['bloque']:
        e.append(Paragraph(l, S['portada_txt']))
    e.append(Spacer(1, 22 * mm))
    e.append(HRFlowable(width='38%', thickness=0.6, color=GRIS_LIN,
                        hAlign='CENTER'))
    e.append(Spacer(1, 5 * mm))
    for l in meta['pie']:
        e.append(Paragraph(l, S['portada_txt']))
    e.append(NextPageTemplate('cuerpo'))
    e.append(PageBreak())
    return e


def componer(ruta_md, ruta_pdf, titulo_encabezado, etiqueta_figura='Figura',
             titulo_llano='EN PALABRAS SENCILLAS'):
    """Compone un Markdown SIN portada ni indice: documentos cortos.

    Es lo que usa el informe para publicacion. Las rutas de las figuras se
    resuelven respecto a la carpeta del Markdown.
    """
    global ETIQUETA_FIGURA, TITULO_LLANO
    with open(ruta_md, encoding='utf-8') as fh:
        md = fh.read()
    base = os.path.dirname(os.path.abspath(ruta_md))
    S = construir_estilos()
    doc = Documento(ruta_pdf, titulo_encabezado)
    # La primera plantilla es la portada del anexo, que aqui no existe.
    doc.pageTemplates = [p for p in doc.pageTemplates if p.id == 'cuerpo']
    doc.titulo_encabezado = titulo_encabezado
    previa = (ETIQUETA_FIGURA, TITULO_LLANO)
    ETIQUETA_FIGURA, TITULO_LLANO = etiqueta_figura, titulo_llano
    try:
        doc.build(parsear(md, base, S))
    finally:
        ETIQUETA_FIGURA, TITULO_LLANO = previa
    return ruta_pdf


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    ruta_md, ruta_pdf = sys.argv[1], sys.argv[2]
    base = os.path.dirname(os.path.abspath(ruta_md))

    with open(ruta_md, encoding='utf-8') as fh:
        md = fh.read()

    # El titulo y el subtitulo se toman del propio Markdown y se retiran del
    # cuerpo, para no repetirlos tras la portada.
    meta = {
        'titulo': 'Anexo de validacion metodologica',
        'subtitulo': 'Verificacion de las metricas morfometricas y de los\n'
                     'algoritmos de optimizacion de AppFinal_v2',
        'bloque': [],
        'pie': [],
    }
    m = re.match(r'^#\s+(.*)', md.split('\n')[0])
    if m:
        meta['titulo'] = m.group(1).strip()
        md = '\n'.join(md.split('\n')[1:])
    m = re.search(r'^\*\*(.+?)\*\*\s*$', md, re.M)
    if m:
        meta['subtitulo'] = m.group(1).strip()
        md = md[:m.start()] + md[m.end():]

    metafile = os.path.join(base, 'portada.txt')
    if os.path.exists(metafile):
        with open(metafile, encoding='utf-8') as fh:
            partes = fh.read().split('---')
        if len(partes) >= 1:
            meta['bloque'] = [l for l in partes[0].strip().split('\n') if l.strip()]
        if len(partes) >= 2:
            meta['pie'] = [l for l in partes[1].strip().split('\n') if l.strip()]

    S = construir_estilos()
    doc = Documento(ruta_pdf, meta['titulo'])

    hist = []
    hist.extend(portada(S, meta))

    toc = TableOfContents()
    toc.levelStyles = [S['toc1'], S['toc2'], S['toc3']]
    hist.append(Paragraph('Indice', S['h1']))
    hist.append(Spacer(1, 4))
    hist.append(toc)
    hist.append(PageBreak())

    hist.extend(parsear(md, base, S))

    doc.multiBuild(hist)
    print('PDF generado: %s (%.1f KB)'
          % (ruta_pdf, os.path.getsize(ruta_pdf) / 1024.0))


if __name__ == '__main__':
    main()
