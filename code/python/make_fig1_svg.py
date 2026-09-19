# -*- coding: utf-8 -*-
"""Emit Figure 1 as a self-contained, structured SVG (true slot geometry)."""
import numpy as np

TAU = 10.719
BS0, HS0, BS1, HS1, BS2, HS2 = 2.0, 0.5, 5.236, 2.5, 8.472, 24.724
BR0, HR0, DR1, DR2, HR2 = 2.0, 1.0, 5.80425613768, 2.038213646116, 16.78981901905
GAP, HS_WIN, HR_WIN, XOFF = 1.0, 4.6, 7.4, 2.2
H = 0.22                                  # height of a tiling column
X0, X1 = -0.7, TAU + 0.7
S = 26.0                                  # px per mm
TITLE_H, LABEL_H = 1.30, 0.95             # mm reserved above the window
PW = (X1 - X0) * S
PH = (HS_WIN + GAP + HR_WIN + TITLE_H + LABEL_H) * S
GUT = 0.9 * S

IRON, TILE = '#d5d8dc', '#aeb6bf'
BLUE, RED, GREEN = '#1b4f72', '#c0392b', '#1e8449'


def X(x): return (x - X0) * S
def Y(y): return (HS_WIN + GAP + TITLE_H + LABEL_H - y) * S
def f(v): return ('%.2f' % v).rstrip('0').rstrip('.')
def pt(x, y): return '%s,%s' % (f(X(x)), f(Y(y)))


def stator_half(hmax):
    w = BS1 / 2 + (BS2 - BS1) / 2 * (hmax - HS0 - HS1) / HS2
    return [(BS0 / 2, 0.0), (BS0 / 2, HS0), (BS1 / 2, HS0 + HS1), (w, hmax)]


def rotor_half(hmax, narc=34):
    r1, r2 = DR1 / 2, DR2 / 2
    c1 = HR0 + r1
    a = np.arcsin((r1 - r2) / HR2)
    t0, t1 = np.arcsin((BR0 / 2) / r1), np.pi / 2 + a
    p = [(BR0 / 2, 0.0), (BR0 / 2, HR0)]
    for t in np.linspace(t0, t1, narc):
        y = c1 - r1 * np.cos(t)
        if y > hmax:
            break
        p.append((r1 * np.sin(t), y))
    ytan, xtan = c1 + r1 * np.sin(a), r1 * np.cos(a)
    if hmax > ytan:
        p.append((xtan - (hmax - ytan) * np.tan(a), hmax))
    return p


def slot_path(half, xc, sign, shift):
    right = [(xc + x, sign * y + shift) for x, y in half]
    left = [(xc - x, sign * y + shift) for x, y in reversed(half)]
    pts = right + left
    return 'M ' + ' L '.join(pt(*p) for p in pts) + ' Z'


def open_walls(half, xc, sign, shift):
    right = [(xc + x, sign * y + shift) for x, y in half]
    left = [(xc - x, sign * y + shift) for x, y in reversed(half)]
    return ('M ' + ' L '.join(pt(*p) for p in right),
            'M ' + ' L '.join(pt(*p) for p in left))


def arrow(p0, p1, rad=0.0):
    """Straight or arc3-like curved arrow, in the style of matplotlib."""
    if rad == 0.0:
        return 'M %s L %s' % (pt(*p0), pt(*p1))
    x0, y0 = X(p0[0]), Y(p0[1]); x1, y1 = X(p1[0]), Y(p1[1])
    cx, cy = (x0 + x1) / 2 - rad * (y1 - y0) / 2, (y0 + y1) / 2 + rad * (x1 - x0) / 2
    return 'M %s,%s Q %s,%s %s,%s' % (f(x0), f(y0), f(cx), f(cy), f(x1), f(y1))


def tiles(xa, xb, y, n, colour, tag):
    out = []
    for i in range(n):
        xx = xa + i * (xb - xa) / n
        out.append('      <rect class="tile %s" x="%s" y="%s" width="%s" height="%s" fill="%s"/>'
                   % (tag, f(X(xx)), f(Y(y + H)), f((xb - xa) / n * S), f(H * S), colour))
    return out


def panel(mode, ox):
    col = RED if mode == 'a' else GREEN
    xs, xr = TAU / 2, TAU / 2 + XOFF
    o = ['  <g id="panel-%s" transform="translate(%s,0)">' % (mode, f(ox))]
    title = ('(a) condition Φ&#8338; = 0 : the opening is a flux barrier' if mode == 'a'
             else '(b) cavity admittance Q : flux enters the slot walls')
    o.append('    <text class="title" x="%s" y="%s" text-anchor="middle">%s</text>' % (f(PW / 2), f(Y(GAP + HS_WIN + 1.55)), title))
    o.append('    <text class="lab" x="%s" y="%s" text-anchor="start">n<tspan class="sub">T</tspan> columns per tooth face</text>'
             % (f(X(X0 + 0.15)), f(Y(GAP + HS_WIN + 0.45))))
    o.append('    <text class="lab" x="%s" y="%s" text-anchor="start" fill="%s">n<tspan class="sub">O</tspan> columns per opening</text>'
             % (f(X(X1 - 5.05)), f(Y(GAP + HS_WIN + 0.45)), col))

    o.append('    <g class="iron">')
    o.append('      <rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="#000" stroke-width="1.1"/>'
             % (f(X(X0)), f(Y(GAP + HS_WIN)), f((X1 - X0) * S), f(HS_WIN * S), IRON))
    o.append('      <rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="#000" stroke-width="1.1"/>'
             % (f(X(X0)), f(Y(0)), f((X1 - X0) * S), f(HR_WIN * S), IRON))
    o.append('    </g>')

    sh, rh = stator_half(HS_WIN), rotor_half(HR_WIN)
    o.append('    <g class="slots">')
    o.append('      <path id="stator-slot-%s" d="%s" fill="#fff"/>' % (mode, slot_path(sh, xs, +1, GAP)))
    for d in open_walls(sh, xs, +1, GAP):
        o.append('      <path class="wall" d="%s"/>' % d)
    o.append('      <path id="rotor-bar-%s" d="%s" fill="#fff"/>' % (mode, slot_path(rh, xr, -1, 0.0)))
    for d in open_walls(rh, xr, -1, 0.0):
        o.append('      <path class="wall" d="%s"/>' % d)
    o.append('    </g>')

    o.append('    <g class="tiling">')
    o += tiles(X0, xs - BS0 / 2, GAP - H, 4, TILE, 'face')
    o += tiles(xs + BS0 / 2, X1, GAP - H, 4, TILE, 'face')
    o += tiles(xs - BS0 / 2, xs + BS0 / 2, GAP - H, 3, col, 'opening')
    o += tiles(X0, xr - BR0 / 2, 0.0, 4, TILE, 'face')
    o += tiles(xr + BR0 / 2, X1, 0.0, 4, TILE, 'face')
    o += tiles(xr - BR0 / 2, xr + BR0 / 2, 0.0, 3, col, 'opening')
    o.append('    </g>')

    o.append('    <g class="flux">')
    for xx in np.linspace(xs - 3.4, xs + 3.4, 9):
        if not (xs - BS0 / 2 - 0.2 < xx < xs + BS0 / 2 + 0.2):
            o.append('      <path class="arrow" d="M %s L %s" stroke="%s" marker-end="url(#head-gap)"/>'
                     % (pt(xx, H + 0.06), pt(xx, GAP - H), BLUE))
    for sgn in (-1, +1):
        xa = xs + sgn * 0.45
        if mode == 'a':                       # the opening is a barrier: flux runs to the tooth faces
            d = 'M %s C %s %s %s' % (pt(xa, H + 0.06), pt(xa, 0.62),
                                     pt(xs + sgn * 1.15, 0.45), pt(xs + sgn * 1.15, GAP - H))
        else:                                 # the flux crosses the opening and lands on the slot wall
            d = 'M %s C %s %s %s' % (pt(xa, H + 0.06), pt(xa, GAP + 0.35),
                                     pt(xs + sgn * 1.15, GAP + 0.95), pt(xs + sgn * 1.62, GAP + HS0 + 1.0))
        o.append('      <path class="arrow" d="%s" stroke="%s" marker-end="url(#head-%s)"/>' % (d, col, mode))
    o.append('    </g>')

    o.append('    <g class="labels">')
    o.append('      <path class="dim" d="M %s L %s M %s L %s M %s L %s"/>'
             % (pt(X0 + 0.10, H), pt(X0 + 0.46, H), pt(X0 + 0.10, GAP - H), pt(X0 + 0.46, GAP - H),
                pt(X0 + 0.28, H), pt(X0 + 0.28, GAP - H)))
    o.append('      <text class="sym" x="%s" y="%s" text-anchor="start" dominant-baseline="middle">g</text>'
             % (f(X(X0 + 0.56)), f(Y(GAP / 2))))
    o.append('      <text class="sym" x="%s" y="%s" text-anchor="middle">b<tspan class="sub">0</tspan></text>' % (f(X(xs)), f(Y(GAP + 0.16))))
    o.append('      <text class="note" x="%s" y="%s" text-anchor="middle">stator slot</text>' % (f(X(xs)), f(Y(GAP + HS_WIN - 1.1))))
    o.append('      <text class="note" x="%s" y="%s" text-anchor="middle">rotor bar</text>' % (f(X(xr)), f(Y(-HR_WIN + 1.2))))
    o.append('      <text class="note" x="%s" y="%s" text-anchor="end">stator</text>'
             % (f(X(X1 - 0.15)), f(Y(GAP + HS_WIN - 0.45))))
    o.append('      <text class="note" x="%s" y="%s" text-anchor="start">rotor</text>'
             % (f(X(X0 + 0.15)), f(Y(-HR_WIN + 0.30))))
    o.append('    </g>')
    o.append('  </g>')
    return o


W = 2 * PW + GUT
head = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" font-family="Inter, Helvetica, Arial, sans-serif">
  <title>Figure 1 — slot-opening closure of the doubly slotted annulus</title>
  <desc>Slot cross-sections to scale: stator opening 2.000 x 0.500 mm, wedge to 5.236 mm over 2.500 mm, trapezoidal body to 8.472 mm; rotor isthmus 2.000 x 1.000 mm, 5.804 mm circle centred 3.902 mm below the bore joined by tangents of half-angle 6.44 deg to a 2.038 mm circle 16.790 mm deeper, 21.711 mm of slot in all. Slot bodies cut by the window edge; air gap drawn wider than scale.</desc>
  <defs>
    <marker id="head-gap" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
      <path d="M 0 1 L 9 5 L 0 9 z" fill="{blue}"/>
    </marker>
    <marker id="head-a" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
      <path d="M 0 1 L 9 5 L 0 9 z" fill="{red}"/>
    </marker>
    <marker id="head-b" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
      <path d="M 0 1 L 9 5 L 0 9 z" fill="{green}"/>
    </marker>
    <style>
      .wall  {{ fill:none; stroke:#000; stroke-width:1.3; stroke-linejoin:round; }}
      .tile  {{ stroke:#000; stroke-width:.5; }}
      .arrow {{ fill:none; stroke-width:1.4; }}
      .dim   {{ fill:none; stroke:#000; stroke-width:1.1; }}
      .title {{ font-size:11.5px; fill:#111; }}
      .lab   {{ font-size:10.5px; }}
      .sub   {{ font-size:.72em; baseline-shift:-25%; }}
      .bold  {{ font-weight:700; }}
      .sym   {{ font-size:13px; font-style:italic; }}
      .note  {{ font-size:11.5px; fill:#111; }}
    </style>
  </defs>
  <rect width="{w}" height="{h}" fill="#fff"/>
""".format(w=f(W), h=f(PH), blue=BLUE, red=RED, green=GREEN)

svg = head + '\n'.join(panel('a', 0) + panel('b', PW + GUT)) + '\n</svg>\n'
open('/tmp/claude-0/-home-claude/7143d4e0-b0a2-5917-b5fc-cf31bb630d6a/scratchpad/figs/fig1_opening_conditions.svg', 'w').write(svg)
print('svg written,', len(svg.splitlines()), 'lines,', len(svg), 'bytes')
