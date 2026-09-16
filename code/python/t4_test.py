# -*- coding: utf-8 -*-
"""Checks of cavity_graded against cavity.py (uniform knots) and of the
rotor bar-current source vector (compatibility, conservation, rectangular
limit)."""
import sys, time, pickle
sys.dont_write_bytecode = True
import numpy as np
from dsop import Machine, MU0
from cavity_graded import cavity_element_knots, cavity_source_rotor, rotor_source_vector
from t4_graded import graded_widths, opening_centres_mm, tile_surface_graded, kc_run
M = Machine(); L = M.L; t0 = time.time()
Qu = pickle.load(open('cavity_Q.pkl', 'rb'))
for nO in [4, 8]:
    xc = -1.0 + (np.arange(nO) + 0.5) * 2.0 / nO
    Qs, _ = cavity_element_knots('stator', xc, L, lc_mouth=0.02); Qr, _ = cavity_element_knots('rotor', xc, L, lc_mouth=0.02)
    print("nO=%d uniform knots: max|Qs-Qs_pkl|/max|Qs| = %.2e ; rotor %.2e" % (nO, np.abs(Qs - Qu[nO][0]).max() / np.abs(Qs).max(), np.abs(Qr - Qu[nO][1]).max() / np.abs(Qr).max()))
# graded widths
for n, q in [(4, 1.5), (5, 1.5), (16, 1.5), (16, 2.0)]:
    w = graded_widths(n, 1.0, q); print("widths n=%d q=%.1f:" % (n, q), np.round(w / w.max(), 4), "sum", w.sum())
print("centres nO=8 q=1.5 (mm):", np.round(opening_centres_mm(8, 1.5), 4))
# rotor source vector
for nO in [4, 16]:
    s, Qr, PhiJ, ramp, xc, info = rotor_source_vector(nO, L, q=1.0)
    print("nO=%d: Abar %.3f mm2 (bar 80.70) | sum f_J %.4f A, sum f_b %.4f A | Phi_J columns %s | Phi_J wall %.3e | sum Phi_J %.2e" %
          (nO, info['Abar_mm2'], info['tot_src'], info['tot_bnd'], np.array2string(PhiJ[:nO], precision=3), PhiJ[nO], PhiJ.sum()))
    print("      Q_r[ramp;0] columns %s wall %.3e" % (np.array2string((Qr[:nO+1, :nO+1] @ ramp)[:nO], precision=3), (Qr[:nO+1, :nO+1] @ ramp)[nO]))
    print("      s (Wb/A) columns %s wall %.3e | sum s %.2e" % (np.array2string(s[:nO], precision=3), s[nO], s.sum()))
    # tangential field check at the mouth: H_x = (1/mu0) dA/dy should be I/b0 = 0.5 A/mm -> in SI 500 A/m
    A = info['A']; xyz = info['xyz']; tri = info['tri']
    # take elements touching the mouth and estimate dA/dy at their centroid from the P1 gradient
    x = xyz[tri, 0]; y = xyz[tri, 1]
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
    area2 = (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
    dAdy = np.sum(c * A[tri], 1) / area2          # per mm
    cen = xyz[tri].mean(1); near = cen[:, 1] < 0.05
    Hx = dAdy[near] / MU0 * 1e3                    # A/m
    print("      H_x near the mouth: mean %.1f A/m (expected 500), std %.1f" % (Hx.mean(), Hx.std()))
# graded scale: reference-free identity of one graded run vs uniform at q=1
r1 = kc_run(17, 4, 2048, 'p1a', 1.0); r2 = kc_run(17, 4, 2048, 'p1a', 1.5)
print("kc (17,4) Nh=2048 p1a uniform %.5f graded q=1.5 %.5f (I1 %.1e/%.1e, I4 %.1e/%.1e) (%.0f s)" % (r1['kC'], r2['kC'], r1['I1'], r2['I1'], r1['I4'], r2['I4'], time.time() - t0))
