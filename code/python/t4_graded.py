# -*- coding: utf-8 -*-
"""
T4 -- tiling graded towards the tooth corners.

Uniform tiling (manuscript, Table 5): tile_surface of dsop.py.  Graded
tiling: on every arc (tooth face and slot opening, both bores) the column
widths follow a geometric progression that shrinks towards the two ends of
the arc, with common ratio q between consecutive columns (q = 1 -> uniform),
at equal column count.  The asymmetric hat 'p1a' (eq. 6 of the manuscript,
partition of unity on any grid) is used for every graded run; the uniform
runs in 'p1a' are those of tiling_sweep_inf_iron.json (Phi_O = 0) and are
recomputed here for the cavity-coupled operator.

Sweeps at phi = 0, infinite iron, N_h = 8192:
  (a) Phi_O = 0 : nT in {9,17,33,65} x nO in {2,4,8,16}, graded q
  (b) cavities  : nT in {17,33,65} x nO in {4,8,16} (+ (33,32)), uniform p1a and graded q
  (c) ratio study at (33,8) and (33,16), q in {1, 1.25, 1.5, 2, 3}
  (d) truncation check of one graded case at N_h = 16384
Results -> t4_graded_results.json, transcript on stdout.
"""
import sys, json, time
sys.dont_write_bytecode = True
import numpy as np
from dsop import Machine, assemble, condense, condense_cavity, solve_linear, fine_potentials, \
    fine_potentials_cavity, Bg1, smooth_surface, tile_surface, invariants
from cavity_graded import cavity_element_knots
import pickle

M = Machine()
KC_FE = 1.2934      # FE (A), real slots, phi = 0, finest mesh (Table 4)
KC_NEU = 1.4598     # FE (B), flux-barrier openings, phi = 0, finest mesh (Table 4)


def graded_widths(n, a, q):
    """n column widths summing to a, geometric with ratio q from each end
    towards the centre (symmetric).  q = 1 -> uniform."""
    if n == 1 or q == 1.0:
        return np.full(n, a / n)
    m = n // 2
    if n % 2 == 0:
        half = q ** np.arange(m)                 # edge -> centre
        w = np.concatenate([half, half[::-1]])
    else:
        half = q ** np.arange(m + 1)             # edge -> centre column
        w = np.concatenate([half[:-1], [half[-1]], half[:-1][::-1]])
    return a * w / w.sum()


def tile_surface_graded(N, R, b0, nT, nO, phi=0.0, q=1.5):
    tau = 2 * np.pi / N
    a_o = b0 / R; a_f = tau - a_o
    wf = graded_widths(nT, a_f, q); wo = graded_widths(nO, a_o, q)
    th, d, isface, tooth = [], [], [], []
    for k in range(N):
        c = tau * k + phi
        e = c - a_f / 2
        for i in range(nT):
            th.append(e + wf[i] / 2); d.append(wf[i]); isface.append(True); tooth.append(k); e += wf[i]
        for i in range(nO):
            th.append(e + wo[i] / 2); d.append(wo[i]); isface.append(False); tooth.append(-1); e += wo[i]
    th = np.array(th); d = np.array(d); isface = np.array(isface); tooth = np.array(tooth)
    thp = np.roll(th, 1); thn = np.roll(th, -1)
    hl = (th - thp) % (2 * np.pi); hr = (thn - th) % (2 * np.pi)
    return dict(th=th, d=d, isface=isface, tooth=tooth, hl=hl, hr=hr, N=N, R=R, wo=wo, wf=wf)


def opening_centres_mm(nO, q):
    """Mouth column centres (mm from the opening centre) for the graded opening, b0 = 2 mm."""
    b0 = 2.0
    w = graded_widths(nO, b0, q)
    e = np.concatenate([[0.0], np.cumsum(w)])
    return -b0 / 2 + 0.5 * (e[:-1] + e[1:])


def kc_run(nT, nO, Nh, basis, q, Q=None):
    """One slotting-ratio evaluation at phi = 0.  q = 1 -> uniform tiling."""
    i3 = np.sqrt(2) * np.array([1.0, -0.5, -0.5]); Us = M.tooth_mmf(i3)
    if q == 1.0:
        ss = tile_surface(M.Ns, M.Rs, M.bs0, nT, nO, 0.0); sr = tile_surface(M.Nr, M.Rr, M.br0, nT, nO, 0.0)
    else:
        ss = tile_surface_graded(M.Ns, M.Rs, M.bs0, nT, nO, 0.0, q); sr = tile_surface_graded(M.Nr, M.Rr, M.br0, nT, nO, 0.0, q)
    A, proj = assemble(M, ss, sr, Nh, basis)
    if Q is None:
        cond = condense(A, ss, sr)
        U, Phi = solve_linear(M, cond, Us); u = fine_potentials(cond, U)
    else:
        cond = condense_cavity(A, ss, sr, Q[0], Q[1], nO)
        U, Phi = solve_linear(M, cond, Us); u = fine_potentials_cavity(cond, U)
    b_slot = Bg1(M, proj, u, cond)
    ss0 = smooth_surface(M.Ns, 0.0); sr0 = smooth_surface(M.Nr, 0.0)
    A0, proj0 = assemble(M, ss0, sr0, Nh, basis); cond0 = condense(A0, ss0, sr0)
    U0, _ = solve_linear(M, cond0, Us); u0 = fine_potentials(cond0, U0)
    b_smooth = Bg1(M, proj0, u0, cond0)
    inv = invariants(A, proj, M)
    return dict(kC=float(b_smooth / b_slot), Bg1_slot=float(b_slot), Bg1_smooth=float(b_smooth), ncol=int(A.shape[0]),
                I1=float(inv['I1']), I4=float(inv['I4']), dmin=float(ss['d'].min()), dmax=float(ss['d'].max()))


def cavity_Q(nO, q, lc=None):
    xc = opening_centres_mm(nO, q)
    if lc is None:
        lc = min(0.02, 0.25 * np.diff(np.concatenate([[-1.0], xc, [1.0]])).min())
        lc = max(lc, 0.003)
    Qs, _ = cavity_element_knots('stator', xc, M.L, lc_mouth=lc)
    Qr, _ = cavity_element_knots('rotor', xc, M.L, lc_mouth=lc)
    return (Qs, Qr), lc


if __name__ == '__main__':
    t0 = time.time(); out = {'q': None, 'neumann_graded': {}, 'neumann_uniform_ref': {}, 'cav_uniform_p1a': {}, 'cav_graded': {},
                             'ratio_study': {}, 'nh_check': {}, 'cavity_lc': {}}
    q = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
    out['q'] = q
    til = json.load(open('tiling_sweep_inf_iron.json'))
    Qu = pickle.load(open('cavity_Q.pkl', 'rb'))
    print("T4 graded tiling sweep, q = %.3f, N_h = 8192, phi = 0, infinite iron, basis p1a" % q, flush=True)
    # --- chain identity check: uniform p1a (33,16) Phi_O=0 must reproduce the JSON
    r = kc_run(33, 16, 8192, 'p1a', 1.0)
    print("CHECK uniform p1a (33,16) Phi_O=0 : kC = %.6f (json %.6f)  I1 %.1e I4 %.1e (%.0f s)" % (r['kC'], til['p1a_33_16']['kC'], r['I1'], r['I4'], time.time() - t0), flush=True)
    out['check_uniform_p1a_33_16'] = r
    # --- (c) ratio study first (cheap, informs the reading)
    for (nT, nO) in [(33, 8), (33, 16)]:
        for qq in [1.0, 1.25, 1.5, 2.0, 3.0]:
            r = kc_run(nT, nO, 8192, 'p1a', qq)
            Q, lc = cavity_Q(nO, qq)
            rc = kc_run(nT, nO, 8192, 'p1a', qq, Q)
            out['ratio_study']["%d_%d_q%.2f" % (nT, nO, qq)] = dict(neumann=r, cav=rc, lc=lc)
            print("RATIO (%d,%d) q=%.2f : Phi_O=0 kC=%.5f | cav kC=%.5f | smallest col %.2e rad, largest %.2e | I1 %.1e I4 %.1e (%.0f s)" %
                  (nT, nO, qq, r['kC'], rc['kC'], r['dmin'], r['dmax'], r['I1'], r['I4'], time.time() - t0), flush=True)
    # --- (a) Phi_O = 0, graded
    for nT in [9, 17, 33, 65]:
        for nO in [2, 4, 8, 16]:
            r = kc_run(nT, nO, 8192, 'p1a', q)
            out['neumann_graded']["%d_%d" % (nT, nO)] = r
            ref = til["p1a_%d_%d" % (nT, nO)]['kC']
            out['neumann_uniform_ref']["%d_%d" % (nT, nO)] = ref
            print("NEU  (%2d,%2d) graded kC=%.5f | uniform p1a (json) %.5f | cols %d (%.0f s)" % (nT, nO, r['kC'], ref, r['ncol'], time.time() - t0), flush=True)
    # --- (b) cavities: uniform p1a and graded
    QG = {}
    for nO in [4, 8, 16, 32]:
        QG[nO], lc = cavity_Q(nO, q); out['cavity_lc'][str(nO)] = lc
    for nT in [17, 33, 65]:
        for nO in [4, 8, 16, 32]:
            if nO == 32 and nT != 33: continue
            ru = kc_run(nT, nO, 8192, 'p1a', 1.0, Qu[nO])
            rg = kc_run(nT, nO, 8192, 'p1a', q, QG[nO])
            out['cav_uniform_p1a']["%d_%d" % (nT, nO)] = ru; out['cav_graded']["%d_%d" % (nT, nO)] = rg
            print("CAV  (%2d,%2d) uniform p1a kC=%.5f | graded kC=%.5f | FE 1.2934 (%.0f s)" % (nT, nO, ru['kC'], rg['kC'], time.time() - t0), flush=True)
            json.dump(out, open('t4_graded_results.json', 'w'), indent=1)
    # --- (d) truncation check on the graded (33,16): N_h = 16384
    r16 = kc_run(33, 16, 16384, 'p1a', q); r16c = kc_run(33, 16, 16384, 'p1a', q, QG[16])
    out['nh_check'] = dict(neumann_8192=out['neumann_graded']['33_16']['kC'], neumann_16384=r16['kC'],
                           cav_8192=out['cav_graded']['33_16']['kC'], cav_16384=r16c['kC'])
    print("NH   graded (33,16): Phi_O=0 %.5f -> %.5f at N_h=16384 ; cav %.5f -> %.5f" % (out['nh_check']['neumann_8192'], r16['kC'], out['nh_check']['cav_8192'], r16c['kC']), flush=True)
    json.dump(out, open('t4_graded_results.json', 'w'), indent=1)
    # --- convergence reading
    print("\nCONVERGENCE along nO (increments and their ratios), Phi_O = 0:")
    for nT in [9, 17, 33, 65]:
        yu = [til["p1a_%d_%d" % (nT, n)]['kC'] for n in [2, 4, 8, 16]]
        yg = [out['neumann_graded']["%d_%d" % (nT, n)]['kC'] for n in [2, 4, 8, 16]]
        du = np.diff(yu); dg = np.diff(yg)
        print("  nT=%2d uniform: %s incr %s ratios %s" % (nT, np.round(yu, 4), np.round(du, 4), np.round(du[1:] / du[:-1], 2)))
        print("  nT=%2d graded : %s incr %s ratios %s" % (nT, np.round(yg, 4), np.round(dg, 4), np.round(dg[1:] / dg[:-1], 2)))
    print("CONVERGENCE along nO, cavities (uniform p1a / graded), FE 1.2934:")
    for nT in [17, 33, 65]:
        ns = [4, 8, 16] + ([32] if nT == 33 else [])
        yu = [out['cav_uniform_p1a']["%d_%d" % (nT, n)]['kC'] for n in ns]
        yg = [out['cav_graded']["%d_%d" % (nT, n)]['kC'] for n in ns]
        du = np.diff(yu); dg = np.diff(yg)
        print("  nT=%2d uniform: %s incr %s ratios %s" % (nT, np.round(yu, 4), np.round(du, 4), np.round(du[1:] / du[:-1], 2) if len(du) > 1 else '-'))
        print("  nT=%2d graded : %s incr %s ratios %s" % (nT, np.round(yg, 4), np.round(dg, 4), np.round(dg[1:] / dg[:-1], 2) if len(dg) > 1 else '-'))
    print("DONE %.0f s" % (time.time() - t0))
