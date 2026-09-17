# -*- coding: utf-8 -*-
"""Production: the three checks of Appendix B / Section 4.1, each with its
element size named.

(1) Current sheet on a slotless bore (formulation A, fem_slots).  The slot
    ampere-turns AT_k are applied as a sheet of zero width at the slot centres
    theta_k = 2 pi (k + 1/2)/N_s: the Neumann boundary term int K phi_i ds with
    K = AT_k delta(theta - theta_k)/R_s, i.e. a nodal load AT_k split linearly
    between the two bore nodes that bracket theta_k.  Compared with the closed
    form of eq. (1) (fem_slots.staircase_Fp) at mid-gap, 8192 samples.  For
    reference the 2 mm x 0.04 mm sheet of prod_fem.py (Jgap) is run on the
    same meshes.  Seven element sizes lc_gap = 0.12 ... 0.028 mm.
(2) Deep slot against a sheet at its mouth (formulation A): stator slotted,
    rotor smooth, phi = 0; (d1) current uniform in the slot body (the archived
    construction, fem_slots.run_case); (d2) same ampere-turns as a 2 mm sheet in
    the gap layer under the opening (R_s - 0.04 mm < r < R_s); (d3) as a sheet
    filling the isthmus (R_s < r < R_s + h_s0).  Fundamentals to seven decimals,
    the decimal at which equality fails, and the whole mid-gap waveform (max
    harmonic difference, rms and pointwise differences).  lc_gap = 0.06,
    0.045, 0.035, 0.028 mm.
(3) Formulation B (fem_annulus) against the staircase potential: Dirichlet on
    the 96 stator arcs (faces at U_j, openings at the mean of the two
    neighbours), fundamental at mid-gap against the analytic Fourier coefficient
    of the same piecewise-constant potential (the construction of
    t_ann_check.py).  lc = 0.08, 0.05, 0.035, 0.025 mm.
Run from code/python/ ; writes outputs/python/prod_appB_checks_results.json and
prod_appB_checks_out.txt (this transcript)."""
import os, time, json
import numpy as np, scipy.sparse.linalg as spla
import fem_slots as F, fem_annulus as FA
from dsop import Machine, MU0
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'python')
LOG = open(os.path.join(W, 'prod_appB_checks_out.txt'), 'w', encoding='utf-8')
def say(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
M = Machine(); t0 = time.time(); i3 = np.sqrt(2) * np.array([1.0, -0.5, -0.5]); AT = M.slot_ampere_turns(i3); Us = M.tooth_mmf(i3)
B_closed = MU0 * (M.p / M.Rm) * abs(F.staircase_Fp(M, 1.0)) * np.cosh(M.p * np.log(M.Rm / M.Rr)) / np.sinh(M.p * M.X)
thk = 2 * np.pi * (np.arange(M.Ns) + 0.5) / M.Ns
out = dict(B_closed=B_closed, sheet_slotless={}, deep_slot={}, formulation_B={})
say("prod_appB_checks: closed-form fundamental of the staircase potential at mid-gap %.9f T/A (I0 = 1 A)" % B_closed)
def solve_load(xyz, tri, fb):
    K, area = F.assemble_P1(xyz, tri, np.full(len(tri), 1.0 / MU0)); n = len(xyz); free = np.arange(1, n)
    A = np.zeros(n); A[free] = spla.spsolve(K[free][:, free].tocsc(), fb[free]); return A
def sheet_elems(xyz, tri, area, halfw_mm, r_lo, r_hi):
    c = xyz[tri].mean(axis=1); r = np.hypot(c[:, 0], c[:, 1]); ang = np.arctan2(c[:, 1], c[:, 0])
    J = np.zeros(len(tri)); m = (r > r_lo) & (r < r_hi); k = np.round((ang / (2 * np.pi)) * M.Ns - 0.5).astype(int) % M.Ns
    dang = ((ang - thk[k]) + np.pi) % (2 * np.pi) - np.pi; m &= np.abs(dang) < (halfw_mm * 1e-3 / M.Rs); sel = np.where(m)[0]
    for kk in range(M.Ns):
        e = sel[k[sel] == kk]; J[e] = AT[kk] / area[e].sum()
    return J
def wave(xyz, tri, A):
    bg1, Bn, a, th = F.midgap_fundamental_A(M, xyz, tri, A)
    Br = np.real(np.fft.irfft(np.fft.rfft(a) * 1j * np.arange(len(a) // 2 + 1), n=len(a))) / M.Rm
    return bg1, Bn, Br
# ------------------------------------------------------------------ (1)
say("\n(1) current sheet on a slotless bore, delta sheet split on the two bracketing bore nodes ; Q = Bg1(FE)/Bg1(closed) - 1")
say("   %-8s %-8s %-12s %-14s %-14s" % ("lc (mm)", "nodes", "bore pitch", "Q delta sheet", "Q 2 mm sheet (prod_fem Jgap)"))
for lc in [0.12, 0.08, 0.06, 0.05, 0.045, 0.035, 0.028]:
    t = time.time(); xyz, tri = F.build_air_mesh(M, 0.0, lc_gap=lc, stator_slots=False, rotor_slots=False)
    r = np.hypot(xyz[:, 0], xyz[:, 1]); bn = np.where(r > M.Rs - 1e-6)[0]; angb = np.arctan2(xyz[bn, 1], xyz[bn, 0]); o = np.argsort(angb); bn = bn[o]; angb = angb[o]
    fb = np.zeros(len(xyz))
    for kk in range(M.Ns):
        d = ((angb - thk[kk]) + np.pi) % (2 * np.pi) - np.pi; j = int(np.argmin(np.abs(d)))
        jl, jr = (j, (j + 1) % len(bn)) if d[j] <= 0 else ((j - 1) % len(bn), j)
        dl = ((thk[kk] - angb[jl]) + np.pi) % (2 * np.pi) - np.pi; dr = ((angb[jr] - thk[kk]) + np.pi) % (2 * np.pi) - np.pi
        wl = dr / (dl + dr); fb[bn[jl]] += wl * AT[kk]; fb[bn[jr]] += (1 - wl) * AT[kk]
    q_delta = wave(xyz, tri, solve_load(xyz, tri, fb))[0] / B_closed - 1
    K, area = F.assemble_P1(xyz, tri, np.ones(len(tri))); J = sheet_elems(xyz, tri, area, 1.0, M.Rs - 0.04e-3, np.inf); A, _ = F.solve_Az(M, xyz, tri, J)
    q_jgap = wave(xyz, tri, A)[0] / B_closed - 1; pitch = 2 * np.pi * M.Rs / len(bn)
    out['sheet_slotless'][str(lc)] = dict(nodes=int(len(xyz)), bore_nodes=int(len(bn)), bore_pitch_mm=float(pitch * 1e3), Q_delta=float(q_delta), Q_jgap_2mm=float(q_jgap))
    say("   %-8.3f %-8d %-12.4f %+13.3e  %+13.3e   (%.0f s)" % (lc, len(xyz), pitch * 1e3, q_delta, q_jgap, time.time() - t))
qd = [v['Q_delta'] for v in out['sheet_slotless'].values()]
say("   delta sheet: max |Q| over the seven sizes %.2e ; monotone decrease in |Q| from lc 0.12 to 0.028: %s" % (max(map(abs, qd)), "yes" if all(abs(qd[i + 1]) < abs(qd[i]) for i in range(len(qd) - 1)) else "no"))
# ------------------------------------------------------------------ (2)
say("\n(2) deep slot (d1) against a sheet at its mouth (d2 gap layer, d3 isthmus) ; stator slotted, rotor smooth, phi = 0")
for lc in [0.06, 0.045, 0.035, 0.028]:
    t = time.time(); xyz, tri = F.build_air_mesh(M, 0.0, lc_gap=lc, stator_slots=True, rotor_slots=False)
    ids = F.slot_body_ids(M, xyz, tri); K, area = F.assemble_P1(xyz, tri, np.ones(len(tri)))
    ba = np.array([area[ids == k].sum() for k in range(M.Ns)]); J1 = np.zeros(len(tri)); m = ids >= 0; J1[m] = AT[ids[m]] / ba[ids[m]]
    Wv = {}
    for lab, J in [('d1', J1), ('d2', sheet_elems(xyz, tri, area, 1.0, M.Rs - 0.04e-3, M.Rs)), ('d3', sheet_elems(xyz, tri, area, 1.0, M.Rs, M.Rs + F.HS0 * 1e-3))]:
        A, _ = F.solve_Az(M, xyz, tri, J); Wv[lab] = wave(xyz, tri, A)
    b1, Bn1, Br1 = Wv['d1']; rec = dict(nodes=int(len(xyz)), Bg1_d1=float(b1), kC_d1=float(B_closed / b1))
    for lab in ['d2', 'd3']:
        b, Bn, Br = Wv[lab]; dBn = np.abs(Bn - Bn1); nmax = int(np.argmax(dBn[1:400]) + 1)
        ndec = next((k for k in range(1, 10) if round(b, k) != round(b1, k)), 10)
        rec[lab] = dict(Bg1=float(b), rel_diff=float(b / b1 - 1), first_decimal_differing=ndec, max_harmonic_diff_rel=float(dBn[nmax] / b1), n_of_max=nmax,
                        rms_rel=float(np.sqrt(np.mean((Br - Br1) ** 2)) / np.sqrt(np.mean(Br1 ** 2))), max_pointwise_rel=float(np.max(np.abs(Br - Br1)) / np.max(np.abs(Br1))))
        say("   lc %.3f nodes %6d : (d1) Bg1 %.7f | (%s) Bg1 %.7f  rel %+.2e  equal to %d decimals, differs at the %dth | max_n |dB_n|/Bg1 %.1e at n = %d | rms(dBr)/rms(Br) %.1e | max|dBr|/max|Br| %.1e   (%.0f s)" % (
            lc, len(xyz), b1, lab, b, b / b1 - 1, ndec - 1, ndec, dBn[nmax] / b1, nmax, rec[lab]['rms_rel'], rec[lab]['max_pointwise_rel'], time.time() - t))
    out['deep_slot'][str(lc)] = rec
# ------------------------------------------------------------------ (3)
say("\n(3) formulation B against the staircase potential (Dirichlet on the 96 stator arcs, t_ann_check.py construction) ; ratio FE/analytic of the fundamental")
tau = 2 * np.pi / M.Ns; a_o = M.bs0 / M.Rs; a_f = tau - a_o; Uop = 0.5 * (Us + np.roll(Us, -1)); nn = np.arange(1, 200)
edges = []; vals = []
for k in range(M.Ns):
    c = tau * k; edges += [c - a_f / 2, c + a_f / 2]; vals += [Us[k], Uop[k]]
cs = FA.piecewise_constant_fourier(edges, vals, nn); Ban = FA.annulus_Bn_from_surface_coeffs(M, cs, np.zeros_like(cs), nn.astype(float))
for lc in [0.08, 0.05, 0.035, 0.025]:
    t = time.time(); xyz, tri, arcs, fs, fr = FA.build_annulus_mesh(M, 0.0, lc)
    phi, R = FA.solve_scalar(M, xyz, tri, arcs, fs, fr, Us, 0.0, opening='dirichlet', Uopen_s=Uop, Uopen_r=np.zeros(M.Nr))
    Bn = FA.midgap_harmonics_scalar(M, xyz, tri, phi); ratio = Bn[M.p] / Ban[M.p - 1]
    out['formulation_B'][str(lc)] = dict(nodes=int(len(xyz)), Bg1_FE=float(Bn[M.p]), Bg1_analytic=float(Ban[M.p - 1]), ratio_minus_1=float(ratio - 1))
    say("   lc %.3f nodes %6d : FE %.7f analytic %.7f ratio - 1 = %+.2e   (%.0f s)" % (lc, len(xyz), Bn[M.p], Ban[M.p - 1], ratio - 1, time.time() - t))
json.dump(out, open(os.path.join(W, 'prod_appB_checks_results.json'), 'w'), indent=1)
say("DONE (%.0f s)" % (time.time() - t0))
