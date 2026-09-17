# -*- coding: utf-8 -*-
"""Production: cavity-coupled operator at the finest uniform tiling (33, 32),
twelve rotor positions, both surface bases, two smooth comparators.

  tiling          n_T = 33, n_O = 32 (5980 bore columns), uniform
  truncation      N_h = 8192
  closure         slot cavities (cavity_Q.pkl of prod_op.py, n_O = 32, lc_mouth 0.005)
  bases           'p1'  symmetric hat, 'p1a' asymmetric hat of eq. (6)
  positions       phi = i/12 tau_r, i = 0..11 (as prod_fem.py / prod_op.py)
  comparators     (a) assembled smooth annulus, same basis, same N_h
                      (dsop.slotting_ratio_cavity -> 'Bg1_smooth', the Table 5 convention);
                  (b) closed form of eq. (1) with the staircase potential
                      (fem_slots.staircase_Fp, the convention of the FE chain)
  reference       FE formulation A, real slots, lc_gap = 0.045, twelve positions
                  (prod_fem_results.json 'pos_real'), mean 1.2686

Run from code/python/ ; writes outputs/python/prod_pos_33_32_results.json and
prod_pos_33_32_out.txt (this transcript)."""
import os, sys, time, json, pickle
import numpy as np
from dsop import Machine, MU0, slotting_ratio_cavity
from fem_slots import staircase_Fp
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'python')
LOG = open(os.path.join(W, 'prod_pos_33_32_out.txt'), 'w', encoding='utf-8')
def say(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
M = Machine(); tau_r = 2 * np.pi / M.Nr; Nh = 8192; nT, nO = 33, 32; t0 = time.time()
Q = pickle.load(open(os.path.join(W, 'cavity_Q.pkl'), 'rb'))[nO]
fem = json.load(open(os.path.join(W, 'prod_fem_results.json')))
kf = {round(float(k), 6): v for k, v in fem['pos_real'].items()}; fr = sorted(kf); kf_mean = float(np.mean([kf[f] for f in fr]))
B_closed = MU0 * (M.p / M.Rm) * abs(staircase_Fp(M, 1.0)) * np.cosh(M.p * np.log(M.Rm / M.Rr)) / np.sinh(M.p * M.X)
say("prod_pos_33_32: tiling (%d, %d), N_h = %d, cavities, bases p1 / p1a, 12 positions" % (nT, nO, Nh))
say("comparator (b) closed form: %.9f T/A ; reference FE(A) lc 0.045: mean %.5f, phi = 0: %.5f" % (B_closed, kf_mean, kf[0.0]))
out = dict(tiling=[nT, nO], Nh=Nh, closure='cavities', positions='i/12 tau_r', reference='prod_fem_results.json pos_real (FE A, lc 0.045)',
           B_closed=B_closed, fe_pos=kf, fe_mean=kf_mean, series={})
for basis in ['p1', 'p1a']:
    s = dict(Bg1_slot={}, Bg1_smooth_asm=None, kC_asm={}, kC_closed={}, ncol=None)
    for f in fr:
        r = slotting_ratio_cavity(M, nT, nO, Nh, Q[0], Q[1], basis, phi=f * tau_r)
        s['Bg1_slot'][f] = float(r['Bg1_slot']); s['Bg1_smooth_asm'] = float(r['Bg1_smooth']); s['ncol'] = int(r['ncol'])
        s['kC_asm'][f] = float(r['kC']); s['kC_closed'][f] = float(B_closed / r['Bg1_slot'])
        say("  %-3s phi/tau_r = %.4f : Bg1_slot %.6f T/A  kC (a) assembled %.5f  kC (b) closed %.5f | FE %.5f  dev (a) %+.3f %%  dev (b) %+.3f %%  (%.0f s)" % (
            basis, f, r['Bg1_slot'], r['kC'], B_closed / r['Bg1_slot'], kf[f], 100 * (r['kC'] / kf[f] - 1), 100 * (B_closed / r['Bg1_slot'] / kf[f] - 1), time.time() - t0))
    for c in ['asm', 'closed']:
        v = np.array([s['kC_' + c][f] for f in fr]); d = 100 * (v / np.array([kf[f] for f in fr]) - 1)
        s['mean_' + c] = float(v.mean()); s['min_' + c] = float(v.min()); s['max_' + c] = float(v.max())
        s['dev_mean_' + c] = float(100 * (v.mean() / kf_mean - 1)); s['dev_pos_min_' + c] = float(d.min()); s['dev_pos_max_' + c] = float(d.max())
        say("  %-3s comparator (%s): kC mean %.5f (min %.5f, max %.5f) ; vs FE mean %.5f: %+.3f %% ; per position %+.3f .. %+.3f %%" % (
            basis, 'a, assembled' if c == 'asm' else 'b, closed form', v.mean(), v.min(), v.max(), kf_mean, 100 * (v.mean() / kf_mean - 1), d.min(), d.max()))
    out['series'][basis] = {k: ({str(kk): vv for kk, vv in val.items()} if isinstance(val, dict) else val) for k, val in s.items()}
for c in ['asm', 'closed']:
    a, b = out['series']['p1']['mean_' + c], out['series']['p1a']['mean_' + c]
    say("basis difference p1a/p1 - 1 on the 12-position mean, comparator (%s): %+.3f %%" % (c, 100 * (b / a - 1)))
json.dump(out, open(os.path.join(W, 'prod_pos_33_32_results.json'), 'w'), indent=1)
say("DONE (%.0f s)" % (time.time() - t0))
