# -*- coding: utf-8 -*-
"""Provenance of the few text values that are arithmetic on archived results:
rotor-position means of the slotting ratio (Table 4 and Section 4.3) from
prod_fem_results.json / prod_op_results.json, and the geometric increments of
Table 5 quoted in Section 4.2.  Run from code/python/ ; writes to stdout."""
import json, os, sys
import numpy as np
sys.dont_write_bytecode = True
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'python')
fem = json.load(open(os.path.join(W, 'prod_fem_results.json')))
op = json.load(open(os.path.join(W, 'prod_op_results.json')))
til = json.load(open(os.path.join(W, 'tiling_sweep_inf_iron.json')))
def stats(d, lab):
    v = np.array(list(d.values()), float)
    print("%-32s n=%2d  min %.4f  max %.4f  mean %.4f  (+-%.2f %% about the mean)" % (lab, len(v), v.min(), v.max(), v.mean(), 100 * (v.max() - v.min()) / 2 / v.mean()))
for k, lab in [('pos_real', 'FE (A), real slots'), ('pos_neu', 'FE (B), flux-barrier openings')]:
    if k in fem: stats(fem[k], lab)
for k, lab in [('pos_cav', 'operator + cavities (33,16)'), ('pos_op_33_16', 'operator Phi_O=0 (33,16)'), ('pos_op_17_4', 'operator Phi_O=0 (17,4)')]:
    if k in op: stats(op[k], lab)
if 'pos_cav' in op and 'pos_real' in fem:
    kc = {float(k): v for k, v in op['pos_cav'].items()}; kf = {float(k): v for k, v in fem['pos_real'].items()}
    common = sorted(set(kc) & set(kf))
    dev = [100 * (kc[x] - kf[x]) / kf[x] for x in common]
    print("cavity (33,16) against FE (A) at the same positions: deviation min %+.2f %% max %+.2f %% (n=%d)" % (min(dev), max(dev), len(common)))
print("Table 5 increments along n_O at n_T = 65, Phi_O = 0, hat: ", np.round(np.diff([til['p1_65_%d' % n]['kC'] for n in [2, 4, 8, 16]]), 4))
print("Table 5 increments along n_T at n_O = 16, Phi_O = 0, hat: ", np.round(np.diff([til['p1_%d_16' % n]['kC'] for n in [9, 17, 33, 65]]), 4))
print("Table 5 asym. hat along n_T at n_O = 16:", [round(til['p1a_%d_16' % n]['kC'], 4) for n in [9, 17, 33, 65]], "increments", np.round(np.diff([til['p1a_%d_16' % n]['kC'] for n in [9, 17, 33, 65]]), 4))
cs = op['cav_sweep']
print("cavity increments along n_O at n_T = 33 (hat):", np.round(np.diff([cs['33_%d' % n] for n in [4, 8, 16, 32]]), 4), "ratios", np.round(np.diff([cs['33_%d' % n] for n in [4, 8, 16, 32]])[1:] / np.diff([cs['33_%d' % n] for n in [4, 8, 16, 32]])[:-1], 2))
