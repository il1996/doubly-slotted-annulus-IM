# -*- coding: utf-8 -*-
"""Production: every derived percentage of Section 4 (slotting ratio) computed
from the archived result files, so that each printed figure has a producing line.

Inputs (outputs/python/): prod_fem_results.json (prod_fem.py), fem_kc_results.json
(t_fem2.py), tiling_sweep_inf_iron.json (t_op3.py), prod_op_results.json (prod_op.py),
t4_graded_results.json (t4_graded.py), prod_pos_33_32_results.json (prod_pos_33_32.py),
prod_limits_results.json (prod_limits.py).  Carter's factors from dsop.carter.
Each line names: tiling, truncation, basis, position convention, reference grid.
Run from code/python/ ; writes outputs/python/prod_section4_pct_out.txt."""
import os, json
import numpy as np
from dsop import Machine, carter
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'python')
LOG = open(os.path.join(W, 'prod_section4_pct_out.txt'), 'w', encoding='utf-8')
def say(s):
    print(s, flush=True); LOG.write(s + '\n')
def J(n): return json.load(open(os.path.join(W, n)))
fem, til, op, gr = J('prod_fem_results.json'), J('tiling_sweep_inf_iron.json'), J('prod_op_results.json'), J('t4_graded_results.json')
pos = J('prod_pos_33_32_results.json') if os.path.isfile(os.path.join(W, 'prod_pos_33_32_results.json')) else None
lim = J('prod_limits_results.json') if os.path.isfile(os.path.join(W, 'prod_limits_results.json')) else None
M = Machine(); C = carter(M); pct = lambda a, b: 100 * (a / b - 1)
A0 = fem['conv_real']['0.0_0.028'][0]; A0_045 = fem['conv_real']['0.0_0.045'][0]; A38_045 = fem['conv_real']['0.375_0.045'][0]; A38_028 = fem['conv_real']['0.375_0.028'][0]
Apos = np.array(list(fem['pos_real'].values())); Bpos = np.array(list(fem['pos_neu'].values())); B0 = fem['conv_neu']['0.0_0.025'][0]
say("=== Section 4 derived percentages (all references: FE formulation A = real slots, B = flux-barrier openings) ===")
say("Carter (dsop.carter): approximate ks %.4f kr %.4f product %.4f ; exact ks %.4f kr %.4f product %.4f" % (C['ks'], C['kr'], C['k'], C['ks_exact'], C['kr_exact'], C['k_exact']))
say("\n[A1] Table 4 / Sec. 4.1: single-surface FE (A, lc 0.045, phi = 0) vs exact conformal Carter: stator %.4f vs %.4f -> %+.2f %% ; rotor %.4f vs %.4f -> %+.2f %%" % (
    fem['single']['0.045'][0], C['ks_exact'], pct(fem['single']['0.045'][0], C['ks_exact']), fem['single']['0.045'][1], C['kr_exact'], pct(fem['single']['0.045'][1], C['kr_exact'])))
say("[A2] FE (A) mesh sequence, phi = 0: " + " / ".join("%.4f (lc %s)" % (fem['conv_real']['0.0_%s' % lc][0], lc) for lc in ['0.08', '0.06', '0.045', '0.035', '0.028']))
say("[A3] FE (A) mesh sequence, phi = 3/8 tau_r: " + " / ".join("%.4f (lc %s)" % (fem['conv_real']['0.375_%s' % lc][0], lc) for lc in ['0.08', '0.06', '0.045', '0.035', '0.028']))
say("[A4] FE (A), lc 0.045, 12 positions: min %.4f max %.4f mean %.4f ; modulation +-%.2f %% about the mean" % (Apos.min(), Apos.max(), Apos.mean(), 100 * (Apos.max() - Apos.min()) / 2 / Apos.mean()))
say("[A5] mesh uncertainty of the position sweep: lc 0.045 vs 0.028 at phi = 0 %+.3f %%, at 3/8 %+.3f %% (FE A, phi = 0 and 3/8)" % (pct(A0_045, A0), pct(A38_045, A38_028)))
say("[A6] FE (A) 12-position mean %.4f (lc 0.045) vs Carter approximate product %.4f: %+.2f %% ; vs exact product %.4f: %+.3f %%" % (Apos.mean(), C['k'], pct(Apos.mean(), C['k']), C['k_exact'], pct(Apos.mean(), C['k_exact'])))
say("[A7] FE (B) phi = 0, lc 0.025: %.4f vs FE (A) phi = 0, lc 0.028: %.4f -> %+.2f %% ('13 %% in the continuous limit')" % (B0, A0, pct(B0, A0)))
say("[A8] FE (B) 12 positions, lc 0.035: min %.4f max %.4f mean %.4f, modulation +-%.2f %% ; mean vs FE (A) mean (lc 0.045) %+.2f %% ('11 %%') ; modulation ratio B/A %.2f" % (
    Bpos.min(), Bpos.max(), Bpos.mean(), 100 * (Bpos.max() - Bpos.min()) / 2 / Bpos.mean(), pct(Bpos.mean(), Apos.mean()), ((Bpos.max() - Bpos.min()) / Bpos.mean()) / ((Apos.max() - Apos.min()) / Apos.mean())))
hat = {k: v['kC'] for k, v in til.items() if k.startswith('p1_')}
kmin, kmax = min(hat, key=hat.get), max(hat, key=hat.get)
say("[A9] Table 5, Phi_O = 0, symmetric hat, N_h 8192, phi = 0, assembled comparator: min %.4f (%s) max %.4f (%s) vs FE (A) phi = 0 lc 0.028 %.4f -> %+.1f %% to %+.1f %% ('16 to 31 %%')" % (
    hat[kmin], kmin[3:], hat[kmax], kmax[3:], A0, pct(hat[kmin], A0), pct(hat[kmax], A0)))
say("[A10] (17,4) Phi_O = 0, hat, phi = 0: %.4f ; FE (A) lc 0.028 / operator - 1 = %+.1f %% ('under-estimates by 19 %%')" % (op['pos_op_17_4']['0.0'], pct(A0, op['pos_op_17_4']['0.0'])))
say("[A11] cavities, hat, phi = 0, assembled comparator: (33,32) %.4f -> %+.2f %% ; (65,16) %.4f -> %+.2f %% vs FE (A) lc 0.028 %.4f ('0.2 to 0.4 %%')" % (
    op['cav_sweep']['33_32'], pct(op['cav_sweep']['33_32'], A0), op['cav_sweep']['65_16'], pct(op['cav_sweep']['65_16'], A0), A0))
for k in ['33_4', '33_8', '33_16', '33_32', '65_4', '65_8', '65_16']:
    say("[A12] graded q = %.1f, asym. hat, cavities, phi = 0, (%s): graded %.4f -> %+.2f %% ; uniform asym. hat %.4f -> %+.2f %% vs FE (A) lc 0.028 %.4f" % (
        gr['q'], k.replace('_', ','), gr['cav_graded'][k]['kC'], pct(gr['cav_graded'][k]['kC'], A0), gr['cav_uniform_p1a'][k]['kC'], pct(gr['cav_uniform_p1a'][k]['kC'], A0), A0))
if 'nh_check' in gr: say("[A12'] N_h check of t4_graded: %s" % json.dumps(gr['nh_check']))
kc = np.array([op['pos_cav'][k] for k in fem['pos_real']]); dev = 100 * (kc / Apos - 1)
say("[A13] cavities (33,16), hat, N_h 8192, 12 positions, assembled comparator: mean %.4f vs FE (A) lc 0.045 mean %.4f -> %+.2f %% ; per position %+.2f .. %+.2f %% ; range %.3f .. %.3f" % (kc.mean(), Apos.mean(), pct(kc.mean(), Apos.mean()), dev.min(), dev.max(), kc.min(), kc.max()))
cs = [op['cav_sweep']['33_%d' % n] for n in [4, 8, 16, 32]]; d = np.diff(cs)
say("[A13'] cavities increments along n_O at n_T = 33: %s ; successive ratios %s (factor 1/r = %s)" % (np.round(d, 4), np.round(d[1:] / d[:-1], 2), np.round(d[:-1] / d[1:], 2)))
say("[A14] Phi_O = 0 at (65,16), phi = 0: symmetric %.4f vs asymmetric %.4f -> %+.2f %% (no FE)" % (til['p1_65_16']['kC'], til['p1a_65_16']['kC'], pct(til['p1a_65_16']['kC'], til['p1_65_16']['kC'])))
say("[A15] Carter approximate product %.4f vs cavities (33,16) phi = 0 %.4f: %+.1f %% ; vs Phi_O = 0 (17,4) %.4f: %+.1f %% (no FE)" % (C['k'], op['cav_sweep']['33_16'], pct(C['k'], op['cav_sweep']['33_16']), op['pos_op_17_4']['0.0'], pct(C['k'], op['pos_op_17_4']['0.0'])))
say("[B3] graded Phi_O = 0, asym. hat, (33,16): %.4f vs FE (B) lc 0.025 %.4f -> %+.2f %% ; uniform asym. hat (33,16) %.4f -> %+.2f %%" % (
    gr['neumann_graded']['33_16']['kC'], B0, pct(gr['neumann_graded']['33_16']['kC'], B0), til['p1a_33_16']['kC'], pct(til['p1a_33_16']['kC'], B0)))
for nT in [33, 65]:
    u = [til['p1a_%d_%d' % (nT, n)]['kC'] for n in [2, 4, 8, 16]]; g = [gr['neumann_graded']['%d_%d' % (nT, n)]['kC'] for n in [2, 4, 8, 16]]
    say("[B3'] n_T = %d increments along n_O: uniform %s ratios %s ; graded %s ratios %s" % (nT, np.round(np.diff(u), 4), np.round(np.diff(u)[1:] / np.diff(u)[:-1], 2), np.round(np.diff(g), 4), np.round(np.diff(g)[1:] / np.diff(g)[:-1], 2)))
if 'ratio_study' in gr: say("[B3''] ratio study q: %s" % json.dumps(gr['ratio_study'])[:600])
h65 = [til['p1_65_%d' % n]['kC'] for n in [2, 4, 8, 16]]; h16 = [til['p1_%d_16' % n]['kC'] for n in [9, 17, 33, 65]]
d65, d16 = np.diff(h65), np.diff(h16); r65, r16 = d65[-1] / d65[-2], d16[-1] / d16[-2]
say("[B4] symmetric hat, Phi_O = 0, phi = 0: along n_O at n_T = 65 increments %s (ratio %.2f, tail A %.4f -> %.4f) ; along n_T at n_O = 16 increments %s (ratio %.2f, tail A %.4f -> %.4f) ; both tails removed from (65,16) %.4f: %.4f ('between 1.44 and 1.47')" % (
    np.round(d65, 4), r65, d65[-1] * r65 / (1 - r65), h65[-1] + d65[-1] * r65 / (1 - r65), np.round(d16, 4), r16, d16[-1] * r16 / (1 - r16), h16[-1] + d16[-1] * r16 / (1 - r16), h65[-1], h65[-1] + d65[-1] * r65 / (1 - r65) + d16[-1] * r16 / (1 - r16)))
if lim:
    e = lim['branches']['p1a']['nT_at_nO_16']
    say("[B4'] asymmetric hat, Phi_O = 0, phi = 0, n_T branch at n_O = 16 (prod_limits.py): terms %s, limit A %.4f, limit B %.4f ('1.455')" % (np.round(e['terms'], 4), e['limit_A'], e['limit_B']))
    for conv in ['A', 'B']:
        e = lim['branches']['p1a'].get('double_from_nO_limits_' + conv)
        if e: say("[B4''] asymmetric hat double limit (n_O-limits by %s): A %s B %s vs FE (B) %.5f" % (conv, ("%.4f (%+.2f %%)" % (e['limit_A'], pct(e['limit_A'], B0))) if e['limit_A'] else "-", "%.4f (%+.2f %%)" % (e['limit_B'], pct(e['limit_B'], B0)), B0))
if pos:
    for b in ['p1', 'p1a']:
        s = pos['series'][b]
        say("[P] (33,32) cavities, N_h 8192, %s, 12 positions, vs FE (A) lc 0.045 (prod_pos_33_32.py): closed-form comparator mean %.5f -> %+.3f %% [%+.3f .. %+.3f] ; assembled comparator mean %.5f -> %+.3f %% [%+.3f .. %+.3f]" % (
            b, s['mean_closed'], s['dev_mean_closed'], s['dev_pos_min_closed'], s['dev_pos_max_closed'], s['mean_asm'], s['dev_mean_asm'], s['dev_pos_min_asm'], s['dev_pos_max_asm']))
say("DONE")
