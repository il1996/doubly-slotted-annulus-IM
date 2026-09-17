# -*- coding: utf-8 -*-
"""Production: refinement branches of the condensed operator with Phi_O = 0 and
their geometric extrapolation, computed and printed by this script.

  closure         Phi_O = 0 (flux-barrier openings), infinite iron
  position        phi = 0 only
  truncation      N_h = 8192
  bases           'p1a' (asymmetric hat, eq. (6)) and 'p1' (symmetric hat)
  comparator      assembled smooth annulus, same basis (dsop.slotting_ratio,
                  the convention of Table 5); the closed-form comparator value
                  (fem_slots.staircase_Fp) is stored alongside
  tilings         the Table 5 grid n_T in {9,17,33,65} x n_O in {2,4,8,16},
                  read from tiling_sweep_inf_iron.json (t_op3.py), extended here
                  by (33,32) (33,64) (65,32) (65,64) (129,16) (129,32) (129,64)
  reference       FE formulation B, phi = 0, finest mesh lc 0.025
                  (prod_fem_results.json conv_neu['0.0_0.025'])

Extrapolation of a monotone sequence v_0..v_n with increments d_k = v_k - v_{k-1}:
  measured ratio      r = d_n / d_{n-1}
  convention A        tail = d_n * r / (1 - r)          (geometric series with ratio r)
  convention B        tail = d_n                        (ratio taken as exactly 1/2)
  limit               v_n + tail
Branches: along n_O at fixed n_T, along n_T at fixed n_O.  Double limit: the
branch limits along n_O, one per n_T, extrapolated along n_T with the same two
conventions.  Run from code/python/ ; writes outputs/python/prod_limits_results.json
and prod_limits_out.txt (this transcript)."""
import os, sys, time, json, gc
import numpy as np
from dsop import Machine, MU0, slotting_ratio
from fem_slots import staircase_Fp
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'python')
LOG = open(os.path.join(W, 'prod_limits_out.txt'), 'w', encoding='utf-8')
def say(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
M = Machine(); Nh = 8192; t0 = time.time()
B_closed = MU0 * (M.p / M.Rm) * abs(staircase_Fp(M, 1.0)) * np.cosh(M.p * np.log(M.Rm / M.Rr)) / np.sinh(M.p * M.X)
til = json.load(open(os.path.join(W, 'tiling_sweep_inf_iron.json')))
feB = json.load(open(os.path.join(W, 'prod_fem_results.json')))['conv_neu']['0.0_0.025'][0]
EXT = [(33, 32), (33, 64), (65, 32), (65, 64), (129, 16), (129, 32), (129, 64)]
res_path = os.path.join(W, 'prod_limits_results.json')
res = json.load(open(res_path)) if os.path.isfile(res_path) else {}
res.setdefault('points', {}); res['reference_feB_phi0_lc0.025'] = feB; res['B_closed'] = B_closed; res['Nh'] = Nh
say("prod_limits: Phi_O = 0, phi = 0, N_h = %d ; reference FE(B) lc 0.025: %.5f ; closed-form comparator %.9f T/A" % (Nh, feB, B_closed))
# ---- points: Table 5 grid from t_op3.py, extension computed here ------------
for basis in ['p1a', 'p1']:
    for nT in [9, 17, 33, 65]:
        for nO in [2, 4, 8, 16]:
            k = "%s_%d_%d" % (basis, nT, nO); res['points'][k] = dict(kC_asm=til[k]['kC'], ncol=til[k]['ncol'], source='tiling_sweep_inf_iron.json (t_op3.py)')
    for (nT, nO) in EXT:
        k = "%s_%d_%d" % (basis, nT, nO)
        if k in res['points'] and res['points'][k].get('source') == 'prod_limits.py':
            say("  %s: already computed, kC %.5f" % (k, res['points'][k]['kC_asm'])); continue
        t = time.time(); r = slotting_ratio(M, nT, nO, Nh, basis)
        res['points'][k] = dict(kC_asm=float(r['kC']), kC_closed=float(B_closed / r['Bg1_slot']), Bg1_slot=float(r['Bg1_slot']), ncol=int(r['ncol']), source='prod_limits.py')
        json.dump(res, open(res_path, 'w'), indent=1); del r; gc.collect()
        say("  (%3d,%2d) %-3s [%5d col.] kC (assembled) %.5f ; (closed form) %.5f  (%.0f s)" % (nT, nO, basis, res['points'][k]['ncol'], res['points'][k]['kC_asm'], res['points'][k]['kC_closed'], time.time() - t))
# ---- extrapolation -----------------------------------------------------------
def val(basis, nT, nO):
    k = "%s_%d_%d" % (basis, nT, nO); return res['points'][k]['kC_asm'] if k in res['points'] else None
def extrap(v):
    d = np.diff(v); out = dict(terms=[float(x) for x in v], increments=[float(x) for x in d])
    if len(d) < 2: out.update(ratios=[], r=None, limit_A=None, limit_B=None, note='fewer than three terms: no extrapolation'); return out
    ratios = [float(d[i] / d[i - 1]) for i in range(1, len(d))]; r = ratios[-1]
    out['ratios'] = ratios; out['r'] = r; out['limit_B'] = float(v[-1] + d[-1])
    out['limit_A'] = float(v[-1] + d[-1] * r / (1 - r)) if 0 < r < 1 else None
    out['note'] = 'geometric (0 < r < 1)' if 0 < r < 1 else 'NOT geometric (r outside (0,1)): convention A undefined'
    return out
def show(lab, e):
    say("  %-38s terms %s | incr %s | ratios %s | r = %s | limit A (r/(1-r)) %s | limit B (tail = last incr) %s | %s" % (
        lab, " ".join("%.4f" % x for x in e['terms']), " ".join("%+.4f" % x for x in e['increments']), " ".join("%.3f" % x for x in e['ratios']),
        ("%.3f" % e['r']) if e['r'] is not None else "-", ("%.4f" % e['limit_A']) if e['limit_A'] is not None else "-", ("%.4f" % e['limit_B']) if e['limit_B'] is not None else "-", e['note']))
res['branches'] = {}
for basis in ['p1a', 'p1']:
    say("\n=== basis %s : branches along n_O (fixed n_T) ===" % basis); B = {}
    for nT in [9, 17, 33, 65, 129]:
        seq = [(nO, val(basis, nT, nO)) for nO in [2, 4, 8, 16, 32, 64] if val(basis, nT, nO) is not None]
        if len(seq) < 3: continue
        e = extrap([x for _, x in seq]); e['nO'] = [n for n, _ in seq]; B['nO_at_nT_%d' % nT] = e; show("n_T = %3d, n_O = %s" % (nT, ",".join(str(n) for n, _ in seq)), e)
    say("=== basis %s : branches along n_T (fixed n_O) ===" % basis)
    for nO in [16, 32, 64]:
        seq = [(nT, val(basis, nT, nO)) for nT in [9, 17, 33, 65, 129] if val(basis, nT, nO) is not None]
        if len(seq) < 3: continue
        e = extrap([x for _, x in seq]); e['nT'] = [n for n, _ in seq]; B['nT_at_nO_%d' % nO] = e; show("n_O = %3d, n_T = %s" % (nO, ",".join(str(n) for n, _ in seq)), e)
    say("=== basis %s : double limit ===" % basis)
    for conv in ['A', 'B']:
        lim = [(nT, B['nO_at_nT_%d' % nT]['limit_' + conv]) for nT in [9, 17, 33, 65] if 'nO_at_nT_%d' % nT in B and B['nO_at_nT_%d' % nT]['limit_' + conv] is not None]
        if len(lim) < 3: say("  n_O-limits (convention %s): fewer than three defined -> no double limit" % conv); continue
        e = extrap([x for _, x in lim]); e['nT'] = [n for n, _ in lim]; B['double_from_nO_limits_' + conv] = e
        show("n_O-limits (conv. %s) at n_T = %s, extrapolated along n_T" % (conv, ",".join(str(n) for n, _ in lim)), e)
        for c2 in ['A', 'B']:
            L = e['limit_' + c2]
            if L is not None: say("     -> double limit (n_O-limits by %s, n_T-extrapolation by %s): %.4f ; vs FE(B) %.5f: %+.2f %%" % (conv, c2, L, feB, 100 * (L / feB - 1)))
    for nO in [16, 32, 64]:
        if 'nT_at_nO_%d' % nO in B:
            e = B['nT_at_nO_%d' % nO]
            say("  n_T-limit at n_O = %2d: A %s, B %s ; vs FE(B): %s / %s" % (nO, ("%.4f" % e['limit_A']) if e['limit_A'] else "-", "%.4f" % e['limit_B'],
                ("%+.2f %%" % (100 * (e['limit_A'] / feB - 1))) if e['limit_A'] else "-", "%+.2f %%" % (100 * (e['limit_B'] / feB - 1))))
    res['branches'][basis] = B
json.dump(res, open(res_path, 'w'), indent=1)
say("DONE (%.0f s)" % (time.time() - t0))
