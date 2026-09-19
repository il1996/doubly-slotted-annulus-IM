# -*- coding: utf-8 -*-
"""Production: sweep of the opening-to-gap ratio b_0/g on the linear chain of
Section IV (infinite iron), to place the single point b_0/g = 8.2 of the paper
on a trend.

What varies    the air gap g alone: R_s fixed, R_r = R_s - g (Machine() with
               R_r, D_r, X = ln(R_s/R_r), R_m overridden; nothing else).
What does not  slot openings b_0 (stator 2.000 mm, rotor isthmus 2.000 mm), the
               slot profiles (fem_slots), the tooth pitches, the tiling
               (n_T, n_O) = (33, 16), the symmetric hat 'p1', N_h = 8192, the
               cavity admittances Q (independent of g; reloaded from
               cavity_Q.pkl of prod_op.py, n_O = 16, lc_mouth 0.01).
Values         b_0/g in {4, 12, 16, 2}, i.e. g = 0.500, 0.167, 0.125, 1.000 mm;
               the nominal point (b_0/g = 8.23, g = 0.24298 mm) is NOT
               recomputed: it is read from prod_fem_results.json ('pos_real')
               and prod_op_results.json ('pos_op_33_16', 'pos_cav').
At each g      the twelve rotor positions phi = i/12 tau_r of pos_real:
  (a) FE formulation A (fem_slots.run_case), real slots, production mesh
      lc_gap = 0.045 mm, ratio (5) with the closed-form smooth comparator of
      eq. (1) recomputed at that g (staircase_Fp, run_case l. 258);
  (b) Carter's product, approximate and exact conformal forms (dsop.carter);
  (c) condensed operator with Phi_O = 0, dsop.slotting_ratio, assembled
      smooth comparator recomputed at that g (definition (a) of Section II-E);
  (d) operator coupled to the slot cavities, dsop.slotting_ratio_cavity,
      same comparator.
Guards         tolerance 0 on the printed digits (k_C to 4 decimals, the
               comparator to 0.1 mT/A); a failing guard stops the script and
               nothing is published (no JSON):
  G1  at nominal g, (a) reproduces pos_real per position and on the mean;
  G2  at nominal g, (c) and (d) reproduce pos_op_33_16 and pos_cav per
      position and on the mean, and the line (33, 16) of Table S3 at phi = 0
      (tiling_sweep_inf_iron.json p1_33_16, prod_op_results.json cav_sweep);
  G3  the assembled smooth comparator at nominal g reproduces 463.8 mT/A;
  G4  Q reloaded from cavity_Q.pkl against Q recomputed (cavity.cavity_element)
      and against code/MEC_IM/cavity_nO16.mat: raw max |dQ| printed; the pass
      criterion is on the printed digits, k_C(d) at phi = 0 with the reloaded
      and with the recomputed Q equal to 4 decimals;
  G5  at the finest g, halving lc moves the twelve-position FE mean by less
      than 0.1 %; otherwise lc is halved again until it does, and the FE row
      of that g is reported at the validated (coarser) lc of the passing pair.
Run from code/python/ :  python prod_sweep_b0g.py [--checkpoint FILE]
Writes outputs/python/prod_sweep_b0g_results.json and prod_sweep_b0g_out.txt."""
import os, sys, time, json, pickle, argparse, datetime, traceback
sys.dont_write_bytecode = True
import numpy as np
import scipy.io as sio
from dsop import Machine, slotting_ratio, slotting_ratio_cavity, carter
import fem_slots as F
from cavity import cavity_element

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.normpath(os.path.join(HERE, '..', '..', 'outputs', 'python'))
MECO = os.path.normpath(os.path.join(HERE, '..', 'MEC_IM'))
NT, NO, NH, BASIS, LC = 33, 16, 8192, 'p1', 0.045
B0 = 2.0e-3                                   # both openings, m
RATIOS = [4.0, 12.0, 16.0, 2.0]               # order of computation; 2 last (largest gap, costliest mesh)
LC_MIN = 0.01                                 # G5 refinement cap (lc below this is not attempted)

ap = argparse.ArgumentParser(); ap.add_argument('--checkpoint', default=os.path.join(W, 'prod_sweep_b0g_checkpoint.json'))
ap.add_argument('--no-b0g2', action='store_true'); args = ap.parse_args()
T0 = time.time(); DUR = {}
def now(): return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
def say(s): print(s, flush=True)
def f4(x): return '%.4f' % x
def machine_at_gap(g):
    M = Machine(); M.g = g; M.Rr = M.Rs - g; M.Dr = 2 * M.Rr; M.X = np.log(M.Rs / M.Rr); M.Rm = 0.5 * (M.Rs + M.Rr); return M

# ---------------------------------------------------------------- archived nominal point
fem = json.load(open(os.path.join(W, 'prod_fem_results.json')))
op = json.load(open(os.path.join(W, 'prod_op_results.json')))
til = json.load(open(os.path.join(W, 'tiling_sweep_inf_iron.json')))
FR = sorted(float(k) for k in fem['pos_real'])            # the twelve positions, fractions of tau_r
def series(d): return {('%.6f' % f): float(d[k]) for f in FR for k in d if abs(float(k) - f) < 1e-9}
def stats(d): v = np.array(list(d.values())); return dict(mean=float(v.mean()), min=float(v.min()), max=float(v.max()))
M0 = Machine(); g0 = M0.g; tau_r = 2 * np.pi / M0.Nr
Q = pickle.load(open(os.path.join(W, 'cavity_Q.pkl'), 'rb'))[NO]

# ---------------------------------------------------------------- one block of the chain at a gap
def fe_positions(M, lc):
    out = {}; nodes = []; bsm = None; t = time.time()
    for f in FR:
        r = F.run_case(M, f * tau_r, lc_gap=lc); out['%.6f' % f] = float(r['kC']); nodes.append(int(r['nnodes'])); bsm = float(r['Bg1_smooth'])
        say('      FE  b0/g=%.3f lc=%.5f phi/tau_r=%.4f : kC %.6f  nodes %d  (%.0f s)' % (B0 / M.g, lc, f, r['kC'], r['nnodes'], time.time() - t))
    return dict(lc=lc, per_position=out, Bg1_smooth_closed=bsm, nodes_min=min(nodes), nodes_max=max(nodes), **stats(out))
def op_positions(M, cav):
    out = {}; bsm = None; ncol = None; t = time.time()
    for f in FR:
        r = slotting_ratio_cavity(M, NT, NO, NH, Q[0], Q[1], BASIS, phi=f * tau_r) if cav else slotting_ratio(M, NT, NO, NH, BASIS, phi=f * tau_r)
        out['%.6f' % f] = float(r['kC']); bsm = float(r['Bg1_smooth']); ncol = int(r['ncol'])
        say('      %s b0/g=%.3f phi/tau_r=%.4f : kC %.6f  (%.0f s)' % ('CAV' if cav else 'NEU', B0 / M.g, f, r['kC'], time.time() - t))
    return dict(per_position=out, Bg1_smooth_assembled=bsm, ncol=ncol, **stats(out))
def carter_at(M): c = carter(M); return {k: float(v) for k, v in c.items()}

# ---------------------------------------------------------------- guards
say('prod_sweep_b0g: start %s' % now()); say('  nominal g = %.6f mm, b0/g = %.4f ; R_s = %.6f mm' % (g0 * 1e3, B0 / g0, M0.Rs * 1e3))
G = {}
t = time.time()
say('  G1: FE formulation A at nominal g, lc %.3f, twelve positions ...' % LC)
fe0 = fe_positions(M0, LC); ref = series(fem['pos_real'])
G['G1'] = dict(recomputed=fe0, archived=dict(per_position=ref, **stats(ref)), Bg1_smooth_closed=fe0['Bg1_smooth_closed'], B_smooth_prod_fem=0.465164,
              per_position_equal_4dec=all(f4(fe0['per_position'][k]) == f4(ref[k]) for k in ref), mean_equal_4dec=f4(fe0['mean']) == f4(stats(ref)['mean']),
              phi0_recomputed=fe0['per_position']['%.6f' % 0.0], phi0_archived=ref['%.6f' % 0.0], phi0_finest_mesh_lc0028=fem['conv_real']['0.0_0.028'][0])
G['G1']['pass'] = G['G1']['per_position_equal_4dec'] and G['G1']['mean_equal_4dec']; DUR['G1'] = time.time() - t
say('  G1 %s : mean %.5f vs archived %.5f ; phi=0 %.5f vs %.5f (finest mesh lc 0.028: %.5f) ; closed comparator %.9f T/A vs prod_fem 0.465164' % (
    'PASS' if G['G1']['pass'] else 'FAIL', fe0['mean'], stats(ref)['mean'], G['G1']['phi0_recomputed'], G['G1']['phi0_archived'], G['G1']['phi0_finest_mesh_lc0028'], fe0['Bg1_smooth_closed']))
t = time.time(); say('  G2: operator Phi_O = 0 and cavities at nominal g, (33,16), p1, N_h 8192, twelve positions ...')
neu0 = op_positions(M0, False); cav0 = op_positions(M0, True); rn = series(op['pos_op_33_16']); rc = series(op['pos_cav'])
G['G2'] = dict(neu_recomputed=neu0, neu_archived=dict(per_position=rn, **stats(rn)), cav_recomputed=cav0, cav_archived=dict(per_position=rc, **stats(rc)),
              tableS3_33_16_neu=til['p1_33_16']['kC'], tableS3_33_16_cav=op['cav_sweep']['33_16'],
              neu_equal_4dec=all(f4(neu0['per_position'][k]) == f4(rn[k]) for k in rn) and f4(neu0['mean']) == f4(stats(rn)['mean']),
              cav_equal_4dec=all(f4(cav0['per_position'][k]) == f4(rc[k]) for k in rc) and f4(cav0['mean']) == f4(stats(rc)['mean']),
              S3_equal_4dec=f4(neu0['per_position']['%.6f' % 0.0]) == f4(til['p1_33_16']['kC']) and f4(cav0['per_position']['%.6f' % 0.0]) == f4(op['cav_sweep']['33_16']))
G['G2']['pass'] = G['G2']['neu_equal_4dec'] and G['G2']['cav_equal_4dec'] and G['G2']['S3_equal_4dec']; DUR['G2'] = time.time() - t
say('  G2 %s : Phi_O=0 mean %.5f vs %.5f ; cavities mean %.5f vs %.5f ; phi=0 %.5f / %.5f vs Table S3 %.5f / %.5f' % (
    'PASS' if G['G2']['pass'] else 'FAIL', neu0['mean'], stats(rn)['mean'], cav0['mean'], stats(rc)['mean'], neu0['per_position']['%.6f' % 0.0], cav0['per_position']['%.6f' % 0.0], til['p1_33_16']['kC'], op['cav_sweep']['33_16']))
bsm_asm = neu0['Bg1_smooth_assembled']
G['G3'] = dict(Bg1_smooth_assembled_T_per_A=bsm_asm, printed_mT_per_A='%.1f' % (bsm_asm * 1e3), expected='463.8', Bg1_smooth_closed_T_per_A=fe0['Bg1_smooth_closed'])
G['G3']['pass'] = G['G3']['printed_mT_per_A'] == '463.8'
say('  G3 %s : assembled smooth comparator (hat) %.6f mT/A -> printed %s (expected 463.8) ; closed form %.6f mT/A' % ('PASS' if G['G3']['pass'] else 'FAIL', bsm_asm * 1e3, G['G3']['printed_mT_per_A'], fe0['Bg1_smooth_closed'] * 1e3))
t = time.time(); say('  G4: cavity matrices n_O = 16 : reload vs recompute vs cavity_nO16.mat ...')
Qs2, _ = cavity_element('stator', NO, M0.L, lc_mouth=0.01); Qr2, _ = cavity_element('rotor', NO, M0.L, lc_mouth=0.01)
mat = sio.loadmat(os.path.join(MECO, 'cavity_nO16.mat'))
k_reload = cav0['per_position']['%.6f' % 0.0]; k_recomp = float(slotting_ratio_cavity(M0, NT, NO, NH, Qs2, Qr2, BASIS, phi=0.0)['kC'])
G['G4'] = dict(dQs_reload_vs_recompute=float(np.max(np.abs(Q[0] - Qs2))), dQr_reload_vs_recompute=float(np.max(np.abs(Q[1] - Qr2))),
              dQs_reload_vs_mat=float(np.max(np.abs(Q[0] - mat['Qs']))), dQr_reload_vs_mat=float(np.max(np.abs(Q[1] - mat['Qr']))),
              maxQs=float(np.max(np.abs(Q[0]))), maxQr=float(np.max(np.abs(Q[1]))), kC_phi0_reloaded=k_reload, kC_phi0_recomputed=k_recomp)
G['G4']['pass'] = f4(k_reload) == f4(k_recomp); DUR['G4'] = time.time() - t
say('  G4 %s : max|dQ| reload vs recompute stator %.2e rotor %.2e (max|Q| %.2e / %.2e) ; reload vs cavity_nO16.mat %.2e / %.2e ; kC(d) phi=0 reloaded %.6f recomputed %.6f' % (
    'PASS' if G['G4']['pass'] else 'FAIL', G['G4']['dQs_reload_vs_recompute'], G['G4']['dQr_reload_vs_recompute'], G['G4']['maxQs'], G['G4']['maxQr'], G['G4']['dQs_reload_vs_mat'], G['G4']['dQr_reload_vs_mat'], k_reload, k_recomp))

def write_transcript(lines, results=None):
    with open(os.path.join(W, 'prod_sweep_b0g_out.txt'), 'w', encoding='utf-8') as fh: fh.write('\n'.join(lines) + '\n')
    if results is not None: json.dump(results, open(os.path.join(W, 'prod_sweep_b0g_results.json'), 'w'), indent=1)
def guard_lines():
    L = ['prod_sweep_b0g.py -- run of %s ; b0/g sweep on the linear chain of Section IV (see the docstring for the settings)' % now(),
         'tiling (n_T, n_O) = (%d, %d) ; basis %s ; N_h = %d ; FE production mesh lc_gap = %.3f mm ; twelve rotor positions i/12 tau_r ; b_0 = %.3f mm on both bores' % (NT, NO, BASIS, NH, LC, B0 * 1e3),
         'nominal point read from the archive: g = %.6f mm, b0/g = %.4f (prod_fem_results.json pos_real ; prod_op_results.json pos_op_33_16, pos_cav)' % (g0 * 1e3, B0 / g0),
         'G1 %s  FE (A) at nominal g, lc %.3f: recomputed mean %.6f (archived %.6f), phi=0 %.6f (archived %.6f ; finest mesh lc 0.028: %.6f), all twelve positions equal to 4 decimals: %s ; closed-form comparator %.9f T/A (prod_fem.py hard-codes 0.465164)' % (
             'PASS' if G['G1']['pass'] else 'FAIL', LC, fe0['mean'], stats(ref)['mean'], G['G1']['phi0_recomputed'], G['G1']['phi0_archived'], G['G1']['phi0_finest_mesh_lc0028'], G['G1']['per_position_equal_4dec'], fe0['Bg1_smooth_closed']),
         'G2 %s  operator at nominal g: Phi_O=0 mean %.6f (archived pos_op_33_16 %.6f), cavities mean %.6f (archived pos_cav %.6f), phi=0 %.6f / %.6f against Table S3 line (33,16) %.6f / %.6f ; all positions equal to 4 decimals: %s / %s' % (
             'PASS' if G['G2']['pass'] else 'FAIL', neu0['mean'], stats(rn)['mean'], cav0['mean'], stats(rc)['mean'], neu0['per_position']['%.6f' % 0.0], cav0['per_position']['%.6f' % 0.0], til['p1_33_16']['kC'], op['cav_sweep']['33_16'], G['G2']['neu_equal_4dec'], G['G2']['cav_equal_4dec']),
         'G3 %s  assembled smooth comparator (hat basis, position-independent) at nominal g = %.6f mT/A, printed %s (expected 463.8) ; closed-form comparator of the FE chain %.6f mT/A' % ('PASS' if G['G3']['pass'] else 'FAIL', bsm_asm * 1e3, G['G3']['printed_mT_per_A'], fe0['Bg1_smooth_closed'] * 1e3),
         'G4 %s  Q (n_O = 16) reloaded from cavity_Q.pkl: max|dQ| against a recompute (cavity.cavity_element, lc_mouth 0.01, gmsh) stator %.3e rotor %.3e of max|Q| %.3e / %.3e ; against code/MEC_IM/cavity_nO16.mat %.3e / %.3e ; k_C (cavities, phi = 0) with the reloaded Q %.6f, with the recomputed Q %.6f -> equal to 4 decimals: %s' % (
             'PASS' if G['G4']['pass'] else 'FAIL', G['G4']['dQs_reload_vs_recompute'], G['G4']['dQr_reload_vs_recompute'], G['G4']['maxQs'], G['G4']['maxQr'], G['G4']['dQs_reload_vs_mat'], G['G4']['dQr_reload_vs_mat'], k_reload, k_recomp, G['G4']['pass'])]
    return L
if not all(G[k]['pass'] for k in ['G1', 'G2', 'G3', 'G4']):
    write_transcript(guard_lines() + ['A GUARD FAILED -- nothing published (no results JSON).'])
    say('GUARD FAILURE -- stop'); sys.exit(1)

# ---------------------------------------------------------------- the sweep
ck = json.load(open(args.checkpoint)) if os.path.isfile(args.checkpoint) else {}
S = ck.get('sweep', {})
for r in ([4.0, 12.0, 16.0] + ([] if args.no_b0g2 else [2.0])):
    key = '%g' % r
    if key in S and S[key].get('complete'): say('  b0/g = %g : from checkpoint' % r); continue
    g = B0 / r; M = machine_at_gap(g); t = time.time()
    say('  b0/g = %g  (g = %.6f mm, R_r = %.6f mm, X = %.6e, 1/X = %.1f) ...' % (r, g * 1e3, M.Rr * 1e3, M.X, 1 / M.X))
    blk = S.get(key, {}); blk.update(b0_over_g=r, g_mm=g * 1e3, Rr_mm=M.Rr * 1e3, X=M.X)
    if 'carter' not in blk: blk['carter'] = carter_at(M)
    if 'fe' not in blk: blk['fe'] = fe_positions(M, LC); S[key] = blk; json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
    if 'op_neu' not in blk: blk['op_neu'] = op_positions(M, False); S[key] = blk; json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
    if 'op_cav' not in blk: blk['op_cav'] = op_positions(M, True); S[key] = blk; json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
    blk['complete'] = True; blk['duration_s'] = blk.get('duration_s', 0) + time.time() - t; S[key] = blk; DUR['b0g_%s' % key] = blk['duration_s']
    json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
    say('  b0/g = %g done : FE mean %.6f [%.6f, %.6f] | Carter %.6f / %.6f | Phi_O=0 %.6f | cav %.6f [%.6f, %.6f]  (%.0f s)' % (
        r, blk['fe']['mean'], blk['fe']['min'], blk['fe']['max'], blk['carter']['k'], blk['carter']['k_exact'], blk['op_neu']['mean'], blk['op_cav']['mean'], blk['op_cav']['min'], blk['op_cav']['max'], blk['duration_s']))

# ---------------------------------------------------------------- G5 at the finest g
rmax = max(float(k) for k in S); key = '%g' % rmax; M = machine_at_gap(B0 / rmax); t = time.time()
G5 = S[key].get('G5', dict(sequence=[]))
if not G5.get('done'):
    seq = G5['sequence'] if G5['sequence'] else [dict(lc=LC, **{k: S[key]['fe'][k] for k in ['mean', 'min', 'max', 'nodes_min', 'nodes_max']}, per_position=S[key]['fe']['per_position'])]
    lc = seq[-1]['lc']; ok = False
    while True:
        if len(seq) >= 2 and abs(100 * (seq[-1]['mean'] / seq[-2]['mean'] - 1)) < 0.1: ok = True; break
        if lc / 2 < LC_MIN: break
        lc = lc / 2; say('  G5: b0/g = %g, lc = %.5f ...' % (rmax, lc)); fe = fe_positions(M, lc)
        seq.append(dict(lc=lc, per_position=fe['per_position'], **{k: fe[k] for k in ['mean', 'min', 'max', 'nodes_min', 'nodes_max']}))
        G5['sequence'] = seq; S[key]['G5'] = G5; json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
    G5.update(done=True, passed=ok, shifts_pct=[100 * (seq[i + 1]['mean'] / seq[i]['mean'] - 1) for i in range(len(seq) - 1)], validated_lc=(seq[-2]['lc'] if ok else None), sequence=seq)
    S[key]['G5'] = G5; json.dump(dict(sweep=S), open(args.checkpoint, 'w'), indent=1)
DUR['G5'] = DUR.get('G5', 0) + time.time() - t
G['G5'] = dict(b0_over_g=rmax, passed=G5['passed'], validated_lc=G5['validated_lc'], shifts_pct=G5['shifts_pct'], lcs=[s['lc'] for s in G5['sequence']], means=[s['mean'] for s in G5['sequence']], nodes=[(s['nodes_min'], s['nodes_max']) for s in G5['sequence']])
say('  G5 %s : lc %s -> means %s ; shifts %s %% ; validated lc %s' % ('PASS' if G5['passed'] else 'FAIL', G['G5']['lcs'], ['%.6f' % m for m in G['G5']['means']], ['%+.4f' % s for s in G['G5']['shifts_pct']], G5['validated_lc']))
if G5['passed'] and G5['validated_lc'] != LC:
    v = [s for s in G5['sequence'] if s['lc'] == G5['validated_lc']][0]
    S[key]['fe_production_lc'] = S[key]['fe']; S[key]['fe'] = dict(lc=v['lc'], per_position=v['per_position'], Bg1_smooth_closed=S[key]['fe_production_lc']['Bg1_smooth_closed'], nodes_min=v['nodes_min'], nodes_max=v['nodes_max'], **stats(v['per_position']))
    say('  G5: FE row of b0/g = %g reported at the validated lc = %.5f (production lc %.3f kept under fe_production_lc)' % (rmax, v['lc'], LC))

# ---------------------------------------------------------------- assemble the results and the transcript
nominal = dict(b0_over_g=B0 / g0, g_mm=g0 * 1e3, Rr_mm=M0.Rr * 1e3, X=M0.X, source='prod_fem_results.json pos_real ; prod_op_results.json pos_op_33_16, pos_cav ; dsop.carter',
               carter=carter_at(M0), fe=dict(lc=LC, per_position=ref, Bg1_smooth_closed=fe0['Bg1_smooth_closed'], **stats(ref)),
               op_neu=dict(per_position=rn, Bg1_smooth_assembled=bsm_asm, **stats(rn)), op_cav=dict(per_position=rc, Bg1_smooth_assembled=bsm_asm, **stats(rc)))
ALL = dict(S); ALL['%g' % (B0 / g0)] = nominal
order = sorted(ALL, key=float)
def dev(a, b): return 100 * (a / b - 1)
L = guard_lines()
if not G5['passed']:
    L.append('G5 FAIL  at b0/g = %g the FE mean still moves by %s %% under halving of lc (sequence lc %s, means %s) ; the refinement cap lc >= %.3f mm was reached -- nothing published.' % (rmax, ['%+.4f' % s for s in G['G5']['shifts_pct']], G['G5']['lcs'], ['%.6f' % m for m in G['G5']['means']], LC_MIN))
    write_transcript(L); say('G5 FAILURE -- stop'); sys.exit(1)
L.append('G5 PASS  at b0/g = %g (g = %.6f mm): FE twelve-position mean at lc %s = %s ; shifts under halving %s %% ; validated lc = %.5f (< 0.1 %% criterion) ; nodes %s ; FE row of that g reported at lc %.5f' % (
    rmax, B0 / rmax * 1e3, ' / '.join('%.5f' % x for x in G['G5']['lcs']), ' / '.join('%.6f' % m for m in G['G5']['means']), ' / '.join('%+.4f' % s for s in G['G5']['shifts_pct']), G5['validated_lc'], G['G5']['nodes'], S[key]['fe']['lc']))
DUR['total'] = time.time() - T0
L.append('durations (s): ' + ', '.join('%s %.0f' % (k, v) for k, v in DUR.items()))
L.append('')
L.append('TABLE 1 -- values per b0/g (twelve-position statistics ; FE at the lc shown ; operator (33,16) p1 N_h 8192 ; Carter from the dimensions ; comparators recomputed at each g)')
L.append('%8s %10s %8s | %10s %10s %10s | %10s %10s | %10s %10s %10s | %10s %10s %10s | %12s %12s' % ('b0/g', 'g [mm]', 'lc [mm]', 'FE mean', 'FE min', 'FE max', 'Carter apx', 'Carter ex', 'PhiO=0 mean', 'min', 'max', 'cav mean', 'cav min', 'cav max', 'comp closed', 'comp asm'))
for k in order:
    b = ALL[k]
    L.append('%8.4f %10.6f %8.5f | %10.6f %10.6f %10.6f | %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %12.6f %12.6f' % (
        b['b0_over_g'], b['g_mm'], b['fe']['lc'], b['fe']['mean'], b['fe']['min'], b['fe']['max'], b['carter']['k'], b['carter']['k_exact'], b['op_neu']['mean'], b['op_neu']['min'], b['op_neu']['max'],
        b['op_cav']['mean'], b['op_cav']['min'], b['op_cav']['max'], b['fe']['Bg1_smooth_closed'] * 1e3, b['op_neu']['Bg1_smooth_assembled'] * 1e3))
L.append('')
L.append('TABLE 2 -- deviations from the FE twelve-position mean (%) ; FE modulation = (min - mean)/mean and (max - mean)/mean over the twelve positions')
L.append('%8s %10s | %12s %12s %12s %12s | %12s %12s | %12s %12s' % ('b0/g', 'g [mm]', 'Carter apx', 'Carter ex', 'PhiO=0', 'cavities', 'FE mod -', 'FE mod +', 'cav mod -', 'cav mod +'))
for k in order:
    b = ALL[k]; m = b['fe']['mean']
    L.append('%8.4f %10.6f | %+12.3f %+12.3f %+12.3f %+12.3f | %+12.3f %+12.3f | %+12.3f %+12.3f' % (
        b['b0_over_g'], b['g_mm'], dev(b['carter']['k'], m), dev(b['carter']['k_exact'], m), dev(b['op_neu']['mean'], m), dev(b['op_cav']['mean'], m),
        dev(b['fe']['min'], m), dev(b['fe']['max'], m), dev(b['op_cav']['min'], b['op_cav']['mean']), dev(b['op_cav']['max'], b['op_cav']['mean'])))
L.append('')
L.append('TABLE 3 -- per position (phi/tau_r) : FE (A) k_C, operator Phi_O = 0, operator + cavities, at each b0/g')
for k in order:
    b = ALL[k]; L.append('  b0/g = %.4f (g = %.6f mm, lc %.5f) :' % (b['b0_over_g'], b['g_mm'], b['fe']['lc']))
    for f in FR:
        kk = '%.6f' % f; L.append('    %.6f  FE %.6f  PhiO=0 %.6f  cav %.6f' % (f, b['fe']['per_position'][kk], b['op_neu']['per_position'][kk], b['op_cav']['per_position'][kk]))
if G5['validated_lc'] != LC:
    L.append(''); L.append('TABLE 4 -- G5 refinement at b0/g = %g : lc / mean / min / max / nodes' % rmax)
    for s in G5['sequence']: L.append('    lc %.5f  mean %.6f  min %.6f  max %.6f  nodes %d-%d' % (s['lc'], s['mean'], s['min'], s['max'], s['nodes_min'], s['nodes_max']))
L.append(''); L.append('end %s' % now())
results = dict(date=now(), settings=dict(nT=NT, nO=NO, Nh=NH, basis=BASIS, lc_gap=LC, b0_m=B0, positions=FR, lc_min_cap=LC_MIN), durations_s=DUR, guards=G, nominal=nominal, sweep=S)
write_transcript(L, results)
say('DONE %s (%.0f s)' % (now(), DUR['total']))
