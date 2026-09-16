"""Production: cavity-augmented operator sweeps (infinite iron) and field comparison."""
import numpy as np, time, json, pickle
from dsop import *
from cavity import cavity_element
import fem_slots as F
M=Machine(); tau_r=2*np.pi/M.Nr
t0=time.time(); out={}
Q={}
for nO in [4,8,16,32]:
    lc=0.02 if nO<=8 else (0.01 if nO==16 else 0.005)
    Q[nO]=(cavity_element('stator',nO,M.L,lc_mouth=lc)[0], cavity_element('rotor',nO,M.L,lc_mouth=lc)[0])
    print("cavity element nO",nO,"done",f"({time.time()-t0:.0f}s)",flush=True)
pickle.dump(Q,open('cavity_Q.pkl','wb'))
# (a) tiling sweep, phi=0
out['cav_sweep']={}
for nT in [17,33,65]:
    for nO in [4,8,16,32]:
        if nT==65 and nO==32: continue
        r=slotting_ratio_cavity(M,nT,nO,8192,Q[nO][0],Q[nO][1],'p1')
        out['cav_sweep'][f"{nT}_{nO}"]=r['kC']
        print(f"CAV sweep nT={nT} nO={nO}: kC={r['kC']:.5f} cols={r['ncol']} ({time.time()-t0:.0f}s)",flush=True)
# p1a basis check at (33,16)
r=slotting_ratio_cavity(M,33,16,8192,Q[16][0],Q[16][1],'p1a'); out['cav_33_16_p1a']=r['kC']
print(f"CAV (33,16) p1a basis: kC={r['kC']:.5f}",flush=True)
# (b) position average, cavity (33,16) and Phi_O=0 (33,16) and (17,4)
out['pos_cav']={}; out['pos_op_17_4']={}; out['pos_op_33_16']={}
for i in range(12):
    frac=i/12; phi=frac*tau_r
    out['pos_cav'][frac]=slotting_ratio_cavity(M,33,16,8192,Q[16][0],Q[16][1],'p1',phi=phi)['kC']
    out['pos_op_17_4'][frac]=slotting_ratio(M,17,4,8192,'p1',phi=phi)['kC']
    out['pos_op_33_16'][frac]=slotting_ratio(M,33,16,8192,'p1',phi=phi)['kC']
    print(f"POS {frac:.4f}: cav(33,16) {out['pos_cav'][frac]:.5f} | op(17,4) {out['pos_op_17_4'][frac]:.5f} | op(33,16) {out['pos_op_33_16'][frac]:.5f} ({time.time()-t0:.0f}s)",flush=True)
for k in ['pos_cav','pos_op_17_4','pos_op_33_16']:
    print(k,"mean",np.mean(list(out[k].values())))
json.dump({k:(v if not isinstance(v,dict) else {str(kk):vv for kk,vv in v.items()}) for k,v in out.items()},open('prod_op_results.json','w'),indent=1)
# (c) mid-gap field waveforms at phi=0 over 2 stator pitches: FEM real, operator Phi_O=0 (33,16), cavity (33,16)
i3=np.sqrt(2)*np.array([1.0,-0.5,-0.5]); AT=M.slot_ampere_turns(i3)
xyz,tri=F.build_air_mesh(M,0.0,lc_gap=0.045)
K,area=F.assemble_P1(xyz,tri,np.ones(len(tri)))
ids=F.slot_body_ids(M,xyz,tri); ba=np.array([area[ids==k].sum() for k in range(M.Ns)])
J=np.zeros(len(tri)); m=ids>=0; J[m]=AT[ids[m]]/ba[ids[m]]
A,area=F.solve_Az(M,xyz,tri,J)
bg1,Bn,a,th=F.midgap_fundamental_A(M,xyz,tri,A,nsamp=16384)
Br_fem=np.gradient(a,th)/M.Rm
def field_op(r):
    proj,u,cond=r['proj'],r['u'],r['cond']
    Brc,Brs=gap_field_harmonics(M,proj,u,cond)
    n=np.arange(1,len(Brc)+1)
    thq=th
    return (np.cos(np.outer(thq,n))@Brc + np.sin(np.outer(thq,n))@Brs)
r1=slotting_ratio(M,33,16,8192,'p1',return_all=True); Br_op=field_op(r1)
r2=slotting_ratio_cavity(M,33,16,8192,Q[16][0],Q[16][1],'p1',return_all=True); Br_cav=field_op(r2)
np.savez('field_waveforms.npz',th=th,Br_fem=Br_fem,Br_op=Br_op,Br_cav=Br_cav,Bn_fem=Bn)
print("waveforms saved; fundamentals: FEM %.5f op %.5f cav %.5f"%(bg1,r1['Bg1_slot'],r2['Bg1_slot']))
print("DONE")
