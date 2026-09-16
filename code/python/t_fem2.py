import numpy as np, time, json
import fem_slots as F
from dsop import Machine, MU0
M=Machine()
i3=np.sqrt(2)*np.array([1.0,-0.5,-0.5]); AT=M.slot_ampere_turns(i3)
B_smooth=0.465164
def Jgap(xyz,tri,area):
    c=xyz[tri].mean(axis=1); r=np.hypot(c[:,0],c[:,1]); ang=np.arctan2(c[:,1],c[:,0])
    J=np.zeros(len(tri)); m=r>M.Rs-0.04e-3
    k=np.round((ang/(2*np.pi))*M.Ns-0.5).astype(int)%M.Ns
    dang=((ang-(2*np.pi*(k+0.5)/M.Ns))+np.pi)%(2*np.pi)-np.pi
    m&= np.abs(dang)<(1.0e-3/M.Rs)
    sel=np.where(m)[0]
    for kk in range(M.Ns):
        e=sel[k[sel]==kk]; J[e]=AT[kk]/area[e].sum()
    return J
def run(label, stator_slots, rotor_slots, phi=0.0, lc=0.06):
    t=time.time()
    xyz,tri=F.build_air_mesh(M,phi,lc_gap=lc,stator_slots=stator_slots,rotor_slots=rotor_slots)
    K,area=F.assemble_P1(xyz,tri,np.ones(len(tri)))
    if stator_slots:
        ids=F.slot_body_ids(M,xyz,tri); ba=np.array([area[ids==k].sum() for k in range(M.Ns)])
        J=np.zeros(len(tri)); m=ids>=0; J[m]=AT[ids[m]]/ba[ids[m]]
    else:
        J=Jgap(xyz,tri,area)
    A,area=F.solve_Az(M,xyz,tri,J)
    bg1,Bn,a,th=F.midgap_fundamental_A(M,xyz,tri,A)
    print(f"{label:40s} lc={lc} phi={phi:.4f} nodes={len(xyz):7d}  Bg1={bg1:.6f}  kC={B_smooth/bg1:.5f}  ({time.time()-t:.0f}s)", flush=True)
    return B_smooth/bg1
res={}
res['stator_only']={lc:run("stator slots only",True,False,lc=lc) for lc in [0.08,0.06,0.045]}
res['rotor_only']={lc:run("rotor slots only",False,True,lc=lc) for lc in [0.08,0.06,0.045]}
res['both']={lc:run("both slotted",True,True,lc=lc) for lc in [0.08,0.06,0.045,0.035]}
tau_r=2*np.pi/M.Nr
res['both_pos']={}
for frac in [0.125,0.25,0.375,0.5,0.625,0.75,0.875]:
    res['both_pos'][frac]=run("both slotted, rotor position",True,True,phi=frac*tau_r,lc=0.06)
json.dump({k:{str(kk):vv for kk,vv in v.items()} for k,v in res.items()},open('fem_kc_results.json','w'),indent=1)
