"""Production runs: FEM (real slots) and Neumann-opening FEM vs rotor position and mesh."""
import numpy as np, time, json
import fem_slots as F, fem_annulus as FA
from dsop import Machine, MU0
M=Machine(); tau_r=2*np.pi/M.Nr
i3=np.sqrt(2)*np.array([1.0,-0.5,-0.5]); AT=M.slot_ampere_turns(i3); Us=M.tooth_mmf(i3)
B_smooth=0.465164   # analytic staircase fundamental at mid-gap, I0=1 A
def Jgap(xyz,tri,area):
    c=xyz[tri].mean(axis=1); r=np.hypot(c[:,0],c[:,1]); ang=np.arctan2(c[:,1],c[:,0])
    J=np.zeros(len(tri)); m=r>M.Rs-0.04e-3
    k=np.round((ang/(2*np.pi))*M.Ns-0.5).astype(int)%M.Ns
    dang=((ang-(2*np.pi*(k+0.5)/M.Ns))+np.pi)%(2*np.pi)-np.pi
    m&= np.abs(dang)<(1.0e-3/M.Rs); sel=np.where(m)[0]
    for kk in range(M.Ns):
        e=sel[k[sel]==kk]; J[e]=AT[kk]/area[e].sum()
    return J
def fem_real(stator_slots, rotor_slots, phi, lc):
    xyz,tri=F.build_air_mesh(M,phi,lc_gap=lc,stator_slots=stator_slots,rotor_slots=rotor_slots)
    K,area=F.assemble_P1(xyz,tri,np.ones(len(tri)))
    if stator_slots:
        ids=F.slot_body_ids(M,xyz,tri); ba=np.array([area[ids==k].sum() for k in range(M.Ns)])
        J=np.zeros(len(tri)); m=ids>=0; J[m]=AT[ids[m]]/ba[ids[m]]
    else: J=Jgap(xyz,tri,area)
    A,area=F.solve_Az(M,xyz,tri,J)
    bg1,Bn,a,th=F.midgap_fundamental_A(M,xyz,tri,A)
    return B_smooth/bg1, len(xyz), Bn
def fem_neumann(phi, lc):
    xyz,tri,arcs,fs,fr=FA.build_annulus_mesh(M,phi,lc)
    phi_,R=FA.solve_scalar(M,xyz,tri,arcs,fs,fr,Us,0.0,opening='neumann')
    Bn=FA.midgap_harmonics_scalar(M,xyz,tri,phi_)
    return B_smooth/Bn[M.p], len(xyz)
out={}
t0=time.time()
# mesh convergence, phi=0 and phi=3/8 tau_r
out['conv_real']={}
for frac in [0.0,0.375]:
    for lc in [0.08,0.06,0.045,0.035,0.028]:
        k,nn,_=fem_real(True,True,frac*tau_r,lc); out['conv_real'][f"{frac}_{lc}"]=(k,nn)
        print(f"REAL conv phi={frac} lc={lc}: kC={k:.5f} nodes={nn} ({time.time()-t0:.0f}s)",flush=True)
out['conv_neu']={}
for frac in [0.0,0.375]:
    for lc in [0.08,0.05,0.035,0.025]:
        k,nn=fem_neumann(frac*tau_r,lc); out['conv_neu'][f"{frac}_{lc}"]=(k,nn)
        print(f"NEUMANN conv phi={frac} lc={lc}: kC={k:.5f} nodes={nn} ({time.time()-t0:.0f}s)",flush=True)
# single surfaces
out['single']={}
for lc in [0.06,0.045,0.035]:
    ks,_,_=fem_real(True,False,0.0,lc); kr,_,_=fem_real(False,True,0.0,lc); out['single'][lc]=(ks,kr)
    print(f"SINGLE lc={lc}: ks={ks:.5f} kr={kr:.5f} product={ks*kr:.5f} ({time.time()-t0:.0f}s)",flush=True)
# rotor position sweep, 12 positions
out['pos_real']={}; out['pos_neu']={}
for i in range(12):
    frac=i/12
    k,nn,Bn=fem_real(True,True,frac*tau_r,0.045); out['pos_real'][frac]=k
    kn,_=fem_neumann(frac*tau_r,0.035); out['pos_neu'][frac]=kn
    print(f"POS {frac:.4f}: real kC={k:.5f}  neumann kC={kn:.5f} ({time.time()-t0:.0f}s)",flush=True)
print("mean real:",np.mean(list(out['pos_real'].values())),"mean neumann:",np.mean(list(out['pos_neu'].values())))
json.dump({k:{str(kk):vv for kk,vv in v.items()} for k,v in out.items()},open('prod_fem_results.json','w'),indent=1)
print("DONE")
