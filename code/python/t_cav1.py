import numpy as np, time
from dsop import *
from cavity import cavity_element
M=Machine(); tau_r=2*np.pi/M.Nr
for nO in [4]:
    Qs,_=cavity_element('stator',nO,M.L,lc_mouth=0.02); Qr,_=cavity_element('rotor',nO,M.L,lc_mouth=0.02)
    for nT in [9,17,33]:
        t=time.time(); r=slotting_ratio_cavity(M,nT,nO,8192,Qs,Qr,'p1'); print(f"cavity-augmented (nT={nT},nO={nO}) phi=0: kC={r['kC']:.5f} Bg1_slot={r['Bg1_slot']:.5f} ({time.time()-t:.0f}s)")
for nO in [2,8]:
    Qs,_=cavity_element('stator',nO,M.L,lc_mouth=0.02); Qr,_=cavity_element('rotor',nO,M.L,lc_mouth=0.02)
    for nT in [17,33]:
        r=slotting_ratio_cavity(M,nT,nO,8192,Qs,Qr,'p1'); print(f"cavity-augmented (nT={nT},nO={nO}) phi=0: kC={r['kC']:.5f}")
print("FEM reference at phi=0 (real slots, lc=0.035): 1.29406")
