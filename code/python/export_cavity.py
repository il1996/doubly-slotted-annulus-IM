import numpy as np, scipy.io as sio
from cavity import cavity_element
from dsop import Machine
M=Machine()
for nO in [2,4,8,16]:
    lc = 0.02 if nO<=8 else 0.01
    Qs,_=cavity_element('stator',nO,M.L,lc_mouth=lc); Qr,_=cavity_element('rotor',nO,M.L,lc_mouth=lc)
    sio.savemat(f'../mec/MEC_IM/cavity_nO{nO}.mat',{'Qs':Qs,'Qr':Qr,'nO':nO})
    print(nO, "saved; Qs diag[:3]", np.diag(Qs)[:3], "Qr wall-wall", Qr[nO,nO])
