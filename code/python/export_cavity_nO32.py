# -*- coding: utf-8 -*-
"""Cavity matrices for n_O = 32, same construction as export_cavity.py (cavity.cavity_element,
lc_mouth = 0.005 as in prod_op.py for n_O = 32), written to code/MEC_IM/cavity_nO32.mat for
RUN_Z9_BASIS_P1A.  Run from code/python/.  Checked against cavity_Q.pkl[32] (prod_op.py):
Qs identical, Qr within 2.7e-10 absolute (3.7e-4 of max|Q|; gmsh rotor mesh)."""
import numpy as np, scipy.io as sio
from cavity import cavity_element
from dsop import Machine
M = Machine(); nO = 32; lc = 0.005
Qs, _ = cavity_element('stator', nO, M.L, lc_mouth=lc); Qr, _ = cavity_element('rotor', nO, M.L, lc_mouth=lc)
sio.savemat('../MEC_IM/cavity_nO32.mat', {'Qs': Qs, 'Qr': Qr, 'nO': nO})
print(nO, "saved; Qs diag[:3]", np.diag(Qs)[:3], "Qr wall-wall", Qr[nO, nO])
