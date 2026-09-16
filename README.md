# Reproduction archive — Paper II

*A Doubly Slotted Annulus Operator for Magnetic Equivalent Circuits:
Trace Conformity, the Slot-Opening Condition, and a Finite-Element
Verification of the Slotting Ratio on an 18.5 kW Cage Induction Machine.*

This archive supports **one paper only**. It is self-contained: the
MATLAB package, the scripts, the Python chain, the finite-element exports
and the transcripts from which every published value originates are here.

## 1. What is here

| directory | contents |
|---|---|
| `code/MEC_IM/+mec/` | the MATLAB package of the performance network (38 functions: geometry, winding, iron, cage, leakage, losses, air-gap operators, Newton solver, equivalent circuit) |
| `code/MEC_IM/` | the `RUN_*.m` scripts that produced every network value, the cavity matrices `cavity_nO{2,4,8,16}.mat` and the bar-current source `cavity_src_nO16.mat` |
| `code/python/` | the independent Python chain: the linear operator rebuilt from the equations (`dsop.py`), the two finite-element references (`fem_slots.py`, `fem_annulus.py`), the slot-cavity element (`cavity.py`, `cavity_graded.py`), the graded tiling (`t4_graded.py`), the production scripts (`prod_fem.py`, `prod_op.py`), the audit of the ANSYS exports (`fea_audit.py`, `sweep_audit.py`, `fea_power.py`, `fea_conv*.py`), the checks (`t_*.py`) and the figures (`make_figures_v2.py`) |
| `code/ansys_noload/` | the COM script that builds and solves the numerical no-load test on a copy of the finite-element project, and its post-processing |
| `code/article/figures/` | the six figures of the paper, PDF and PNG |
| `code/SET_REFERENCE_PATH.m` | repoints the older scripts, archived unchanged, to `reference/` |
| `outputs/MEC_IM/` | MATLAB transcripts (`diary`) at full precision, and the `.mat` files the figures read |
| `outputs/python/` | results of the Python chain (`*.json`, `field_waveforms.npz`, `cavity_Q.pkl`, `t4_graded_results.json`) |
| `outputs/ansys_noload/` | results of the numerical no-load test |
| `reference/ANSYS_18_5kW/` | the finite-element exports used as the reference, including the no-load voltage sweep |
| `notes/` | the manifest binding each published quantity to the chain that produced it, the delivery notes of the verification blocks and the independent audit |
| `tools/octave_compat/` | shims used only to run the chain under GNU Octave; not needed under MATLAB |

## 2. The one-minute test

Do this before anything else. The block below reproduces Table 8 of the
paper, the five closures of the air gap at the no-load and rated points,
and prints the reference values it reads from `reference/`. Under MATLAB
R2024a it takes about 45 s.

```matlab
cd code/MEC_IM
RUN_Z1_CAVITY
```

Expected, to all printed digits: `Xm0 71.990 | I0 8.612 | T 116.07`
(Carter), `61.015 | 10.323 | 114.31` (operator, Φ_O = 0), `67.713 |
9.542 | 115.15`, `69.209 | 9.407 | 115.29` and `70.029 | 9.153 | 115.53`
(operator with slot cavities at (17, 4), (33, 8) and (33, 16)).

If it passes, the archive works on your machine. The transcript it writes,
`Z1_cavity_out.txt`, is the one archived in `outputs/MEC_IM/`.

## 3. The scripts of the paper

| script | what it produces |
|---|---|
| `RUN_Z1_CAVITY.m` | Table 8 (five closures, no-load and rated points) |
| `RUN_Z6_COUPLED_CONV_CAV.m` | Table 10 (tiling sweep of the coupled solution, both closures) |
| `RUN_Z2_FIELDS.m`, `RUN_Z4_FIELDS_CAV.m` | mid-gap waveforms of Fig. 6, both closures |
| `RUN_M11_FIELD_ERR.m`, `RUN_M11B_ROTORPOS.m`, `RUN_M11C_MATCHED.m`, `RUN_Z4B_FIELD_ERR_CAV.m` | Table 11 |
| `RUN_B10_B1_SKEWOFF.m`, `RUN_Z5_SWEEP_CAV.m` | slip characteristics of Fig. 5, both closures; Table 9; standstill |
| `RUN_Z3_LEAKAGE.m` | leakage sensitivity of Section 5.3 |
| `RUN_Z7_NOLOAD_NET.m` | network side of the numerical no-load test (Section 5.3) |
| `RUN_Z8_CAVITY_LOAD.m` | rotor cavity map with the bar current (Section 6.3) |
| `RUN_R8_TABLE2.m`, `RUN_INVARIANTS.m` | Table 2 |
| `RUN_B2_KC.m` | Table 3 |
| `code/python/prod_fem.py` | Table 4, Fig. 2 (finite-element references) |
| `code/python/prod_op.py`, `t_op*.py` | Table 5 (uniform tilings), Fig. 2, Fig. 4 |
| `code/python/t4_graded.py` | Table 5 (graded tiling), Fig. 3 |
| `code/ansys_noload/noload_sweep_com.py`, `noload_postprocess.py` | the numerical no-load test |

Every script that reads the finite-element reference does so through the
relative path `../../reference/ANSYS_18_5kW`; the older scripts archived
unchanged carry the absolute path of the machine of origin and are
repointed by `code/SET_REFERENCE_PATH.m`.

## 4. How to read a published number

Every value in the paper comes from a transcript in `outputs/`. The rule
the project follows is **one quantity, one chain**: when two scripts
produced the same quantity, one of them was archived rather than kept.
`notes/MANIFEST.md` binds each published quantity to the chain that
produced it.

Each transcript opens with the configuration that produced it: tiling
$n_T$ and $n_O$, truncation $N_h$, surface basis, slip, skew setting,
and the path of the finite-element reference used.

## 5. Environment

MATLAB R2021b or later, no toolbox beyond the base product (verified under
R2024a). Python 3 with numpy, scipy, matplotlib and gmsh for
`code/python/`; pywin32 and ANSYS Electronics Desktop 2023 R1 for
`code/ansys_noload/` only. The finite-element exports in `reference/` are
plain text, so everything except the no-load sweep itself can be checked
without an ANSYS licence. The shims in `tools/octave_compat/` are needed
only under GNU Octave.

## 6. Licence

Copyright (c) 2026 Idris Laouar, Ahcene Boukadoum, Nabil Mezhoud,
Electrotechnical Laboratory Skikda (LES), University 20 August 1955 Skikda.

The archive is released under **CC BY 4.0**, with one exception. The
finite-element exports in `reference/ANSYS_18_5kW/` are carved out: the
right to redistribute them has not yet been confirmed against the ANSYS
licence agreement, so no licence is granted over that directory until it
is. See `LICENSE.txt`.

## 7. How to cite

    I. Laouar, A. Boukadoum and N. Mezhoud, "Reproduction archive for a
    doubly slotted annulus operator for magnetic equivalent circuits
    (Paper II)", release v1.1.0, GitHub, 2026.
    https://github.com/il1996/doubly-slotted-annulus-IM/releases/tag/v1.1.0

Cite the **tagged release**, not the branch. A branch moves; a tag does
not, and the paper refers to the state of the archive at `v1.1.0`.

`CITATION.cff` holds the same metadata in machine-readable form.
