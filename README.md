# Reproduction archive — Paper II

*A Doubly Slotted Annulus Operator for Magnetic Equivalent Circuits:
Carter's Factor Is Low in the Thin-Gap, Wide-Opening Regime, on an
18.5 kW Cage Induction Machine.*

This archive supports **one paper only**. The companion paper on trace
conformity has its own archive, and neither depends on the other to be
read or re-run.

## 1. What is here

| directory | contents |
|---|---|
| `code/MEC_IM/` | the MATLAB chain for the 18.5 kW cage induction machine |
| `code/article/` | the manuscript source, its bibliography and the 5 figures it calls |
| `outputs/MEC_IM/` | execution transcripts, produced by `diary`, at full precision |
| `reference/ANSYS_18_5kW/` | the finite-element exports used as the reference |
| `notes/` | delivery notes for each verification block, and the independent audit |

## 2. The two-minute test

Do this before anything else. The block below reproduces a published
quantity and checks itself. It takes about 95 seconds.

```matlab
cd code/MEC_IM
RUN_R8_TABLE2
```

Expected: the converged row of Table 2, no-load current **+21.6 %** and
saturated magnetising reactance **−10.6 %**, and `GARDE PASSEE`.

If it passes, the archive works on your machine.

## 3. How to read a published number

Every value in the paper comes from a transcript in `outputs/`. The rule
the project follows is **one quantity, one chain**: when two scripts
produced the same quantity, one of them was archived rather than kept.
`notes/MANIFEST.md` binds each published quantity to the chain that
produced it.

Each transcript opens with the configuration that produced it: tiling
$n_T$ and $n_O$, truncation $N_h$, surface basis, slip, skew setting,
and the path of the finite-element reference used. Two of those settings
change published numbers by more than a point, so none of them is
optional reading.

## 4. Verification blocks

- `R1_NOTE.md` reconciles the two implementations that returned
  different values for the single-slice ripple.
- `R2_NOTE.md` sweeps the tiling over a factor of four in each
  direction. It is the block that establishes the paper's central
  reservation: the slotting ratio has not converged in tiling, and its
  dispersion is of the same order as the correction it reports.
- `R3_NOTE.md` states the finite-element verification of the ratio,
  which requires two magnetostatic solves and is not delivered.
- `R8_NOTE.md` regenerates Table 2 in the configuration of the rest of
  the paper and removes the note that had been absorbing the difference.

Each block carries a guard: a check that would contradict its own result
if the result were wrong. `R2_NOTE.md` is the case where the guard
failed and the failure is reported: the ratio is not better conditioned
than its denominator.

## 5. What the archive establishes, and what it does not

An independent audit, recorded in `notes/AUDIT_INDEPENDANT.md`, searched
every value published at three or more significant digits against the
numbers held in these transcripts, matching numerically rather than by
string. The match rate is 93 %. The remainder are declared machine
inputs, such as the friction and windage loss, or arithmetic carried out
in the text, such as the ratio of two reactances.

That establishes that no published result stands detached from an
execution. It does not establish reproducibility: a number present in a
transcript proves that some run produced it once, not that this code
reproduces it today on another machine. Only the test of Section 2
proves that.

Two limitations are declared in the paper itself and are visible in the
transcripts. The tiling has not converged. And the magnetising reactance
of the reference is recorded as not re-measurable from the available
exports, so the deviations that depend on it inherit that status.

## 6. Environment

MATLAB R2021b or later, no toolbox beyond the base product. The
finite-element exports in `reference/` are plain text, so the archive
can be checked without an ANSYS licence. The verification of the
smooth-to-slotted permeance ratio described in `R3_NOTE.md` does require
the ANSYS project and is not part of this archive.

## 7. Licence

Copyright (c) 2026 Idris Laouar, Ahcene Boukadoum, Nabil Mezhoud,
Electrotechnical Laboratory Skikda (LES), University 20 August 1955 Skikda.

The archive is released under **CC BY 4.0**, with one exception. The
finite-element exports in `reference/ANSYS_18_5kW/` are carved out: the
right to redistribute them has not yet been confirmed against the ANSYS
licence agreement, so no licence is granted over that directory until it
is. See `LICENSE.txt`.

## 8. How to cite

    I. Laouar, A. Boukadoum and N. Mezhoud, "Reproduction archive for a
    doubly slotted annulus operator for magnetic equivalent circuits
    (Paper II)", release v1.0.0, GitHub, 2026.
    https://github.com/il1996/doubly-slotted-annulus-IM/releases/tag/v1.0.0

Cite the **tagged release**, not the branch. A branch moves; a tag does
not, and the paper refers to the state of the archive at `v1.0.0`. The
release page carries a downloadable snapshot of exactly that state.

`CITATION.cff` holds the same metadata in machine-readable form, so
GitHub's "Cite this repository" button and reference managers pick it up
without retyping.

A `.zenodo.json` is present but inert: this archive is distributed
through GitHub only. Enabling the Zenodo integration on the repository
would mint a DOI from the next release without any further editing.
