# MANIFEST — une grandeur, une chaîne

> Bloc **C1** de `SPEC_CLAUDE_CODE_v3.md` §8. Ce fichier déclare, pour
> chaque grandeur publiée, **la chaîne qui fait foi**, sa sortie datée, et
> la configuration qui la définit. Toute cellule d'un tableau des deux
> articles doit se retrouver dans une ligne de ce manifeste.
>
> **Établi le 6 août 2026.**

---

## 0. Comment lire une sortie `.txt`

**Règle** : c'est le **dernier bloc complet** qui fait foi, pas le second.

La formulation d'origine — « c'est le second bloc » — est fausse sur au
moins un fichier. `A1_table7_out.txt` contient **quatre** blocs, et les
**trois premiers** donnent une FEM de maillage nulle : c'est le bug du
démarrage à chaud `U0`, qui n'est licite qu'avec le solveur Newton. Un
lecteur appliquant l'ancienne règle prendrait un bloc cassé.

Vérifier la **cohérence interne** du bloc retenu, pas seulement son rang.

---

## 1. PMSM 15/14, 750 W — Article I

| grandeur publiée | chaîne qui fait foi | sortie | date |
|---|---|---|---|
| Table 7, toutes lignes | `RUN_A1_TABLE7.m` → `mesh_bldc` + `solve_bldc_mesh` (maillage) ; `cogging_mec` + `inductance_mec` (localisé) | `A1_table7_out.txt` **dernier bloc** | 4 août |
| Table 5(a), $n_{sh}=1$ | `RUN_A2_TABLE5.m` | `A2_table5_out.txt` 2ᵉ bloc | 4 août |
| Table 5(b) **et** encadrement | `RUN_X1_TABLE5B_RECONCILE.m`, $n_{sh}=4$ | `X1_table5b_reconcile_out.txt` | 5 août |
| bases P0/P1 sur le bore | `RUN_V1_PMSM_BASE.m` | `V1_pmsm_basis_out.txt` 2ᵉ bloc | 4 août |
| pertes aimants | `RUN_A3_PMLOSS.m` → `pm_loss`, $N_p=241$ | `A3_pmloss_out.txt` | 4 août |
| pertes fer sur maillage | `RUN_X2_IRONLOSS.m`, intégration **par cellule** | `X2_ironloss_out.txt` 2ᵉ bloc | 5 août |
| pertes fer, décomposition et dérive en $n_{sh}$ | `RUN_A4_IRONLOSS.m` | `A4_ironloss_out.txt` | 5 août |
| trois cellules *Lumped* | `RUN_A5_LUMPED.m` → `cogging_mec` + `krkl` | `A5_lumped_out.txt` | 5 août |
| mutuelle $M$, compensation $L_d$ | `RUN_A5BIS_MUTUAL.m` | `A5bis_mutual_out.txt` | 5 août |
| comparaison de fermetures N1/N2/N2b/N3 | `RUN_CARTER_CMP.m` | `A6_n2b_out.txt` | 5 août |
| temps de calcul | `RUN_T9_TEMPS.m` | `T9_temps_out.txt` | 4 août |
| **figures** | `BLDC_MEC_COMPLET.m` | `BLDC_FIG*.pdf` (vectoriel) | 6 août |

## 2. MAS 48/44, 18,5 kW — Article II

| grandeur publiée | chaîne qui fait foi | sortie | date |
|---|---|---|---|
| schéma équivalent, point nominal, trois essais, champs, caractéristique | `RUN_B1_IM_P1.m` | `B1_im_p1_out.txt` 2ᵉ bloc | 4 août |
| $k_C(N_h)$, deux bases | `RUN_B2_KC.m` | `B2_kC_out.txt` | 4 août |
| le prix (Carter / P0 / P1) | `RUN_B3_PRIX.m` | `B3_prix_out.txt` | 4 août |
| $R'_r(s)$, effet de peau | `RUN_B4_BARSKIN.m` | `B4_barskin_out.txt` | 5 août |
| vrillage | `RUN_B5_SKEW.m` | `B5_skew_out.txt` | 5 août |
| anneau **et** audit de `ansys_ref` | `RUN_B6_ANNEAU.m` | `B6_anneau_out.txt` | 5 août |
| sondes locales Table 17 | `RUN_B7_PROBES.m` | `B7_probes_out.txt` | 5 août |
| conformité de la trace (Art. I §3) | `RUN_X3_TRACE_TABLE.m` | `X3_trace_table_out.txt` | 5 août |

---

## 3. Chaînes de sondes ≠ chaîne de performance

**À déclarer dans toute légende de la Table 17.** Les deux chaînes de la
machine asynchrone n'emploient pas le même opérateur d'entrefer :

| | réseau de performance | chaîne de sondes |
|---|---|---|
| opérateur | `mec.airgap_dtn_tooth` | `mec.airgap_fourier` |
| base | **P1** | **P0**, imposée (`airgap_fourier.m:63` ; `mesh_refined.m:271` appelle sans argument) |
| troncature | $N_h = 8192$ | $N_h = 100$ (`mesh_refined.m:258`) |

Rapport de troncature **82**. Les sondes dérivent de jusqu'à **5,55 %**
sur un facteur 8 de troncature, le maximum sur la dent rotor.

---

## 4. Variantes qui coexistent — et pourquoi elles restent

La spécification supposait des variantes orphelines à déplacer dans
`archive/`. **Vérification faite, six des sept sont vivantes** : elles
alimentent des scripts de diagnostic distincts. Les archiver casserait
la chaîne.

| fichier | référencé par | verdict |
|---|---|---|
| `pm_loss.m` | `RUN_A3_PMLOSS`, `RUN_T18_CONV` | **chaîne des pertes aimants à vide** |
| `pm_loss_load.m` | `BLDC_MEC_COMPLET` + 6 scripts | **chaîne des pertes aimants en charge** |
| `pm_loss_R.m` | `RUN_SLOT2D_VALID:52` | diagnostic — conserver |
| `cogging_mec.m` | `BLDC_MEC_COMPLET`, `RUN_A1`, `RUN_A5`, `RUN_CARTER_CMP` | **chaîne du champ denté** |
| `cogging_mec2.m` | `RUN_A8:70` | diagnostic (bouche 2D maillée) — conserver |
| `inductance_mec.m` | `RUN_A1_TABLE7:96` | **chaîne des inductances localisées** |
| `inductance_mec2.m` | **personne** | ⚠ seule variante réellement morte |
| `subdomain_mec.m` | 7 scripts de diagnostic | conserver |
| `subdomain_mec2.m` | `RUN_SHOE` | conserver |

**Décision.** Aucune variante n'est déplacée. Ce qui manquait n'était pas
un ménage mais une **déclaration** : c'est l'objet des §1 et §2 ci-dessus.
`inductance_mec2.m` est laissée en place avec cette mention — elle est
sans appelant, mais la supprimer n'apporterait rien et le dossier vient
de perdre deux fichiers sans cause identifiée.

---

## 5. Utilitaires promus — fin des duplications forcées

Quatre fonctions étaient **locales à un script**, donc non appelables :
toute reprise devait les dupliquer, ce qui créait mécaniquement une
seconde chaîne pour la même grandeur. Elles sont promues :

| fonction promue | était locale à | bloc qui avait dû la dupliquer |
|---|---|---|
| `mec.inst_currents` | `RUN_ARTICLE.m:508` | B7 |
| `mec.surf_potentials` (ex-`surfU`) | `RUN_ARTICLE.m:514` | B7 |
| `mec.regional_max` | `RUN_ARTICLE.m:531` | B7 |
| `krkl` (MEC_BLDC) | `RUN_A1_TABLE7.m:139` | A5 |

Les scripts hôtes conservent leur copie locale, qui les masque
localement : **leur comportement est inchangé**. Tout script nouveau doit
appeler la version promue.

---

## 6. Figures — export et exceptions

`BLDC_MEC_COMPLET.m`, fonction locale `savefigure` :

- **export vectoriel `.pdf` par défaut** (`ContentType`, `vector`) — c'est
  ce que le `.tex` appelle ;
- `.png` 130 dpi conservé pour la relecture rapide ;
- `.fig` conservé pour l'édition ;
- **une seule exception raster** : `BLDC_FIG8_cartes2D`, carte de champ
  dense (des dizaines de milliers de patches) → PNG **600 dpi**.

`BLDC_FIG5_charge` a été **scindée en trois** — neuf sous-graphiques
tombent sous 45 mm de large en simple colonne :

| figure | contenu |
|---|---|
| `BLDC_FIG5_1_transient` | vitesse, courant de phase, régime établi, courant de bus |
| `BLDC_FIG5_2_conversion` | couple, pertes Joule, rendement, répartition |
| `BLDC_FIG5_3_saturation` | $L_{\ell\ell}(i)$ et $\psi_{ab}(i)$ |

`BLDC_FIG5b_charge_local` est inchangée.

---

## 7. `RUN_IND_MESH.m` — bloc C2, tranché

La spécification demandait : « en-tête *NE PAS UTILISER* sur un script qui
**alimente la Table 7** — à corriger ; lignes 78-79 : les deux termes se
terminent par `*0`, le diagnostic vaut zéro par construction, alors que la
conclusion imprimée juste en dessous porte l'argument central de la §3.5.
Rétablir le calcul ou retirer l'affirmation. »

**La prémisse est fausse sur un point et juste sur l'autre.**

1. **`RUN_IND_MESH.m` n'alimente pas la Table 7.** `RUN_A1_TABLE7.m`
   calcule $L_a$ et $M$ directement depuis `mesh_bldc` (l. 75-79) et la
   colonne localisée depuis `inductance_mec` (l. 96). `RUN_IND_MESH`
   n'est appelé nulle part. Son propre en-tête le dit d'ailleurs
   (l. 10) : « les inductances du scorecard restent celles
   d'`inductance_mec` ». **L'en-tête d'avertissement est correct et doit
   rester.**

2. **Le `*0` est réel.** Lignes 78-79, les deux termes de la part
   traversant la surface d'entrefer se terminent par `*0` : la valeur
   imprimée est identiquement nulle. Et la phrase imprimée juste en
   dessous — « la fuite d'encoche n'est plus une formule ajoutée » — n'est
   donc soutenue par aucun chiffre.

**Décision — corrigée le 7 août 2026.** La rédaction précédente de ce
paragraphe annonçait « $L_a = 12{,}25$ mH contre 50,21 mH mesurés, soit
**−75,6 %** ». **Ce chiffre n'est pas dans `ind_mesh_out.txt`** et a été
retiré : la vérification faite directement sur la sortie donne
$L_{a,\text{EF}} = 50{,}209$ mH, et le maillage la reproduit à
**+0,3 à +0,6 %** ($M_s = 360$ et $540$, $n_{sh} = 1$ et $2$). Le seul
écart de cet ordre dans le fichier est celui du **réseau à une dent**,
$35{,}351$ mH, soit **−29,6 %** — c'est-à-dire d'une *autre* chaîne.
La leçon est celle de la règle 7 : ce paragraphe avait été écrit depuis
un souvenir de run et non depuis la sortie.

Ce qui reste vrai, et qui suffit : **le diagnostic du script vaut zéro
par construction**. Les deux parts imprimées — « part portée par les
branches d'AIR » et « part traversant la SURFACE d'entrefer » —
s'affichent l'une et l'autre à `0.0 %`, et la phrase imprimée juste en
dessous n'est donc étayée par aucun nombre de ce fichier. **Aucune
affirmation de la §3.5 ne peut reposer sur lui**, non parce que ses
inductances seraient fausses, mais parce que le chiffre qu'il devait
produire n'existe pas.

L'affirmation reste néanmoins vraie et publiable — mais par une autre
chaîne : `mesh_bldc` résout la fuite d'encoche dans ses cellules d'air et
donne $L_a = 50{,}379$ mH contre $50{,}209$ mH, soit **+0,34 %**
($n_{sh} = 1$, `A5bis_mutual_out.txt`). C'est cette chaîne, et elle
seule, qui étaye la §3.5. Deux chaînes calculent donc $L_a$ ; une seule
fait foi, et c'est A5bis.

---

## 8. Vocabulaire — décision

Les figures écrivent « MEC », le manuscrit écrit « mesh model » et
« lumped model ». **Dans chaque figure prise isolément, « MEC » n'est pas
ambigu** : une seule chaîne y apparaît. L'ambiguïté est dans les
**tableaux du manuscrit**, où les deux variantes sont côte à côte.

Décision : les figures gardent « MEC » ; les tableaux écrivent
**« MEC (mesh) »** et **« MEC (lumped) »**. Aucun label de figure n'est
touché — les modifier romprait des légendes déjà serrées sans lever
d'ambiguïté réelle.

---

## 9. Article II, version « Trace Conformity, the Slot-Opening Condition » — grandeurs nouvelles (15–16 septembre 2026)

Même règle : une grandeur, une chaîne. Les chaînes MATLAB tournent depuis
`code/MEC_IM/` et lisent la référence par le chemin relatif
`../../reference/ANSYS_18_5kW` ; les chaînes Python tournent depuis
`code/python/` et écrivent dans `outputs/python/`.

| grandeur publiée | chaîne qui fait foi | sortie | date |
|---|---|---|---|
| Table 1 (données machine), Table 6 (références relues des `.tab`, fenêtre t ∈ [1, 2] s) | `fea_audit.py`, `fea_power.py`, `sweep_audit.py`, `fea_conv*.py` ; relecture MATLAB dans `RUN_Z1_CAVITY.m` | `Z1_cavity_out.txt` (en-tête), `outputs/python/` | 15 sept. |
| Table 2 (invariants, chapeau asymétrique) | `RUN_INVARIANTS.m`, `RUN_R8_TABLE2.m` ; `t_op1.py`, `t_op2.py` (I-1, I-4 en base p1a) | `INV_tests_out.txt`, `R8_table2_out.txt` | 12 août / 15 sept. |
| Table 3 (troncature, deux bases) | `RUN_B2_KC.m` | `B2_kC_out.txt` | 4 août |
| Table 4 (référence EF du rapport d'encochage, formulations A et B, maillages, 12 positions) | `prod_fem.py`, `t_fem2.py`, `fem_slots.py`, `fem_annulus.py` | `prod_fem_results.json`, `fem_kc_results.json` | 15 sept. |
| Table 5, colonnes uniformes (Φ_O = 0, deux bases ; cavités, base chapeau) | `t_op3.py` ; `prod_op.py` | `tiling_sweep_inf_iron.json` ; `prod_op_results.json` (`cav_sweep`, `pos_cav`, `cav_33_16_p1a`) | 15 sept. |
| Table 5, dernière colonne (pavage gradué, q = 1,5, chapeau asymétrique) et §4.2 (étude de q, contrôle N_h = 16 384) | `t4_graded.py` (+ `cavity_graded.py`) | `t4_graded_results.json`, `t4_graded_log.txt` | 16 sept. |
| Fig. 2, Fig. 4 | `prod_op.py`, `prod_fem.py` → `make_figures_v2.py` | `prod_op_results.json`, `field_waveforms.npz` | 15–16 sept. |
| Fig. 3 | `t_op3.py`, `prod_op.py`, `t4_graded.py` → `make_figures_v2.py` | idem | 16 sept. |
| Table 7 (constantes identifiées) | `+mec/leakage.m`, `machine_18_5kW.m`, `stator_resistance.m`, `mech_losses.m`, `stray_losses.m` (valeurs déclarées dans le code) ; sensibilité §5.3 : `RUN_Z3_LEAKAGE.m` | `Z3_leakage_out.txt` | 15 sept. |
| Table 8 (cinq fermetures, à vide et nominal) | `RUN_Z1_CAVITY.m` | `Z1_cavity_out.txt` | 15 sept. (reproduit sous MATLAB R2024a le 16 sept., tous chiffres identiques) |
| Table 9 (schéma équivalent Φ_O = 0, calage) | `RUN_B10_B1_SKEWOFF.m` | `B10_b1_skewoff_out.txt` | 6 août |
| Table 10 (convergence en pavage, Φ_O = 0 et cavités, huit pavages) | `RUN_Z6_COUPLED_CONV_CAV.m` (garde : reproduit `M5_coupled_conv_out.txt` et la ligne (33,16) de Z1) | `Z6_coupled_conv_cav_out.txt` | 16 sept. |
| Table 11, colonnes Φ_O = 0 | `RUN_M11_FIELD_ERR.m`, `RUN_M11B_ROTORPOS.m`, `RUN_M11C_MATCHED.m` | `M11_field_err_out.txt`, `M11B_rotorpos_out.txt`, `M11C_matched_out.txt` | 12 août |
| Table 11, colonnes cavités (33,16) | `RUN_Z4B_FIELD_ERR_CAV.m` | `Z4B_field_err_cav_out.txt` | 16 sept. |
| Fig. 5 (caractéristiques, deux fermetures ; décrochage) | `RUN_B10_B1_SKEWOFF.m` (`B10_b1_skewoff.mat`), `RUN_Z5_SWEEP_CAV.m` (`Z5_sweep_cav.mat`) → `make_figures_v2.py` | `B10_b1_skewoff_out.txt`, `Z5_sweep_cav_out.txt` | 16 sept. |
| §5.5, calage avec cavités (74,66 N m, 107,5 A) | `RUN_Z5_SWEEP_CAV.m` | `Z5_sweep_cav_out.txt` | 16 sept. |
| Fig. 6 (champs au mi-entrefer, deux fermetures) | `RUN_Z2_FIELDS.m`, `RUN_Z4_FIELDS_CAV.m` → `make_figures_v2.py` | `Z2_fields_{avide,charge}.txt`, `Z4_fields_cav_{avide,charge}.txt` | 15–16 sept. |
| §6.3, carte de cavité en charge (variantes A/B/C) | `RUN_Z8_CAVITY_LOAD.m` avec `cavity_src_nO16.mat` (`export_cavity_src.py` ← `cavity_graded.cavity_source_rotor`) ; `+mec/airgap_dtn_tooth_cav.m` (`set_source`), `+mec/solve_network.m` (`AG.f`) | `Z8_cavity_load_out.txt` | 16 sept. |
| §5.7 et Table 12, essai à vide numérique, côté réseau (mêmes tensions, même construction) | `RUN_Z7_NOLOAD_NET.m` | `Z7_noload_net_out.txt` | 16 sept. |
| §5.7 et Table 12, essai à vide numérique, côté EF (cinq tensions, rotor entraîné à 1500 tr/min) | `code/ansys_noload/noload_sweep_com.py` (copie du projet, designs `NL_V*`), `solve_one_com.py`, `noload_export_com.py`, `run_rest_com.ps1`, puis `noload_postprocess.py` (phaseurs sur [1, 2) s, E_2D = V − R_s I − jωL_ext I, contrôle contre jωΨ) et `rms_check.py` (THD des courants de référence) | `reference/ANSYS_18_5kW/noload_sweep/NL_V*_NL_{wave,misc}.tab`, `outputs/ansys_noload/noload_results.{txt,json}`, `sweep_log.txt` | 16 sept. |
| §4.3, opérateur à cavités au pavage (33,32), deux bases (p1, p1a), douze positions, deux comparateurs (assemblé = convention de Table 5 ; forme fermée = convention de la chaîne EF), référence `pos_real` (formulation A, lc 0,045) — moyennes, étendues position par position, écart entre bases | `prod_pos_33_32.py` (fonctions de production `dsop.slotting_ratio_cavity`, `fem_slots.staircase_Fp`, matrices `cavity_Q.pkl` de `prod_op.py`) | `prod_pos_33_32_results.json`, `prod_pos_33_32_out.txt` | 17 sept. |
| §4.2, branches de raffinement de la condition Φ_O = 0 à φ = 0, N_h = 8192, deux bases : grille de Table 5 (`t_op3.py`) étendue à (33,32), (33,64), (65,32), (65,64), (129,16), (129,32), (129,64) ; incréments, rapports, extrapolations géométriques (conventions A : queue = d·r/(1−r) ; B : queue = dernier incrément), limites de branche et limite double contre la formulation B (lc 0,025) — le « 1,455 » et la limite double sortent de ce script | `prod_limits.py` | `prod_limits_results.json`, `prod_limits_out.txt` | 17 sept. |
| §4, tous les pourcentages dérivés (facteurs de Carter contre EF, incertitude de maillage, 13 % / 11 % / 16–31 % / 19 % / 0,2–0,4 % / 0,3–0,7 % / 1,4–2,2 % / 0,5 % / 0,4 % / 6 % / 0,3 % / −2,8 % / −20,4 %, « between 1.44 and 1.47 »), chacun avec pavage, troncature, base, position et grille de référence | `prod_section4_pct.py` (lit les `.json` ci-dessus et `dsop.carter`) | `prod_section4_pct_out.txt` | 17 sept. |
| §4.1 et Annexe B, les trois contrôles avec leur taille d'élément : (1) nappe de courant de largeur nulle sur alésage lisse contre la forme fermée, sept lc (0,12–0,028) ; (2) encoche profonde contre nappe à l'embouchure, fondamental à sept décimales et forme d'onde complète, lc 0,06–0,028 ; (3) formulation B contre le potentiel en escalier (construction de `t_ann_check.py`), lc 0,08–0,025 | `prod_appB_checks.py` | `prod_appB_checks_results.json`, `prod_appB_checks_out.txt` | 17 sept. |
| Scripts producteurs de Table 2 (invariants, `INV_tests_out.txt`) et de Table 11, colonnes Φ_O = 0 (`M11_*_out.txt`), absents de l'archive jusqu'au 17 sept. : copies **verbatim** des originaux du 12 août 2026 (`MEC\MEC_IM\`, sha256 dans `RAPPORT_MESURES_R25_R26.md`), ré-exécutés le 17 sept. contre le paquet `+mec` de cette archive : `INV_tests_out.txt`, `M11_field_err_out.txt`, `M11C_matched_out.txt` reproduits ligne à ligne (durées exceptées) ; `M11B_rotorpos_out.txt` reproduit à tous ses nombres, le verdict de la garde 2 différant parce que le script porte un seuil corrigé après coup et déclaré (l. 124–134 : 1e-6 → 1e-4 T). Ces trois scripts `RUN_M11*` lisent la référence par le chemin ABSOLU `C:\Users\hp\Desktop\ANSYS résultat 18.5KW` (fichiers identiques à `reference/ANSYS_18_5kW`, sha256 vérifiés), non par le chemin relatif | `RUN_INVARIANTS.m`, `RUN_M11_FIELD_ERR.m`, `RUN_M11B_ROTORPOS.m`, `RUN_M11C_MATCHED.m` | sorties déjà archivées (12 août) | 17 sept. |
| Sensibilité de base au niveau machine (R30, 17 sept.) : fermeture cavités, (33,16) et (33,32), N_h = 8192, configuration de RUN_Z1_CAVITY, chapeau symétrique p1 contre chapeau asymétrique p1a — X_m0, I_0, E_1, X_m(s), I_1, couple, P_fe, P_cu,r, cos φ, η, écarts à la référence Maxwell. **Le paquet de production n'accepte pas 'p1a' et n'a pas été modifié** : la base p1a vient de `variant_p1a/+mec/airgap_fourier.m` (copie de la production + cas 'p1a', identique à la production pour p0/p1 et à `dsop.assemble` pour p1a, MES_R3_VALIDATE), rendue prioritaire comme dossier courant par le script ; garde 1 : points p1 identiques à la production seule (`RUN_Z9_BASIS_P1A_PROD.m`, session séparée) ; garde 2 : (33,16) p1 identique à Z1. `cavity_nO32.mat` par `export_cavity_nO32.py` (construction de `export_cavity.py`, lc_mouth 0,005) | `RUN_Z9_BASIS_P1A_PROD.m` puis `RUN_Z9_BASIS_P1A.m` ; `variant_p1a/+mec/airgap_fourier.m` ; `export_cavity_nO32.py` | `Z9_basis_p1a_out.txt`, `Z9_basis_p1a.mat`, `Z9_basis_p1a_prod.mat` | 17 sept. |
| Robustesse du gain des cavités au changement de base (R31, 18 sept.) : les trois fermetures de RUN_Z1_CAVITY (Carter ; Φ_O = 0 ; cavités) aux pavages (17,4), (33,16), (33,32), N_h = 8192, bases p1 et p1a, une seule exécution — X_m0, I_0, E_1, X_m(s), I_1, couple, P_fe, P_cu,r, cos φ, η, écarts à Maxwell, écart entre bases en points, gain Φ_O = 0 → cavités par base (la paire du résumé « +21.5 % to +7.7 % » = Φ_O = 0 (17,4) → cavités (33,16)). Même variante `variant_p1a` et mêmes gardes que Z9 (points p1 et Carter contre la production seule, session séparée ; contre Z1) | `RUN_Z10_CLOSURES_BASIS_PROD.m` puis `RUN_Z10_CLOSURES_BASIS.m` | `Z10_closures_basis_out.txt`, `Z10_closures_basis.mat`, `Z10_closures_basis_prod.mat` | 18 sept. |

Deux précisions de provenance :

- `cavity_nO{2,4,8,16}.mat` (matrices **Q** des cavités) sortent de
  `export_cavity.py` (élément `cavity.py`, gmsh) ; `cavity_Q.pkl` porte les
  mêmes matrices pour n_O = 4, 8, 16, 32 (`prod_op.py`).
- Sous MATLAB, `+mec/equivalent_circuit.m` porte la seule modification
  faite au paquet depuis le 6 août (ligne 82, scission en deux instructions,
  sans effet sur aucun résultat) ; `+mec/airgap_dtn_tooth_cav.m` et
  `+mec/solve_network.m` ont reçu le 16 septembre un chemin optionnel de
  source de flux (`cav.sr`, `AG.f`) inactif tant qu'aucun script n'appelle
  `set_source` — la garde de `RUN_Z8_CAVITY_LOAD.m` vérifie que la chaîne
  sans source reproduit la Table 8 au 1e-4 près.
