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
