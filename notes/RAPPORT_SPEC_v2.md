# Exécution de `SPEC_CLAUDE_CODE_v2.md` — rapport

Suite de `RAPPORT_SPEC_CLAUDE_CODE.md`.

| Bloc | État |
|---|---|
| **V1** verrou | **TRANCHÉ — issue n° 2** |
| **A1** Table 7 | **FAIT — 17 lignes, une exécution** |
| **A2** Tables 5(a)/5(b) | **FAIT — 2 cellules comblées, colonne peak B ajoutée** |
| **A3** pertes aimants | **FAIT — la dichotomie du manuscrit s'effondre** |
| A4–A7 (Article I) | débloqués, non exécutés |
| **B2** k_C(N_h) publiable | **FAIT — §6.4 réfutée par lecture directe** |
| **B3** le prix chiffré | **FAIT — P1 perd sur les 5 colonnes** |
| **B1** chaîne MAS en P1 | **FAIT — sondes Table 17 exceptées** |
| **A4** pertes fer | **FAIT — mon −14,8 % était faux : +25,3 %** |
| **A5** cellules *Lumped* | **FAIT — RÉSOLUES, φ = 39,429°** |
| **A6** ligne N2b | **FAIT — γ_m = 0 à a_r = a_t** |
| **B4** `bar_skin` | **FAIT — prémisse de T12 invalide, signe inversé** |
| **B6** anneau | **FAIT — 105,5 % de la bascule vient de la barre** |
| **A7** figures | **traité dans `MANIFEST.md` §6** |
| **B5** vrillage | **FAIT — la §6.7 est fausse CÔTÉ MEC** |
| **C1** chaînes figées | **déjà fait — `MANIFEST.md`, 8 sections** |
| **C2** `RUN_IND_MESH` | **FAIT — 2W = λ à 1,3×10⁻¹²** |
| **C3** multiples de σ | **FAIT — deux essais, non deux normalisations** |
| **C4** colonnes d'écart | **FAIT — Tables 8, 10, 16, 17 régénérées** |
| **C5** ondulation / vrillage | **FAIT — vrillage sans objet ; 101,4 périmé** ⚠ *le « P1 dégrade » est réfuté : voir A-1 puis R1* |

### Série R — réponse au rapport d'évaluation (`SPEC_CLAUDE_CODE_v4`)

| bloc | état | résultat |
|---|---|---|
| **R1** écart 108,2 / 128 N·m | ✅ **TRANCHÉ** | **108,224433 N·m publiable** ; l'écart de +18,27 % vient de la **rotation des courants statoriques omise** (`ph = p·θ`). ⚠ **La garde échoue** : couple moyen de la carte à **−15,36 %** du schéma équivalent — défaut de modèle distinct, signalé |
| **R2** sensibilité $k_C(n_T,n_O)$ | ✅ **TRANCHÉ — défavorable** | $k_C$ **disperse de 6,82 %** sur neuf pavages (1,2978 → 1,3889). **La dispersion excède la correction annoncée** : +2,5 % à +9,7 % au-dessus de Carter contre +5,2 % publié. $X_m$ encoché **croît encore** dans les deux directions : pavage **non convergé**. Garde : le rapport n'est **pas** plus stable que ses composantes (6,8155 % contre 6,7523 %) — mais le numérateur étant constant par construction, **aucune compensation n'est possible sur cet axe**, ce qui répond à l'objection anticipée |
| **R3** vérification EF du rapport | ⛔ **BLOQUÉ** | Exige **deux résolutions ANSYS** ; le projet ne contient **aucune géométrie lisse** (un seul `Setup1`). **Livré** : garde analytique **4,3779802798e-04 H**, convention de rayon tranchée (0,1486 %, sous tolérance), spécification des deux runs, traitement prêt. **Aucun chiffre EF de substitution** — la règle 5 l'interdit |
| **R4** déverrouillage troncature–pavage | ✅ **TRANCHÉ — favorable** | Verrou confirmé : crochet = **ln π à six décimales**, +0,0005 % sur un facteur 16 en $M_s$. Déverrouillé à pavage fixe : crochet **+182 %**, pente **2,302585/décade** contre $\ln 10 = 2{,}302585$ prédit (**−0,0000 %**), $B_{g1}$ dérive de **+0,29 %** sans se stabiliser. **La prédiction falsifiable du §3.5 devient une démonstration.** Garde : identité **au bit près** (0,000e+00 T) au point de recouvrement. Point $N_h = 8M_s$ écarté (perte de rang, mécanisme nommé) |
| **R5** trancher $n_{sh}$ | ✅ **TRANCHÉ — $n_{sh} = 2$** | Deux gardes concordantes. Compte d'inconnues : `mesh_bldc` produit **5400 à $n_{sh}=1$ et 6480 à $n_{sh}=2$**, exactement les deux comptes déclarés par la Table 12 ($n_{sh}=4$ donnerait 8640). $B_{g1}$ : colonne 1 ← $n_{sh}=1$ (1,75e−06), colonne 2 ← $n_{sh}=2$ (4,96e−06), $n_{sh}=4$ écarté (3,94e−04, **80× la tolérance**). **Les trois tableaux sont justes** ; défaut unique **ligne 2358** du `.tex`, correction fournie, rien à régénérer. Aucun `.tex` modifié |
| **R6** carte de saturation | ✅ **TRANCHÉ — défavorable, autre cause** | Quatre gardes passées (deux routes EF à **3,78 %**, calage à 9,7e−17, Fig. 12 reproduite à 0,02 %). **La spécification se trompe de signe** : le réseau **sous**-sature, +39 / +21 / +139 % à 2 / 5 / 10 A. **Deux hypothèses réfutées** : la grille (2 % de l'écart, mauvais sens, anti-corrélée) et la calibration des aimants (réseau −0,5 %, EF −31 %). **Ce qui reste** : l'effondrement publié 104,3 → 12,1 mH n'est produit par aucun réseau (brut 133,7 → 38,9, facteur 3,4 contre 11,0) — il vient d'une **soustraction constante de 29,36 mH dont la prémisse est réfutée** (l'écart change de signe, +29,9 → −19,3 mH). 15 et 20 A **non mesurables**, rien n'est extrapolé |
| **R7** recompter les 29 grandeurs | ✅ **TRANCHÉ — l'annonce est exacte** | **23/29 vérifié**, identique au compteur du programme (garde passée). Le comptage visuel du rapporteur s'explique par trois grandeurs à moins d'un point du seuil, **toutes à l'intérieur** — dont le couple à l'arrêt à **+4,9520 %**. Hors grandeurs que le manuscrit déclare non validées : **21/26**. **Défaut découvert** : l'entrée #20 régénère les +2,45 % retirés (0,341619 W), donc la phrase du §6.6 « *not produced by any chain that can be reconstructed* » est **fausse** — la chaîne est `BLDC_MEC_COMPLET.m:338-340` (`pm_loss_load`), différente de celle que le §6.5 gèle (facteur 1,7) |
| **R8** régénérer la Table 2 (Art. II) | ✅ **FAIT — garde passée** | Table 2 = `tab:price`, régénérée à **vrillage harmonique neutralisé**. Les chiffres coïncident avec ceux cités ailleurs **sans note** : $I_0$ +21,6 % (0,014 pt), $X_m$ −10,6 % (0,010 pt). La note des lignes 808–818 **n'a plus d'objet**. **Conséquence à propager** : le contraste de tête passe de « +1,1 % → +21,6 % » à « **+1,5 % → +21,6 %** », le +1,1 % étant une valeur vrillage actif — quatre occurrences listées. Panneau LaTeX fourni, aucun `.tex` modifié |
| R9 | ⬜ | — |

---

## V1 — La base P1 sur le bore du PMSM : VERDICT

Script : `MEC_BLDC/RUN_V1_PMSM_BASE.m`, sortie `V1_pmsm_basis_out.txt`.
Machine 15/14 750 W, kfringe = 0,325, solveur linéaire.
Portage de la base P1 dans `airgap_magnet.m` et `cogging_mec.m`.

### Écart à la spécification, et sa raison

**Il est impossible de faire varier N seul à pavage fixe sur cette chaîne.**
`cogging_mec.m:27` impose `numax = Nsurf/2` : en dessous, l'opérateur de
couronne Y est de **rang déficient** ; au-dessus, les harmoniques replient sur
une base à Nsurf colonnes.

C'est une **différence structurelle** avec la MAS, où `airgap_dtn_tooth`
condense une grille fine par complément de Schur et laisse N libre — ce qui
avait permis à T1/T16 de balayer N jusqu'à 8192 à pavage constant.

Le balayage porte donc sur `Nsurf`, avec N = Nsurf/2 qui suit : il mesure
l'effet **conjoint** pavage + troncature. Déclaré en en-tête de sortie.

### Résultats

| Nsurf | N | B_g1 P0 | B_g1 P1 | rang 8 P0 | rang 8 P1 |
|---|---|---|---|---|---|
| 540 | 270 | 1.07390 | 1.07267 | 0.021761 | 0.023000 |
| 1080 | 540 | 1.07865 | 1.07801 | 0.016944 | 0.017596 |
| 1260 | 630 | 1.07755 | 1.07698 | 0.018062 | 0.018637 |
| 2160 | 1080 | 1.07793 | 1.07759 | 0.017672 | 0.018020 |
| 4320 | 2160 | 1.07759 | 1.07741 | 0.018014 | 0.018195 |
| 8640 | 4320 | **1.07743** | **1.07734** | 0.018180 | 0.018271 |

**Dispersions sur les six discrétisations**

| Grandeur | P0 | P1 |
|---|---|---|
| B_g1 | **0,440 %** | 0,496 % |
| \|B_r\| moyen | 0,461 % | 0,519 % |
| B_r crête | 1,856 % | 1,613 % |
| B_t rms | 5,592 % | 5,921 % |
| rang 8 | **26,1 %** | **28,5 %** |

**Écart P0 → P1 à la plus fine discrétisation** : B_g1 −0,01 %, |B_r| −0,01 %,
B_r crête +0,03 %, B_t rms +0,10 %, rang 8 +0,50 %.

### Verdict : issue n° 2, et plus nettement qu'attendu

**1. P0 ne dérive pas sur le PMSM.** Jusqu'à N = 4320 — bien au-delà des 3088
où la MAS dérivait de 9,3 % — B_g1 reste stable à 0,44 %. Les incréments
**alternent de signe** (rapports −0,23 / −0,35 / −0,89 / +0,48) : convergence
**oscillante**, non dérive monotone. Rien à voir avec les rapports tendant
vers 1 de T1.

**2. P1 ne change rien** : 0,01 % sur B_g1. **L'issue n° 3 est exclue —
la Table 7 n'est PAS à refaire, ni les chiffres du résumé.**

**3. Explication structurelle, à écrire dans l'Article I.**
- La couronne du PMSM inclut l'aimant de 3,5 mm : X = ln(R_s/r_mi) = **0,139**
  contre **0,00297** pour l'entrefer mince de la MAS — un facteur **47**.
- La surface y est discrétisée en 1260 colonnes uniformes, **non condensée**
  sur 92 nœuds dentaires.

Le mécanisme de la divergence — élimination des colonnes d'ouverture sur un
entrefer très mince — n'existe pas ici.

**4. Confirmation qui renforce §8.1.** Le rang 8 disperse de 26–28 % **dans
les deux bases**. Aucune base ne le fait converger : signature de la
singularité de coin de pointe de dent. La bande de denture reste une grandeur
**bornée, non prédite**, indépendamment du choix de projection.

### Conséquences pour la rédaction

- **Article I** conserve ses chiffres et **gagne un argument** : la
  formulation est robuste au choix de base sur le PMSM, et l'on sait dire
  pourquoi.
- **Article II** garde l'intégralité du problème de convergence, qui devient
  sa **spécificité** plutôt qu'un défaut commun aux deux machines.
- Le contraste entre les deux machines — X différant d'un facteur 47,
  condensation présente ou absente — devient un résultat en soi, et justifie
  la séparation en deux articles.

### Limite à déclarer

La question « la troncature masque-t-elle le défaut ? » ne reçoit qu'une
réponse **partielle** : le balayage ne sépare pas pavage et troncature. Le
test pur exigerait de porter la condensation de Schur au PMSM —
développement, non réglage.

---

## A1 — Table 7 régénérée en une exécution

Script : `MEC_BLDC/RUN_A1_TABLE7.m`, sortie `A1_table7_out.txt`.
Machine 15/14, **M_s = 540**, **N_p = 61**, **base P0**, **solveur linéaire
(μ_i = 3000)**, kfringe = 0,325. Références EF **relues du projet**, aucune
transcription. Inconnues : n_sh=1 → 5400, n_sh=2 → 6480, 1 dent → 1260.

### Valeurs

| Grandeur | Mesh n_sh=1 | Mesh n_sh=2 | Lumped | FEA |
|---|---|---|---|---|
| B_g1 (T) | 1.07787 | 1.07908 | 1.07755 | 1.07455 |
| \|B_r\| moyen (T) | 0.76182 | 0.76273 | 0.76147 | 0.76091 |
| B_r crête (T) | 0.98051 | 0.97536 | 0.98600 | 0.99073 |
| B_t rms (T) | 0.13755 | 0.13478 | 0.14317 | 0.13326 |
| bande ν=8 (T) | 0.01704 | 0.01581 | 0.01806 | 0.01966 |
| ligne de flux A (Wb/m) | 0.00586 | 0.00587 | 0.00586 | 0.00585 |
| flux totalisé (Wb) | 0.25894 | 0.25908 | 0.26121 | 0.25833 |
| FMM /aimant (A) | 761.674 | 759.329 | 762.912 | 769.445 |
| FMM /entrefer (A) | 711.386 | 711.880 | 725.064 | 721.504 |
| dispersion FMM (%) | 6.369 | 5.959 | 6.810 | 6.585 |
| FEM de phase (V) | 238.813 | 238.818 | 240.877 | 233.031 |
| enveloppe six-step (V) | 464.556 | 464.645 | 468.719 | 459.918 |
| L_a (mH) | 50.379 | 50.484 | 50.773 | 50.209 |
| M (mH) | −2.016 | −2.021 | −2.220 | −2.133 |
| L_d (mH) | 52.395 | 52.505 | 52.993 | 52.342 |
| k_r | 1.07069 | 1.06665 | 1.05220 | 1.08693 |
| k_l | 0.85326 | 0.85368 | 0.85371 | 0.86723 |

### Écarts (valeurs pleines)

| Grandeur | Mesh n_sh=1 | Mesh n_sh=2 | Lumped |
|---|---|---|---|
| B_g1 | +0,309 % | +0,422 % | +0,279 % |
| \|B_r\| moyen | +0,119 % | +0,240 % | +0,074 % |
| B_r crête | −1,031 % | −1,551 % | −0,477 % |
| **B_t rms** | +3,218 % | **+1,143 %** | +7,432 % |
| **bande ν=8** | **−13,340 %** | **−19,584 %** | −8,119 % |
| ligne de flux A | +0,115 % | +0,221 % | +0,087 % |
| flux totalisé | +0,236 % | +0,292 % | +1,114 % |
| FMM /aimant | −1,010 % | −1,315 % | −0,849 % |
| FMM /entrefer | −1,402 % | −1,334 % | +0,494 % |
| dispersion FMM | −3,270 % | −9,498 % | +3,416 % |
| FEM de phase | +2,481 % | +2,484 % | +3,367 % |
| enveloppe six-step | +1,008 % | +1,028 % | +1,914 % |
| L_a | +0,339 % | +0,548 % | +1,123 % |
| M | **−5,482 %** | **−5,241 %** | **+4,089 %** |
| **L_d** | **+0,102 %** | +0,312 % | +1,244 % |
| k_r | −1,494 % | −1,865 % | −3,195 % |
| k_l | −1,611 % | −1,562 % | −1,559 % |

### Trois constats

**1. Une configuration VRAIMENT unique est impossible sans changer ce qui est
comparé.** La FEM de phase sort à **+2,48 %** ici, contre +0,9 % publié.
L'écart vient du **solveur** : la valeur publiée utilise le Newton (note (a)
du tableau), A1 impose le linéaire partout. Le manuscrit mélange les deux pour
une raison défendable — les inductances de référence viennent d'un essai
magnétostatique à 1 A où le fer n'est pas saturé, et les confronter à une
solution Newton opposerait des grandeurs de nature différente. **Cet arbitrage
doit être déclaré ligne par ligne**, ce que la note (a) ne fait qu'à moitié.

**2. Aucun n_sh ne domine.** n_sh=2 gagne sur B_t, k_l, la FEM ; n_sh=1 sur
B_g1, B_r crête, ν=8, L_d, k_r. Présenter la colonne *Mesh* comme « le »
modèle raffiné sans déclarer lequel n'est pas soutenable.

**3. La mutuelle M encadre la référence sans l'atteindre** : −5,48 % côté
maillage, +4,09 % côté localisé. Inhabituel pour une grandeur intégrale, et
commenté nulle part dans le manuscrit. À examiner.

### Bug corrigé et documenté

Le **démarrage à chaud n'est pas licite avec le solveur linéaire** : passer
`U0` fait accepter la solution précédente inchangée, le flux totalisé devient
constant sur le balayage et la FEM sort à zéro. Il n'est valide qu'avec le
Newton, où `U0` sert d'initialisation. Erreur née de la composition de deux
recettes correctes séparément (démarrage à chaud de `RUN_SCORE_MESH`, solveur
linéaire de `RUN_KRKL_MESH`). Diagnostic d'étendue de `lam_a` laissé en place.

---

## A2 — Tables 5(a) et 5(b) reconstruites

Script : `MEC_BLDC/RUN_A2_TABLE5.m`, sortie `A2_table5_out.txt`.
Références EF **relues du projet** : B_g1 = 1,07455 T, ν=8 = 0,019659 T,
**L_d = 52,3415 mH** (le tableau publié utilisait 52,24 mH — périmé).

### Table 5(a) — entièrement à n_sh = 1, solveur linéaire

| M_s | inconnues | B_g1 | L_d | ν = 8 | peak iron B (T) |
|---|---|---|---|---|---|
| 180 | **1800** | +0,9 % | **+0,5 %** | −42,3 % | 2.2087 |
| 360 | **3600** | +0,4 % | +0,2 % | −19,8 % | 2.5696 |
| 540 | **5400** | +0,3 % | +0,1 % | −13,3 % | 2.6986 |
| 900 | **9000** | **+0,2 %** | +2,5 % | −9,7 % | 2.7989 |

**Les deux cellules manquantes sont comblées** : L_d(180) = +0,5 %,
B_g1(900) = +0,2 %. La colonne des inconnues est corrigée. B_g1 et ν=8
**reproduisent exactement** les valeurs publiées — la reconstruction est donc
validée par les colonnes qui existaient déjà.

### Table 5(b) — n_sh = 2, les deux solveurs

| M_s | inconnues | ν=8 lin. | ν=8 Newton | peak B lin. | peak B NL |
|---|---|---|---|---|---|
| 360 | 4320 | −25,2 % | +11,8 % | 2.7651 | 1.9271 |
| 540 | 6480 | −19,6 % | +16,8 % | 2.9736 | 1.9778 |
| 900 | 10800 | −16,0 % | **+26,2 %** | **3.1522** | 2.0095 |

### Ce que « peak iron B » démontre

**En linéaire (μ = 3000), le champ crête atteint 2,80 T à n_sh = 1 et
3,15 T à n_sh = 2** — physiquement impossible pour du M350-50A, qui sature
vers 2 T. Le modèle laisse le champ croître **sans borne** au coin de bec,
exactement comme un maillage EF raffiné sur une singularité.

La colonne Newton plafonne à 2,01 T parce que la saturation l'y contraint, et
reproduit les 1,928 / 1,979 / 2,013 publiées **à trois décimales**.

### L'argument de la §8.1, dans sa forme la plus nette

Le tableau montre un **encadrement, jamais une convergence** :

- le **linéaire** sous-estime ν=8 (−25,2 → −16,0 %) et laisse **diverger** B ;
- le **Newton** surestime ν=8 (+11,8 → **+26,2 %**) et **borne** B.

Les deux solveurs s'éloignent de la référence **en sens opposés et de plus en
plus** à mesure qu'on raffine. Et les deux directions de raffinement agissent
elles aussi en sens contraires : angulaire (5a) fait remonter ν=8 de −42 % à
−10 %, radial (5b) le fait redescendre. **Une grandeur dont la valeur raffinée
dépend de la direction du raffinement n'est pas sous-résolue : elle n'est pas
définie.**

### Réserve à lever avant substitution

Mes ν=8 de la Table 5(b) diffèrent d'environ **1,5 point** des publiées
(−25,2 vs −26,8 ; +11,8 vs +10,4 ; −16,0 vs −18,8). L'écart est
**systématique**, donc il désigne un paramètre — probablement `kfringe` ou la
référence EF employée. La colonne « peak iron B » coïncide, elle. À trancher.

---

## A3 — Pertes aimants : chaîne gelée, et la dichotomie s'effondre

Script : `MEC_BLDC/RUN_A3_PMLOSS.m`, sortie `A3_pmloss_out.txt`.
**Chaîne gelée : `pm_loss`, N_p = 241** — la seule dont la convergence soit
établie (0,2034 à N_p=181 → 0,2038 à N_p=481, soit 0,2 % sur un facteur 2,7).

### Résultats

| Point | Source du champ | P (W) | FEA (W) | écart |
|---|---|---|---|---|
| à vide | réseau à une dent | 0.2036 | 0.334 | **−39,1 %** |
| à vide | maillage polaire | 0.2516 | 0.334 | **−24,7 %** |
| en charge | brut | 1.593 | 3.183 | −50,0 % |
| en charge | **normalisé** | — | — | **−32,0 %** |

Normalisation T11 : (1257/1340)² = 0,8800 ; (1,419/1,552)² = 0,8360 ;
produit **0,7356**.

**Périmètre déclaré** : la valeur en charge vient de `BLDC_MEC_COMPLET`
section 5b. `pm_loss_load` exige la chaîne d'entraînement complète (instants,
courants triphasés, position, carte de potentiel d'induit par ampère) et ne
s'appelle pas isolément. Seule la normalisation est calculée ici.

### LA CONCLUSION : la dichotomie du manuscrit n'existe pas

Le manuscrit oppose un point à vide **validant** la chaîne (+2,4 %) à un point
en charge **ouvert** (−50 %). **Les deux chiffres sont des artefacts, et ils se
compensaient.**

- Le +2,4 % repose sur 0,342 W, valeur que **six chaînes testées** ne
  reproduisent pas (§T18). La chaîne gelée donne −39,1 % ou −24,7 % selon la
  source du champ.
- Le −50 % compare **deux machines à des points différents** : 1257 vs
  1340 tr/min, 1,419 vs 1,552 A. Normalisé, il vaut −32,0 %.

**Les trois écarts sont du même ordre : −25 % à −39 %.** Il s'agit d'un écart
**unique**, distribué différemment selon le point de fonctionnement, non de
deux situations qualitativement distinctes.

### Réécriture requise

**§5.7** ne peut plus annoncer que la formulation des pertes aimants est
validée à vide. **§8.4** ne peut plus présenter la charge comme la seule
« open discrepancy » : les quatre hypothèses testées et rejetées portaient sur
une cible de −50 % qui n'en fait que −32, et le point à vide présente le même
ordre d'écart.

**Ce qui reste vrai et publiable** : la source du champ améliore le résultat de
**14 points** (−39,1 → −24,7 % à vide), ce qui est l'apport mesurable du
maillage sur cette grandeur.

---

## B2 — Tableau k_C(N_h) publiable, P0 et P1 côte à côte

Script : `MEC_IM/RUN_B2_KC.m`, sortie `B2_kC_out.txt`.
Pavage fixe nT=17, nO=4 (quasi uniforme, condition du chapeau symétrique).

| N_h | X_m lisse P0 | X_m dent P0 | k_C P0 | X_m lisse P1 | X_m dent P1 | k_C P1 |
|---|---|---|---|---|---|---|
| 512 | 81.912 | 61.900 | 1.3233 | 81.270 | 60.639 | 1.3402 |
| 1024 | 81.944 | 63.094 | 1.2988 | 81.270 | 60.954 | 1.3333 |
| 2048 | 81.973 | 64.038 | 1.2801 | 81.270 | 60.998 | 1.3323 |
| 4096 | 82.003 | 64.863 | 1.2642 | 81.270 | 61.011 | 1.3321 |
| 8192 | 82.032 | 65.607 | 1.2503 | 81.270 | **61.015** | **1.3320** |

### Le contraste qui réfute la §6.4

| colonne | dérive sur le balayage |
|---|---|
| X_m **lisse** P0 | **+0,15 %** |
| X_m **encoché** P0 | **+5,99 %** |
| k_C P0 | −5,51 % |
| X_m **lisse** P1 | **+0,00 %** (rigoureusement constant) |
| X_m **encoché** P1 | +0,62 % |
| k_C P1 | −0,62 % |

**Le numérateur ne bouge pas, le dénominateur dérive, le rapport hérite
intégralement de la dérive.** Rapport des deux dérives P0 : **40**. En P1 le
cas lisse est *rigoureusement* constant — une surface sans saut de potentiel
est représentée exactement par la base chapeau dès les premiers harmoniques.

**L'argument de compensation de la §6.4 est réfuté par lecture directe.**

### Incréments et rapports

| N_h | ΔP0 | rap. P0 | ΔP1 | rap. P1 |
|---|---|---|---|---|
| 1024 | 1.1935 | 0.791 | 0.3155 | 0.138 |
| 2048 | 0.9444 | 0.874 | 0.0436 | 0.298 |
| 4096 | 0.8252 | 0.902 | 0.0130 | 0.304 |
| 8192 | 0.7443 | — | 0.0040 | — |

P0 : rapports **→ 1**, queue logarithmique (T17 la démontre en forme fermée).
P1 : rapports **≈ 0,3**, convergence géométrique.

**Dispersion de k_C : 5,68 % (P0) → 0,62 % (P1).**
**Convergé : k_C = 1,3320, X_m = 61,015 Ω, soit +5,2 % au-dessus de Carter.**
Cohérent au millième avec T16, obtenu sur une autre liste de troncatures —
contrôle croisé.

---

## B3 — Le prix, chiffré

Script : `MEC_IM/RUN_B3_PRIX.m`, sortie `B3_prix_out.txt`.
**Trois** colonnes, pour isoler le prix du changement de base seul : P0 est
construit dans la configuration exacte du manuscrit (pavage 6/2, N_h = 3088).

| modèle | X_m0 (Ω) | X_m charge | I_0 (A) | E_1 (V) | couple (N·m) |
|---|---|---|---|---|---|
| Carter | 71.990 | 50.648 | 8.584 | 378.53 | 116.530 |
| P0 (publié) | 64.781 | 43.214 | 9.893 | 375.80 | 115.115 |
| **P1 (convergé)** | **61.015** | **41.249** | **10.296** | **374.95** | **114.655** |
| Référence EF | — | 46.0 | 8.49 | 382.1 | 121.63 |

### Écarts à la référence

| modèle | X_m charge | I_0 | E_1 | couple |
|---|---|---|---|---|
| Carter | +10,1 % | **+1,1 %** | −0,9 % | −4,2 % |
| P0 (publié) | **−6,1 %** | +16,5 % | −1,7 % | −5,4 % |
| P1 (convergé) | −10,3 % | +21,3 % | −1,9 % | −5,7 % |

### Le prix isolé : P0 → P1

| grandeur | P0 | P1 | variation |
|---|---|---|---|
| X_m0 non saturé | 64.781 | 61.015 | **−5,81 %** |
| X_m saturé charge | 43.214 | 41.249 | −4,55 % |
| courant à vide | 9.893 | 10.296 | **+4,08 %** |
| f.e.m. entrefer | 375.795 | 374.953 | −0,22 % |
| couple | 115.115 | 114.655 | −0,40 % |

### Correction de mes estimations

Mon rapport annonçait **+23 %** sur I_0 et **≈ −11,5 %** sur X_m. Les valeurs
calculées sont **+21,3 %** et **−10,3 %**. J'extrapolais en 1/X_m au lieu de
résoudre le schéma équivalent avec la réactance saturée au point nominal.

### Le fait à regarder en face

**P1 est moins bon que P0 sur les cinq colonnes.** Aucune grandeur n'y gagne.

**Le seul argument pour P1 est l'existence, pas la précision** — et il tient :
T17 démontre en forme fermée qu'en P0 le terme propre vaut `−ln|2 sin 0|`, la
série diverge, et les 64,78 Ω ne sont qu'un point sur une courbe sans limite.
**Un nombre qui varie avec la troncature n'est pas une prédiction, même s'il
tombe plus près de la référence.**

**Et aucun des trois modèles ne domine** : Carter est le meilleur sur I_0
(+1,1 %) et le pire sur X_m saturé (+10,1 %). L'article doit dire lequel il
retient et pourquoi ; le choix ne peut plus rester implicite.

---

## B1 — Chaîne MAS complète en base P1

Script : `MEC_IM/RUN_B1_IM_P1.m`, sortie `B1_im_p1_out.txt`.
Configuration : pavage nT=17, nO=4, **N_h = 8192, base P1**, opérateur 92×92,
**X_m0 = 61,015 Ω** (cohérent au millième avec B2 et T16, chemins différents).

### 1. Schéma équivalent (s = 0,0188)

| paramètre | P1 | EF |
|---|---|---|
| R_s (Ω) | 0.4302 | 0.4450 |
| R'_r (Ω) | 0.4234 | 0.4400 |
| X_σs (Ω) | 0.8673 (`Lk.Xs_slot`) | — |
| X_σr (Ω) | 0.7325 (`r.Xr`) | — |
| X_m saturé (Ω) | **41.249** | 46.0 |
| R_fe (Ω) | 1787.2 | 1740 |

### 2. Point nominal

| grandeur | P1 | EF | écart |
|---|---|---|---|
| couple (N·m) | 114.656 | 121.53 | −5,7 % |
| I₁ (A) | 19.084 | 19.70 | −3,1 % |
| rendement | 0.9190 | 0.9160 | +0,3 % |
| cos φ | 0.8205 | 0.8640 | −5,0 % |
| pertes fer (W) | 227.37 | 232.6 | −2,2 % |

### 3. Trois conditions

| | P1 | EF | écart |
|---|---|---|---|
| **en charge** couple | 114.656 | 121.63 | −5,7 % |
| I₁ | 19.084 | 19.73 | −3,3 % |
| I barre | 296.61 | 324.7 | −8,6 % |
| pertes fer | 227.37 | 232.6 | −2,2 % |
| **à vide** I_m | **10.296** | 8.49 | **+21,3 %** |
| f.e.m. | 374.95 | 382.1 | −1,9 % |
| pertes fer | 235.67 | 249.3 | −5,5 % |
| **rotor bloqué** couple | 102.60 | 104.31 | **−1,6 %** |
| I₁ | 109.53 | 108.21 | +1,2 % |
| I barre | 1966.7 | 1929 | +2,0 % |

### 4. Champs à mi-entrefer

| | P1 | EF | écart |
|---|---|---|---|
| B_g1 à vide (T) | 0.9452 | 0.942 | **+0,3 %** |
| B_t rms à vide (T) | 0.1201 | 0.128 | −6,1 % |
| B_g1 en charge (T) | 0.8983 | 0.920 | −2,4 % |
| B_t rms en charge (T) | 0.1188 | 0.131 | **−9,3 %** |

### 5. Caractéristique

**Bilan de puissance B6 : 1,34×10⁻⁵** sur 30 glissements — l'assemblage en
base P1 est validé. Décrochage **362,1 N·m à s = 0,115** (EF 324,9 à s = 0,105).

---

## CORRECTION de l'affirmation de B3

**B3 concluait que « P1 est moins bon que P0 sur les cinq colonnes ». C'est
vrai pour les cinq grandeurs de B3 — toutes liées à la branche magnétisante —
mais FAUX en général.**

| grandeur | P0 | P1 |
|---|---|---|
| rotor bloqué, couple | −4,6 % | **−1,6 %** |
| B_t rms en charge (chemin dentaire) | −12,6 % | **−9,3 %** |
| B_g1 à vide | — | **+0,3 %** |

**Formulation correcte** : la base P1 **concentre** l'erreur sur la branche
magnétisante — courant à vide et X_m saturé — et **améliore** tout ce qui en
dépend peu : le calage, le champ fondamental, la composante tangentielle.

C'est cohérent avec T12, qui établit que l'erreur au calage vient de
`bar_skin` et non du couplage : à s = 1, I_m ne vaut que 1 % de I₂, donc une
réactance magnétisante plus faible y pèse peu.

**C'est une affirmation plus précise et plus défendable que « P1 coûte en
précision ».**

### Reste ouvert dans B1

Les **sondes locales (Table 17)** exigent le maillage `mesh_refined`, chaîne
distincte du réseau de performance. Non produites ici.

---

## PRIORITÉ 4 — A4, A5, A6, B4, B6

**Trois de ces cinq tâches corrigent des affirmations de mes rapports
précédents.** Elles sont consignées comme telles.

---

### A5 — RÉSOLU : les trois cellules NE SONT PAS irréproductibles

Script `RUN_A5_LUMPED.m`, diary `A5_lumped_out.txt`.

**Ma conclusion « provenance perdue, aucune chaîne ne les produit » était
FAUSSE.**

L'agent a posé un **test d'identité interne** : dans `krkl.m`,
`k_r = mean(mm)/mean(mg)` par construction. Or les FMM publiées donnent
k_r = 1,05445 alors que le manuscrit publie 1,05220 — le quotient des FMM du
programme **à φ = 0**. Le système est donc surdéterminé : aucune configuration
unique ne peut satisfaire les deux. Cela oriente vers le seul paramètre qui
déplace les FMM sans être un paramètre de modèle : **la position rotorique**.

`outF.txt` déclare N_surf = 1260, N_p = 721, et §2 de `BLDC_MEC_COMPLET`
évalue les FMM **à la position identifiée**, non à φ = 0. Rejouée :

| grandeur | réexécuté | cible `outF` | écart |
|---|---|---|---|
| B_r crête (T) | 0.98458 | 0.9846 | −0,0024 % |
| FMM aimant (A) | 763.24273 | 763.2427 | **0,0000 %** |
| FMM entrefer (A) | 723.82915 | 723.8292 | **−0,0000 %** |
| k_l | 0.85377 | 0.8538 | −0,0040 % |

**Correspondance parfaite à φ = 39,429° mécaniques** (indice 553/721).

**Le défaut réel, plus précis** : la colonne *Lumped* évalue les **FMM** à la
position identifiée mais **k_r** à φ = 0. Elle mélange deux positions.

---

### A4 — mon −14,8 % était FAUX ; correct : +25,3 %

Script `RUN_A4_IRONLOSS.m`, diary `A4_ironloss_out.txt`.

**Trois défauts compensatoires** dans mon calcul T18, isolés par l'agent :
sommation sur les **branches** ×1,814, k_Fe ×0,900, fréquence ×0,830 (1500 au
lieu de 1688 tr/min).

| n_sh | inconnues | stator (W) | éc./EF | +rotor (W) | éc./EF | peak B |
|---|---|---|---|---|---|---|
| 1 | 5400 | 22.576 | +14,4 % | 24.436 | +23,8 % | 1.9742 |
| 2 | 6480 | 22.854 | +15,8 % | **24.713** | **+25,3 %** | 1.9779 |
| 4 | 8640 | 22.996 | +16,6 % | 24.856 | +26,0 % | 1.9791 |

**Décomposition** : hystérésis 40,8 %, Foucault 59,2 %, excès 0,0 % — ce
dernier **par donnée** (k_e = 0 dans `BLDC.aedt`), non par calcul.
**Répartition** : dents 76,6 %, culasse 23,4 % ; dont **bec 10,0 %**.
Le modèle localisé **n'a pas de région de bec** : sa dent est un bloc.

#### Conclusion §8.5 — plus nuancée que mon hypothèse

| | dérive n_sh 1→4 |
|---|---|
| **total** des pertes fer | **+1,86 %** |
| **part de bec** seule | **+20,28 %** |
| bande ν=8 (rappel A2) | doublement |

**Publier « les pertes dérivent comme ν=8 » serait faux.** Les pertes sont la
somme d'une part volumique **convergée** (90 %) et d'une part de coin qui **ne
converge pas** (10 %). Ici la part de coin est petite, donc le total paraît
robuste — **mais c'est une propriété de la répartition des pertes de cette
machine, non une garantie.** Sur une machine où le bec porte davantage, le
même défaut dominerait.

**Et la saturation masque la singularité sans la supprimer** : en Newton le
champ crête reste borné (1,9742 → 1,9791 T) ; en linéaire le même coin diverge
(2,7989 → 3,4709 T, cf. A2). Elle borne l'amplitude, pas la dépendance au
maillage.

---

### B4 — la prémisse de T12 est INVALIDE, et le signe est inversé

Script `RUN_B4_BARSKIN.m`, diary `B4_barskin_out.txt`.

**« R'_r 0,4234 contre 0,4400, soit −3,8 % » n'est pas une comparaison
valide** : 0,4234 est la valeur au **point nominal** (s = 0,0188), 0,4400
appartient à **s = 1**. Deux glissements différents.

| | valeur |
|---|---|
| R'_r chaîne, s = 0,0188 | 0.423392 |
| R'_r chaîne, s = 1 | **0.467080** |
| référence implicite, s = 1 | 0.452834 |
| **écart au calage** | **+3,15 %** |

**Le MEC SURESTIME R'_r au calage. `bar_skin` ne le sous-estime pas.**

Ma conclusion de T12 — « l'erreur au calage vient de `bar_skin`, qui
sous-estime R'_r de 3,8 %, ce qui explique la moitié du −4,6 % » — est fausse
sur le signe et sur l'amplitude.

---

### B6 — décomposition exacte, et il résout la question de B4

Script `RUN_B6_ANNEAU.m`, diary `B6_anneau_out.txt`.

**Contrôle fondateur** : le rapport fondamental du réseau reproduit le modèle
idéal `1/(2 sin(πp/N_r)) = 3.513337` à **2,22×10⁻¹⁶ près**. L'écart de 0,2 %
sur le RMS vient des branches de cage **harmoniques**, dont le report d'anneau
est `1/(2 sin(πνp/N_r))`.

| source | condition | rapport | écart/idéal |
|---|---|---|---|
| réseau | en charge | 3.5053 | −0,23 % |
| réseau | calage | 3.5042 | −0,26 % |
| **référence** | en charge | 3.3589 | **−4,39 %** |
| **référence** | calage | 3.3842 | **−3,68 %** |

#### Décomposition exacte : e_anneau = e_barre + e_rapport + croisé

| condition | e_barre | e_rapport | e_anneau calculé | publié |
|---|---|---|---|---|
| en charge | −8,10 % | +4,33 % | −4,12 % | −4,1 % |
| calage | +0,78 % | +3,48 % | +4,29 % | +4,3 % |

**Bascule** : ligne barre **+8,88 pts**, terme d'anneau **−0,84 pts** →
**105,5 % de l'inversion est portée par la barre.** Le terme d'anneau garde le
même signe dans les deux conditions.

**À qui appartient e_rapport** : écart du réseau au modèle idéal −0,23 % et
−0,29 % ; écart de la **référence** −4,36 % et −3,65 %. **C'est un écart de la
référence, pas une erreur du réseau.**

Conclusion tenue avec **et sans vrillage** (103,6 % et 106,2 %) — elle ne
dépend ni du millésime du modèle ni du vrillage.

#### B6 résout la question de provenance de B4

Les deux couples EF au calage ne sont **pas deux références concurrentes** :

| essai | T (N·m) |
|---|---|
| transitoire rotor bloqué, moyenne régime établi | 104.309 |
| caractéristique T(s) au point s = 1 | 98.982 |
| MEC au calage (vrillage ON) | 102.602 |

**Deux essais différents, à ne pas confondre.** Le « problème de provenance »
de B4 est une confusion d'essais dans les scripts, non une ambiguïté de
référence.

---

### A6 — ligne N2b, mise en forme finale

Script `RUN_CARTER_CMP.m` **enrichi en place** (et non dupliqué : deux chaînes
pour une même grandeur reproduiraient le défaut que C1 vise), diary
`A6_n2b_out.txt`.

| variante | 2π/a_r | γ_m (deg) | γ_z (deg) | γ_m/γ_z | forme |
|---|---|---|---|---|---|
| a_r = 0,25 a_t | 69.58 | 7.7608 | 12.9347 | 60,0 % | plateau + gaussienne |
| a_r = 0,50 a_t | 34.79 | 5.1739 | 15.5217 | 33,3 % | plateau + gaussienne |
| a_r = 1,00 a_t | 17.40 | **0.0000** | 20.6956 | **0,0 %** | **gaussienne seule** |

**À a_r = a_t, γ_m = 0 exactement : le plateau disparaît et l'éq. (14) change
de NATURE sous l'effet du maillage seul.** Argument plus fort que la
dispersion de 219 %, car il porte sur l'identité du modèle, non sa précision.

**Le nombre d'éléments n'est pas entier** (69,58 / 34,79 / 17,40 par tour) :
preuve qu'il s'agit d'un choix de discrétisation et non d'un découpage de la
machine, qui a 15 dents.

**⟨λ⟩ = 1,000000** pour les quatre modulations (écart max 6,66×10⁻¹⁵) ;
B_g1 modulés identiques à B_g1(N1) à 2,31×10⁻⁷ T.

**Limite déclarée** : la colonne B_t des trois N2b **n'est pas calculée
indépendamment** — elle reprend B_t(N1) par construction, une perméance
scalaire ne pouvant générer aucune composante tangentielle. Les trois −19,7 %
sont une propriété du modèle, pas trois mesures.

---

## PRIORITÉ 5 (C1–C4) et PRIORITÉ 6 (A7, B5)

### C1 ✅ — déjà fait avant mon intervention

`MANIFEST.md` existe (6 août, 8 sections) et va au-delà de v2 : règle de
lecture des `.txt` (« le **dernier** bloc complet fait foi »), utilitaires
promus, décision de vocabulaire. `inductance_mec2.m` y est identifiée comme
**seule variante réellement morte**.

### C2 ✅ — corrigé, et le diagnostic ferme enfin

**En-tête réécrit** avec son historique plutôt que supprimé.

**Le calcul mort (`*0`) est remplacé par une décomposition ÉNERGÉTIQUE :**

```
AIR : branches (encoche + becs)    60.0 %
AIR : couronne d'entrefer          38.5 %
FER : dents + culasses              1.5 %
CONTROLE 2W vs lambda : 50.4837 = 50.4837 mH   (residu -1.3e-12)
```

**Le maillage porte 60,0 % de l'énergie dans l'air d'encoche ; le modèle
localisé y ajoute une perméance analytique valant 60 %.** Deux routes
indépendantes, le même nombre. **La §3.5 est enfin étayée par un calcul actif.**

**Quatre pièges traversés, chacun attrapé par le contrôle `2W = λ`** : λ = ΣF·Φ
est une forme **duale** (0 % d'air par construction) ; terme source omis dans
la chute de réluctance ; **l'entrefer n'est pas une branche** mais un bloc
dense (19,45 mH manquants) ; signe de la couronne, assemblée en `−Y`.

La part « surface d'entrefer » est **retirée, non réparée** : les branches de
surface portent E = 0, il n'y avait pas de calcul à rétablir.

### C3 ✅ ⚠ — les deux σ viennent de deux essais

Définition trouvée, `RUN_COG_MATCH.m:49` : `sig = std(res)*sqrt(2/numel(th))`
— **amplitude spectrale** de bruit blanc, non écart-type temporel.

**⚠ Ma prémisse de T15 était fausse.** Les deux valeurs ne partagent pas une
normalisation parce qu'elles viennent de **deux essais différents** :
77,0 mN·m en magnétostatique **remaillée**, 0,88 mN·m en transitoire à maillage
**préservé**. Chacune est correctement normalisée à **son propre** plancher.
`RUN_COG_MATCH` donne 2σ = **±1,401 mN·m** pour la référence transitoire.

**Défaut de rédaction, non de calcul** : les deux valeurs sont rétablissables
à condition d'énoncer σ par essai.

### C4 ✅ — LES QUATRE TABLES RÉGÉNÉRÉES

Chacune a désormais son `.mat` à pleine précision.

#### Table 8 — spectre spatial de B_r (`C4c_table8.mat`)

| ordre | réseau | EF | écart plein | publié |
|---|---|---|---|---|
| p = 7 | 1.0775474 | 1.0745514 | **+0,2788 %** | +0,3 % |
| \|N_s−p\| = 8 | 0.0180619 | 0.0196579 | **−8,1190 %** | −8,1 % |
| 3p = 21 | 0.2593437 | 0.2656534 | **−2,3751 %** | −2,4 % |
| N_s+p = 22 | 0.0352377 | 0.0327597 | **+7,5644 %** | +7,6 % |
| \|2N_s−p\| = 23 | 0.0325496 | 0.0202107 | **+61,0513 %** | +61,1 % |
| 5p = 35 | 0.1047221 | 0.1066984 | **−1,8523 %** | −1,9 % |
| 2N_s+p = 37 | 0.0426914 | 0.0270805 | **+57,6465 %** | +57,6 % |

**Les sept ordres reproduisent le publié.** Cela règle une des huit cellules
signalées par T14 : l'ordre 22 y figurait parce que le calcul naïf donnait
+7,32 %. La pleine précision donne **+7,5644 %**, qui s'arrondit bien à +7,6.
**La valeur publiée était juste ; c'est mon test qui manquait de résolution.**

Hiérarchie confirmée : harmoniques de travail ≤ 2,4 %, premières bandes ≈ 8 %,
**secondes bandes ≈ 61 %** — elles échantillonnent le coin à une longueur
d'onde valant la moitié de l'ouverture d'encoche.

#### Table 10 — transitoire en charge PMSM (`T10_full.mat`)

Obtenue en **instrumentant `BLDC_MEC_COMPLET`** (section 5b) : ces grandeurs
n'existent que dans la chaîne d'entraînement et n'étaient imprimées qu'à la
précision d'affichage.

| grandeur | MEC | FEA | écart |
|---|---|---|---|
| vitesse établie (tr/min) | 1257.029231 | 1339.826255 | **−6,1797 %** |
| puissance d'arbre (W) | 651.466845 | 682.820131 | **−4,5917 %** |
| couple électromagnétique (N·m) | 4.949010 | 4.866870 | **+1,6877 %** |
| pertes cuivre (W) | 69.431939 | 71.563994 | **−2,9792 %** |
| courant de phase rms (A) | 1.418566 | 1.551620 | **−8,5752 %** |
| courant déduit des pertes (A) | 1.504682 | 1.527609 | **−1,5009 %** |
| pertes fer (W) | 19.083041 | 18.464992 | **+3,3471 %** |
| pertes aimant (W) | 1.592913 | 3.182682 | **−49,9506 %** |
| rendement (%) | 86.332272 | 86.514779 | **−0,2110 pt** |

**Deux notes inscrites dans la structure sauvegardée** :
- Le **courant de phase a deux mesures qui ne sont pas la même grandeur**. Le
  rms direct est biaisé bas si la fenêtre ne couvre pas un nombre entier de
  périodes ; le courant déduit ne l'est pas — mais **il n'est pas un contrôle
  indépendant** de la ligne « pertes cuivre », il en est la racine carrée.
- Le rendement s'écarte de **0,2110 point**, non de 0,2 % — points de
  pourcentage, non écart relatif.

**Bug de mon instrumentation, corrigé** : j'avais pris `dOL(:,3)` pour le
couple EF alors que c'est la perte fer en mW (18465 N·m). Le couple est formé
depuis la puissance d'arbre et la vitesse, comme le MEC forme `P_out = T·ω`.

#### Tables 16 et 17 — MAS (`C4b_tables_16_17.mat`)

**Table 16**, base P1, valeurs pleines : couple en charge 114.655478
(−5,7342 %), I₁ 19.084155 (−3,2734 %), I barre 296.614417 (−8,6497 %),
I anneau 1042.106434 (−4,4815 %), P entrée 18713.545689 (−8,1273 %),
I_m à vide 10.296202 (**+21,2745 %**), E₁ 374.953104 (−1,8704 %),
couple au calage 102.602253 (**−1,6372 %**).

**Contrôle des références EF relues** : I₁ en charge **19,7253** contre 19,73
publié ; I₀ à vide **8,4886** contre 8,49 ; E₁ **382,1209** contre 382,1.

**Table 17** — la seule que B1 ne pouvait pas produire, car elle relève de la
chaîne **`mesh_refined`**, distincte du réseau de performance :

| région | MEC vide | EF vide | écart | MEC charge | EF charge | écart |
|---|---|---|---|---|---|---|
| culasse stator | 1.879314 | 1.840 | +2,1366 % | 1.798220 | 1.900 | −5,3569 % |
| **dent stator** | 1.968712 | 1.630 | **+20,7798 %** | 1.902159 | 1.680 | +13,2238 % |
| **dent rotor** | 2.151882 | 1.930 | **+11,4965 %** | 2.087811 | 1.870 | +11,6476 % |
| culasse rotor | 1.710736 | 1.670 | +2,4393 % | 1.585654 | 1.580 | +0,3578 % |

**En base P1 la Table 17 diffère du publié** : dent stator +20,78 % contre
+19,0 %, dent rotor +11,50 % contre +10,1 %. La base P1 relève les champs
locaux de dent. À reporter si P1 est retenue.

---

### A7 ⬜ — traité dans `MANIFEST.md` §6

Scission de `BLDC_FIG5` en `_1_transient`, `_2_conversion`, `_3_saturation`.
Libellés anglais des FIG 4, 5, 6 faits.

### B5 ✅ ⚠ — la §6.7 est fausse CÔTÉ MEC, non côté référence

Déjà tranché le 6 août (`B5_skew_out.txt`, `+mec/ansys_ref.m:13-18`).

```
Cote REFERENCE  : VRAI. NumberOfSlices = 1 ne produit aucun
  moyennage axial ; la section est droite.
Cote MEC        : FAUX. mec.cage applique ksq = 0.997147 au
  fondamental ET ksq_nu aux harmoniques, sans que
  M.opt.skew_harm ne soit jamais mis a 0 dans la chaine.
```

**⚠ Mon inférence aurait été inverse.** J'allais conclure du couple au calage
(MEC vrillé 102,6 contre non vrillé 85,7, EF 104,31) que la **référence** était
vrillée. C'est le **MEC** qui l'est.

**Données de la voie (a)**, produites par `RUN_B10_B1_SKEWOFF` :

| grandeur | vrillage ON | vrillage OFF | EF |
|---|---|---|---|
| I_m à vide | +21,3 % | +21,59 % | 8,49 A |
| B_g1 à vide | +0,3 % | +0,37 % | 0,942 T |
| B_t rms en charge | −9,3 % | −9,23 % | 0,131 T |
| cos φ | −5,0 % | −5,05 % | 0,864 |
| **décrochage** | 362,1 (+11,4 %) | **344,7 (+6,1 %)** | 324,9 N·m |
| **couple au calage** | −1,6 % | **−17,84 %** | 104,31 N·m |

**Sur le fondamental, neutraliser le vrillage ne change rien** (< 0,1 point) —
confirmation directe des 0,57 % de B5 et des 0,3 % que T10 avait estimés.

**Le décrochage s'améliore** : +11,4 % → **+6,1 %**. Argument en faveur de la
voie (a) que le diagnostic n'avait pas anticipé.

**Le calage s'effondre** (−1,6 % → −17,84 %) : à s = 1 tous les glissements
harmoniques valent 1 et la cage harmonique gouverne. Ce point est **déjà exclu
du tableau de comparaison** (§5.6 Article II) — ce n'est donc pas un obstacle,
à condition que l'exclusion soit maintenue et déclarée.

**Décision à prendre** : (a) neutraliser le vrillage MEC et republier
l'ondulation — comparaison loyale, mais recalibration C5 à refaire ; ou
(b) le conserver et déclarer que la référence ne l'a pas — honnête, mais
l'accord sur l'ondulation perd sa valeur de validation.

---

## C5 — calibration de l'ondulation : le vrillage est SANS OBJET, et deux découvertes

Scripts : `RUN_C5_RIPPLE_SKEWOFF.m` (réseau — abandonné), puis
`RUN_C5B_RIPPLE_CHAIN.m` (maillage — la bonne chaîne). Diary
`C5b_ripple_chain_out.txt`.

### La question posée : neutraliser le vrillage MEC change-t-il l'ondulation ?

**Non.** Mesuré dans la configuration exacte de `RUN_ARTICLE` — glissement
déduit de la vitesse mesurée (s = 0,020151), 19 positions sur deux pas
d'encoche, charge **mesurée** `[trq(:,1),trq(:,2)]`, fenêtre `t ≥ 1,90` :

| | carte brute | après `dq_startup` |
|---|---|---|
| effet de `skew_harm = 0` | **+0,015 %** | **+0,46 %** |

**Poser `skew_harm = 0` dans la chaîne de production est licite et sans
conséquence.** La raison est structurelle : la carte est **magnétostatique et
n'utilise que le fondamental** (`inst_currents` ne lit que `r.I1c`, `r.I2c`),
alors que `skew_harm` agit sur les branches de cage **harmoniques**, qui ne
l'alimentent pas.

**La §6.7 est donc défendable telle quelle.** Le facteur « plusieurs
centaines » que B5 mesure porte sur le **schéma équivalent** et ne se
transporte pas à la carte d'ondulation. **La voie (a) recommandée par B5 est
sans objet pour cette grandeur** — il n'y a rien à recalibrer.

### Découverte 1 — les 101,4 N·m sont PÉRIMÉS

```
RUN_ARTICLE, chaine actuelle (base P0 par defaut) :
  carte en charge   : 112.1 N.m brut
  apres dq_startup  : 108.2 N.m   vs FEM 105.5  ->  +2.6 %
```

**Le manuscrit publie 101,4 N·m, soit −3,9 %. La chaîne actuelle donne
108,2 N·m, soit +2,6 %.** Le signe de l'écart change.

**Cinquième cellule non régénérable**, après les trois cellules *Lumped*
(résolues par A5) et les pertes aimants à 0,342 W (non résolues).

### Découverte 2 — la base P1 dégrade l'ondulation de +2,6 % à +20,9 %

| configuration | carte brute | après dq | vs FEM 105,5 |
|---|---|---|---|
| **P0** (défaut, `RUN_ARTICLE`) | 112.1 | 108.2 | **+2,6 %** |
| **P1** (nT=17 nO=4 N_h=8192) | 132.8 | 127.5 | **+20,9 %** |

Cause identifiée : P1 relève X_m0 de 43,2 à **61,0 Ω**, ce qui augmente les
courants et donc la modulation de denture. **+18,5 % sur la carte.**

**Cela ajoute une SIXIÈME grandeur au tableau du prix de B3**, qui ne comptait
que cinq quantités de la branche magnétisante.

**Bilan corrigé de la base P1** :
- **améliore** : couple au calage (−4,6 → −1,6 %), B_t rms (−12,6 → −9,3 %),
  B_g1 à vide (+0,3 %) ;
- **dégrade** : courant à vide (+16,5 → +21,3 %), X_m saturé (−6,1 → −10,3 %),
  **ondulation de denture (+2,6 → +20,9 %)**.

### Quatre hypothèses éliminées en chemin

| hypothèse | effet mesuré | verdict |
|---|---|---|
| glissement 0,0188 → 0,0202 | T_pp 126,8 → 132,8 | **aggrave** |
| base P0 → P1 | +18,5 % | **aggrave** |
| filtrage inertiel (J = 0,17) | atténuation **0,960** seulement | insuffisant |
| charge mesurée vs constante | 127,7 → 127,5 | négligeable |

Aucune n'expliquait l'écart aux 101,4 — **parce que ces 101,4 n'existent
plus**.

### Erreurs de ma séquence C5, consignées

| # | erreur | détectée par |
|---|---|---|
| 1 | ondulation reconstruite sur le **réseau dentaire** au lieu du maillage (146–158 N·m) | écart aux 101,4 |
| 2 | FMM rotorique **fondamentale seule** → `skew_harm` sans effet (+0,14 %) | insensibilité au test |
| 3 | `Xmn = Xm0/ν²` au lieu de `Xm·H.sig` (saturée) | couple moyen −38,9 N·m |
| 4 | **signe** de la FMM rotorique omis (`torque_vw:46` : « le rotor s'oppose ») | garde couple moyen |
| 5 | carte **brute** comparée à une référence **transitoire** | attenuation mesurée 0,96 |

**La garde « couple moyen contre schéma équivalent » a attrapé 3 et 4**, que
l'inspection du code n'avait pas révélés. Elle est en place dans le script.

**La leçon de T18, que j'ai réenfreinte** : piloter la chaîne existante, ne pas
la refaire. Cinq itérations ont été nécessaires pour revenir à
`RUN_ARTICLE:324-347`.

---

## Modifications de code apportées

| Fichier | Modification |
|---|---|
| `MEC_BLDC/airgap_magnet.m` | argument `basis` ('p0'/'p1'), noyau chapeau `(4/(pi*nu^2*d))*sin^2(nu*d/2)` |
| `MEC_BLDC/cogging_mec.m` | argument `basis` transmis à `airgap_magnet` ; note sur `numax = Nsurf/2` |
| `MEC_IM/+mec/airgap_fourier.m` | argument `basis` (T16) |
| `MEC_IM/+mec/airgap_dtn_tooth.m` | `basis` transmis (T16) |
| `MEC_BLDC/pm_loss.m` | argument `ext` : source externe de potentiels (T18) |

Les deux machines partagent désormais la **même** option de base, condition
pour que toute comparaison mesure la base et non une différence
d'implémentation.

---

## Reste à faire

**Priorité 2** — A1 (Table 7 en configuration unique), A2 (Tables 5a/5b),
A3 (pertes aimants, chaîne gelée `pm_loss` Np ≥ 181).

**Priorité 3** — B1 (chaîne MAS en base P1, le bloc le plus lourd),
B2 (tableau k_C(N_h) publiable), B3 (le prix chiffré).

**Priorité 4** — A4, A5, A6, B4, B6.

**Priorité 5** — C1 à C4, obligatoire avant diffusion du code.

**Priorité 6** — A7 (figures, dont la scission de `BLDC_FIG5_charge`),
B5 (vrillage).
