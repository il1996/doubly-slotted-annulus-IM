# R9 — note de livraison

> Sortie : `article/R9_archive_out.txt` (recopiée dans l'archive)
> Script : `article/RUN_R9_ARCHIVE.m`
> Archive : `article/zenodo/` — **359 fichiers, 6,84 Mo**
> Exécuté le 10 août 2026.

## L'archive est montée, et elle fonctionne

**359 fichiers, 6,84 Mo** — assez petit pour un dépôt Zenodo ordinaire, assez
complet pour être vérifiable sans licence ANSYS.

| dossier | contenu | poids |
|---|---|---|
| `code/` | 201 fichiers — MEC_BLDC (107), MEC_IM (87 dont 37 dans `+mec/`), article (4 + 2 `.tex`) | 1,47 Mo |
| `outputs/` | 61 sorties `diary` `.txt`, telles que produites | 0,33 Mo |
| `reference/` | **84 exports ANSYS `.tab`** | 4,89 Mo |
| `notes/` | manifeste daté + 8 notes R + 2 rapports de synthèse | 0,13 Mo |
| racine | `README.md`, `LICENSE.txt` | |

**La référence EF est embarquée**, et c'est le point qui compte : elle pèse
4,89 Mo, donc rien n'obligeait à la laisser dehors. C'est ce qui sépare « le
lecteur peut lire le code » de « le lecteur peut refaire le calcul ».

**Versions mesurées, non supposées** : MATLAB R2024a (24.1.0.2537033),
**aucune toolbox** — établi par `matlab.codetools.requiredFilesAndProducts`
sur les deux programmes maîtres. ANSYS Electronics Desktop **2023 R1**, lu
dans les deux `.aedt`.

## Le code est archivé inchangé

Octet pour octet celui qui a produit les nombres publiés. Les chemins absolus
de la machine d'origine sont donc encore là, et un utilitaire séparé les
adapte : `code/SET_REFERENCE_PATH.m`, qui liste chaque fichier et chaque ligne
qu'il touche, écrit un `.bak` à côté, et accepte `('dryrun')`.

**Testé : 16 fichiers, 19 lignes.** L'arborescence est un miroir de
l'arborescence de travail — `code/` joue le rôle de `MEC/` — parce que **cinq
gardes relisent leur manuscrit par un chemin relatif**. Une arborescence
« propre » aurait cassé ces gardes.

## Ce que j'ai vérifié moi-même

**Les deux tests du README ont tourné dans l'archive**, sur ses propres
données, après `SET_REFERENCE_PATH` :

| test | article | résultat |
|---|---|---|
| `RUN_R5_NSH` | I | 5400 / 6480 inconnues, $B_{g1}$ 1,07787 / 1,07908 T — **GARDE PASSÉE** |
| `RUN_R8_TABLE2` | II | $I_0$ +21,6 %, $X_m$ −10,6 % — **GARDE PASSÉE** (0,014 et 0,010 pt) |

**La garde de complétude du manifeste passe** aussi : les 20 sorties et les
21 scripts qu'il nomme sont tous présents.

*Elle n'est pas passée du premier coup* : la première exécution a signalé
`RUN_T9_TEMPS.m` et sa sortie manquants — ils vivent dans `article/` et non
dans un dossier machine. C'est exactement ce qu'une garde doit faire.

## Ce que je n'ai pas pu vérifier

**La garde d'acceptation de la v4 n'est pas passée.** Elle exige qu'**un tiers
qui n'a pas écrit le code** régénère une grandeur de chaque article avec le
seul README. Je ne peux pas m'auto-administrer ce test.

L'archive est donc **montée et fonctionnelle, mais non validée comme
reproductible**. Le §4 du README est le test à lui remettre : moins de deux
minutes, et chaque script imprime sa propre garde.

## Huit scripts ne peuvent pas s'exécuter dans l'archive

Cette liste est **produite par le script de montage**, pas écrite à la main :

| dépendance absente | scripts |
|---|---|
| projet ANSYS (`.aedt`, `.aedtresults`) | `RUN_A4_IRONLOSS`, `RUN_B5_SKEW`, `RUN_R3_FE_RATIO`, `RUN_T9_TEMPS` |
| dossier externe de code analytique | `RUN_ARTICLE`, `RUN_COMPARE_ANALYTIC` |
| manuscrit hérité `MEC_DtN_paper_v2.tex` | `RUN_A5_LUMPED`, `RUN_B6_ANNEAU` |

**`RUN_ARTICLE` est le programme maître de l'Article II** : son `addpath` de
la ligne 64 pointe vers un dossier de comparaison analytique hors périmètre.
Commenter cette ligne suffit à le faire tourner, au prix de la colonne
analytique. C'est signalé dans le README plutôt que laissé à découvrir.

## Deux décisions qui vous appartiennent

**Les manuscrits sont dans l'archive**, et c'est un choix qui doit être
confirmé. Ils y sont parce que cinq gardes les relisent pour se vérifier — les
retirer rend l'archive incapable de se contrôler elle-même, y compris pour le
test du §4. **Mais l'accord éditeur Springer prime**, et il n'a pas été
vérifié. S'il l'interdit, supprimer `code/article/*.tex` et le déclarer.

**La licence est un modèle.** La ligne de copyright est vide. Une licence est
un acte juridique qui appartient au titulaire des droits, et je ne l'ai pas
choisie à sa place. Le fichier propose MIT pour `code/` et CC BY 4.0 pour les
données, et signale une question non tranchée : **la redistribution des
exports ANSYS au regard de la licence du solveur.**

## Le DOI

Les deux manuscrits portent déjà `10.5281/zenodo.XXXXXXX` —
`ArticleI_DtN_PMSM.tex:3607` et `ArticleII_Carter_IM.tex:1591` — et la formule
« available from the corresponding author on reasonable request » n'apparaît
plus dans aucun des deux. **Il ne reste qu'à substituer le DOI réel une fois
le dépôt effectué.**

**Aucun dépôt n'a été effectué.** Téléverser publie ces fichiers de façon
durable et indexable ; c'est une action sortante, et elle vous revient.

*Détail de comptage : le journal de montage annonce 359 fichiers, sa propre
copie ayant été ajoutée à `outputs/article/` après coup.*
