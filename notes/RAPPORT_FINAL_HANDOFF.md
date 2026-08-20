# Rapport final — transmission pour la poursuite des deux articles

> **Objet.** Bilan méthodologique complet de la campagne de vérification menée
> sur les programmes MEC des deux machines (PMSM 15/14 750 W, MAS 48/44
> 18,5 kW) et sur les manuscrits qui en découlent. Destiné à qui poursuivra la
> rédaction : il porte les résultats, les corrections, les pièges, et les
> contributions scientifiques que ce travail rend possibles.
>
> **À lire avec** `RAPPORT_SPEC_v2.md` (détail tâche par tâche),
> `MANIFEST.md` (une chaîne par grandeur) et `SPEC_CLAUDE_CODE_v3.md`.

---

# 1. Ce que la campagne a établi

Trente tâches exécutées, chacune avec son script, son `diary` et son `.mat`.
**Aucune valeur transcrite à la main.** Les résultats sont classés ici par
importance **scientifique**, non par ordre d'exécution.

## 1.1 Le résultat central : la base de projection décide de l'existence

L'opérateur de Dirichlet-to-Neumann projette le potentiel de surface sur une
base. La campagne établit que **ce choix, jamais discuté dans la littérature
MEC, décide si la grandeur calculée existe**.

**Base P0** (constante par colonne, celle de tous les travaux publiés) :
la projection décroît en $1/n$, donc le terme d'énergie $K_n|W|^2 \sim 1/n$ —
**série harmonique, divergente**. T17 le démontre en forme fermée : en
décomposant le produit de trois fonctions circulaires en quatre cosinus et en
sommant par $\sum \cos(n\alpha)/n = -\ln|2\sin(\alpha/2)|$, le terme propre
($i=j$, $\Delta\theta = 0$) vaut $-\ln|2\sin 0| = +\infty$. **L'auto-énergie
d'un potentiel discontinu est infinie.** Vérifié numériquement : la croissance
en $\tfrac12\ln N$ est prédite au quatrième chiffre, constante d'Euler
comprise.

**Base P1** (chapeau, convolution de deux créneaux) : projection en $1/n^2$,
énergie en $1/n^3$, **absolument convergente**. La queue résiduelle décroît en
$1/N^2$, vérifié sur trois décades.

**Conséquence chiffrée (MAS)** : $k_C$ dérive de **5,68 %** en P0 et se
stabilise à **0,62 %** en P1, avec 0,02 % sur les trois derniers points.
Valeurs convergées : $k_C = 1{,}3320$, $X_m = 61{,}015\ \Omega$.

**Et l'accord publié avec Carter était une coïncidence de troncature.** La
courbe $k_C(N_h)$ traverse 1,2665 vers $N_h = 3088$ — précisément la troncature
retenue — puis continue de descendre. La valeur convergée est **+5,2 %
au-dessus de Carter**.

> **Ce point transforme le résultat plutôt qu'il ne le détruit.** « L'opérateur
> reproduit Carter à trois décimales » devient « la solution exacte de la
> couronne converge 5,2 % au-dessus de l'approximation de Carter ». C'est une
> **correction quantifiée d'un résultat de 1926**, plus fort que la coïncidence
> qu'elle remplace : Carter suppose une encoche isolée dans un plan infini, et
> l'écart mesure ce que cette approximation néglige — interaction entre
> encoches, et seconde surface encochée.

## 1.2 La réfutation de l'argument de compensation

La §6.4 du manuscrit soutenait que la dérive de $X_m$ se compense dans le
rapport $k_C = X_m(\text{lisse})/X_m(\text{encoché})$. **C'est faux, et le
tableau B2 le montre en deux colonnes** :

| colonne | dérive sur $N_h \in [512, 8192]$ |
|---|---|
| $X_m$ **lisse**, P0 | **+0,15 %** |
| $X_m$ **encoché**, P0 | **+5,99 %** |
| $X_m$ **lisse**, P1 | **+0,00 %** (rigoureusement constant) |

**Rapport des deux dérives : 40.** Le numérateur ne bouge pas, le dénominateur
dérive, le rapport hérite intégralement. La raison est structurelle et
prédictible : une surface **lisse n'a aucun saut de potentiel**, une surface
**encochée en a un par dent** — les deux queues harmoniques ne sont pas de même
nature et n'ont aucune raison de se compenser.

En P1 le cas lisse est **rigoureusement constant** sur les cinq troncatures :
une surface sans discontinuité est représentée exactement par la base chapeau
dès les premiers harmoniques.

## 1.3 Les deux machines ne sont pas dans le même régime — et cela justifie deux articles

**Le PMSM ne souffre pas de la divergence.** Jusqu'à $N = 4320$ — au-delà des
3088 où la MAS dérivait de 9,3 % — $B_{g1}$ reste stable à **0,44 %**, et les
incréments **alternent de signe** (−0,23, −0,35, −0,89, +0,48) : convergence
oscillante, non queue logarithmique. Passer en P1 déplace $B_{g1}$ de
**0,01 %**.

**Deux différences structurelles l'expliquent** :

| | PMSM | MAS |
|---|---|---|
| $X = \ln(R_s/R_{\text{int}})$ | **0,139** | **0,00297** |
| surface | 1260 colonnes uniformes | 736 colonnes **condensées** sur 92 nœuds |

La couronne du PMSM inclut l'aimant de 3,5 mm — un facteur **47** sur $X$ — et
sa surface n'est pas condensée par complément de Schur. Le mécanisme de la
divergence n'y existe pas.

> **Recommandation éditoriale.** Ce contraste justifie la séparation en deux
> articles mieux que ne le faisait l'architecture initiale : ce n'est pas une
> commodité, c'est un **résultat**. L'Article I montre une formulation robuste
> au choix de base ; l'Article II montre où et pourquoi elle cesse de l'être.

**Limite à déclarer** : sur la chaîne PMSM, `cogging_mec` impose
$n_{\max} = N_{\text{surf}}/2$ — la troncature **n'est pas un paramètre libre**.
Le balayage varie pavage et troncature ensemble. Le test pur exigerait de
porter la condensation de Schur au PMSM.

> ### ⚠ CORRECTION — l'explication ci-dessus est FAUSSE
>
> *Ajoutée le 8 août 2026 après relecture croisée
> (`CONSIGNES_COWORK_07AOUT.md` §0). Le texte d'origine est conservé,
> conformément à la pratique du dossier : il documente l'erreur.*
>
> **L'immunité du PMSM ne vient ni du facteur 47 sur $X$, ni de l'absence
> de condensation.** L'Article I §3.5 (`sub:v1`) la réfute déjà, et son
> argument est le bon :
>
> La queue s'installe dès $n \gtrsim 1/X$, soit $n \approx 7$ sur le PMSM
> contre $n \approx 337$ sur la MAS. Elle commence donc **plus tôt**, pas
> plus tard. **Si l'épaisseur était le mécanisme, cette machine dériverait
> davantage, pas moins.** L'explication par le facteur 47 prédit l'inverse
> de ce qu'on observe.
>
> **Le mécanisme réel est le verrouillage troncature–pavage.** La chaîne
> impose $N = M_s/2$ sur un pavage uniforme $d = 2\pi/M_s$, donc le terme
> divergent de la forme fermée de T17 devient
>
> $$\ln N + \ln|2\sin(d/2)| \longrightarrow \ln\frac{M_s}{2} + \ln\frac{2\pi}{M_s} = \ln\pi$$
>
> **indépendant de $M_s$** : raffiner la surface raffine la troncature
> d'autant, et les deux logarithmes s'annulent. Le terme est **stationnaire
> par construction**. Vérifié : le crochet tient à
> $1{,}721940 \to 1{,}721946$ ($\ln\pi + \gamma$) sur un facteur 16 en
> $M_s$, contre **+166 %** sur la chaîne MAS où le pavage est fixe.
>
> **Ce que cela change.** La stabilité observée sur le PMSM n'est pas une
> robustesse de la formulation : c'est un **accident de chaîne**. La
> « limite à déclarer » ci-dessus n'était pas une limite du test — elle en
> **était le résultat**, et je ne l'ai pas reconnue comme telle.
> L'Article I l'énonce correctement : *« locking, not immunity »*, et
> *« reporting the observed stability as robustness would be reporting an
> accident »*.
>
> **Consigne.** Ne pas promouvoir le contraste inter-machines fondé sur le
> facteur 47 ; renvoyer à §3.5 de l'Article I. Le contraste est réel, sa
> cause n'est pas celle-là. La **recommandation éditoriale** ci-dessus
> reste valable — deux articles se justifient — mais **pour ce motif-ci**,
> non pour celui qu'elle invoque.

## 1.4 La singularité de coin : un encadrement, jamais une convergence

Trois faits indépendants, tous mesurés.

**(a) Les deux solveurs encadrent sans atteindre.** Sur la bande de denture
$\nu = 8$ à $n_{sh} = 2$ :

| $M_s$ | linéaire | Newton | peak $B$ lin. | peak $B$ NL |
|---|---|---|---|---|
| 360 | −25,2 % | +11,8 % | 2,7651 | 1,9271 |
| 900 | **−16,0 %** | **+26,2 %** | **3,1522** | 2,0095 |

Le linéaire sous-estime $\nu=8$ et laisse **diverger** $B$ (3,15 T pour une
tôle qui sature vers 2 T) ; le Newton surestime $\nu=8$ et **borne** $B$. Les
deux s'écartent **en sens opposés et de plus en plus** sous raffinement.

**(b) Les deux directions de raffinement agissent en sens contraires.** Le
raffinement **angulaire** fait remonter $\nu = 8$ de −42,3 % à −9,7 % ; le
raffinement **radial** du bec le fait redescendre.

> **Formulation à retenir** : *une grandeur dont la valeur raffinée dépend de
> la direction du raffinement n'est pas sous-résolue — elle n'est pas définie.*

**(c) La hiérarchie spectrale confirme le mécanisme.** Table 8 régénérée :
harmoniques de travail ≤ 2,4 %, premières bandes ≈ 8 %, **secondes bandes
≈ 61 %**. Ces dernières échantillonnent la région de coin à une longueur d'onde
valant la moitié de l'ouverture d'encoche.

**La saturation masque la singularité sans la supprimer** : le Newton borne
l'amplitude (1,9742 → 1,9791 T) mais pas la dépendance au maillage.

## 1.5 Les pertes ne sont pas une grandeur intégrale — mais pas pour la raison attendue

L'hypothèse de travail était que les pertes, dépendant de $B^{1{,}5}$ à $B^2$
**locaux**, dérivent comme la bande de denture. **C'est faux, et le résultat
est plus fin.**

| | dérive $n_{sh}: 1 \to 4$ |
|---|---|
| **total** des pertes fer | **+1,86 %** |
| **part de bec** seule | **+20,28 %** |

Les pertes sont la **somme** d'une part volumique **convergée** (90 %) et d'une
part de coin qui **ne converge pas** (10 %). Ici la part de coin est petite,
donc le total paraît robuste — **mais c'est une propriété de la répartition des
pertes de cette machine, non une garantie.** Sur une machine où le bec porte
davantage (fréquence plus haute, bec plus mince, saturation plus forte), le
même défaut dominerait.

**Décomposition** : hystérésis 40,8 %, Foucault 59,2 %, excès **0,0 % par
donnée** ($k_e = 0$ dans le projet, non par calcul). Dents 76,6 %, culasse
23,4 %, dont **bec 10,0 %** — région que le modèle localisé **ne possède pas**.

## 1.6 La dichotomie « à vide validé / en charge ouvert » n'existe pas

Le manuscrit oppose un point à vide **validant** la chaîne des pertes aimants
(+2,4 %) à un point en charge **ouvert** (−50 %). **Les deux chiffres sont des
artefacts, et ils se compensaient.**

| point | source | écart |
|---|---|---|
| à vide | réseau à une dent | **−39,1 %** |
| à vide | maillage | **−24,7 %** |
| en charge | brut | −50,0 % |
| en charge | **normalisé** $(1257/1340)^2(1{,}419/1{,}552)^2$ | **−32,0 %** |

Le +2,4 % reposait sur 0,342 W, que **six chaînes testées ne reproduisent
pas**. Le −50 % comparait **deux machines à des points de fonctionnement
différents**. Les trois écarts sont du même ordre : **−25 % à −39 %**. Écart
**unique**, distribué différemment.

**§5.7 et §8.4 sont à réécrire** : les quatre hypothèses testées et rejetées en
§8.4 visaient une cible de −50 % qui n'en fait que −32.

## 1.7 Deux lignes de tableau ne mesurent pas ce qu'elles annoncent

**Les courants d'anneau (Table 16).** Le rapport fondamental du réseau
reproduit le modèle géométrique idéal $1/(2\sin(p\pi/N_r)) = 3{,}513337$ à
**$2{,}22\times10^{-16}$ près**. La décomposition exacte
$e_{\text{anneau}} = e_{\text{barre}} + e_{\text{rapport}} + \text{croisé}$
donne :

| condition | $e_{\text{barre}}$ | $e_{\text{rapport}}$ | calculé | publié |
|---|---|---|---|---|
| en charge | −8,10 % | +4,33 % | −4,12 % | −4,1 % |
| calage | +0,78 % | +3,48 % | +4,29 % | +4,3 % |

**105,5 % de l'inversion de signe est portée par la ligne de barre.** Et
$e_{\text{rapport}}$ appartient à la **référence** (−4,36 % / −3,65 %), non au
réseau (−0,23 % / −0,29 %). **La référence n'emploie pas le modèle d'anneau
géométrique idéal** — segmentation, résistance de contact, effets 3-D.

**La résistance rotorique au calage.** La comparaison « 0,4234 contre 0,4400,
soit −3,8 % » **n'est pas valide** : 0,4234 est au **point nominal**, 0,4400 à
$s = 1$. Au calage la chaîne donne 0,4671 contre 0,4528 : le MEC **surestime de
+3,15 %**, il ne sous-estime pas. Le signe est inversé.

## 1.8 Le vrillage : une objection levée, à un prix nul

Le projet EF déclare `UseSkewModel=true`, `SkewAngle=7.5deg`,
`NumberOfSlices=1`. **La référence est donc une section droite** (une tranche
ne produit aucun moyennage axial). Le **réseau**, lui, applique bien un facteur
de vrillage à ses branches de cage — la §6.7 était fausse **côté MEC**, non
côté référence.

**Mais la conséquence est nulle sur la grandeur en jeu** : la carte
d'ondulation est **magnétostatique et pilotée par le fondamental seul**, alors
que `skew_harm` agit sur les branches de cage **harmoniques**. Mesuré dans la
configuration exacte de `RUN_ARTICLE` (charge mesurée, fenêtre $t \ge 1{,}90$) :
**+0,015 %** sur la carte brute, **+0,46 %** après injection mécanique.

**La §6.7 est défendable — mais elle doit énoncer la raison au lieu de
supposer.** La rédaction corrigée est en place dans `MEC_DtN_paper_v2.tex`.

## 1.9 Le coût de la méthode, chiffré

**Temps de calcul (T9).** Sources : profils `*.results\DV*.profile`, hôte
`DESKTOP-3TEJ3FH`, Maxwell2D 2023.1, 2 cœurs.

| base | rapport |
|---|---|
| totaux bruts | 268× (PMSM), 137× (MAS) |
| **par point de la caractéristique** | **30×** |

Les totaux comparent des périmètres différents et flattent le résultat. **La
base homogène est le point de fonctionnement** : 649 s/point en EF contre
21,6 s/point pour le MEC. **C'est ce chiffre qu'il faut publier**, avec ses
trois réserves — matériel non apparié, aucune cible de précision commune,
totaux minorants.

**Le prix de la base P1 (MAS).** Le bilan est **symétrique**, non uniforme :

| P1 **améliore** | P1 **dégrade** |
|---|---|
| couple au calage (−4,6 → **−1,6 %**) | courant à vide (+16,5 → **+21,3 %**) |
| $B_t$ rms (−12,6 → **−9,3 %**) | $X_m$ saturé (−6,1 → **−10,3 %**) |
| $B_{g1}$ à vide (**+0,3 %**) | ~~**ondulation** (+2,6 → **+20,9 %**)~~ ⚠ |

> ### ⚠ CORRECTION — l'ondulation ne figure PAS parmi les dégradations
>
> *Ajoutée le 8 août 2026 (item A-1). La ligne barrée est conservée.*
>
> **Le changement de base n'a aucun effet sur l'ondulation : −0,38 %.**
> Mesuré en pilotant **la même chaîne** avec l'option `basis`
> (`RUN_A1_RIPPLE_BASIS.m`, sortie `A1_ripple_basis_out.txt`) :
>
> | base | carte brute | post-dq | vs EF | $B_t$ rms |
> |---|---|---|---|---|
> | P0 | 132,805 | 127,994 | +21,34 % | 0,1138 |
> | P1 | 132,791 | 127,511 | +20,88 % | 0,1227 |
>
> Le « +2,6 → +20,9 % » comparait **deux chaînes** — le 108,2 de
> `RUN_ARTICLE` en P0, le 127,5 d'une reconstruction en P1 — et non deux
> bases.
>
> **La garde exigée a tranché.** $B_t$ s'améliore nettement sous P1
> (−13,13 → −6,36 %) tandis que l'ondulation ne bouge pas (+21,34 →
> +20,88 %). Les deux ne divergent donc pas : **l'ondulation ne répond
> simplement pas à la base.**
>
> **Mécanisme.** La carte est calculée par le traitement d'entrefer
> **propre à `mesh_refined`** (`me.gapF`), non par l'opérateur dentaire
> `ctx.AG`. Changer la base de `ctx.AG` n'agit que par les **courants**, et
> le couple du schéma ne bouge que de 0,40 %. **L'ondulation n'est pas une
> grandeur de gradient du potentiel de surface** — elle relève d'une autre
> chaîne, et c'est pourquoi elle échappe au raisonnement par analogie
> avec $B_t$.
>
> **Conséquence rédactionnelle.** Publier la valeur P0 n'est **pas** une
> incohérence silencieuse : elle est **insensible à la base**. Mais il faut
> le **déclarer**, avec sa raison. C'est une quatrième issue, préférable
> aux trois envisagées.
>
> **Réserve maintenue** : cette chaîne donne 128 N·m là où `RUN_ARTICLE` en
> donne 108,2. L'écart n'est **pas** expliqué et ne relève pas de la base.
> À élucider avant de citer un chiffre d'ondulation autre que celui de
> `RUN_ARTICLE`.

> **L'argument pour P1 n'est pas la précision, c'est l'existence.** Un nombre
> qui varie avec la troncature n'est pas une prédiction, même s'il tombe plus
> près de la référence. Cette formulation est défendable ; conserver une valeur
> non convergée parce qu'elle flatte ne l'est pas.

## 1.10 Un argument enfin étayé : la fuite résolue plutôt qu'ajoutée

La §3.5 affirme que le maillage **résout** la fuite d'encoche là où le modèle
localisé l'**ajoute** par formule. Cette affirmation reposait sur un
diagnostic **désactivé** (`*0` en fin d'expression). Rétabli par décomposition
énergétique :

```
AIR : branches (encoche + becs)    60.0 %
AIR : couronne d'entrefer          38.5 %
FER : dents + culasses              1.5 %
CONTROLE  2W = lambda : 50.4837 = 50.4837 mH   (residu -1.3e-12)
```

**Le maillage porte 60,0 % de l'énergie dans l'air d'encoche ; le modèle
localisé y ajoute une perméance analytique valant 60 %.** Deux routes
indépendantes, le même nombre — l'une par formule, l'autre par résolution.

---

# 2. Contributions scientifiques possibles

Ce que la campagne rend **publiable au-delà de ce que les manuscrits portent
aujourd'hui**. Classé par rapport valeur/effort.

## C-1. La base de surface comme objet de la formulation *(fort, faisable)*

**Aucun travail MEC ne discute le choix de base de projection du potentiel de
surface.** Tous emploient implicitement P0. La campagne établit que ce choix
décide de la **convergence même** de la réactance magnétisante sur un entrefer
mince, avec une démonstration en forme fermée du mécanisme.

*Matériel disponible* : T16 (mesure), T17 (démonstration), B2 (tableau
publiable), V1 (contre-exemple sur une seconde topologie).

*Ce qui manque* : une base P1 **asymétrique** adaptée aux largeurs réelles de
colonne — le chapeau symétrique plafonne à 0,62 % à cause d'une non-uniformité
résiduelle de 2,5 % du pavage.

## C-2. La régularisation analytique de la queue *(fort, moyen)*

T17 fournit la forme fermée de la somme infinie. Elle permet de **remplacer la
troncature par une régularisation explicite** : sommer analytiquement tout ce
qui dépasse $N_h$ au lieu de le négliger. Combinée à C-1, elle donnerait une
formulation **sans paramètre numérique libre** — ce qui est l'aboutissement
naturel de la thèse « sans constante ajustée ».

## C-3. Le critère de définition d'une grandeur locale *(fort, faible effort)*

Les trois faits de §1.4 constituent, ensemble, un **critère opérationnel** :
une grandeur dont la valeur raffinée dépend de la **direction** du raffinement,
et que deux solveurs encadrent en s'en écartant tous deux davantage, n'est pas
sous-résolue — elle n'est pas définie. **Ce critère est transposable** à toute
comparaison MEC/EF sur des grandeurs de coin (bandes de denture, détente, bruit
acoustique).

*Matériel* : A2 (les deux solveurs, les deux directions), Table 8 (hiérarchie
spectrale), A4 (part de bec).

## C-4. L'audit de la référence comme méthode *(déjà amorcé, à généraliser)*

La procédure en cinq tests de la §7 est, de l'avis même de la campagne, aussi
réutilisable que le résultat de formulation. Elle a démontré son rendement :
elle a établi que la référence de détente était **inexploitable** (raie à
l'ordre 42, impossible pour une machine 15/14) et que son modèle d'anneau
**diffère du modèle géométrique idéal de 4 %**, de façon reproductible.

*Enrichissement possible* : ajouter un sixième test — **cohérence interne des
rapports géométriques** (anneau/barre, ici $2{,}22\times10^{-16}$ côté modèle
contre 4 % côté référence). C'est un test qui ne coûte rien et qui distingue un
écart de modèle d'une erreur de précision.

## C-5. La décomposition énergétique comme validation de « résolu vs ajouté »
*(moyen, faible effort)*

Le contrôle $2W = \lambda$ de §1.10 est une **preuve** que la fuite est portée
par des cellules d'air et non par une formule. Il est généralisable à toute
affirmation du type « telle constante ajustée a été supprimée » : il suffit de
montrer que l'énergie correspondante est portée par des degrés de liberté
résolus.

## C-6. Le contraste inter-machines comme résultat *(faible effort, fort rendement éditorial)*

$X$ différant d'un facteur 47, condensation présente ou absente : ces deux
paramètres **prédisent** lequel des deux régimes s'applique. Énoncés comme
critère — *à quelle condition un opérateur de couronne condensé converge-t-il ?*
— ils donnent au lecteur un moyen de savoir, avant de calculer, dans quel cas
il se trouve.

---

# 3. Ce qui reste à faire

## 3.1 Réécritures imposées par les résultats

| section | ce qui change |
|---|---|
| Résumé, §6.2, §6.4, §9 | $k_C$ convergé 1,3320, **+5,2 %** au-dessus de Carter ; divergence P0 **exacte** |
| §5.7 et §8.4 | la dichotomie à vide / en charge **n'existe pas** |
| §8.5 | pertes = part volumique convergée **+** part de coin non convergée |
| §8.1 | encadrement linéaire/Newton, deux directions de raffinement |
| §6.7 | vrillage : énoncer la raison (fait dans `v2`) |
| Table 16 | lignes d'anneau = **écart de modèle**, non précision |
| Résumé | « 0,1 % sur $L_d$ » → **0,3 %** à configuration unifiée |

## 3.2 Décisions qui appartiennent aux auteurs

1. **Quel modèle l'article retient** — aucun des trois (Carter, P0, P1) ne
   domine. Carter est le meilleur sur $I_0$ (+1,1 %) et le pire sur $X_m$
   saturé (+10,1 %).
2. **Quel $n_{sh}$** pour la colonne *Mesh* de la Table 7 — aucun ne domine.
3. **Déclarer le solveur ligne par ligne** : une configuration vraiment unique
   est **impossible** sans changer ce qui est comparé (les inductances de
   référence viennent d'un essai à 1 A non saturé).
4. **Les 0,342 W de pertes aimants** : adopter 0,204 W (`pm_loss`, $N_p \ge 181$).
5. **Archiver ou maintenir** `MEC_DtN_paper*.tex` : trois copies héritées
   coexistent avec `ArticleI_DtN_PMSM.tex` et `ArticleII_Carter_IM.tex`. Tant
   qu'elles coexistent, **chaque correction doit être faite quatre fois** —
   c'est ainsi que les 101,4 N·m ont survécu à leur péremption.

## 3.3 Calculs restants

- **Sondes locales de la Table 17** en base P1 : produites (§C4b), à reporter
  si P1 est retenue — dent stator **+20,78 %** contre +19,0 % publié.
- **Base P1 asymétrique** (C-1) : lever le plafond de 0,62 %.
- **Condensation de Schur au PMSM** : rendrait V1 concluant.
- **Tables 8, 10, 16, 17** : régénérées, `.mat` disponibles.

---

# 4. Pièges méthodologiques du dossier

**À lire avant toute reprise.** Chacun a coûté au moins un run faux, et
plusieurs ont survécu plusieurs mois.

| # | piège | garde qui le détecte |
|---|---|---|
| 1 | démarrage à chaud `U0` **illicite** avec un solveur linéaire | étendue nulle du flux sur un balayage |
| 2 | Bertotti **volumique** (W/m³), non massique — facteur 7650 | ordre de grandeur |
| 3 | sommation des pertes sur les **cellules**, non les branches (×1,81) | — |
| 4 | fonctions locales d'un **script** : pas de partage d'espace de travail | erreur d'exécution |
| 5 | `gradient(X,h)` sur une matrice dérive selon la **dim 2** | — |
| 6 | l'entrefer **n'est pas une branche** : bloc dense, signe opposé, énergie $-\tfrac12 U'YU$ | **$2W = \lambda$** |
| 7 | $\lambda = \sum F\Phi$ est une forme **duale** : ne décompose pas « où passe le flux » | 0 % d'air absurde |
| 8 | FMM rotorique : **signe négatif** (« le rotor s'oppose ») | couple moyen négatif |
| 9 | $X_{m\nu} = X_m\,H.\text{sig}$ (**saturée**), non $X_{m0}/\nu^2$ | couple moyen |
| 10 | position rotorique : la colonne *Lumped* évalue les FMM à **φ = 39,429°**, $k_r$ à φ = 0 | test d'identité $k_r = \overline{mm}/\overline{mg}$ |
| 11 | **deux essais EF au calage** : transitoire (104,309) et T(s) à s=1 (98,982) | — |
| 12 | comparer une **carte brute** à une grandeur **transitoire** | facteur d'atténuation mesuré 0,96 |
| 13 | `ws1`/`ws2` ont été ×1000 pendant une partie du projet | affichage sans point-virgule |
| 14 | `Set-Content -Encoding UTF8` écrit un **BOM** en PowerShell 5.1 | premiers octets du fichier |

**Deux règles qui résument le reste.**

> **Piloter la chaîne existante, ne pas la refaire.** Réécrire `pm_loss` a coûté
> un facteur 336 ; reconstruire la carte d'ondulation a coûté cinq itérations et
> quatre chiffres différents. Les fonctions du dossier portent dans leurs
> commentaires des pièges qu'une réécriture réintroduit.

> **Toute affirmation doit avoir une garde qui la contredit si elle est
> fausse.** $2W = \lambda$ a attrapé quatre erreurs successives qu'aucune
> relecture n'avait vues.

---

# 5. Corrections consignées

## 5.1 Erreurs de la campagne — conservées, non effacées

| # | affirmation | correction | source |
|---|---|---|---|
| 1 | pertes fer maillage **−14,8 %** | **+25,3 %** | A4 |
| 2 | trois cellules *Lumped* **irréproductibles** | **résolues**, φ = 39,429° | A5 |
| 3 | `bar_skin` **sous-estime** $R'_r$ de 3,8 % | **surestime** de 3,15 % | B4 |
| 4 | « P1 moins bon sur **les cinq colonnes** » | faux hors branche magnétisante | B1 |
| 5 | $I_0$ **+23 %**, $X_m$ **−11,5 %** (estimés) | **+21,3 %**, **−10,3 %** | B3 |
| 6 | les deux σ **mélangent** deux normalisations | **deux essais différents** | C3 |
| 7 | réseau à une dent **−29,6 %** sur $L_a$ | **+1,1 %** (artefact `ws1`) | — |
| 8 | T14 : **14 cellules** hors intervalle | **0** (arrondi de l'écart négligé) | T14 |
| 9 | ANSYS **43,03 h** | **24,60 h** | T9 |
| 10 | « **9 421 s** pour un point » | itération de conception, non le run publié | T9 |
| 11 | l'inertie explique l'écart aux 101,4 | atténuation **0,96** seulement | C5 |
| 12 | immunité du PMSM expliquée par le **facteur 47** sur $X$ et l'absence de condensation (§1.3) | **verrouillage troncature–pavage** : $N = M_s/2$ rend $\ln N + \ln\lvert 2\sin(d/2)\rvert \to \ln\pi$ stationnaire. L'explication réfutée prédisait l'**inverse** de l'observation — la queue s'installe dès $n\approx 7$ sur le PMSM contre $n\approx 337$ sur la MAS, donc **plus tôt** | Art. I §3.5 ; relecture croisée 7 août |
| 13 | P1 **dégrade** l'ondulation de +2,6 % à +20,9 % (§1.9) | **aucun effet** : −0,38 % à chaîne constante. Les deux chiffres venaient de **deux chaînes différentes**. La carte relève du traitement d'entrefer de `mesh_refined`, non de `ctx.AG` | A-1 |

## 5.2 Valeurs publiées que le programme ne régénère pas

| cellule | état |
|---|---|
| $B_r$ crête, FMM aimant, FMM entrefer (*Lumped*) | **RÉSOLUES** — φ = 39,429°, écarts 0,0000 % |
| $k_r$ *Lumped* 1,0522 | **incohérent** — évalué à φ = 0 alors que les FMM le sont à 39,429° |
| Pertes aimants 0,342 W | **NON RÉSOLU** — six chaînes testées. Adopter 0,204 W |
| Ondulation 101,4 N·m | **PÉRIMÉ** — la chaîne donne **108,2 N·m** (+2,6 % au lieu de −3,9 %). **Corrigé dans les quatre `.tex`** |

---

# 6. Inventaire

**Scripts produits** : `RUN_V1_PMSM_BASE`, `RUN_A1_TABLE7`, `RUN_A2_TABLE5`,
`RUN_A3_PMLOSS`, `RUN_A4_IRONLOSS`, `RUN_A5_LUMPED`, `RUN_B1_IM_P1`,
`RUN_B2_KC`, `RUN_B3_PRIX`, `RUN_B4_BARSKIN`, `RUN_B6_ANNEAU`,
`RUN_B10_B1_SKEWOFF`, `RUN_C4_ECARTS`, `RUN_C4B_TABLES_16_17`,
`RUN_C4C_TABLE8`, `RUN_C5B_RIPPLE_CHAIN`, `RUN_T1_T2_KC_TRUNCATION`,
`RUN_T16_BASE_P1`, `RUN_T17_QUEUE`, `RUN_T9_TEMPS`, `RUN_T12_CALAGE`,
`RUN_T14_AUDIT_ECARTS`, `RUN_T18_*`, `RUN_VERIFY_OPERATOR`.

**Modifications de fonctions** : option `basis` ('p0'/'p1') dans
`+mec/airgap_fourier`, `+mec/airgap_dtn_tooth`, `airgap_magnet`,
`cogging_mec` ; argument `ext` (source externe de potentiels) dans `pm_loss` ;
`R.Bg1_field` dans `+mec/magnetizing` ; en-tête et diagnostic énergétique dans
`RUN_IND_MESH` ; sauvegarde `T10_full.mat` dans `BLDC_MEC_COMPLET`.

**`.mat` à pleine précision** : `A1_table7`, `B1_im_p1`, `B2_kC`, `B3_prix`,
`C4b_tables_16_17`, `C4c_table8`, `T10_full`, `C5b_ripple`, `T1_T2_kc`,
`T16_p1`, `V1_pmsm`, `A4_ironloss`.

**Documents** : `RAPPORT_SPEC_CLAUDE_CODE.md` (v1),
`RAPPORT_SPEC_v2.md` (v2, détail par tâche),
`SPEC_CLAUDE_CODE_v2_MAJ.md` (bilan de la spécification),
`MANIFEST.md` (chaînes), le présent rapport.
