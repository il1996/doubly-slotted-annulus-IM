# R8 — note de livraison

> Sortie : `MEC_IM/R8_table2_out.txt`
> Script : `MEC_IM/RUN_R8_TABLE2.m`
> Exécuté le 9 août 2026, 93 s. MAS 48/44, 18,5 kW, 690 V. Pavage $n_T = 17$,
> $n_O = 4$ ; $N_h = 8192$ ; $s = 0{,}0188$. **Vrillage harmonique neutralisé**
> (`M.opt.skew_harm = 0`), et aussi actif, pour vérifier la note plutôt que
> l'admettre.

## La garde passe : la note n'a plus d'objet

| | cité ailleurs dans le `.tex` | Table 2 régénérée | écart |
|---|---|---|---|
| $I_0$ | +21,6 % | **+21,6 %** | 0,014 pt |
| $X_m$ saturée | −10,6 % | **−10,6 %** | 0,010 pt |

Les chiffres du tableau coïncident avec ceux de la conclusion **sans note
explicative**. La note des lignes 808–818 de `ArticleII_Carter_IM.tex` —
celle qui annonçait *« +21,3 % et −10,3 % here »* contre les valeurs du
réglage loyal — **peut être supprimée**.

Les deux valeurs de référence ont été **relues du `.tex` par expression
régulière**, non transcrites.

## La Table 2 régénérée

| fermeture | $X_{m0}$ (Ω) | $X_m$ charge (Ω) | $I_0$ (A) | $E_1$ (V) | couple (N·m) |
|---|---|---|---|---|---|
| Carter | 71,990 | 50,446 | 8,612 | 378,68 | 116,075 |
| piecewise const. | 64,781 | 43,075 | 9,921 | 375,92 | 114,751 |
| **hat (converged)** | **61,015** | **41,129** | **10,323** | **375,07** | **114,314** |
| référence EF | — | 46,0 | 8,4886 | 382,12 | 121,63 |

| écarts | $X_m$ charge | $I_0$ | $E_1$ | couple |
|---|---|---|---|---|
| Carter | +9,7 % | **+1,5 %** | −0,9 % | −4,6 % |
| piecewise const. | −6,4 % | +16,9 % | −1,6 % | −5,7 % |
| hat (converged) | **−10,6 %** | **+21,6 %** | −1,8 % | −6,0 % |

Le panneau LaTeX prêt à coller est en fin de sortie.

## Ce que la régénération change ailleurs, et il faut le savoir

**Le contraste de tête passe de « +1,1 % → +21,6 % » à « +1,5 % → +21,6 % ».**

C'est la conséquence la moins visible et la plus importante. Le +1,1 % de
Carter était une valeur **vrillage actif** ; en configuration loyale elle vaut
+1,5 %. Or ce chiffre porte l'argument central et apparaît au moins quatre
fois :

| ligne | phrase |
|---|---|
| 217–218 | *« moves the computed no-load current … from +1.1 % to +21.6 % »* |
| 800 | ligne « Carter » de la Table 2 |
| 1076 | *« against 1.1 % for the smooth-gap closure it replaces »* |
| 1523 | *« matched to 1.1 % by construction »* |

**La ligne 217–218 est précisément le défaut que le rapporteur signale** :
elle compare un +1,1 % vrillage actif à un +21,6 % vrillage neutralisé. La
régénération supprime ce mélange — les deux deviennent loyaux.

**L'argument ne bouge pas en substance** : le facteur passe d'environ 20 à
environ 14, et la thèse — la fermeture classique tient le courant à vide *par
construction*, l'opérateur non — est intacte.

## La note disait-elle vrai ? Presque

Elle affirmait que l'écart entre les deux réglages est *« le cage harmonique,
non la fermeture »*. Testé : le décalage doit alors être le même pour les trois
fermetures.

| fermeture | décalage sur $I_0$ | décalage sur $X_m$ |
|---|---|---|
| Carter | +0,339 pt | −0,439 pt |
| piecewise const. | +0,326 pt | −0,301 pt |
| hat (converged) | +0,319 pt | −0,262 pt |

**Sur $I_0$ l'affirmation tient bien** : les trois décalages s'étalent sur
0,02 point. **Sur $X_m$ elle tient moins** : de −0,262 à −0,439 pt, soit un
étalement de 0,18 point — le réglage de vrillage n'y est pas rigoureusement
séparable de la fermeture. L'affirmation était donc défendable pour la
grandeur qu'elle servait à défendre, et approximative pour l'autre. Elle
devient sans objet de toute façon.

## Deux points de provenance à déclarer

**Le couple de référence.** Le tableau publié utilisait 121,53 N·m — le point
nominal obtenu en résolvant le schéma équivalent — avec une note [a]
expliquant qu'il coexiste avec 121,63 N·m, moyenne temporelle du régime
établi. La régénération emploie **121,63**, qui est la valeur du reste du
manuscrit (Tables 3 et 4). Les écarts de couple bougent donc de 0,1 point.
*La note [a] documente une provenance réelle et n'absorbe rien : sa
suppression est un choix éditorial, pas une correction.*

**Le courant à vide de référence.** Relu à **8,4886 A** sur
`transitoire\a vide\Winding Plot 4.tab`, fenêtre $t > 1{,}0$ s — c'est bien le
8,49 employé par le manuscrit. En revanche l'en-tête de `+mec/ansys_ref.m`
annonce 8,4986 A pour le même fichier, soit **+0,12 %**. La différence tient
vraisemblablement à la fenêtre de moyennage, qui n'y est pas déclarée. Sans
effet ici — la garde passe à 0,014 point — mais l'en-tête et la mesure ne
disent pas la même chose.

## Ce que je n'ai pas fait

**Aucun `.tex` n'est modifié.** La note à supprimer est localisée
(lignes 808–818), le panneau de remplacement est produit, et les quatre
occurrences du +1,1 % sont listées avec leur ligne. L'application vous revient.

**Je n'ai pas relu $X_m = 46$ Ω de la référence** : `mec.ansys_ref` la déclare
*« non re-mesurable depuis les .tab disponibles, statut NON VÉRIFIÉ »*. Elle
est employée telle quelle et son statut est déclaré plutôt que masqué. Les
écarts sur $X_m$ — dont le −10,6 % qui porte l'argument — héritent de cette
réserve.
