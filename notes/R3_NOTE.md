# R3 — note de livraison

> Sortie : `MEC_IM/R3_fe_ratio_out.txt`
> Script : `MEC_IM/RUN_R3_FE_RATIO.m`
> Exécuté le 9 août 2026.

## État : BLOQUÉ par une dépendance externe

**R3 exige deux résolutions magnétostatiques ANSYS que je ne peux pas
lancer**, et dont l'une n'existe pas.

Vérification faite sur `IM_18kW_690V.aedt` : **un seul design, un seul
`Setup1`, aucune géométrie lisse**. Les quatre occurrences de « smooth » dans
le fichier sont incidentes. La géométrie lisse équivalente doit donc être
**créée**, ce qui suppose l'interface Maxwell ou un script VBS, puis résolue.

**Ce n'est pas un contournement possible.** Le rapport que R3 demande est
précisément celui que la méthode de référence produit sur *sa* géométrie ; le
dériver du modèle reviendrait à valider le modèle par lui-même.

## Ce qui est livré et exploitable

**La garde analytique, calculée.** Le run lisse devra redonner :

```
Lambda_lisse = mu0*(2*pi*R*L)/g = 4,3779802798e-04 H
  avec R = 0,08176725 m (mi-entrefer), L = 0,16478245 m, g = 0,2430 mm
```

**La convention de rayon est sans enjeu** : au rayon d'alésage on obtient
4,3844856301e-04 H, soit **0,1486 % d'écart** — inférieur à la tolérance de
« quelques dixièmes de pour cent » que la v4 fixe. Ce point méritait d'être
tranché avant les runs plutôt qu'après, car il aurait pu être invoqué pour
absorber un désaccord.

**La valeur MEC à confronter, relue de R2** — et R2 change ce que R3 doit
tester :

| | valeur |
|---|---|
| $k_C$ au pavage de production (17, 4) | 1,331971 |
| **intervalle sur les neuf pavages** | **1,297821 à 1,388919** |
| dispersion | 6,8155 % |
| Carter classique | 1,266458 |

**La confrontation EF doit porter sur l'intervalle, non sur la valeur de
production.** Un chiffre EF tombant à 1,31 confirmerait le modèle au pavage
fin et l'infirmerait au pavage grossier — ce qui ne se lit pas si l'on compare
à 1,332 seul.

**La spécification des deux runs**, écrite dans la sortie : ce qui doit être
tenu identique — entrefer mécanique 0,2430 mm, rayons, bobine, FMM, courant,
perméabilité *déclarée*, maillage, critère de convergence, vrillage, et le
**même projet** — et la seule différence admise : encoches supprimées **des
deux côtés**.

**Le traitement est en place** : dès que `MEC_IM/R3_fe/slotted_flux.tab` et
`smooth_flux.tab` existeront, le script forme le rapport et vérifie la garde
sans autre intervention.

## Ce que je n'ai pas fait, et pourquoi

**Je n'ai pas produit de chiffre EF de substitution.** La v4 interdit de
réconcilier un écart en changeant la référence, et fournir une estimation
analytique du rapport EF reviendrait à cela.

**Je n'ai pas créé le design lisse.** Modifier un projet ANSYS de référence
dont dépend toute la validation n'est pas une opération que je ferai sans
demande explicite — une erreur de géométrie y invaliderait silencieusement les
comparaisons antérieures.

## Conséquence pour la soumission

R3 est **bloquant** au sens de la v4, et il le reste. L'objection centrale de
l'Article II — *« le modèle se trompe de 21,6 % sur la grandeur même que sa
correction est censée corriger »* — n'a toujours pas de réponse.

**R2 l'aggrave** : la correction n'est plus un chiffre mais un intervalle de
6,8 %, non convergé. Un rapporteur peut désormais objecter que la correction
et son incertitude sont du même ordre.

**Deux runs ANSYS restent la seule voie.** Ils sont courts — magnétostatique,
une position, même maillage — et ils trancheraient une objection que rien
d'autre ne peut lever.
