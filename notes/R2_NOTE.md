# R2 — note de livraison

> Sortie : `MEC_IM/R2_kc_tiling_out.txt`
> Script : `MEC_IM/RUN_R2_KC_TILING.m`
> Exécuté le 9 août 2026. Base chapeau, $N_h = 8192$ fixe sur les neuf points.

## Ce qui n'est pas publiable en l'état

**Le +5,2 % au-dessus de Carter n'est pas robuste au pavage.**

| $n_T$ | $n_O$ | colonnes | $X_m$ encoché | $k_C$ | écart / production |
|---|---|---|---|---|---|
| 9 | 2 | 1012 | 58,512985 | **1,388919** | +4,2755 % |
| 9 | 4 | 1196 | 59,866654 | 1,357513 | +1,9177 % |
| 9 | 8 | 1564 | 60,835972 | 1,335884 | +0,2938 % |
| 17 | 2 | 1748 | 59,795480 | 1,359129 | +2,0390 % |
| **17** | **4** | **1932** | **61,014695** | **1,331971** | **0** |
| 17 | 8 | 2300 | 61,937595 | 1,312124 | −1,4900 % |
| 33 | 2 | 3220 | 60,901734 | 1,334441 | +0,1855 % |
| 33 | 4 | 3404 | 61,953174 | 1,311794 | −1,5148 % |
| 33 | 8 | 3772 | 62,620182 | **1,297821** | −2,5638 % |

$X_m$ lisse, définition publiée : **81,269790 Ω**, invariant du pavage.

**Dispersion de $k_C$ : 6,8155 %**, de 1,297821 à 1,388919.

**La dispersion excède la correction annoncée.** L'écart de $k_C$ à Carter
(1,266458) vaut **+9,67 %** au pavage le plus grossier et **+2,48 %** au plus
fin — contre +5,2 % annoncé au pavage de production. **La correction
centenaire est du même ordre que l'incertitude de discrétisation.**

## Le point aggravant : le pavage n'est pas convergé

$X_m$ encoché **croît de façon monotone dans les deux directions de
raffinement**, sans se stabiliser :

- à $n_T = 33$ : 60,90 → 61,95 → 62,62 quand $n_O$ passe de 2 à 8 ;
- à $n_O = 8$ : 60,84 → 61,94 → 62,62 quand $n_T$ passe de 9 à 33.

Le point le plus fin des neuf n'est donc **pas** une valeur convergée mais un
minorant, et $k_C$ y décroît encore. **Rien dans ces données n'exclut que
$k_C$ rejoigne, voire franchisse, la valeur de Carter sous raffinement
supplémentaire.**

C'est le même défaut structurel que T17 a démontré en forme fermée sur l'axe
des troncatures, transposé à l'axe du pavage. Il n'est pas établi qu'une
limite existe.

## Ce que dit la GARDE

| grandeur | dispersion |
|---|---|
| $X_m$ lisse (a), définition publiée | **0,0000 %** — constant par construction |
| $X_m$ encoché | 6,7523 % |
| $k_C$ (a) = lisse(a)/encoché | **6,8155 %** |
| $X_m$ lisse (b), pavage apparié | 13,6451 % |
| $k_C$ (b) = lisse(b)/encoché | 9,4858 % |

**Le rapport n'est PAS plus stable que ses composantes** : 6,8155 % contre
6,7523 %. C'est attendu et cela répond à l'objection anticipée par la v4.

Avec la définition publiée, le **numérateur est constant sur l'axe du pavage**
— sa dispersion est nulle *par construction*, non par compensation. La
dispersion de $k_C$ est donc **exactement** celle de $X_m$ encoché, à 0,06
point près. **Il n'existe aucun degré de liberté pour une compensation
numérateur/dénominateur sur cet axe.**

C'est ce qui distingue cet axe de celui des **troncatures**, où le §4.3 a
réfuté un argument de compensation : là, les *deux* termes variaient — le
numérateur de 0,15 %, le dénominateur de 5,99 %. Ici un seul varie.
**L'analogie qu'un rapporteur relèverait ne tient pas**, et le tableau
l'établit par les chiffres.

La colonne de contrôle (b) confirme : $k_C$ (b) disperse de 9,49 %, entre les
6,75 % et 13,65 % de ses composantes. Aucune compensation là non plus.

## Ce qui est publiable

**Le tableau des neuf points, avec la barre de dispersion.** La v4 le prévoit :
*« le résultat central devra alors être annoncé avec sa barre de dispersion,
ce qui reste publiable »*.

Formulation défendable : *l'opérateur exact, en base chapeau, donne un rapport
d'encochage compris entre 1,298 et 1,389 selon le pavage, soit +2,5 % à +9,7 %
au-dessus de la valeur de Carter ; la valeur au pavage de production est
1,332 (+5,2 %)*.

**Ce qui n'est plus défendable** : citer 1,3320 comme *la* valeur, ni annoncer
+5,2 % sans sa dispersion. Signalé pour correction — le résumé de l'Article II
et son §4 sont concernés. **La chaîne n'a pas été ajustée vers l'affirmation.**

## Ce qui reste à faire avant soumission

Le pavage n'étant pas convergé, deux voies :

1. **Étendre la grille** vers $n_T = 65$, $n_O = 16$ et montrer une limite. Le
   coût croît en $O(N_h \cdot M^2)$ — le point (33, 8) demande déjà 3 772
   colonnes et 10 s ; (65, 16) en demanderait environ 7 500 et
   quatre fois plus de temps. Faisable.
2. **Déclarer la non-convergence** comme T17 l'a fait pour la troncature, et
   publier l'intervalle. Cohérent avec la posture du dossier.

La voie 1 est préférable : si une limite existe, elle change la nature du
résultat. Si elle n'existe pas, la voie 2 devient la seule honnête — et il
faudra dire que **le rapport d'encochage n'est pas une propriété de la
géométrie mais de la discrétisation**, ce qui affaiblit considérablement
l'Article II.

## Critère d'acceptation v4

*« Si $k_C$ tient à quelques dixièmes de pour cent sur les neuf points, le
résultat devient une propriété de la géométrie »* — **non satisfait** :
6,82 %.

*« S'il dérive, il vaut infiniment mieux le découvrir soi-même »* — c'est le
cas, et le tableau est produit.
