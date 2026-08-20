# R1 — note de livraison

> Sortie : `MEC_IM/R1_ripple_reconcile_out.txt`
> Script : `MEC_IM/RUN_R1_RIPPLE_RECONCILE.m`
> Exécuté le 9 août 2026.

## Ce qui est publiable

**La valeur publiable est 108,224433 N·m — la chaîne (A), celle de
`RUN_ARTICLE.m:334-342`.** Elle donne +2,5955 % contre la référence EF
105,486548 N·m, recalculée de `Plot 1.tab` colonne 3 sur `t ≥ 1,90 s`.

**Le mécanisme de l'écart est nommé et certain.** Les deux chaînes ne
diffèrent que par une ligne : le déphasage `ph = p·θ` appliqué aux courants
statoriques.

```
(A) i3k = sqrt(2)*I1*[cos(psi1+ph); cos(psi1+ph-2pi/3); cos(psi1+ph+2pi/3)]
(B) i3k = sqrt(2)*I1*[cos(psi1);    cos(psi1-2pi/3);    cos(psi1+2pi/3)]
```

Sur les **30,0000 degrés électriques** que balaient les deux pas d'encoche,
figer les courants fait varier l'**angle de charge** d'autant. La variation
du couple **fondamental** s'ajoute alors à l'ondulation de denture :
**+18,5082 %** sur la carte brute, **+18,2676 %** après injection mécanique.
La chaîne (B) ne mesure donc pas l'ondulation de denture seule, et ses
127,994395 N·m ne sont pas publiables.

**Deux causes formellement exclues, comme exigé.**

*Piège n° 12* — l'atténuation par l'inertie vaut **0,965737** pour (A) et
**0,963776** pour (B). Identique aux deux ; elle ne rend compte d'aucune part
des 18 %.

*Amont du calcul de couple* — à position rotor identique (θ = 0, où les deux
variantes ont `ph = 0`), le champ tangentiel le long du bore est identique au
**bit près** : `B_t` rms = 0,10302311 T des deux côtés, écart maximal point à
point **0,000e+00 T**. Le maillage, le traitement d'entrefer interne
(`me.gapF`) et le solveur sont donc hors de cause. **L'écart est entièrement
dans l'excitation.**

## Ce qui n'est pas publiable — la garde a sauté

**La garde exigée par la v4 §R1 n'est pas satisfaite, et son échec est un
résultat.** Les deux chaînes devaient redonner le même couple moyen. Elles ne
le font pas :

| chaîne | couple moyen | vs schéma équivalent (122,897523 N·m) |
|---|---|---|
| (A) courants tournants | 104,025647 | **−15,3558 %** |
| (B) courants figés | 120,496227 | −1,9539 % |

La v4 pose que dans ce cas *« ce n'est pas un problème d'ondulation mais de
modèle, et le §7 n'est pas le seul concerné »*. C'est le cas.

**Le couple moyen de la carte publiée est 15,4 % en dessous du couple du
schéma équivalent au même glissement.** Ce défaut est indépendant de
l'ondulation — celle-ci est correcte — et il n'était pas connu. Il touche
toute grandeur tirée de la carte magnétostatique en charge, pas seulement §7.

**Signalé pour correction, non corrigé.** Conformément à la règle 4, la chaîne
n'est pas ajustée vers l'affirmation. Deux pistes, à trancher par l'auteur :

1. **Recalage angulaire.** Si l'origine angulaire du bobinage dans
   `mesh_refined` ne coïncide pas avec celle que suppose `ph = p·θ`, la
   chaîne (A) tourne les courants à partir d'un angle de charge décalé. Le
   couple moyen serait alors systématiquement biaisé sans que l'ondulation le
   soit, ce qui correspond exactement à l'observation.
2. **Différence de route.** Le couple de la carte vient du tenseur de Maxwell
   harmonique ; celui du schéma de `m·I₂′²·R_r′/(s·ω_s)`. Un écart est
   attendu, mais 15,4 % contre 2,0 % pour la variante fautive demande une
   explication que le dossier ne fournit pas.

**Le §7 peut donc citer 108,2 N·m avec son mécanisme d'écart nommé**, ce qui
lève le motif de rejet immédiat. Mais l'affirmation implicite selon laquelle
la carte reproduit le point de fonctionnement doit être vérifiée avant que
d'autres grandeurs n'en soient tirées.

## Critère d'acceptation v4

*« Une seule valeur survit, avec le mécanisme de l'écart nommé. »* — satisfait :
**108,224433 N·m**, mécanisme = rotation des courants statoriques omise dans
(B), inflation mesurée +18,27 %.

*Garde exécutée* — oui, et **elle échoue** ; son échec est rapporté ci-dessus
plutôt que masqué.
