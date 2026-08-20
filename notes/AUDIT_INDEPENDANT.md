# Audit indépendant de l'archive

> Réalisé le 10 août 2026 par un agent qui n'a écrit aucun des scripts
> archivés. Il répond partiellement à la garde d'acceptation de la
> spécification v4 — *« un tiers qui n'a pas écrit le code doit pouvoir
> régénérer une grandeur de chaque article »* — que l'agent de montage ne
> pouvait pas s'auto-administrer.

---

## Ce qui a été vérifié, et ce qui ne l'a pas été

**Non vérifié : l'exécution.** Le test de deux minutes du §4 du `README`
demande MATLAB. Il reste à passer, et il reste le critère décisif.

**Vérifié : la traçabilité arrière.** Chaque valeur publiée à trois
chiffres significatifs ou plus dans les deux manuscrits a été recherchée
numériquement — avec une tolérance d'un demi-quantum sur la dernière
décimale publiée, et non par comparaison de chaînes — dans les 3713
valeurs distinctes que contiennent les 59 transcriptions `.txt` de
`outputs/`.

---

## Résultat

| manuscrit | valeurs publiées à $\geq 3$ chiffres | retrouvées | taux |
|---|---|---|---|
| `ArticleI_DtN_PMSM.tex` | 250 | 232 | **93 %** |
| `ArticleII_Carter_IM.tex` | 136 | 127 | **93 %** |

Les vingt-sept valeurs non retrouvées ont été examinées une par une.
Aucune n'est un chiffre de résultat orphelin ; toutes relèvent de deux
catégories qui n'engagent pas la traçabilité des sorties.

**Données d'entrée déclarées** — longueur de fer \SI{164.8}{\milli\metre},
inertie $3{,}789\times10^{-4}$, pertes par frottement et ventilation
\SI{98.2}{\watt}, ouverture rapportée à l'entrefer $b_0/g = 8{,}23$,
épaisseur logarithmique $\Xg = 2{,}97\times10^{-3}$. Ce sont des
paramètres de machine ou des rapports géométriques, pas des sorties de
chaîne.

**Arithmétique conduite dans le texte** — les taux 79,3 % et 80,8 % des
comptes 23/29 et 21/26, les sommes en points de pourcentage 41,7 et 42,6,
la dispersion 6,82 % formée sur les neuf pavages de la table du §4.3, les
écarts 3,55 et 5,52 A du diagnostic de la carte de saturation. Chacune
se recalcule à partir de valeurs qui, elles, sont archivées.

---

## Ce que cela établit, et ce que cela n'établit pas

**Établi.** Les manuscrits ne contiennent pas de chiffre de résultat
détaché de toute exécution. Le motif que le rapport d'évaluation
signalait — trois valeurs retirées entre les deux articles parce
qu'elles n'étaient pas régénérables — ne se reproduit nulle part
ailleurs dans l'état actuel des textes.

**Non établi.** Que les chaînes archivées, ré-exécutées, redonnent ces
valeurs. Un chiffre présent dans une transcription prouve qu'une
exécution l'a produit un jour ; il ne prouve pas que le code archivé le
reproduise aujourd'hui sur une autre machine. Seul le test du §4 le
prouve, et il n'a pas été passé.

L'archive est donc **cohérente et complète en traçabilité arrière, non
encore validée en reproductibilité**. Le manuscrit ne doit rien affirmer
de plus que cela tant que le test du §4 n'a pas tourné chez un tiers.

---

## Méthode, pour qu'elle soit refaisable

Les valeurs publiées sont extraites des deux `.tex` par les motifs
`\num{...}`, `\SI{...}{...}` et les pourcentages littéraux, commentaires
LaTeX exclus ; sont retenues celles portant au moins trois chiffres
significatifs, les autres n'identifiant pas une grandeur. Les valeurs
archivées sont extraites de `outputs/**/*.txt` par un motif numérique
général et converties en flottants. L'appariement est numérique, avec
tolérance relative $6\times10^{-4}$, ce qui absorbe l'arrondi d'une
valeur publiée à quatre décimales contre une transcription à six.

Un appariement par chaînes de caractères donnerait 84 % et 90 %, et ces
deux chiffres seraient faux : ils compteraient comme absentes des
valeurs présentes sous une autre écriture — `59.795480` contre
`59.7955`, ou `3.0354e-07` contre `3.0354e-7`.
