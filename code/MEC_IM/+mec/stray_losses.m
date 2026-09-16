function Padd = stray_losses(M, I2)
%STRAY_LOSSES  Pertes supplementaires en charge (additional / stray load losses).
%
%   Padd = mec.stray_losses(M, I2) renvoie les PERTES SUPPLEMENTAIRES en
%   charge de la machine, pour un courant rotorique fondamental rapporte au
%   stator I2 [A rms]. Ce sont les pertes de charge que le reseau de
%   reluctances FONDAMENTAL ne capture PAS :
%     * pertes par pulsation de denture (harmoniques de permeance) dans le
%       fer stator et rotor ;
%     * pertes de surface (rotor/stator) dues aux harmoniques d'espace ;
%     * pertes additionnelles de cage (courants inter-barres, harmoniques de
%       rang eleve a la frequence de denture) ;
%     * pertes d'extremite 3D (flux de fuite des tetes de bobines et des
%       anneaux dans les pieces massives, hors d'un modele 2D).
%
%   JUSTIFICATION ET CALAGE. Sans ce terme le rendement MEC atteint 94,1 % au
%   nominal (irrealiste). Le modele de champ pas-a-pas (RUN_STEPPING) ne
%   reproduit PAR LA PHYSIQUE de denture que ~17,5 W (4,6 %) de ces pertes ;
%   les mecanismes dominants (surface, inter-barres, extremite 3D) sont hors
%   de portee d'un modele 2D fondamental. Elles sont donc introduites comme
%   une CONSTANTE DE PERTES SUPPLEMENTAIRES = ALLOCATION IEC 60034-2-1,
%   Pn_ref = 381,5 W (~2 % Pin). Cette valeur est COMMUNE aux trois methodes
%   de l'etude comparative (analytique PLL = 381,5 W ; residu non identifie de
%   l'EF ANSYS 2D = 381,5 W) : les rendements se comparent alors a hypothese
%   de pertes supplementaires IDENTIQUE.
%   NB : a 381,5 W le rendement MEC vaut ~92,1 %, AU-DESSUS de la plage
%   experimentale 90,5-91,3 %, parce que les pertes IDENTIFIEES du MEC sont
%   ~215 W sous celles de l'EF (surtout la cage). Pour placer le rendement au
%   MILIEU de la plage (90,9 %) il faudrait Pn_ref ~ 606 W.
%
%   LOI DE CHARGE (norme IEEE 112 / IEC 60034-2-1). Les pertes
%   supplementaires EN CHARGE sont proportionnelles au CARRE du courant de
%   charge (elles s'annulent a vide, preservant la validation a vide) :
%       Padd = Pn_ref * (I2 / I2_ref)^2
%   Pn_ref est cale au point nominal (courant I2_ref). La loi quadratique
%   n'etant caracterisee que jusqu'a ~150 % de charge (plage des essais
%   normalises), I2 est plafonne a load_cap*I2_ref au-dela (extrapolation non
%   caracterisee : la fuite de fuite sature -> plateau, pas de divergence en
%   I2^2 au fort glissement).
%
%   Voir aussi : mec.equivalent_circuit, mec.power_balance, mec.mech_losses,
%   RUN_STEPPING, RUN_VALIDATION.

st = M.stray;
I2eff = min(abs(I2), st.load_cap * st.I2_ref);   % plafond 150 % de charge
Padd  = st.Pn_ref * (I2eff / st.I2_ref).^2;
end
