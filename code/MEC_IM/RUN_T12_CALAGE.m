%% RUN_T12_CALAGE - ou est l'erreur de couple au rotor bloque ?
%
%  T12. Au calage : courant stator -0,3 %, courant barre +0,7 %, mais couple
%  -4,6 %. Deux amplitudes justes a moins de 1 % et un couple faux de 4,6 % :
%  l'erreur ne peut pas etre dans les amplitudes. Elle est dans la PHASE ou
%  dans le mutuel. Hypothese de la specification : cela pointe vers la meme
%  branche magnetisante que le +16,5 % a vide.
%
%  METHODE. Le couple au calage vaut T = m*I2'^2*Rr'/w_s. On decompose donc
%  l'ecart en ses facteurs mesurables, puis on ECHANGE la fermeture
%  d'entrefer (Carter <-> operateur) a tout le reste egal : si la branche
%  magnetisante est en cause, le couple doit bouger alors que les amplitudes
%  ne bougent pas.
clear; clc;
diary('T12_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M);
ctxC=ctx; ctxC.AG=ctx.AGcarter; ctxC.Xm0=mec.magnetizing(ctxC,0.2).Xm;

%  References EF au rotor bloque (Table 16)
TF=104.31; I1F=108.21; IbF=1929; IrF=6530;

fprintf('=== T12 : localisation de l''erreur de couple au calage (s = 1) ===\n\n');
fprintf('  %-22s %10s %10s %10s\n','grandeur','CARTER','DtN','EF');
rC=mec.equivalent_circuit(ctxC,1,ctxC.Xm0);
rD=mec.equivalent_circuit(ctx ,1,ctx.Xm0);
pr=@(l,a,b,f)fprintf('  %-22s %10.2f %10.2f %10.2f   DtN %+6.1f %%\n', ...
    l,a,b,f,100*(b-f)/f);
pr('couple [N.m]',rC.Tem,rD.Tem,TF);
pr('I1 [A]',rC.I1,rD.I1,I1F);
pr('I2'' rapporte [A]',rC.I2,rD.I2,NaN);
pr('Im magnetisant [A]',rC.Im,rD.Im,NaN);
fprintf('  %-22s %10.4f %10.4f %10s\n','cos(phi)',rC.cosphi,rD.cosphi,'--');

%% ---- phases : c'est la ou la specification pointe -------------------
angD=@(z)angle(z)*180/pi;
fprintf('\n  --- phases (Uph reference reelle) ---\n');
fprintf('  %-22s %10s %10s\n','angle [deg]','CARTER','DtN');
fprintf('  %-22s %10.3f %10.3f\n','arg(I1)',angD(rC.I1c),angD(rD.I1c));
fprintf('  %-22s %10.3f %10.3f\n','arg(I2)',angD(rC.I2c),angD(rD.I2c));
fprintf('  %-22s %10.3f %10.3f\n','arg(E1)',angD(rC.E1c),angD(rD.E1c));
dC=angD(rC.I2c)-angD(rC.I1c); dD=angD(rD.I2c)-angD(rD.I1c);
fprintf('  %-22s %10.3f %10.3f  <- dephasage barre/stator\n','arg(I2)-arg(I1)',dC,dD);

%% ---- decomposition de l'ecart de couple ------------------------------
%  T = m*I2'^2*Rr'/w_s : l'ecart se factorise en amplitude et resistance.
fprintf('\n  --- decomposition de l''ecart de couple (DtN vs EF) ---\n');
fA=(IbF*(1+0.007)/IbF)^2;                 % contribution amplitude barre (+0,7 %)
fprintf('  facteur amplitude  (I_barre +0.7 %%)^2      : %+6.2f %%\n',100*(fA-1));
fprintf('  facteur resistance (Rr'' 0.4234 vs 0.4400)  : %+6.2f %%\n', ...
    100*(0.4234/0.4400-1));
fprintf('  produit des deux                           : %+6.2f %%\n', ...
    100*(fA*0.4234/0.4400-1));
fprintf('  ecart OBSERVE sur le couple                : %+6.2f %%\n', ...
    100*(rD.Tem-TF)/TF);
fprintf('  RESIDU non explique par les amplitudes     : %+6.2f %%\n', ...
    100*((rD.Tem/TF)/(fA*0.4234/0.4400)-1));

%% ---- test decisif : la branche magnetisante est-elle en cause ? -----
fprintf('\n  --- test : echange de la fermeture, tout le reste egal ---\n');
fprintf('  Xm0 Carter %.2f ohm -> DtN %.2f ohm  (%+.1f %%)\n', ...
    ctxC.Xm0,ctx.Xm0,100*(ctx.Xm0-ctxC.Xm0)/ctxC.Xm0);
fprintf('  couple    %.2f -> %.2f N.m  (%+.2f %%)\n', ...
    rC.Tem,rD.Tem,100*(rD.Tem-rC.Tem)/rC.Tem);
fprintf('  I1        %.2f -> %.2f A    (%+.2f %%)\n', ...
    rC.I1,rD.I1,100*(rD.I1-rC.I1)/rC.I1);
fprintf('  I2        %.2f -> %.2f A    (%+.2f %%)\n', ...
    rC.I2,rD.I2,100*(rD.I2-rC.I2)/rC.I2);
%  Le couple varie comme I2^2 : on teste si l'echange le respecte.
fI2=(rD.I2/rC.I2)^2; fT=rD.Tem/rC.Tem;
fprintf('\n  (I2_DtN/I2_Carter)^2 = %.4f   T_DtN/T_Carter = %.4f\n',fI2,fT);
fprintf('  residu apres correction par I2^2 : %+.2f %%\n',100*(fT/fI2-1));
fprintf(['\n  LECTURE. Si le rapport des couples suit (I2)^2, l''echange de\n' ...
         '  fermeture n''agit QUE par l''amplitude : la branche magnetisante\n' ...
         '  n''est alors PAS en cause, contrairement a l''hypothese de T12.\n' ...
         '  Au calage Im vaut ~1 %% de I2 : la branche magnetisante ne peut\n' ...
         '  structurellement pas porter 4,6 %% de couple.\n']);
%  Reste alors la resistance rotorique, seule grandeur libre : T = m*I2''^2*Rr''.
fprintf('\n  --- ce qui reste : la resistance rotorique au calage ---\n');
fprintf('  Rr'' MEC 0.4234 vs EF 0.4400 ohm : %+.1f %%\n',100*(0.4234/0.4400-1));
fprintf('  couple attendu de (I_barre +0.7 %%)^2 x Rr'' : %+.2f %%\n', ...
    100*(1.007^2*0.4234/0.4400-1));
fprintf('  couple observe                             : %+.2f %%\n', ...
    100*(rD.Tem-TF)/TF);
fprintf(['  A s = 1 la frequence rotorique vaut 50 Hz : mec.bar_skin est\n' ...
         '  sollicite a son extreme. C''est la, et non dans le couplage, que\n' ...
         '  se situe l''erreur de couple au calage.\n']);
fprintf('\n=== T12 termine ===\n');
diary off;
