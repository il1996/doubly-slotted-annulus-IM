%% RUN_C5_RIPPLE_SKEWOFF - calibration C5 refaite, vrillage MEC neutralise
%
%  MOTIF. B5 etablit que la §6.7 est fausse COTE MEC : la reference est
%  NON vrillee (NumberOfSlices = 1 ne produit aucun moyennage axial), mais
%  mec.cage applique ksq au fondamental ET ksq_nu aux harmoniques sans que
%  M.opt.skew_harm ne soit jamais mis a 0. L'ondulation, portee par les
%  harmoniques de denture, comparait donc un MEC vrille a une reference
%  droite -- les branches de denture y etant attenuees d'un facteur
%  plusieurs centaines.
%
%  VOIE (a) : neutraliser le vrillage MEC et republier l'ondulation. C'est
%  la seule qui preserve la valeur de la comparaison.
%
%  METHODE. Couple par TENSEUR DE MAXWELL harmonique de l'operateur
%  (AT.torque), a excitation complete -- courant statorique TOTAL et onde de
%  FMM de barres, comme mec.torque_vw. La carte est construite sur DEUX PAS
%  D'ENCOCHE statoriques, ou le motif de denture est periodique.
clear; clc; t0=tic;
diary('C5_ripple_out.txt'); diary on;
nT=17; nO=4; Nh=8192; s_ch=0.0188; Npos=19;
TppF=105.5;                         % reference EF, ondulation c-c [N.m]

fprintf('=== C5 : ondulation, vrillage MEC neutralise ===\n');
fprintf('  base P1 | nT=%d nO=%d N_h=%d | s = %.4f\n',nT,nO,Nh,s_ch);
fprintf('  %d positions sur DEUX pas d''encoche statoriques\n',Npos);
fprintf('  reference EF : %.1f N.m c-c (section droite, NON vrillee)\n\n',TppF);

R=struct(); nm={'vrillage ON (etat publie)','vrillage OFF (voie a)'};
for k=1:2
    M=mec.machine_18_5kW();
    if k==2, M.opt.skew_harm=0; end          % <-- neutralisation
    ctx=mec.build_context(M); G=ctx.G;
    A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
    ctx.AG=A; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
    r=mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);

    %  excitation statorique complete
    p1=angle(r.I1c);
    i3=sqrt(2)*abs(r.I1c)*[cos(p1);cos(p1-2*pi/3);cos(p1+2*pi/3)];
    Fs=ctx.W.slotMMF(i3);

    %  ---- FMM ROTORIQUE : fondamental + SPECTRE HARMONIQUE ----
    %  res.I2n ne donne que des MODULES : sommer des ondes sans leurs
    %  phases n'a pas de sens. On reconstruit donc les impedances
    %  harmoniques par mec.cage -- ce qui rend AUSSI la chaine sensible a
    %  M.opt.skew_harm, dont c'est precisement le point d'action (ksq_nu).
    tb=2*pi*(0:M.Nr-1)/M.Nr;
    H=ctx.H; nH=numel(H.nu);
    ib=sqrt(2)*abs(r.I2c)*(2*M.m*ctx.W.kw1*ctx.W.Nph/M.Nr) ...
       *cos(M.p*tb-angle(r.I2c));                    % fondamental
    nharm=0;
    for kk=1:nH
        if H.dir(kk)>0, sn=1-H.nu(kk)*(1-s_ch); else, sn=1+H.nu(kk)*(1-s_ch); end
        if abs(sn)<1e-6, continue; end
        Ck=mec.cage(M,G,ctx.W,ctx.Lk,max(abs(sn),1e-3),H.nu(kk),H.kw(kk));
        %  Xmn EXACTEMENT comme mec.equivalent_circuit:66 -- reactance
        %  SATUREE au point de fonctionnement multipliee par le coefficient
        %  de fuite differentielle H.sig, et NON Xm0/nu^2. Prendre Xm0
        %  (61,0 ohm au lieu de 41,2) gonfle chaque courant harmonique et
        %  retourne le couple moyen a -38,9 N.m. Le controle "couple moyen"
        %  est la garde qui l'attrape.
        Xmn=r.Xm*H.sig(kk);
        if Xmn<=0, continue; end
        Zh=1/( 1/(1i*Xmn) + 1/(Ck.Rr/sn + 1i*Ck.Xr) );
        I2c_n=r.I1c*Zh/(Ck.Rr/sn + 1i*Ck.Xr);        % COMPLEXE : phase gardee
        Ibn=abs(I2c_n)*2*M.m*abs(H.kw(kk))*ctx.W.Nph/M.Nr;
        if ~isfinite(Ibn)||Ibn<=0, continue; end
        ib=ib+sqrt(2)*Ibn*cos(H.nu(kk)*M.p*tb-angle(I2c_n));
        nharm=nharm+1;
    end
    %  FMM = somme cumulee des courants de barre, moyenne retiree.
    %  SIGNE NEGATIF : mec.torque_vw:46 porte "Fr = -Fu*(...)" avec le
    %  commentaire "signe : rotor s'oppose (charge)". L'omettre laisse le
    %  couple moyen a -39 N.m au lieu de +114,7 -- c'est la garde "couple
    %  moyen" qui l'attrape, et non l'inspection du code.
    Fr=-(cumsum(ib(:))-mean(cumsum(ib(:))));
    if k==1
        fprintf('  [%d harmoniques de cage retenus sur %d]\n',nharm,nH);
    end

    %  carte d'ondulation : l'operateur est REBATI a chaque position, la
    %  couronne bi-encochee n'etant PAS invariante en position (C2 restreint)
    taus=2*pi/M.Ns; phis=linspace(0,2*taus,Npos);
    T=nan(1,Npos);
    for q=1:Npos
        Aq=mec.airgap_dtn_tooth(M,G,phis(q),nT,nO,Nh,'p1');
        S=mec.solve_network(ctx.net,G,ctx.BH,Aq,Fs,Fr,M.opt);
        T(q)=Aq.torque(S.Usurf);
    end
    R(k).T=T; R(k).Tm=mean(T); R(k).Tpp=max(T)-min(T);
    R(k).Xm=ctx.Xm0; R(k).Tem=r.Tem; R(k).I2=abs(r.I2c);
    %  GARDE. Le couple moyen de la carte doit retomber sur celui du schema
    %  equivalent. S'il s'en ecarte, l'excitation harmonique est fausse et
    %  les valeurs absolues d'ondulation ne veulent rien dire.
    R(k).ecT=100*(R(k).Tm-r.Tem)/r.Tem;
    fprintf('  %-26s Tmoy %9.3f (schema %7.3f, %+6.2f %%) | Tpp %9.3f N.m\n', ...
        nm{k},R(k).Tm,r.Tem,R(k).ecT,R(k).Tpp);
end

%% ---- comparaison -----------------------------------------------------
fprintf('\n  ---- ondulation crete-a-crete ----\n');
fprintf('  %-26s %12s %12s %11s\n','etat','Tpp (N.m)','EF (N.m)','ecart');
for k=1:2
    fprintf('  %-26s %12.3f %12.1f %10.2f %%\n', ...
        nm{k},R(k).Tpp,TppF,100*(R(k).Tpp-TppF)/TppF);
end
fprintf('\n  effet de la neutralisation : Tpp %+.2f %% | Tmoy %+.2f %%\n', ...
    100*(R(2).Tpp-R(1).Tpp)/R(1).Tpp,100*(R(2).Tm-R(1).Tm)/R(1).Tm);
fprintf('  courant rotorique I2 : %.4f -> %.4f A (%+.2f %%)\n', ...
    R(1).I2,R(2).I2,100*(R(2).I2-R(1).I2)/R(1).I2);

fprintf(['\n  LECTURE. La reference est une SECTION DROITE. Seule la ligne\n' ...
         '  "vrillage OFF" la compare a un modele dans le MEME etat. La\n' ...
         '  ligne "ON" est conservee pour montrer ce que valait l''accord\n' ...
         '  publie -- elle ne doit pas etre citee comme validation.\n']);
save('C5_ripple.mat','R','TppF','Npos');
fprintf('\n  duree %.0f s\n=== C5 termine ===\n',toc(t0));
diary off;
