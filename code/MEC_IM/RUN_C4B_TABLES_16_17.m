%% RUN_C4B_TABLES_16_17 - Tables 16 et 17 a PLEINE PRECISION
%
%  C4 (reste). Ces deux tables n'avaient pas de .mat : leurs colonnes
%  d'ecart ne pouvaient pas etre regenerees. On relance donc leurs chaines
%  AVEC SAUVEGARDE des valeurs pleines.
%
%  Table 16 : trois conditions d'essai (en charge, a vide, rotor bloque),
%             chaine du RESEAU DE PERFORMANCE, base P1 convergee.
%  Table 17 : densites de flux LOCALES, chaine du MAILLAGE RAFFINE
%             (mec.mesh_refined) -- c'est une chaine DISTINCTE, et c'est
%             pourquoi B1 ne l'avait pas produite.
%
%  Les references EF sont RELUES des .tab du projet, jamais transcrites.
clear; clc; t0=tic;
diary('C4b_tables_16_17_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G; W=ctx.W; BH=ctx.BH;
nT=17; nO=4; Nh=8192; s_ch=0.0188;
A1=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
ctx.AG=A1; ctx.Xm0=mec.magnetizing(ctx,0.2).Xm;
%  M.FEA est un champ du PMSM ; la MAS porte sa racine ailleurs.
ROOT='C:\Users\hp\Desktop\ANSYS résultat 18.5KW';
rd=@(f)readmatrix(f,'FileType','text','NumHeaderLines',1,'Delimiter','\t');
rmsw=@(A,c,t0_)sqrt(mean(A(A(:,1)>=t0_,c).^2,'omitnan'));
avgw=@(A,c,t0_)mean(A(A(:,1)>=t0_,c),'omitnan');

fprintf('=== C4b : Tables 16 et 17, pleine precision ===\n');
fprintf('  base P1 | nT=%d nO=%d N_h=%d | X_m0 = %.6f ohm\n\n',nT,nO,Nh,ctx.Xm0);

%% ======================= TABLE 16 ======================================
r  = mec.equivalent_circuit(ctx,s_ch,ctx.Xm0);
r0 = mec.equivalent_circuit(ctx,1e-4,ctx.Xm0);
rb = mec.equivalent_circuit(ctx,1.0 ,ctx.Xm0);
kb = (2*M.m*W.kw1*W.Nph)/M.Nr;              % report barre
kr_ring = 1/(2*sin(M.p*pi/M.Nr));           % rapport anneau/barre ideal

%  --- references EF, relues ---
Dc=fullfile(ROOT,'transitoire','en charge');
D0=fullfile(ROOT,'transitoire','a vide');
Db=fullfile(ROOT,'transitoire','rotor bloque');
tc=1.0;
try, trq =rd(fullfile(Dc,'Torque Plot 1.tab'));            Tc_F=avgw(trq,3,tc);  catch, Tc_F=NaN; end
try, cur =rd(fullfile(Dc,'Winding Plot 4.tab'));           I1c_F=rmsw(cur,2,tc); catch, I1c_F=NaN; end
try, cur0=rd(fullfile(D0,'Winding Plot 4.tab'));           I10_F=rmsw(cur0,2,tc);catch, I10_F=NaN; end
try, vi0 =rd(fullfile(D0,'Winding Plot 2.tab'));           E10_F=rmsw(vi0,2,tc); catch, E10_F=NaN; end

fprintf('---- TABLE 16 : trois conditions, valeurs PLEINES ----\n');
fprintf('  %-26s %14s %14s %12s\n','grandeur','reseau P1','EF','ecart');
p6=@(l,a,f)fprintf('  %-26s %14.6f %14.6f %11.4f %%\n',l,a,f,100*(a-f)/f);
fprintf('  -- EN CHARGE (s = %.4f) --\n',s_ch);
p6('couple [N.m]',r.Tem,121.63);
p6('I1 [A]',r.I1,19.73);
p6('I barre [A]',r.I2*kb,324.7);
p6('I anneau [A]',r.I2*kb*kr_ring,1091);
p6('P entree [W]',r.Pin,20369);
p6('pertes fer [W]',r.Pfe,232.6);
fprintf('  -- A VIDE --\n');
p6('I magnetisant [A]',r0.I1,8.49);
p6('f.e.m. entrefer [V]',r0.E1,382.1);
p6('pertes fer [W]',r0.Pfe,249.3);
fprintf('  -- ROTOR BLOQUE (s = 1) --\n');
p6('couple [N.m]',rb.Tem,104.31);
p6('I1 [A]',rb.I1,108.21);
p6('I barre [A]',rb.I2*kb,1929);
p6('I anneau [A]',rb.I2*kb*kr_ring,6530);
fprintf('  (rapport anneau/barre ideal = %.6f, applique au reseau)\n',kr_ring);
fprintf('  controle EF relu : couple charge %.4f | I1 charge %.4f\n',Tc_F,I1c_F);
fprintf('                     I0 a vide %.4f | E1 a vide %.4f\n',I10_F,E10_F);

%% ======================= TABLE 17 ======================================
%  Chaine DISTINCTE : maillage raffine mec.mesh_refined + solve_mesh.
fprintf('\n---- TABLE 17 : densites locales, chaine du MAILLAGE ----\n');
n1=1469.772776; s1=(M.ns*60-n1)/(M.ns*60);
r1f=mec.equivalent_circuit(ctx,s1,ctx.Xm0);
mk=@(rr,ph)deal_mesh(M,G,BH,W,rr,ph);
try
    [me0,Se0]=mk(r0,0);            % a vide
    [me1,Se1]=mk(r1f,0);           % en charge
    reg={'culasse stator','dent stator','dent rotor','culasse rotor'};
    F0=[1.84 1.63 1.93 1.67];      % sondes EF a vide
    F1=[1.90 1.68 1.87 1.58];      % sondes EF en charge
    %  regional_max renvoie [Bts,Bys,Btr,Byr] ; l'ordre du tableau est
    %  culasse stator, dent stator, dent rotor, culasse rotor.
    [ts0,ys0,tr0,yr0]=mec.regional_max(me0,Se0); B0=[ys0 ts0 tr0 yr0];
    [ts1,ys1,tr1,yr1]=mec.regional_max(me1,Se1); B1v=[ys1 ts1 tr1 yr1];
    fprintf('  %-22s %12s %12s %11s | %12s %12s %11s\n', ...
        'region','MEC vide','EF vide','ecart','MEC charge','EF charge','ecart');
    for k=1:4
        fprintf('  %-22s %12.6f %12.6f %10.4f %% | %12.6f %12.6f %10.4f %%\n', ...
            reg{k},B0(k),F0(k),100*(B0(k)-F0(k))/F0(k), ...
            B1v(k),F1(k),100*(B1v(k)-F1(k))/F1(k));
    end
catch ME
    fprintf('  NON PRODUITE : %s\n',ME.message);
    fprintf('  (la chaine mesh_refined n''a pas pu etre appelee ici ;\n');
    fprintf('   la Table 17 reste a regenerer depuis RUN_ARTICLE section IV.)\n');
end

save('C4b_tables_16_17.mat','r','r0','rb','kb','kr_ring');
fprintf('\n  duree %.0f s\n=== C4b termine ===\n',toc(t0));
diary off;

% ======================================================================
function [me,S]=deal_mesh(M,G,BH,W,rr,ph)
    [i3,ib]=deal_currents(M,W,rr);
    me=mec.mesh_refined(M,G,6,i3,ph,3,ib);
    S=mec.solve_mesh(me,BH,M.opt);
end
function [i3,ib]=deal_currents(M,W,r)
    psi1=angle(r.I1c); psi2=angle(r.I2c);
    i3=sqrt(2)*r.I1*[cos(psi1);cos(psi1-2*pi/3);cos(psi1+2*pi/3)];
    Ibar=r.I2*(2*M.m*W.kw1*W.Nph)/M.Nr;
    ib=sqrt(2)*Ibar*cos(M.p*(2*pi*(0:M.Nr-1)'/M.Nr)-psi2);
end
