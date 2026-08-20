%% RUN_R2_KC_TILING - sensibilite de k_C au pavage, en base CHAPEAU
%
%  BLOC R2 de SPEC_CLAUDE_CODE_v4. Le rapporteur : "publier une correction
%  de +5,2 % a un resultat centenaire sans un tableau k_C(n_T, n_O) en base
%  chapeau est insuffisant". Le 0,52 % d'invariance actuellement cite vient
%  de la base CONSTANTE PAR MORCEAUX au pavage 6/2 -- ni la base ni le
%  pavage de production.
%
%  DEUX DEFINITIONS DU X_m LISSE, et il faut les distinguer.
%
%   (a) DEFINITION PUBLIEE (B2) : arcs pavant 100 % de l'alesage a UN NOEUD
%       PAR DENT. Elle NE DEPEND PAS de (n_T, n_O). C'est la reference du
%       k_C = 1,3320 publie ; la regle 5 interdit de la changer.
%   (b) DEFINITION A PAVAGE APPARIE (diagnostic) : (n_T + n_O) colonnes
%       uniformes par pas, TOUTES traitees comme des faces. Elle suit le
%       pavage et rend la comparaison symetrique.
%
%  Les deux sont produites. (a) fait foi ; (b) sert la GARDE.
%
%  GARDE (v4 §R2). Le rapport k_C doit etre PLUS STABLE que chacune de ses
%  composantes. Point decisif a etablir : avec la definition (a), le
%  NUMERATEUR est CONSTANT par construction sur l'axe du pavage. Il ne peut
%  donc y avoir AUCUNE compensation numerateur/denominateur sur cet axe --
%  contrairement a l'axe des TRONCATURES, ou le §4.3 a refute un argument
%  de compensation. L'analogie que releverait un rapporteur ne tient pas,
%  et ce script le montre par les chiffres plutot que par l'argument.
clear; clc; t0=tic;
diary('R2_kc_tiling_out.txt'); diary on;
M=mec.machine_18_5kW(); ctx=mec.build_context(M); G=ctx.G;
Rs=G.Rs; Rr=G.Rr; L=M.L; Nh=8192;
taus=2*pi/M.Ns; taur=2*pi/M.Nr;
thsT=2*pi*(0:M.Ns-1)/M.Ns; thrT=2*pi*(0:M.Nr-1)/M.Nr;
nTl=[9 17 33]; nOl=[2 4 8];
nTp=17; nOp=4;                      % pavage de PRODUCTION

fprintf('=== R2 : k_C(n_T, n_O) en base chapeau ===\n\n');
fprintf('  CONFIGURATION\n');
fprintf('    machine       : MAS 48/44, 18,5 kW | Ns=%d Nr=%d p=%d\n',M.Ns,M.Nr,M.p);
fprintf('    base          : P1 (chapeau)\n');
fprintf('    troncature    : N_h = %d, FIXE sur les neuf points\n',Nh);
fprintf('    operateur     : mec.airgap_dtn_tooth (condense par Schur)\n');
fprintf('    X_m lisse (a) : mec.airgap_fourier, UN noeud par dent, arcs\n');
fprintf('                    pavant 100 %% -- INDEPENDANT du pavage\n');
fprintf('    X_m lisse (b) : (n_T+n_O) colonnes uniformes, toutes faces\n');
fprintf('    courant       : I_m = 0,2 A (reference non saturee)\n');
fprintf('    pavage de production : n_T = %d, n_O = %d\n',nTp,nOp);
fprintf('    Carter classique : %.6f\n',ctx.AGcarter.kC);

%% ---- X_m lisse, definition (a) : calcule UNE fois ---------------------
AL=mec.airgap_fourier(thsT,repmat(taus,1,M.Ns),thrT,repmat(taur,1,M.Nr), ...
                      Rs,Rr,L,Nh,'p1');
cl=ctx; cl.AG=AL; XmS_a=mec.magnetizing(cl,0.2).Xm;
fprintf('\n  X_m lisse (a), invariant du pavage : %.6f ohm\n',XmS_a);

%% ---- la grille des neuf points ---------------------------------------
nP=numel(nTl)*numel(nOl);
R=nan(nP,7);   % nT nO Ncol XmS_b XmD kC_a kC_b
i=0;
fprintf('\n  %4s %4s %8s %12s %12s %12s %10s %10s\n', ...
    'n_T','n_O','colonnes','Xm lisse(b)','Xm encoche','k_C (a)','k_C (b)','t [s]');
for nT=nTl
  for nO=nOl
    i=i+1; tk=tic;
    A=mec.airgap_dtn_tooth(M,G,0,nT,nO,Nh,'p1');
    cd_=ctx; cd_.AG=A; XmD=mec.magnetizing(cd_,0.2).Xm;
    %  (b) lisse a pavage APPARIE : (nT+nO) colonnes uniformes par pas
    nc=nT+nO;
    ds=taus/nc; dr=taur/nc;
    ths=reshape((0:M.Ns-1).'*taus + ((0:nc-1)+0.5)*ds,1,[]);
    thr=reshape((0:M.Nr-1).'*taur + ((0:nc-1)+0.5)*dr,1,[]);
    ALb=mec.airgap_fourier(ths,repmat(ds,1,numel(ths)), ...
                           thr,repmat(dr,1,numel(thr)),Rs,Rr,L,Nh,'p1');
    %  CONDENSATION. Une surface LISSE n'a aucune colonne d'ouverture :
    %  il n'y a donc RIEN a eliminer par complement de Schur. Les colonnes
    %  d'une meme dent sont simplement equipotentielles, et la condensation
    %  se reduit a la prolongation Y_c = P' Y P.
    %  NB : en base P1 cette subdivision N'EST PAS neutre -- la somme de
    %  chapeaux sur une subdivision ne represente pas le meme profil qu'un
    %  chapeau unique sur l'arc entier. C'est ce qui rend (b) informatif.
    tS=repelem((1:M.Ns).',nc); tR=repelem((1:M.Nr).',nc)+M.Ns;
    tt=[tS;tR]; nAll=numel(tt);
    P=sparse(1:nAll,tt,1,nAll,M.Ns+M.Nr);
    Yb=full(0.5*(P.'*ALb.Y*P + (P.'*ALb.Y*P).'));
    %  On NE pose PAS de champ 'expand' : sa presence ferait appeler
    %  AG.field a TROIS arguments dans mec.magnetizing (signature de
    %  l'operateur condense), alors que celui herite d'airgap_fourier en
    %  attend QUATRE. Seul AG.Y est necessaire ici.
    ALc=rmfield(ALb,intersect(fieldnames(ALb),{'expand'}));
    ALc.Y=Yb;
    clb=ctx; clb.AG=ALc; XmS_b=mec.magnetizing(clb,0.2).Xm;
    R(i,:)=[nT nO A.Msf+A.Mrf XmS_b XmD XmS_a/XmD XmS_b/XmD];
    fprintf('  %4d %4d %8d %12.6f %12.6f %12.6f %10.6f %10.1f\n', ...
        nT,nO,R(i,3),XmS_b,XmD,R(i,6),R(i,7),toc(tk));
  end
end

%% ---- ecart a la valeur de production ---------------------------------
ip=find(R(:,1)==nTp & R(:,2)==nOp,1);
kCp=R(ip,6);
fprintf('\n  ---- ecart relatif a la valeur de production (n_T=%d, n_O=%d) ----\n',nTp,nOp);
fprintf('  %4s %4s %14s %12s %14s\n','n_T','n_O','k_C (a)','ecart','X_m encoche');
for i=1:nP
    fprintf('  %4d %4d %14.6f %11.4f %% %14.6f\n', ...
        R(i,1),R(i,2),R(i,6),100*(R(i,6)-kCp)/kCp,R(i,5));
end

%% ================= GARDE ==============================================
d=@(x)100*(max(x)-min(x))/mean(x);
fprintf('\n  ---- GARDE : le rapport est-il plus stable que ses composantes ? ----\n');
fprintf('  %-34s %12s\n','grandeur','dispersion');
fprintf('  %-34s %11.4f %%   <- CONSTANT par construction\n', ...
    'X_m lisse (a), def. publiee',0);
fprintf('  %-34s %11.4f %%\n','X_m encoche',d(R(:,5)));
fprintf('  %-34s %11.4f %%\n','k_C (a) = lisse(a)/encoche',d(R(:,6)));
fprintf('  %-34s %11.4f %%\n','X_m lisse (b), pavage apparie',d(R(:,4)));
fprintf('  %-34s %11.4f %%\n','k_C (b) = lisse(b)/encoche',d(R(:,7)));

fprintf(['\n  LECTURE DE LA GARDE, ET REPONSE A L''OBJECTION ANTICIPEE.\n' ...
  '  Avec la definition publiee (a), le NUMERATEUR est CONSTANT sur\n' ...
  '  l''axe du pavage : sa dispersion est nulle PAR CONSTRUCTION, non\n' ...
  '  par compensation. La dispersion de k_C (a) est donc EXACTEMENT\n' ...
  '  celle de X_m encoche -- il n''existe aucun degre de liberte pour\n' ...
  '  une compensation numerateur/denominateur.\n' ...
  '  C''est ce qui distingue cet axe de celui des TRONCATURES, ou le\n' ...
  '  §4.3 a refute un argument de compensation : la, les DEUX termes\n' ...
  '  variaient. Ici un seul varie. L''analogie qu''un rapporteur\n' ...
  '  releverait ne tient pas, et le tableau le montre.\n' ...
  '  La colonne (b) est le controle : si k_C (b) etait nettement plus\n' ...
  '  stable que ses deux composantes, il y aurait compensation SUR CET\n' ...
  '  AXE et il faudrait le dire.\n']);

%% ---- tableau LaTeX pret a inserer -------------------------------------
fprintf('\n  ---- TABLEAU LaTeX ----\n\n');
fprintf('\\begin{tabular}{rrrrrr}\n\\toprule\n');
fprintf('$n_T$ & $n_O$ & columns & $X_m$ smooth & $X_m$ slotted & $k_C$ \\\\\n');
fprintf('\\midrule\n');
for i=1:nP
    st=''; if i==ip, st='\\textbf'; end %#ok<NASGU>
    if i==ip
        fprintf('\\textbf{%d} & \\textbf{%d} & \\textbf{%d} & \\textbf{%.3f} & \\textbf{%.3f} & \\textbf{%.4f} \\\\\n', ...
            R(i,1),R(i,2),R(i,3),XmS_a,R(i,5),R(i,6));
    else
        fprintf('%d & %d & %d & %.3f & %.3f & %.4f \\\\\n', ...
            R(i,1),R(i,2),R(i,3),XmS_a,R(i,5),R(i,6));
    end
end
fprintf('\\bottomrule\n\\end{tabular}\n');
fprintf('\n  (ligne en gras = pavage de production ; X_m smooth est\n');
fprintf('   invariant du pavage par definition)\n');

save('R2_kc_tiling.mat','R','XmS_a','Nh','nTp','nOp');
fprintf('\n  duree totale %.0f s\n=== R2 termine ===\n',toc(t0));
diary off;
