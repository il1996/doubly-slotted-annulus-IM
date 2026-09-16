function AT = airgap_dtn_tooth_cav(M, G, phi_rot, nT, nO, Nh, basis, cav)
%AIRGAP_DTN_TOOTH_CAV  Variante de airgap_dtn_tooth avec l'admittance des
%   cavites d'encoche ajoutee sur les colonnes d'ouverture avant le
%   complement de Schur (cav.Qs, cav.Qr : matrices (nO+2)^2 calculees par
%   elements finis dans le repere local d'une encoche ; dofs = colonnes de
%   bouche 1..nO, dent gauche, dent droite (rotor : paroi unique en nO+1)).
%   Sans 'cav', identique a airgap_dtn_tooth.
%
%   COURANT DE BARRE (optionnel). Si cav.sr est present (vecteur (nO+1)x1,
%   flux par ampere de courant de barre sur les nO colonnes de bouche puis
%   sur la paroi, calcule par cavity_graded.cavity_source_rotor), l'element
%   accepte une source : AT2 = AT.set_source(Iz) avec Iz (Nr x 1) le courant
%   de la barre logee dans l'encoche qui SUIT la dent rotorique i (sens +z,
%   convention U_droite - U_gauche = -Iz). AT2.f est le vecteur de flux
%   source condense sur les noeuds dentaires (a ajouter au bilan de flux,
%   voir mec.solve_network), et AT2.expand / AT2.field tiennent compte du
%   potentiel d'ouverture induit par la source. Sans appel a set_source,
%   rien ne change.
%AIRGAP_DTN_TOOTH  Operateur DtN sur grille de surface FINE, condense sur
%                  les noeuds dentaires du reseau de performance.
%
%   POURQUOI UNE GRILLE FINE. Construit directement sur un noeud par dent,
%   l'operateur ne peut pas produire l'effet d'encoche : les arcs de face
%   couvrent 81 % de l'alesage et il n'existe AUCUN degre de liberte au
%   droit de l'ouverture. La projection de Fourier y impose implicitement
%   phi = 0, ce qui n'est pas la condition physique, et le coefficient de
%   Carter implicite tend vers 1 quand on raffine en harmoniques (mesure :
%   1.053 a Nh=288, 1.024 a Nh=1400, contre 1.266 pour Carter). Le PMSM
%   n'a pas ce probleme parce qu'il porte 84 noeuds de surface par pas
%   dentaire, dont plusieurs au-dessus de chaque ouverture.
%
%   CONSTRUCTION. La surface est PAVEE : nT colonnes sur la face de dent,
%   nO colonnes sur l'ouverture, de chaque cote. L'operateur exact est
%   assemble sur cette grille, puis condense sur les noeuds dentaires par
%   complement de Schur, avec deux conditions PHYSIQUES et non ajustees :
%
%     (a) les colonnes de FACE d'une meme dent sont equipotentielles --
%         le bec est du fer, il porte le potentiel du noeud dentaire ;
%         U_T = P * U_dent, P matrice de prolongement 0/1 ;
%     (b) les colonnes d'OUVERTURE ne debitent aucun flux dans le fer --
%         c'est de l'air : le flux qui les atteint se redistribue dans la
%         couronne. Phi_O = 0.
%
%   D'ou, en partitionnant Y = [Y_TT Y_TO ; Y_OT Y_OO] :
%       U_O   = -Y_OO \ (Y_OT * P * U_dent)
%       Y_cond = P' * ( Y_TT - Y_TO * (Y_OO \ Y_OT) ) * P
%   symetrique par construction, et sans aucun parametre de frange : la
%   redistribution sous l'ouverture est resolue, pas modelisee.
%
%   AT.Y             : operateur condense, ordre [stator ; rotor]
%   AT.expand(Ut)    : potentiel de surface FIN reconstruit
%   AT.field(Ut,r,thq), AT.torque(Ut) : champs et couple, sur la grille fine
%   AT.AF            : operateur fin (diagnostic)

if nargin<3||isempty(phi_rot), phi_rot=0; end
if nargin<4||isempty(nT), nT=6; end
if nargin<5||isempty(nO), nO=2; end
Ns=M.Ns; Nr=M.Nr; Rs=G.Rs; Rr=G.Rr; L=M.L;
taus=2*pi/Ns; taur=2*pi/Nr;
aS=(G.taus-M.bs0)/Rs;  oS=taus-aS;        % face / ouverture stator [rad]
aR=(G.taur-M.br0)/Rr;  oR=taur-aR;        % face / ouverture rotor
if nargin<6||isempty(Nh)
    %  L'entrefer est mince : les harmoniques restent couples jusqu'a
    %  n ~ 1/X. Il faut aussi resoudre la plus fine colonne.
    X=log(Rs/Rr);
    Nh=ceil(max(4/X, 6*2*pi/min([aS/nT,oS/nO,aR/nT,oR/nO])));
    Nh=min(Nh,4000);
end

%  ---- grille pavante, cote stator ----
[ths,dths,isFaceS,idT_S]=tile(Ns,taus,aS,oS,nT,nO,0);
%  ---- grille pavante, cote rotor (tournee de phi_rot) ----
[thr,dthr,isFaceR,idT_R]=tile(Nr,taur,aR,oR,nT,nO,phi_rot);

if nargin<7||isempty(basis), basis='p0'; end
AF=mec.airgap_fourier(ths,dths,thr,dthr,Rs,Rr,L,Nh,basis);
Msf=numel(ths); Mrf=numel(thr);

%  ---- partition face (T) / ouverture (O) ----
isFace=[isFaceS(:); isFaceR(:)];
iT=find(isFace); iO=find(~isFace);
%  prolongement : colonne de face -> noeud dentaire
tooth=[idT_S(:); idT_R(:)+Ns];            % dents stator 1..Ns, rotor Ns+1..Ns+Nr
P=sparse(1:numel(iT),tooth(iT),1,numel(iT),Ns+Nr);

Y=AF.Y;
Ytt=Y(iT,iT); Yto=Y(iT,iO); Yoo=Y(iO,iO); Yot=Y(iO,iT);
nt=Ns+Nr; nOt=numel(iO);
AT.src_ok=false;
if nargin>=8 && ~isempty(cav)
    %  ---- admittance des cavites d'encoche sur les colonnes d'ouverture ----
    pos=zeros(Msf+Mrf,1); pos(iO)=1:nOt;
    Qoo=zeros(nOt,nOt); Qto=zeros(nt,nOt); Qtt=zeros(nt,nt);
    Qs=cav.Qs; Qr=cav.Qr;
    for i=1:Ns                                   % ouverture apres la dent i
        cols=find(~isFaceS(:) & idT_S(:)==i); io=pos(cols);
        L2=i; R2=mod(i,Ns)+1;
        Qoo(io,io)=Qoo(io,io)+Qs(1:nO,1:nO);
        Qto(L2,io)=Qto(L2,io)+Qs(nO+1,1:nO); Qto(R2,io)=Qto(R2,io)+Qs(nO+2,1:nO);
        Qtt([L2 R2],[L2 R2])=Qtt([L2 R2],[L2 R2])+Qs(nO+1:nO+2,nO+1:nO+2);
    end
    rotor_io=zeros(Nr,nO);
    for i=1:Nr
        cols=find(~isFaceR(:) & idT_R(:)==i)+Msf; io=pos(cols);
        t=Ns+i;
        Qoo(io,io)=Qoo(io,io)+Qr(1:nO,1:nO);
        Qto(t,io)=Qto(t,io)+Qr(nO+1,1:nO);
        Qtt(t,t)=Qtt(t,t)+Qr(nO+1,nO+1);
        rotor_io(i,:)=io(:).';
    end
    B  = Yot*P + Qto.';                          % (nO x nt)
    YoQ = Yoo+Qoo;
    Wt = YoQ\B;                                  % U_O = -Wt*Ut
    Yc = P.'*Ytt*P + Qtt - B.'*Wt;
    Yc = full(0.5*(Yc+Yc.'));
    W  = [];
    AT.Wt=Wt; AT.cav=true;
    if isfield(cav,'sr') && ~isempty(cav.sr)
        AT.src_ok=true; AT.sr=cav.sr(:); AT.rotor_io=rotor_io; AT.YoQ=YoQ; AT.B=B; AT.nOt=nOt; AT.nO=nO;
    end
else
    %  Y_OO est definie positive (bloc principal d'une matrice SDP de rang
    %  plein) : la factorisation est licite.
    W  = Yoo\Yot;                              % (nO x nT)
    Sc = Ytt - Yto*W;                          % complement de Schur
    Yc = P.'*Sc*P;
    Yc = full(0.5*(Yc+Yc.'));                  % symetrisation numerique
    AT.Wt=W*P; AT.cav=false;
end

AT.Y=Yc; AT.AF=AF; AT.P=P; AT.iT=iT; AT.iO=iO; AT.W=[];
AT.Ns=Ns; AT.Nr=Nr; AT.Msf=Msf; AT.Mrf=Mrf; AT.Nh=Nh;
AT.ths=ths; AT.dths=dths; AT.thr=thr; AT.dthr=dthr;
AT.isFace=isFace; AT.phi=phi_rot; AT.X=AF.X;
AT.f=[]; AT.uO=[];
AT=bind_handles(AT);
end

% ======================================================================
function AT=bind_handles(AT)
    AT.expand=@(Ut) expand_fine(AT,Ut);
    AT.field =@(Ut,r,thq) fieldT(AT,Ut,r,thq);
    AT.torque=@(Ut) torqueT(AT,Ut);
    if AT.src_ok
        AT.set_source=@(Iz) set_source(AT,Iz);
    end
end

function AT2=set_source(AT,Iz)
%  Iz (Nr x 1) : courant (+z) de la barre logee dans l'encoche qui suit la
%  dent rotorique i. Flux source sur les colonnes de bouche et les parois,
%  condense sur les noeuds dentaires : f = s_t - B'*(YoQ \ s_O) ;
%  potentiel d'ouverture induit : uO = -(YoQ \ s_O).
    Iz=Iz(:); nO=AT.nO; nt=AT.Ns+AT.Nr;
    sO=zeros(AT.nOt,1); st=zeros(nt,1);
    for i=1:AT.Nr
        io=AT.rotor_io(i,:);
        sO(io)=sO(io)+AT.sr(1:nO)*Iz(i);
        st(AT.Ns+i)=st(AT.Ns+i)+AT.sr(nO+1)*Iz(i);
    end
    y=AT.YoQ\sO;
    AT2=AT; AT2.f=st-AT.B.'*y; AT2.uO=-y; AT2.Iz=Iz;
    AT2=bind_handles(AT2);
end

function [th,dth,isFace,idT]=tile(N,tau,aFace,aOpen,nT,nO,shift)
%  Pave un pas dentaire : nT colonnes de face centrees sur la dent, puis
%  nO colonnes d'ouverture. Les colonnes sont contigues par construction.
K=N*(nT+nO); th=zeros(1,K); dth=zeros(1,K); isFace=false(1,K); idT=zeros(1,K);
k=0; edge=-aFace/2;
for i=1:N
    for c=1:nT
        k=k+1; dth(k)=aFace/nT; th(k)=edge+dth(k)/2; isFace(k)=true; idT(k)=i;
        edge=edge+dth(k);
    end
    for c=1:nO
        k=k+1; dth(k)=aOpen/nO; th(k)=edge+dth(k)/2; isFace(k)=false; idT(k)=i;
        edge=edge+dth(k);
    end
end
th=th+shift;
end

function Uf=expand_fine(AT,Ut)
    Ut=Ut(:);
    Uf=zeros(AT.Msf+AT.Mrf,1);
    UT=AT.P*Ut;
    Uf(AT.iT)=UT;
    Uf(AT.iO)=-AT.Wt*Ut;                   % Phi_O = 0 (ou cavite)
    if ~isempty(AT.uO), Uf(AT.iO)=Uf(AT.iO)+AT.uO; end   % source de courant de barre
end

function [Br,Bt]=fieldT(AT,Ut,r,thq)
    Uf=expand_fine(AT,Ut);
    [Br,Bt]=AT.AF.field(Uf(1:AT.Msf),Uf(AT.Msf+1:end),r,thq);
end

function T=torqueT(AT,Ut)
    Uf=expand_fine(AT,Ut);
    T=AT.AF.torque(Uf(1:AT.Msf),Uf(AT.Msf+1:end));
end
