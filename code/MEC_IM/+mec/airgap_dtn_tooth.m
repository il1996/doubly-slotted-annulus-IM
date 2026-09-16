function AT = airgap_dtn_tooth(M, G, phi_rot, nT, nO, Nh, basis)
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
%  Y_OO est definie positive (bloc principal d'une matrice SDP de rang
%  plein) : la factorisation est licite.
W  = Yoo\Yot;                              % (nO x nT)
Sc = Ytt - Yto*W;                          % complement de Schur
Yc = P.'*Sc*P;
Yc = full(0.5*(Yc+Yc.'));                  % symetrisation numerique

AT.Y=Yc; AT.AF=AF; AT.P=P; AT.iT=iT; AT.iO=iO; AT.W=W;
AT.Ns=Ns; AT.Nr=Nr; AT.Msf=Msf; AT.Mrf=Mrf; AT.Nh=Nh;
AT.ths=ths; AT.dths=dths; AT.thr=thr; AT.dthr=dthr;
AT.isFace=isFace; AT.phi=phi_rot; AT.X=AF.X;
AT.expand=@(Ut) expand_fine(AT,Ut);
AT.field =@(Ut,r,thq) fieldT(AT,Ut,r,thq);
AT.torque=@(Ut) torqueT(AT,Ut);
end

% ======================================================================
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
    Uf(AT.iO)=-AT.W*UT;                    % condition Phi_O = 0
end

function [Br,Bt]=fieldT(AT,Ut,r,thq)
    Uf=expand_fine(AT,Ut);
    [Br,Bt]=AT.AF.field(Uf(1:AT.Msf),Uf(AT.Msf+1:end),r,thq);
end

function T=torqueT(AT,Ut)
    Uf=expand_fine(AT,Ut);
    T=AT.AF.torque(Uf(1:AT.Msf),Uf(AT.Msf+1:end));
end
