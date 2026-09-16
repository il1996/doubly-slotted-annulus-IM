function kC = carter(tau, b0, g)
%CARTER  Coefficient de Carter d'une denture.
%
%   kC = mec.carter(tau, b0, g) renvoie le coefficient de Carter pour un
%   pas d'encoche tau, une ouverture d'encoche b0 et un entrefer g :
%
%       gamma = (b0/g)^2 / (5 + b0/g)
%       kC    = tau / (tau - gamma*g)
%
%   Le coefficient de Carter traduit l'augmentation apparente de
%   l'entrefer due à l'ouverture des encoches (perte de perméance
%   d'entrefer). Grandeur PHYSIQUE (dépend de b0, g), indépendante de
%   toute discrétisation.

if b0 <= 0
    kC = 1;
    return;
end
ratio = b0/g;
gamma = ratio^2 / (5 + ratio);
kC = tau / (tau - gamma*g);
end
