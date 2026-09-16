function sat = saturation_curve(ctx, im_list)
%SATURATION_CURVE  Courbe de magnétisation saturée Lm(Im) par le réseau MEC.
%
%   sat = mec.saturation_curve(ctx, im_list) évalue, pour une liste de
%   courants magnétisants efficaces im_list [A], l'inductance magnétisante
%   SATURÉE Lm issue du réseau de réluctances (mec.magnetizing). Renvoie un
%   interpolant Lm(im) réutilisable (ex. dans le modèle dq transitoire, où
%   la saturation du flux principal évolue pendant le démarrage).
%
%   sat.im   : courants magnétisants [A]
%   sat.Lm   : inductances magnétisantes saturées [H]
%   sat.Lmof : fonction interpolée Lm(im) (extrapolation plate aux bornes)
%
%   Voir aussi : mec.magnetizing, mec.dq_startup.

if nargin < 2 || isempty(im_list)
    im_list = [0.2 1 2 4 6 8 10 13 16 20 25 32];
end
im = im_list(:);
Lm = zeros(size(im));
for k = 1:numel(im)
    R = mec.magnetizing(ctx, im(k));
    Lm(k) = R.Lm;
end
sat.im = im;  sat.Lm = Lm;
sat.Lmof = @(x) interp1(im, Lm, min(max(abs(x),im(1)),im(end)), 'linear');
end
