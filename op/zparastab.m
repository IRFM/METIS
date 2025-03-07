% parametre pour la stabilite
function [alpha,ate,ati,ane,ani,gE,shear] = zparastab(param,data)

% le scisaillement magnetique
shear = data.prof.shear;
% le scisaillement de rotation
gE    = data.neo.gammae;
% parametre alpha
alpha = data.neo.alpha;

% les alpha
Epsilonx = data.equi.rhomax.* param.gene.x ./ data.equi.rmoy;
ate      = data.prof.gte ./ data.prof.te .* data.equi.rmoy .* (1 + Epsilonx);
ati      = data.prof.gti ./ data.prof.ti .* data.equi.rmoy .* (1 + Epsilonx);
ane      = data.prof.gne ./ data.prof.ne .* data.equi.rmoy .* (1 + Epsilonx);
ani      = data.prof.gni ./ data.prof.ni .* data.equi.rmoy .* (1 + Epsilonx);
