% calcul la puissance de fusion du a l'interaction faisceau-plasma (idn)
% salpha est le nombre de neutron dd
function [salpha,emean] = znbi_dd(nD,ti,pnbi_th,taus_nbi,e0nbi,ecritnbi)

persistent svt
if isempty(svt)
   svt = load('sigmav_nbi_d_plasma_d');
end

if length(e0nbi) == 1
    e0nbi = ones(size(nD)) .* e0nbi;
end
e0nbi(e0nbi == 0) = 1e6;

% fonction de distribution des ions D issu de l'injection de neutre
mp          = 1.6726485e-27;            % masse au repos du proton (kg)
e           = 1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
sn0         = pnbi_th ./ (e .* e0nbi);
vc          = sqrt(2 .*e .*  ecritnbi ./ mp ./ 2);
v0          = sqrt(2 .*e .*  e0nbi ./ mp ./ 2);
inde        = min(find(max(e0nbi) <= svt.enbi));
if isempty(inde)
   inde = length(svt.enbi);
end
vi          = sqrt(2 .*e .*  svt.enbi(1:inde)' ./ mp ./ 2);
ut          = ones(size(nD));
ue          = ones(size(vi));
fn0         = ((sn0 .* taus_nbi) * ue)./ 4 ./ pi .* (((v0 * ue) - ut * vi) >0) ./ ((ut*vi) .^ 3 + (vc*ue) .^ 3);
emean       = trapz(vi,fn0.*(ut* svt.enbi(1:inde)'),2) ./ max(trapz(vi,fn0,2),eps);
 
% integration sur la fonction des ions rapides
svnbidd     = svt.svnbidd(1:inde,:);
svm         = tsplinet(ones(size(svnbidd,1),1)*svt.eti,svnbidd,ones(size(svnbidd,1),1)*min(ti',max(svt.eti)))';
salpha      = 4 .* pi .* nD .* trapz(vi,fn0 .* svm .* (svm >0) .* (ut *vi) .^ 2,2);
% securite
salpha(pnbi_th <=0) = 0;

