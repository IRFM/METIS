% calcul la puissance de fusion du a l'interaction faisceau-plasma (idn)
function [salpha,emean] = zbpi0_pb(nB,ti,pnbi_th,taus_nbi,e0nbi,ecritnbi,ftnbi,flag)

persistent svs
if isempty(svs)
   svs = load('sigmav_nbi_p_plasma_b');
end

% constantes
phys = cphys;

% fonction de distribution des ions boron issu de l'injection de neutre
sn0         = ftnbi .* pnbi_th ./ (phys.e .* e0nbi);
vc          = sqrt(2 .* phys.e .*  ecritnbi ./ phys.mp);
v0          = sqrt(2 .* phys.e .*  e0nbi ./ phys.mp);
inde        = min(find(e0nbi <= svs.enbi));
if isempty(inde)
   inde = length(svs.enbi);
end
vi          = sqrt(2 .* phys.e .*  svs.enbi(1:inde)' ./ phys.mp);
ut          = ones(size(nB));
ue          = ones(size(vi));
fn0         = ((sn0 .* taus_nbi) * ue)./ 4 ./ pi .* ((v0 .* (ut* ue) - ut * vi) >0) ./ ((ut*vi) .^ 3 + (vc*ue) .^ 3);
emean       = trapz(vi,fn0.*(ut* svs.enbi(1:inde)'),2) ./ max(trapz(vi,fn0,2),eps);
 
% integration sur la fonction des ions rapides
svnbipb     = svs.svnbipb(1:inde,:);
svm         = tsplinet(ones(size(svnbipb,1),1)*svs.eti(:)',svnbipb,ones(size(svnbipb,1),1)*min(ti',max(svs.eti)))';
salpha      = 4 .* pi .* nB .* trapz(vi,fn0 .* svm .* (ut *vi) .^ 2,2) * 3; % reaction with 3 alpha particles
% securite
salpha(pnbi_th <=0) = 0;

if nargin > 7
 	disp('in zbpi0')
	keyboard
end
%disp(salpha(end-2)' .* 3.56e6 .* 1.602176462e-19)






function cs = csdt(ecm)
%ref : improved formulas for fusion cross-sctions and thermal reactivities,
% H-S Bosch and G.M. Hale, NF vol 32, p 611- , 1992.
% ecm en keV et cs en 1e-3 barn , soit 1e-3 * 1e-4 * 1e-24 m ^ 2 = 1e-31 m ^ 2

% ecm  <550 keV
a1 = 6.927e4; a2 = 7.454e8;  a3 = 2.05e6;   a4 = 5.2002e4; a5 = 0;
b1 = 6.38e1;  b2 = -9.95e-1; b3 = 6.981e-5; b4 = 1.728e-4;
sse1 = (a1 + ecm .* (a2 + ecm .* (a3 + ecm .* (a4 + ecm .* a5)))) ./ ...
       (1 + ecm .* (b1 + ecm .* (b2 + ecm .* (b3 + ecm .* b4))));
% de 550 a 4700
sse2 = -1.4714e6 ./(1 + ecm .* (-8.4127e-3 + ecm .* (4.7983e-6 + ecm .* (-1.0748e-9 +ecm .* 8.5184e-14))));
see = sse1 .* (ecm <= 550) + sse2 .* (ecm >550);
cs = sse ./ ecm ./exp(34.3827./sqrt(ecm));
