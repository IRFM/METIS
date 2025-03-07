% ZPELLETCOOR cette fonction calcule sur la droite d'ablation les grandeurs plasmas
%---------------------------------------------------------------------------------------
% fichier zpelletcoor.m ->  zpelletcoor
%
%
% fonction Matlab 7 :
%
% Cette fonction calcule pour l'injection de glacon, sur la droite d'ablation les grandeurs plasmas 
%  
% syntaxe  :
%  
%         [r,z]=zthetapsi2rz(equi,theta,psi)
%
% entrees :
%
%    equi   =  structure  equi.
%    theta  =  angle 0..2*pi 
%    psi    =  flux poloidal normalise (0 au centre et 1 au bord)
%
% sorties :
% 
%     r,z            = coordonnees (R,Z) des points correspondants  (theta,psi)
%     
% fonction ecrite par J-F Artaud
% version 4.0, du 27/11/2007.
% 
% 
% liste des modifications :  cf cvs.
% 
%---------------------------------------------------------------------------------
%



%
% test : 
%    datak = zget1t(data,k);
%    psi = rand(30000,1);  
%    theta  = 2.*pi.*rand(size(psi));
%    [r,z] = zthetapsi2rz(datak.equi,theta,psi);
%
%
%
function [r,z] = zthetapsi2rz(equi,theta,psi)

% mise en ligne
sout  = size(theta);
theta = theta(:);
psi   = psi(:);
ind0 = find(psi  == 0);
psi(psi == 0) = eps;

% donnees de l'equilibre
req = double(squeeze(equi.R(1,:,2:(end-1)))); 
zeq = double(squeeze(equi.Z(1,:,2:(end-1)))); 
vth = ones(1,size(req,2));
psi_eq = double(squeeze(equi.psiRZ))'; 
% normalisation (psi = 0 au centre , psi = 1 au bord);
psi_eq = (psi_eq - psi_eq(1)) ./ (psi_eq(end) - psi_eq(1));

% interpolation sur psi
ru = interp1(psi_eq,req,psi,'pchip','extrap');
zu = interp1(psi_eq,zeq,psi,'pchip','extrap');

% calcul de l'angle
cl   = (ru - req(1,1)) + sqrt(-1) .* (zu - zeq(1,1));
thl  = unwrap(angle(cl),[],2);
thl(thl<0) = thl(thl<0) + 2 .* pi;
thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
rl  = cat(2,ru,ru,ru);
zl  = cat(2,zu,zu,zu);
% interpolation sur tetha
r = tsplinet(thl,rl,theta);
z = tsplinet(thl,zl,theta);
% correction 0
r(ind0) = req(1,1);
z(ind0) = zeq(1,1);








