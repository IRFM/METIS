function [fn,inte] = fluxneutron(td,nd,ro,Rm,a)

%function fn = fluxneutron(td,nd,ro,R,a[,d])
%
% Input: 
%  td : profil temperature Deuterium (keV)
%        (autant de colonne que de profil)
%  nd : profil densite de Deuterium (M-3) 
%        (autant d'elements que de lignes dans td)
%  ro : rayon normalise
%        (autant d'elements que de lignes dans td)
%  Rm : rayon des surfaces de flux (M)
%  a  : petit rayon plasma (M)
%
% Output: 
%  fn : flux de neutron (S-1)
%
% Warning:
%  Taux de reaction de fusion DD d'apres these G. Martin
%  (fonction sigvdd) (La moitie de ces fusions donne un neutron...)

%Last Update: 11/12/96 11:26:13
%
%by Thierry HUTTER			CEA-DRFC 10/12/96
%Tel: 4840				e-mail: hutter@cea.fr

% fn=4 pi^2 a^2 Somme[0->1](1/2 nd^2 sigmav (R+ delta(ro)) ro dro

[l,c] = size(td);
dro   = gradient(ro);
pia2  = (pi .* a) .^ 2;

ppp   = sigvddn(td) .* nd .* nd;
ppx   = ppp .* (ones(size(td,1),1) * (ro .* dro)) .* Rm;

% 4 x 1/2 x 1/2 = 1 !

fn    = pia2 .* (sum(ppx(:,1:c-1)') + sum(ppx(:,2:c)'))';
inte  = pia2(1) .* (ones(size(td,1),1) * (ro .* dro)) .* Rm .* nd .* nd;
inte  = inte ./ (max(nd' .* nd')' * ones(1,size(nd,2)));
