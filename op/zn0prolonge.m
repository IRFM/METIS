% ZN0PROLONGE prolonge le profil de densite de neutre froid vers le centre
%---------------------------------------------------------------
% fichier zn0prolonge.m ->  zn0prolonge
%
%
% fonction Matlab 5 :
%
% Cette fonction prolonge le profil de densite de neutre froid 
% vers le centre. En effet les fonctions qui calcul la densite
% de neutres froid en provenance du bord effectue le calcul sur 
% une tranche de faible eppaisseur (typiqument 0.2 m). Certain
% module ont besoin de la densite jusqu'au centre (NeoClass). 
% Le profil est prolonge en utilisant une fonction de la forme :
%    n = n0 exp( - sqrt( 1 - x^2)/la)
% 
% n0 et la sont calculees a partir de la densite donnee en entree
% les points donnees en entree ne sont pas modifies en sortie.
%  
% syntaxe  :
%  
%     np = zn0prolonge(x,n);
%    
% entree :
%
%     x  = coordonnee normalisee (param.gene.x)
%     n  = profil de densite de neutre a prolonger
%
% sorties :
% 
%     np = profil de densite de neutre prolonger
%     
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 25/10/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function np = zn0prolonge(x,n)

% point de depart du prolongement
ind = min(find(n >0));
if isempty(ind)
	np =n;
	return;
elseif ind ==1
	np =n;
	return;
end

% raccordement
xc = x(ind);
nc = n(ind);
n0 = n(end);

% longueur de decroissance
la = - sqrt(1 - xc .^ 2) ./ log(nc ./ n0);

% nouveau profil
nn = n0 .* exp( - sqrt(1 - x .^ 2) ./ la);

% profil de sortie
np            = n;
np(1:(ind-1)) = nn(1:(ind-1));
