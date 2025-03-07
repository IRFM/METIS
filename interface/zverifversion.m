% ZVERIFVERSION  verifie la coherence entre le numero de version de Cronos et celle des donnees
%------------------------------------------------------------------------------- 
% fichier : zverifversion.m  ->  zverifversion
% 
% fonction Matlab 5 : 
%	verifie la coherence entre le numero de version courante
%	de cronos et le numero de version des donnees
%
% syntaxe :  
% 	zverifversion
%  
% entrees :  
%  
% sorties :  
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
%  
% liste des modifications :  
%  
%-------------------------------------------------------------------------------  
%  
function zverifversion

% version courante
[zver,zdate] = zinebversion;

% version des donnees
version = evalin('base','param.gene.version_zineb');

if zver ~= version
	% on affiche un warning en modal
	h = warndlg(sprintf('La version courante de Cronos est la %g \net la version des donnees est la %g :\n%s', ...
	             zver,version,'Vous devez mettre a jour les donnees avant toute execution de Cronos'), ...
		     'Difference de numero de version');
end
