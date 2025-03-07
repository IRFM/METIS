% ZRM_RAPSAUVE supprime les fichiers intermediaires de sauvegarde rapide
%------------------------------------------------------------------------
% fichier zrm_rapsauve.m ->  zrm_rapsauve
%
%
% fonction Matlab 5 :
%
% Cette fonction supprime les fichiers intermediaires de sauvegarde rapide
%
% syntaxe  :
%  
%       cr =zrm_rapsauve(file,kmin,kmax)
%       
%   ou   
%      
%       cr =zrm_rapsauve(file);
%       
% entrees :
%
%       file = chemin + racine des fichiers de sauvegarde rapide 
%              (par exemple : /usr/drfc/toto/zineb/data/rapsauve/choc007)
%       kmin = premier indice a supprimer
%       kmax = dernier indice a supprimer       
%
%  sorties :
% 
%       cr   = compte rendu d'execution
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 17/01/2005.
% 
% 
% liste des modifications : 
%
%    * 22/05/2002 -> modification du cas nargin = 1 pour qu'il marche en reprise
%    * 17/01/2005 -> remplace rm par rm -f
%
%--------------------------------------------------------------
%
function cr =zrm_rapsauve(file,kmin,kmax)

cr =0;
if nargin  == 3
	for k=kmin:kmax
		nom=strcat(file,'_',int2str(k),'.mat');
		[s,t] = unix(['rm -f ',nom]);
		if s ~=0
			cr =s;
			disp('Probleme lors de la suppression des fichiers intermediaires');
			disp(t)
			return
		end
	end
elseif nargin == 2
	s = kmin;
	k = 1;
	while s == 0
		nom=strcat(file,'_',int2str(k),'.mat');
		[s,t] = unix(['rm -f ',nom]);
		k = k + 1;
	end
else
	%s = 0;
	%k = 1;
	%while s == 0
	%	nom=strcat(file,'_',int2str(k),'.mat');
	%	[s,t] = unix(['rm -f ',nom]);
	%	k = k + 1;
	%end
	[p,f] =fileparts(file);
	[s,t] = unix(sprintf('find %s -name "%s_*.mat" -exec rm -f {} \\;',p,f));
	if s ~=0
			disp('Probleme lors de la suppression des fichiers intermediaires');
			disp(t)
	end
end

