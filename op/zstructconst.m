% ZSTRUCTCONST met a une valeur constante les donnees d'une structure
%----------------------------------------------------------
% fichier zstructconst.m ->  zstructconst
%
%
% fonction Matlab 5 :
%
% Cette fonction met a une valeur constante les donnees d'une structure.
% La valeur par defaut est 0
% 
%
% syntaxe  :
%  
%       data = zstructconst(data,{valeur});
%       
% entrees :
%
%       data   = structure a modifier
%       valeur = valeur a mettre                    
% sorties :
% 
%       data   = structure  modifiee
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 21/06/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function data = zstructconst(data,valeur)

% test des entrees
if nargin <2
	valeur = 0;
elseif isempty(valeur)
	valeur = 0;
end

% liste des champs la structures
champ = fieldnames(data);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
    		champnew{k}=strcat(champc,'.',champnew{k});
    	end
    	% ajout a la liste des champs
    	if isempty(champ)
    		champ =champnew;
    	else
    		champ=cat(1,champ,champnew);
    	end
    else
    	eval(strcat(champc,' = valeur .* ones(size(',champc,'));'));
    end
end
