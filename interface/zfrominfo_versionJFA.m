% ZSHOTINFO   permet d'afficher LA LISTE DES info de param.from.shot
%----------------------------------------------------------------
% fichier zshotinfo_versionJFA.m 
%
% fonction Matlab 5 :
%	Cette fonction permet d'afficher les donnees "param.from.shot" de zineb.
%	Les donnees doivent etre dans l'espace de travail de base.
%	Elles sont contenue dans la structureparam
% 
% syntaxe  :
%     zshotinfo;
%
% entrees :
%
% sorties
%
% fonction écrite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zshotinfo

% liste des donnees
try 
	param.from.shot.info = evalin('base','param.from.shot.info') ;
catch
	param.from.shot.info=[] ;
end

if isempty(param)  
	herror = errordlg('Pas de données "param" dans la base','ATTENTION') ;
	zwaitfor(herror) ;
	zuicloseone(hform) ;
	return
end

if ~isstruct(param)
	return
end
	 
% Ouverture d'un fichier temporaire
% file = tempname
file = 'liste2-shot-info' ;
[fid,mess] = fopen(file,'w') ;
if mess~=''
	herror = errordlg('a l''ouverture du fichier','Probleme') ;
	zwaitfor(herror) ;
	zuicloseone(hform) ;
end

fprintf(fid,' param.from.shot.info \n -------------------- \n') ;

% liste des champs la structures
champ = fieldnames(param.from.shot.info) ;
nom   = fieldnames(param.from.shot.info) ; 

for k=1:length(nom)
	% nom complet pour acceder a la variable
	champ{k} = strcat('param.from.shot.info.',champ{k}) ;
end
 
test = 0 ;
while (~isempty(champ))
	%premier champ de la liste
	champc=champ{1} ;
	%fprintf(fid,'%s \n',champ{1}) ;
	champ(1)=[] ;
	eval(strcat('test=isstruct(',champc,');')) ;
	if test
	   	% cas d'une sous structure -> ajout de champs
	   	eval(strcat('champnew=fieldnames(',champc,');')) ;
		
	   	for k=1:length(champnew)
			nom2 = champnew{k} ;
	   		% nom complet pour acceder a la variable
	   		champnew{k}=strcat(champc,'.',champnew{k}) ;
		end
					
	   	% ajout a la liste des champs
	   	if isempty(champ)
	   		champ =champnew;
	   	else
	   		champ=cat(1,champ,champnew);
	   	end
	else
		fprintf(fid,'%s : %s \n',champc,eval(champc)) ;
	end
end
         
