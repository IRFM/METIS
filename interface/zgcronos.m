% boite de dialogue pour zuigcronos
function [param,data,post,cr] = zgcronos(option)

% initialisation des sorties
param =[];
data =[];
post =[];
cr = -1;

if nargin < 1
 
  	valeur.numchoc  	= 0;
  	valeur.occurrence = 0;
	
  	type.numchoc  	   = 'integer';
  	type.occurrence  	= 'integer';
	
  	borne.numchoc  	= [1500,inf];
  	borne.occurrence 	= [0,9];

  	defaut.numchoc  	= 0;
  	defaut.occurrence = 0;

  	info.numchoc    	= 'Numero du choc TS';
  	info.occurrence   = 'Numero d''occurence de tprof';

	interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	param.valeur     = valeur;
	param.type       = type;
	param.borne      = borne;
	param.defaut     = defaut;
	param.info       = info;
	param.interface  = interface;
	
	param.description = 'Appel de gcronos pour la creation d''un jeu de donnees TS';   % description (une ligne) de la fonction
	
	param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	param.gui      = '';                            % nom de l'interface graphique specifique si elle existe
	param.controle = '';                            % nom de la fonction de controle des valeurs si elle existe
	
	return
end

numchoc = option.numchoc;
occurrence = option.occurrence;
evalin('base','cr = -314;');
[param,data,post,cr] = gcronos(numchoc,occurrence,'write');
