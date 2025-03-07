% ZCONNEXION lit les parametres en provenance des modules externes a l'initialisation
%----------------------------------------------------------------------------------------
% fichier zconnexion.m ->  zconnexion
%
%
% fonction Matlab 5 :
%
% Cette fonction lit les parametres en provenance des modules externes lors
% de l'initialisation des structures de donnees. Les modules externes doivent etre
% auto-declarant
%
% syntaxe  :
%  
%     [cr,data,param] = zconnexion(data,param);
%     
% entrees :
%
%     data       =  structure de donnees contenant les variables 
%                   dependant du temps.
%                  
%     param      = structure de donnees contenant les variables 
%                  independantes du temps. 
%                       
% sorties :
% 
%     cr         =  compte rendu d'execution
% 
%     data       =  structure de donnees contenant les variables 
%                   dependant du temps.
%                  
%     param      = structure de donnees contenant les variables 
%                  independantes du temps. 
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 30/08/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [crs,data,param] = zconnexion(data,param)

% initialisation
crs =0;
cr = NaN;


% liste des fonction
liste = fieldnames(param.fonction);
contenu = struct2cell(param.fonction);

while ~isempty(liste)
	fonction = liste{1};
	valeur   = contenu{1};
	if isstruct(valeur)
		liste    = cat(1,liste,strcat(strcat(fonction,'.'),fieldnames(valeur)));
		contenu  = cat(1,contenu,struct2cell(valeur));
	elseif ~isempty(valeur)
		% recherche du parametre de nombre
		if isfield(param.nombre,fonction)
			nb = getfield(param.nombre,fonction);
                        if (nb ==0)|~isfinite(nb)
 	                      sortie.valeur=[];
	                      sortie.type=[];
	                      sortie.borne=[];
	                      sortie.defaut=[];
	                      sortie.info='';
	                      sortie.interface.jet='';
	                      sortie.interface.ts='';
             	              sortie.description = 'vide';   % description (une ligne) de la fonction
	                      sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	                      sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	                      sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	                      sortie.resume ='';                           % nom de la fonction qui ecrit le resume du parametrage
                        else
			     eval(['sortie = ',valeur,'(',num2str(nb),');cr =0;'],'cr =1;');
			end
                        if cr ~=0
				disp('-----------------')
				fprintf('La fonction %s est inconnue ou a produit une erreur :\n%s\n',valeur,lasterr);	
			end
			
		else
			eval(['sortie = ',valeur,';cr =0;'],'cr =1;');
			if cr ~=0
				disp('-----------------')
				fprintf('La fonction %s est inconnue ou a produit une erreur :\n%s\n',valeur,lasterr);	
			end
		end
		if cr ==0
			% recopie des parametres par default
			eval(['param.cons.',fonction,' = sortie.valeur;']);
			if strcmp(fonction,'equi')
				param.gene.nbrhorz   = sortie.nbrho;        % nombre de points en rho de la grille RZ des surfaces magnetiques
				param.gene.nbthetarz = sortie.nbtheta;      % nombre de points en theta de la grille RZ des surfaces magnetiques
				param.gene.nbmoderz = sortie.nbmode;      % nombre de points dans la transformee de Fourrier de la DSMF
				% la grille des surface magnetique
				data.equi.R         = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));       % R des surfaces magentiques (m)
				data.equi.Z         = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));       % Z des surfaces magentiques (m)
				data.equi.BR        = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));       % Br des surfaces magentiques (T)
				data.equi.BZ        = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));       % Bz des surfaces magentiques (T)
				data.equi.BPHI      = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));       % Bphi des surfaces magentiques (T)
				data.equi.rhoRZ     = single(NaN .* ones(param.gene.nbt,sortie.nbrho));                      % rho correspondant aux surfaces magentiques (m)
				data.equi.psiRZ     = single(NaN .* ones(param.gene.nbt,sortie.nbrho)); 
				data.equi.df2RZ     = single(NaN .* ones(param.gene.nbt,sortie.nbrho)); 
				data.equi.dprRZ     = single(NaN .* ones(param.gene.nbt,sortie.nbrho)); 
				data.equi.frmode    = single(NaN .* ones(param.gene.nbt,sortie.nbmode)); 
			elseif strcmp(fonction,'machine')
				data.cons.asser.pfcur = NaN .* ones(param.gene.nbt,sortie.valeur.nbpfcoils)
				param.gene.nbpfcoils  = sortie.nbpfcoils;
			end	
		else	
		   crs = crs+cr;
		end
		
	end
	liste(1)=[];
	contenu(1)=[];
end

if crs ~= 0
	disp('  ')
	disp('----------------------------------------------------');
        disp('Probleme dans les modules externes. Fin du programme'); 
	disp('----------------------------------------------------');
end

    
