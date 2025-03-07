% ZIMPORT importe une donneee dans la structure data de cronos
%--------------------------------------------------------------
% fichier zimport.m ->  zimport
%
%
% fonction Matlab 5 :
% 
% Cette fonction importe une donneee dans la structure data de cronos. 
% Elle permet d'importer soit une consigne soit un profil
% La source des donnees peut etre :
%    * le workspace matlab
%    * un fichier matlab
%    * un fichier ascii
%    * la base de donnees arcad
%    
% La structure des fichier ascii est :
%    * pour les consignes, un tableau, dont la premiere colonne est le temps et un des autres colonne
%      la donnee 
%    
%    * pour les profils, un tableau, dont la premiere colonne a partir de la deuxieme ligne est le temps, 
%      dont la premiere ligne a partir de la deuxieme colonne est la coordonnee d'espace et et le reste du tableau
%      sauf l'element (1,1) est la donnee. Chaque ligne representant un temps.
%      
% Dans cette version, la coordonnee d'espace est obligatoirement sqrt(phi/pi/B0) ou sa valeur normalisee a 1.    
% 
%
% syntaxe  :
%  
%   cr = zimport(nom,mode,source,filename,type,nom_temps,nom_espace,nom_data,coordonnee,valeur_defaut,positif)
%  
% entrees :
%  
%   nom           =   nom de la donnee a modifier dans le workspace (ex: data.cons.hyb)
%   mode          =   format de la donnee : 0 -> consigne, 1 -> profil
%   source        =   source des donnees a importer: 0 -> workspace, 1 -> fichier, -1 -> arcad
%   filename      =   chemin et nom complet du fichier (comprimer ou non) ou <objet>@<numchoc>.<occurence> pour les donnees
%                     lues dans la base (ex : tprof@28011.1).
%   type          =   type du fichier: 0 -> matlab, 1-> ascii
%   nom_temps     =   nom de la variable de temps dans le fichier ou provenant de la base (ex: times) [times]
%   nom_espace    =   pour les consignes, numero de la voie de la donnee importe.
%                     pour les profils, nom de la variable d'espace ans le fichier ou provenant de la base (ex: rhofit) [rhofit]
%   nom_data      =   nom de la veriable donnee dans le fichier ou provenant de la base (amin pour samin)
%   coordonnee    =   usage reserve pour plustard
%   valeur_defaut =   valeur prise par la donnee, apres reechantillonage, hors de l'intervalle de temps ou la donnee importee est definie sur t
%   positif       =   1 -> force la donnee reechntillone a etre positive si positif, 0 -> rien
%   
% sortie
% 
%   cr = comterendu d'execution, 0 si ok 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 17/01/2005.
% 
% 
% liste des modifications : 
%
%   * 01/10/2001 -> securite anti espace dans les noms
%   * 17/10/2002 -> support specifique pour Te,Pe, Ti et pion
%   * 03/06/2003 -> correction bug fin d'intervalle de temps
%   * 18/08/2004 -> charge les donnees en provenance de la base meme si nom_data est vide
%   * 17/01/2005 -> remplace rm par rm -f
%   * 21/07/2005 -> bug ligne 279
%
%--------------------------------------------------------------
%
function cr = zimport(nom,mode,source,filename,type,nom_temps,nom_espace,nom_data,coordonnee,valeur_defaut,positif)

% test des entrees
% juste le nombre d'argument
if nargin < 11
	disp('Nombre d''argument incorrect ...')
	cr = -1;
end
voie = [];
cr   = 0;

nom_temps  = strrep(nom_temps,' ','');

nom_espace = strrep(nom_espace,' ','');

nom_data   = strrep(nom_data,' ','');

% chargement des donnees
if source == 0
	% source matlab
	% temps
	try
	    temps = evalin('base',nom_temps);
	catch
	    errordlg(sprintf('Impossible d''importer la variable de temps (%s) depuis le workspace', ...
	                      nom_temps),'Erreur d''importation depuis le workspace','non-modal');
	    disp('Probleme dans zimport :')
	    disp(lasterr)  
	    cr = 1;
	    return
	end
	% espace   
	if mode  == 1   
		try
		    espace = evalin('base',nom_espace);
		catch
		    errordlg(sprintf('Impossible d''importer la variable d''espace (%s) depuis le workspace', ...
	                     nom_espace),'Erreur d''importation depuis le workspace','non-modal');
	            disp('Probleme dans zimport :')
	            disp(lasterr)                  
	            cr = 2;
	            return
	        end
	else
		espace =[];
	end
	% donnees
	try
	    data = evalin('base',nom_data);
	catch
	    errordlg(sprintf('Impossible d''importer la donnee (%s) depuis le workspace', ...
	             nom_data),'Erreur d''importation depuis le workspace','non-modal');
	    disp('Probleme dans zimport :')
	    disp(lasterr)                  
	    cr = 3;
	    return
	end
	if mode ==0   
	    	if ischar(nom_espace)
	    		voie = str2num(nom_espace);
	    	else
	    		voie = nom_espace;
	    	end
	    	data  = data(:,voie);
	end
elseif source == -1 
	% acces a la base
	try 
	    % separation des champs
	    [producteur,numchoc] = strtok(filename,'@');
	    [donnees,sout]       = cgcgettrait(str2num(numchoc(2:end)),producteur);
	    if isempty(nom_temps)
	    	   if isempty(nom_data)
                      nom_data = producteur(2:end);
                   end
                   donnees = getfield(sout,nom_data);
	    	   nom_temps = 'temps';
	    	   nom_data  = 'data';
	    	   nom_espace = 'espace';
	    end
	    if mode  == 1
	    	% mode profil
	    	temps   = getfield(donnees,nom_temps);
	    	espace  = getfield(donnees,nom_espace);
	    	data    = getfield(donnees,nom_data);
                clear donnees
	    else
	   	% mode consigne
	    	temps   = getfield(donnees,nom_temps);
	    	espace  = [];
	    	data    = getfield(donnees,nom_data);
	    	if ischar(nom_espace)
	    		voie = str2num(nom_espace);
	    	else
	    		voie = nom_espace;
	    	end
	    	data  = data(:,voie);
                clear donnees
	    end
	catch
	    errordlg(sprintf('Impossible d''importer la donnee (%s) \n depuis ce producteur (%s)', ...
	             nom_data,filename),'Erreur d''importation depuis la base de donnees','non-modal');
	    disp('Probleme dans zimport :')
	    disp(lasterr)                  
	    cr = 4;
	    return
        end

elseif type == 0
	% lecture dans un fichier matlab
	rmfile ='';
	try
	    % decompression si se finit par .gz
	    hout =warndlg('Chargement en cours ...',sprintf('Importation de %s',nom));
	    delete(findobj(hout,'string','OK'));
	    pause(1)
	    [lfile,rmfile,cr]=zgzip(filename,'uncompress');
	
	    if mode  == 1
	    	% mode profil
	    	donnees = load(lfile);
	    	temps   = getfield(donnees,nom_temps);
	    	espace  = getfield(donnees,nom_espace);
	    	data    = getfield(donnees,nom_data);
                clear donnees
	    else
		% mode consigne
		donnees = load(lfile);
	    	temps   = getfield(donnees,nom_temps);
	    	espace  = [];
	    	data    = getfield(donnees,nom_data);
	    	if ischar(nom_espace)
	    		voie = str2num(nom_espace);
	    	else
	    		voie = nom_espace;
	    	end
	    	data  = data(:,voie);
                clear donnees
	    end
	catch
	    errordlg(sprintf('Impossible d''importer la donnee (%s) \n depuis le fichier(%s)', ...
	             nom_data,filename),'Erreur d''importation depuis un fichier matlab','non-modal');
	    disp('Probleme dans zimport :')
	    disp(lasterr)                  
	    delete(hout);
	    cr = 5;
            return
        end
        % supression du fichier temporaire s'il existe
        if ~isempty(rmfile)
		[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
	end
        delete(hout);
else
	% lecture dans un fichier ascii
	try
	    donnees = load(filename,'-ascii');
	
	    if mode  == 1
	    	% mode profil
	    	temps   = donnees(2:end,1);
	    	espace  = donnees(1,2:end);
	    	data    = donnees(2:end,2:end);

	    else
		% mode consigne
	    	temps   = donnees(:,1);
	    	espace  = [];
	    	if ischar(nom_espace)
	    		voie = str2num(nom_espace);
	    	else
	    		voie = nom_espace;
	    	end
	    	data  = donnees(:,voie+1);
	    end
	    clear donnees
	catch
	    errordlg(sprintf('Impossible d''importer la donnee (%s) \n depuis le fichier(%s)', ...
	             nom_data,filename),'Erreur d''importation depuis un fichier matlab','non-modal');
	    disp('Probleme dans zimport :')
	    disp(lasterr)                  
	    cr = 6;
	    return
        end
	
end


% reechantillonage + importation dans le workspace
if mode == 1
	% recupeartion de la parametrisation du reechantillonage
	try 
	   groupe = evalin('base','param.from.sample.groupe');
       catch
           disp('info parametrage zsample indisponible => utiliation des valeurs par defauts');
           groupe.ondelette         = 0;
           groupe.energie           = 1;
           groupe.defaut.temps     = NaN;
           groupe.defaut.espace    = 0;
           groupe.defaut.inf       = [];
           groupe.plus             = 0;
       end
       groupe.plus =positif;
       if ~isempty(valeur_defaut)
       	   groupe.defaut.temps     = valeur_defaut;
       end
       % lecture de x et  de t dans le workspace
       tt = evalin('base','data.gene.temps');
       xx = evalin('base','param.gene.x');
       
       % il n'y a pas de choix de la coordonnee pour le moment
       % la coordonee d'entree est normalisee a 1
       if size(espace,1) >1
          	espace = espace ./ (max(espace,[],2) * ones(1,size(espace,2)));
       else
                espace = espace ./ max(espace);	
       end
       
       % suppression des NaN
       ind =find(~isfinite(data));
       if ~isempty(ind)
	   disp('Presence de NaN ou Inf dans le signal importe');
	   data(ind) = valeur_defaut;    
       end
       
       % securite bornes
       tt = ztsecure(tt,temps);
       % reechantillonage
       data  = zsample(data,temps,espace,tt,xx,groupe);
       
       % mise a jour du workspace
       zassignin('base',nom,data);
       zuisavenonok;
       
else
	% recupeartion de la parametrisation du reechantillonage
	try 
	   signal = evalin('base','param.from.sample.groupe');
       catch
           disp('info parametrage zsample indisponible => utiliation des valeurs par defauts');
           signal.ondelette        = 0;  
           signal.defaut.temps     = NaN;
           signal.defaut.espace    = 0;
           signal.defaut.inf       = [];
           signal.plus             = 0;
       end
       signal.plus =positif;
       if ~isempty(valeur_defaut)
       	   signal.defaut.temps     = valeur_defaut;
       end
       % lecture de t dans le workspace
       tt = evalin('base','data.gene.temps');
       
       % securite bornes
       tt = ztsecure(tt,temps);
       
       % reechantillonage
       data  = zsample(data,temps,tt,signal);
       
       % mise a jour du workspace
       data_new = evalin('base',nom);
       
       % test de dimension
       if voie > size(data_new,2)
       	     errordlg(sprintf('la donnee a %d voies dans le workspac, or, la voie choisie est la %d', ...
	             size(data_new,2),voie),'Erreur d''importation -> #voie non valide','non-modal');
	     cr = 10;
             return
       end	
       data_new(:,voie) = data;
       zassignin('base',nom,data_new);
       zuisavenonok;

end

switch nom
case 'data.prof.te'
	evalin('base','param.edit.tepe =''te'';');
	evalin('base','zrecalcultemppres(''el'');');
case 'data.prof.pe'
	evalin('base','param.edit.tepe =''pe'';');
	evalin('base','zrecalcultemppres(''el'');');
case 'data.prof.ti'
	evalin('base','param.edit.tipion =''ti'';');
	evalin('base','zrecalcultemppres(''ion'');');
case 'data.prof.pion'
	evalin('base','param.edit.tipion =''pion'';');
	evalin('base','zrecalcultemppres(''ion'');');
end

% securite si borne ==	
function tt = ztsecure(tt,temps)

   if tt(1) == min(temps(:))
      tt(1) = tt(1) + 1e-6;
   end
   if tt(end) == max(temps(:))
      tt(end) = tt(end) - 1e-6;
   end
