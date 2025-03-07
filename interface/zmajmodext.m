% ZMAJMODEXT   connecte un module externe dans l'interface
%----------------------------------------------------------------------------------------
% fichier zmajmodext.m 
%
% fonction Matlab 5 :
%	Cette fonction connecte un module extrene dans l'interfacel.
%	Les modules externes doivent etre auto-declarant
%
% syntaxe  :
%	zmajmodext(fonction,module);
%     
% entrees :
%	fonction = nom de la fonction (param.fonction.xxx)
%	module   = nom du module externe
% 
% sorties
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.7, du 19/10/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zmajmodext(fonction,module)

% initialisation
crs =0;
cr = NaN;

param =evalin('base','param') ;
%data = evalin('base','data') ;

% type de la fonction (param.fonction.xxx)
%fonction   = get(getappdata(hfig,'hfct'),'String') ;
% nom du module externe
%module = get(getappdata(hfig,'hmodule'),'String') ;

if isempty(fonction)
   error('type de la fonction vide');
end
   

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
		sortie.help = '';			     % nom du fichier d'aide s'il existe, sinon aide de la fonction
		sortie.gui  ='';			     % nom de l'interface graphique specifique si elle existe
		sortie.controle = '';			     % nom de la fonction de controle des valeurs si elle existe
		sortie.resume ='';			     % nom de la fonction qui ecrit le resume du parametrage
	else
		eval(['sortie = ',module,'(',num2str(nb),');cr =0;'],'cr =1;');
	end
        if cr ~=0
		warndlg(sprintf('La fonction %s est inconnue ou a produit une erreur :\n%s\n',module,lasterr), ...
		        'Erreur de connexion a un module externe');
	        return	
	end
			
else
	eval(['sortie = ',module,';cr =0;'],'cr =1;');
	if cr ~=0
		warndlg(sprintf('La fonction %s est inconnue ou a produit une erreur :\n%s\n',module,lasterr), ...
		        'Erreur de connexion a un module externe');
	        return	
	end
end
% sauvegarde les anciennes valeurs
evalin('base',['memcons.',fonction,' = param.cons.',fonction,';']);
% recopie des parametres par default
eval(['param.cons.',fonction,' = sortie.valeur;']);
if strcmp(fonction,'equi')
	param.gene.nbrhorz   = sortie.nbrho;        % nombre de points en rho de la grille RZ des surfaces magnetiques
        param.gene.nbthetarz = sortie.nbtheta;      % nombre de points en theta de la grille RZ des surfaces magnetiques
        % la grille des surface magnetique
	m3 = single(NaN .* ones(param.gene.nbt,sortie.nbrho,sortie.nbtheta));
	m2 = single(NaN .* ones(param.gene.nbt,sortie.nbrho));
        zassignin('base','data.equi.R',m3);        % R des surfaces magentiques (m)
        zassignin('base','data.equi.Z',m3);        % Z des surfaces magentiques (m)
        zassignin('base','data.equi.BR',m3);       % Br des surfaces magentiques (T)
        zassignin('base','data.equi.BZ',m3);       % Bz des surfaces magentiques (T)
        zassignin('base','data.equi.BPHI',m3);     % Bphi des surfaces magentiques (T)
        zassignin('base','data.equi.rhoRZ',m2);    % rho correspondant aux surfaces magentiques (m)
end
if strcmp(fonction,'machine')
		m11 = NaN .* ones(param.gene.nbt,sortie.valeur.nbpfcoils);
		param.gene.nbpfcoils  = sortie.valeur.nbpfcoils;
		zassignin('base','data.cons.asser.pfcur',m11);        % declaration du nombre de pfcoils
end

zassignin('base','param',param);

