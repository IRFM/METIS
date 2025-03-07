% ZPROFILE valeurs par defaut du formulaire 'configuration profiler' sous Edition
%-----------------------------------------------------------------

function d_out = zprofile

% valeurs des donnees
    valeur.varlog1    = evalin('base','param.profile.varlog1') ;              
	 if isempty(valeur.varlog1) valeur.varlog1 = '                       ' ;end
    valeur.varlog2    = evalin('base','param.profile.varlog2') ;              
	 if isempty(valeur.varlog2) valeur.varlog2 = '                       ' ;end
    valeur.varlog3    = evalin('base','param.profile.varlog3') ;              
	 if isempty(valeur.varlog3) valeur.varlog3 = '                       ' ;end
    valeur.onoff	    = evalin('base','param.profile.onoff') ; 
    valeur.rapport 	 = evalin('base','param.profile.rapport') ; 
    valeur.memoire	 = evalin('base','param.profile.memoire') ;              
    valeur.erreurmode = evalin('base','param.profile.erreurmode') ;              
       
% type des donnees
    type.varlog1    = 'texte' ;           	
    type.varlog2    = 'texte' ;             	
    type.varlog3    = 'texte' ;             	
    type.onoff	     = 'integer' ;
    type.rapport	  = 'integer' ;	 
    type.memoire	  = 'integer' ;             	
    type.erreurmode = 'integer' ;             	
         
% bornes
    borne.varlog1   	= '   ' ;        
    borne.varlog2   	= '   ' ;         
    borne.varlog3   	= '   ' ;         
    borne.onoff	   = {0,1} ;
    borne.rapport   	= {0,1} ;
    borne.memoire	   = {0,1} ;         
    borne.erreurmode	= {0,1} ;         

%valeur par defaut
    defaut = valeur;
              
    info = zinfo ;

    info.onoff	     = info.param.profile.onoff ; 
    info.rapport 	  = info.param.profile.rapport ; 
    info.memoire	  = info.param.profile.memoire ;
    info.erreurmode = info.param.profile.erreurmode ;
    info.varlog1    = info.param.profile.varlog1 ;
    info.varlog2    = info.param.profile.varlog2 ;
    info.varlog3    = info.param.profile.varlog3 ;
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
	
