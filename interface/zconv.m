% ZCONV  valeurs par defaut du formulaire 'convergence' sous Edition
%---------------------------------------------------------------------
% fichier : zconv.m
% 
% fonction Matlab 5 : 
%	fonction qui definit les valeurs par defaut du formulaire 'convergence' 
%	sous Edition
%  
% syntaxe :  
%   d_out = zconv
%  
% entrees :  
%  
% sorties :  
%	d_out  = structure contenant les valeurs par defaut
%  
% fonction M-icrite par C. Passeron , poste 61-19
% version  2.2  du  23/09/2003  
%  
% liste des modifications :  
%  * 11/10/2001 -> ajout de la variable creux
%  * 23/09/2003 -> ajout de factti et nb force
%  
%----------------------------------------------------------------------
%  
function d_out = zconv

info = zinfo ;

% valeurs des donnees
    valeur.delta_adia	= evalin('base','param.gene.delta_adia') ;              
    valeur.nmax		= evalin('base','param.gene.nmax') ;              
    valeur.nequi_ini	= evalin('base','param.gene.nequi_ini') ;               
    valeur.nequi	= evalin('base','param.gene.nequi') ;              
    valeur.dpsi_ini	= evalin('base','param.gene.dpsi_ini') ;               
    valeur.djmoy	= evalin('base','param.gene.djmoy') ;               
    valeur.creux	= evalin('base','param.gene.creux') ;               
    valeur.amorti	= evalin('base','param.gene.amorti') ;               
    valeur.mjmoy	= evalin('base','param.gene.mjmoy') ;               
    valeur.factti	= evalin('base','param.gene.factti') ;               
    valeur.nbforce	= evalin('base','param.gene.nbforce') ;               
    valeur.critere.ne  	= evalin('base','param.gene.critere.ne') ;		
    valeur.critere.pe  	= evalin('base','param.gene.critere.pe') ;		
    valeur.critere.pion	= evalin('base','param.gene.critere.pion') ;		  
    valeur.critere.psi 	= evalin('base','param.gene.critere.psi') ;		 
    valeur.critere.fluce= evalin('base','param.gene.critere.fluce') ;		    
    valeur.critere.flucion= evalin('base','param.gene.critere.flucion') ;		
    valeur.critere.rot	= evalin('base','param.gene.critere.rot') ;		 
       
% type des donnees
    type.delta_adia	= 'reel' ;
    type.nmax		= 'integer' ;
    type.nequi_ini	= 'integer' ;
    type.nequi 		= 'integer' ;
    type.dpsi_ini 	= 'reel' ;
    type.djmoy 		= 'reel' ;
    type.creux 		= 'reel' ;
    type.amorti 	= 'reel' ;
    type.mjmoy 		= 'reel' ;
    type.factti 	= 'reel' ;
    type.nbforce 	= 'integer' ;
    type.critere.ne	= 'reel' ;
    type.critere.pe	= 'reel' ;
    type.critere.pion	= 'reel' ;
    type.critere.psi	= 'reel' ;
    type.critere.fluce	= 'reel' ;
    type.critere.flucion= 'reel' ;
    type.critere.rot	= 'reel' ;  
         
% bornes
    borne.delta_adia	= [1.e-3,1.e-1] ;         
    borne.nmax		= [0,100] ;         
    borne.nequi_ini	= [0,100] ;         
    borne.nequi 	= [0,100] ;         
    borne.dpsi_ini 	= [1.e-4,1] ;         
    borne.djmoy 	= [1.e-4,1] ;         
    borne.creux 	= [0.01,1] ;         
    borne.amorti 	= [0,1] ;         
    borne.mjmoy 	= [0,1] ;     
    borne.factti 	= [1,1000] ;         
    borne.nbforce 	= [0,1000] ;     
    borne.critere.ne	= [1.e-6,1.e-1] ;
    borne.critere.pe	= [1.e-6,1.e-1] ;
    borne.critere.pion	= [1.e-6,1.e-1] ;
    borne.critere.psi	= [1.e-8,1.e-3] ;
    borne.critere.fluce	= [1.e-6,1.e-1] ;
    borne.critere.flucion= [1.e-8,1.e-3] ;
    borne.critere.rot	= [1.e-8,1.e-3] ;           

%valeur par defaut
    defaut.delta_adia	= evalin('base','param.gene.delta_adia') ;              
    defaut.nmax		= evalin('base','param.gene.nmax') ;              
    defaut.nequi_ini	= evalin('base','param.gene.nequi_ini') ;               
    defaut.nequi	= evalin('base','param.gene.nequi') ;              
    defaut.dpsi_ini	= evalin('base','param.gene.dpsi_ini') ;               
    defaut.djmoy	= evalin('base','param.gene.djmoy') ;               
    defaut.creux	= evalin('base','param.gene.creux') ;               
    defaut.amorti	= evalin('base','param.gene.amorti') ;               
    defaut.mjmoy	= evalin('base','param.gene.mjmoy') ;               
    defaut.factti	= evalin('base','param.gene.factti') ;               
    defaut.nbforce	= evalin('base','param.gene.nbforce') ;               
    defaut.critere.ne  	= evalin('base','param.gene.critere.ne') ;		
    defaut.critere.pe  	= evalin('base','param.gene.critere.pe') ;		
    defaut.critere.pion	= evalin('base','param.gene.critere.pion') ;		  
    defaut.critere.psi 	= evalin('base','param.gene.critere.psi') ;		 
    defaut.critere.fluce= evalin('base','param.gene.critere.fluce') ;		    
    defaut.critere.flucion= evalin('base','param.gene.critere.flucion') ;		
    defaut.critere.rot	= evalin('base','param.gene.critere.rot') ;		 

    info.delta_adia	= info.param.gene.delta_adia ;
    info.nmax		= info.param.gene.nmax ;
    info.nequi_ini	= info.param.gene.nequi_ini ;
    info.nequi		= info.param.gene.nequi ;
    info.dpsi_ini	= info.param.gene.dpsi_ini ;
    info.djmoy		= info.param.gene.djmoy ;
    info.creux		= info.param.gene.creux ;
    info.amorti		= info.param.gene.amorti ;
    info.mjmoy		= info.param.gene.mjmoy ;
    info.factti		= info.param.gene.factti ;
    info.nbforce	= info.param.gene.nbforce ;
    info.critere.ne	= info.param.gene.critere.ne ;
    info.critere.pe	= info.param.gene.critere.pe ;
    info.critere.pion	= info.param.gene.critere.pion ;
    info.critere.psi	= info.param.gene.critere.psi ;
    info.critere.fluce  = info.param.gene.critere.fluce ; 
    info.critere.flucion= info.param.gene.critere.flucion ;
    info.critere.rot	= info.param.gene.critere.rot ;     
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
	
