% ZMULTI valeurs par defaut du formulaire 'multiplicateur' sous Edition
%-----------------------------------------------------------------
% fichier : zmulti.m 
% 
% fonction Matlab 5 : 
%	fonction qui definit les valeurs par defaut du formulaire "multiplicateur"
%  
% syntaxe :  
%	d_out = zmulti
%  
% entrees :  
%  
% sorties :  
%	d_out  = structure contenant les valeurs par defaut
%  
% fonction écrite par C. Passeron , poste 61-19
% version  2.2  du  02/09/2003  
%  
% liste des modifications :  
%	* 02/09/2003 -> ajout des multiplicateur pour la rotation
%  
%----------------------------------------------------------------------
%
function d_out = zmulti

info = zinfo ;

% valeurs des donnees
    valeur.ee = evalin('base','param.cons.neomulti.ee') ;              
    valeur.ve = evalin('base','param.cons.neomulti.ve') ;              
    valeur.ii = evalin('base','param.cons.neomulti.ii') ;               
    valeur.vi = evalin('base','param.cons.neomulti.vi') ;              
    valeur.nn = evalin('base','param.cons.neomulti.nn') ;               
    valeur.vn = evalin('base','param.cons.neomulti.vn') ;               
    valeur.rot = evalin('base','param.cons.neomulti.rot') ;               
    valeur.rotv = evalin('base','param.cons.neomulti.rotv') ;               
       
% type des donnees
    type.ee= 'reel' ;
    type.ve= 'reel' ;
    type.ii= 'reel' ;
    type.vi= 'reel' ;
    type.nn= 'reel' ;
    type.vn= 'reel' ;
    type.rot = 'reel' ;
    type.rotv = 'reel' ;
         
% bornes
    borne.ee = [0,1.e3] ;	    
    borne.ve = [0,1.e3] ;	      
    borne.ii = [0,1.e3] ;	      
    borne.vi = [0,1.e3] ;	      
    borne.nn = [0,1.e3] ; 	
    borne.vn = [0,1.e3]; 	
    borne.rot = [0,1.e3] ; 	
    borne.rotv = [0,1.e3]; 	

%valeur par defaut
    defaut.ee = evalin('base','param.cons.neomulti.ee') ;		  
    defaut.ve = evalin('base','param.cons.neomulti.ve') ;  	    
    defaut.ii = evalin('base','param.cons.neomulti.ii') ;		  
    defaut.vi = evalin('base','param.cons.neomulti.vi') ; 	     
    defaut.nn = evalin('base','param.cons.neomulti.nn') ;		 
    defaut.vn = evalin('base','param.cons.neomulti.vn') ; 	      
    defaut.rot = evalin('base','param.cons.neomulti.rot') ;		 
    defaut.rotv = evalin('base','param.cons.neomulti.rotv') ; 	      

    info.ee = info.param.cons.neomulti.ee ;
    info.ve = info.param.cons.neomulti.ve ;
    info.ii = info.param.cons.neomulti.ii ;
    info.vi = info.param.cons.neomulti.vi ;
    info.nn = info.param.cons.neomulti.nn ;
    info.vn = info.param.cons.neomulti.vn ;
    info.rot = info.param.cons.neomulti.rot ;
    info.rotv = info.param.cons.neomulti.rotv ;
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
	
