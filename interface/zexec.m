% ZEXEC valeurs par defaut du formulaire 'execution' sous Edition
%------------------------------------------------------------------------------- 
% fichier : zxec.m 
% 
% fonction Matlab 5 : 
% 
% fonction qui definit les valeurs par defaut du formulaire "execution"
%  
% syntaxe :  
%   d_out = zexec
%  
% entrees :  
%  
% sorties :  
%   d_out  = structure contenant les valeurs par defaut
%  
% fonction écrite par xxxxxxx , poste XX-XX  
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function d_out = zexec

% valeurs des donnees
    valeur.tdeb		= evalin('base','param.gene.tdeb') ;              
    valeur.tfin		= evalin('base','param.gene.tfin') ;               
    valeur.t		= evalin('base','param.gene.t') ;              
    valeur.kmin		= evalin('base','param.gene.kmin') ;		 
    valeur.kmax		= evalin('base','param.gene.kmax') ;		  
    valeur.k	       	= evalin('base','param.gene.k') ;	      
    valeur.verbose     	= evalin('base','param.gene.verbose') ;  	     
       
% type des donnees
    type.tdeb		= 'reel' ;             	
    type.tfin		= 'reel' ;
    type.t		= 'reel' ;
    type.kmin		= 'integer' ;
    type.kmax	       	= 'integer' ;
    type.k	       	= 'integer' ;
    type.verbose       	= '' ;
         
% bornes
    tdeb = min(evalin('base','data.gene.temps')) ;
    tfin = max(evalin('base','data.gene.temps'));
    nbt  =length(evalin('base','data.gene.temps'));
    borne.tdeb		= [tdeb,tfin];                 	
    borne.tfin		= [tdeb,tfin] ;         
    borne.t		= [tdeb,tfin] ;         
    borne.kmin         	= [1, nbt-1] ;  
    borne.kmax         	= [2, nbt] ;	
    borne.k	       	= [1, nbt] ; 
    borne.verbose      	= {0,1} ;    

%valeur par defaut
    defaut.tdeb		= evalin('base','param.gene.tdeb') ;                  
    defaut.tfin		= evalin('base','param.gene.tfin') ;
    defaut.t		= evalin('base','param.gene.t') ; 
    defaut.kmin		= evalin('base','param.gene.kmin') ;
    defaut.kmax		= evalin('base','param.gene.kmax') ; 
    defaut.k	       	= evalin('base','param.gene.k') ;    
    defaut.verbose     	= evalin('base','param.gene.verbose') ;  
       
info = zinfo ;

    info.tdeb		= info.param.gene.tdeb ;
    info.tfin		= info.param.gene.tfin ;
    info.t		= info.param.gene.t;
    info.kmin	       	= info.param.gene.kmin ;
    info.kmax	       	= info.param.gene.kmax ;
    info.k	       	= info.param.gene.k ;
    info.verbose       	= info.param.gene.verbose ;
	%
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
