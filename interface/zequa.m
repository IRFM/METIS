% ZEQUA valeurs par defaut du formulaire 'configuration equations' sous Edition
%-----------------------------------------------------------------
% fichier : zxec.m 
% 
% fonction Matlab 5 : 
% fonction qui definit les valeurs par defaut du formulaire "execution"
%  
% syntaxe :  
%   d_out = zequa
%
% entrees :  
%  
% sorties :  
%   d_out  = structure contenant les valeurs par defaut
%  
% fonction �rite par C. Passeron , poste 61-19
% version  1.7  du  29/09/2001  
%  
% liste des modifications :  
%	* 11/10/2001 -> ajout de la variable guido
%  
%-------------------------------------------------------------------------------  

function d_out = zequa

% valeurs des donnees
    valeur.nbeq_mode	= evalin('base','param.gene.nbeq_mode') ; 
    valeur.source_bord 	= evalin('base','param.gene.source_bord') ; 
    valeur.ti_invar 	= evalin('base','param.gene.ti_invar') ; 
    valeur.psiequi	= evalin('base','param.gene.psiequi') ;              
    valeur.guido	= evalin('base','param.gene.guido','0') ;
    valeur.qdds		= evalin('base','param.gene.qdds','0') ;
    valeur.xidds	= evalin('base','param.gene.xidds','0') ;
    valeur.adiabatic	= evalin('base','param.gene.adiabatic') ;
    valeur.cn		= evalin('base','param.gene.cn') ;               
    valeur.signe.ip	= evalin('base','param.gene.signe.ip') ;              
    valeur.signe.b0	= evalin('base','param.gene.signe.b0') ;               
    valeur.lambda	= evalin('base','param.gene.lambda') ;              
    valeur.modecoef	= evalin('base','param.gene.modecoef') ;              
    valeur.self		= evalin('base','param.gene.self') ;               
    valeur.fast		= evalin('base','param.gene.fast') ;              
    valeur.force	= evalin('base','param.gene.force') ;              
    valeur.nonneg	= evalin('base','param.gene.nonneg') ;              
    valeur.corrae	= evalin('base','param.gene.corrae') ;              
    valeur.coefplat	= evalin('base','param.gene.coefplat') ;              
    valeur.coefbord	= evalin('base','param.gene.coefbord') ;              
    valeur.evx_inter	= evalin('base','param.gene.evx_inter') ;               
       
% type des donnees
    type.nbeq_mode	= '' ;
    type.source_bord	= '' ;	 
    type.ti_invar	= '' ;	 
    type.psiequi	= '' ;             	
    type.guido		= '' ;
    type.qdds		= '' ;
    type.xidds		= '' ;
    type.adiabatic	= '' ;
    type.cn 		= '' ;
    type.signe.ip	= '' ;             	
    type.signe.b0	= '' ;
    type.lambda		= '' ;
    type.modecoef	= '' ;
    type.self		= '' ;
    type.fast 		= '' ;
    type.force		= '' ;
    type.nonneg 	= '' ;
    type.corrae  	= '' ;
    type.coefplat  	= '' ;
    type.coefbord  	= '' ;
    type.evx_inter 	= '' ;
         
% bornes
    borne.nbeq_mode	= {0,1,2,3,4,5} ;
    borne.source_bord	= {0,1} ;
    borne.ti_invar	= {0,1} ;
    borne.psiequi	= {0,1,2} ;                 	
    borne.guido		= {0,1} ;
    borne.qdds		= [0,3] ;
    borne.xidds		= [0,10] ;
    borne.adiabatic	= {0,1,2} ;
    borne.cn 		= {-1,0,0.3,0.5} ;         
    borne.signe.ip	= {-1,1} ;                 	
    borne.signe.b0	= {-1,1} ;         
    borne.lambda	= {0,3/2,5/2} ;         
    borne.modecoef	= {0,1} ;  
    borne.self		= {-1,0,1} ;   
    borne.fast  	= {0,1,2} ; 
    borne.force		= {0,1,2} ;   
    borne.nonneg	= {0,1} ;   
    borne.corrae	= {0,1} ;   
    borne.coefplat	= {0,3,5,7} ;   
    borne.coefbord	= {0,1} ;   
    borne.evx_inter  	= {0,1} ;    

%valeur par defaut
    defaut.nbeq_mode	= evalin('base','param.gene.nbeq_mode') ; 
    defaut.source_bord 	= evalin('base','param.gene.source_bord') ; 
    defaut.ti_invar 	= evalin('base','param.gene.ti_invar') ; 
    defaut.psiequi	= evalin('base','param.gene.psiequi') ;              
    defaut.guido	= evalin('base','param.gene.guido','0') ;
    defaut.qdds		= evalin('base','param.gene.qdds','0') ;
    defaut.xidds	= evalin('base','param.gene.xidds','0') ;
    defaut.adiabatic	= evalin('base','param.gene.adiabatic') ;
    defaut.cn		= evalin('base','param.gene.cn') ;               
    defaut.signe.ip	= evalin('base','param.gene.signe.ip') ;              
    defaut.signe.b0	= evalin('base','param.gene.signe.b0') ;               
    defaut.lambda	= evalin('base','param.gene.lambda') ;              
    defaut.modecoef	= evalin('base','param.gene.modecoef') ;              
    defaut.self		= evalin('base','param.gene.self') ;               
    defaut.fast		= evalin('base','param.gene.fast') ;              
    defaut.force	= evalin('base','param.gene.force') ;              
    defaut.nonneg	= evalin('base','param.gene.nonneg') ;              
    defaut.corrae	= evalin('base','param.gene.corrae') ;              
    defaut.coefplat	= evalin('base','param.gene.coefplat') ;              
    defaut.coefbord	= evalin('base','param.gene.coefbord') ;              
    defaut.evx_inter	= evalin('base','param.gene.evx_inter') ;               
              
info = zinfo ;

    info.nbeq_mode	= info.param.gene.nbeq_mode ; 
    info.source_bord 	= info.param.gene.source_bord ; 
    info.ti_invar 	= info.param.gene.ti_invar ; 
    info.psiequi	= info.param.gene.psiequi ;
    info.guido		= info.param.gene.guido ;
    info.qdds		= info.param.gene.qdds ;
    info.xidds		= info.param.gene.xidds ;
    info.adiabatic	= info.param.gene.adiabatic ;
    %info.cn		= info.param.gene.cn ; implicite=0, C-N=0.5'
    info.cn		= info.param.gene.cn;
    info.signe.ip	= info.param.gene.signe.ip ;
    info.signe.b0	= info.param.gene.signe.b0 ;
    info.lambda		= info.param.gene.lambda ;
    info.modecoef	= info.param.gene.modecoef ;
    info.self		= info.param.gene.self ;
    info.fast		= info.param.gene.fast ;
    info.force		= info.param.gene.force ;
    info.nonneg		= info.param.gene.nonneg ;
    info.corrae 	= info.param.gene.corrae ;
    info.coefplat 	= info.param.gene.coefplat ;
    info.coefbord 	= info.param.gene.coefbord ;
    info.evx_inter	= info.param.gene.evx_inter ;
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
	
