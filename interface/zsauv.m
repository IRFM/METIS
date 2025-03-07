% ZSAUV valeurs par defaut du formulaire 'sauvegarde' sous Edition
%------------------------------------------------------------------------------- 
% fichier : zsauv.m 
% 
% fonction Matlab 5 : 
% 
% fonction qui definit les valeurs par defaut du formulaire "sauvegarde"
%  
% syntaxe :  
%   d_out = zsauv
%  
% entrees :  
%  
% sorties :  
%   d_out  = structure contenant les valeurs par defaut
%  
% fonction �rite par xxxxxxx , poste XX-XX  
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
function d_out = zsauv

% valeurs des donnees
    valeur.file     = evalin('base','param.gene.file') ;   
    valeur.rapsauve = evalin('base','param.gene.rapsauve') ;   
    valeur.nbsauve  = evalin('base','param.gene.nbsauve');               
       
% type des donnees
    type.file		= '' ;             	
    type.rapsauve	= '' ;
    type.nbsauve	= '' ;
         
% bornes
    borne.file		= '' ;                 	
    borne.rapsauve	= '' ;         
    borne.nbsauve	= {'at the end','every time ','1/2','1/5','1/10'} ;

%valeur par defaut
    defaut.file		= evalin('base','param.gene.file') ;                
    defaut.rapsauve	= evalin('base','param.gene.rapsauve');  
    defaut.nbsauve	= evalin('base','param.gene.nbsauve') ;
       
info = zinfo ;

    info.file		= info.param.gene.file ;
    info.rapsauve	= info.param.gene.rapsauve ;
    info.nbsauve	= strcat(info.param.gene.nbsauve,'1=all the time, 2=1/2, 5=1/5, 10=1/10') ;
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
		
    sortie.help = '' ;			% nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='' ;			% nom de l'interface graphique specifique si elle existe
    sortie.controle = '' ;		% nom de la fonction de controle des valeurs si elle existe
	
    d_out =sortie ;
	
