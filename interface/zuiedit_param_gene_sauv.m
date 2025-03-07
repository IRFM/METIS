% ZUIEDIT_PARAM_GENE_SAUV   formulaire 'informations' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_gene_sauv.m ->  
%	zuicreeform : creation du formulaire
%	zuiedit_param_gene_sauv_crtl	: fonction de control des donnees
%
% fonction Matlab 5 :
%	fonction de creation du formulaire 
%	de sauvegarde , commandes de parametres  
% 	sous le mode edition du formulaire principal
%
% syntaxe  :
%	zuiedit_param_gene_sauv;
%
% entrees
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------

function zuiedit_param_gene_sauv

[hfig,h] = zuiformhandle('ed_param_sauv');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;

%--------------------------------------------------------------------
% le formulaire

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Global simulation parameters ',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','Saving',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Informations : param.plot  
% -------------------------tous les tegroupemps
colj = {'void','jump','void',[],''} ;

col1 = {'libel','text','file',8,info.param.gene.file,''};
file = evalin('base','param.gene.file','[]') ;
col2 = {'edit_file','edit',file,30,info.param.gene.file,[],'param.gene.file'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','rapsauve',8,info.param.gene.rapsauve,''};
file = evalin('base','param.gene.rapsauve','[]') ;
col2 = {'edit_rapsauve','edit',file,30,info.param.gene.rapsauve,[],'param.gene.rapsauve'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','nbsauve',8,info.param.gene.nbsauve,''};
val  = evalin('base','param.gene.nbsauve') ;
liste = {'at the end','each time','1/2','1/5','1/10'} ;
borne = [0,1,2,5,10] ;
ind = findstr(borne,val) ;
col2 = {'pop_nbsauve','popup',liste,ind,info.param.gene.nbsauve,borne,'param.gene.nbsauve'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','rebuilt',8,info.param.gene.rebuilt,''};
val = evalin('base','param.gene.rebuilt','[]') ;
col2 = {'pop_rebuilt','popup',{'no','yes'},val+1,info.param.gene.rebuilt,{0,1},'param.gene.rebuilt'};
form{length(form)+1} = {colj,col1,col2};

%col1 = {'libel','text','post-treatment',8,info.param.gene.post,''};
%val = evalin('base','param.gene.post','[]') ;
%col2 = {'pop_psot','popup',{'no','yes'},val+1,info.param.gene.post,{0,1},'param.gene.post'};
%form{length(form)+1} = {colj,col1,col2};


hout=zuicreeform('Edition','ed_param_sauv','zuiedit_param_gene_sauv_fct','zuiedit_param_gene_sauv_ctrl',form);
