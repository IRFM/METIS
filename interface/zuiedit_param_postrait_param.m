function zuiedit_param_postrait_param

[hfig,h] = zuiformhandle('ed_postrait_param');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;

%--------------------------------------------------------------------
% le formulaire
colj = {'void','jump','void',[],''} ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'libel_loadfile','text@full','Post-treatment',[],''};
form{length(form)+1} = {col1};
col1 = {'libel_loadfile','text@full','calculation parameter',[],''};
form{length(form)+1} = {col1};
form{length(form)+1} = {sepa};

col1 = {'libel','text','post-treatment',8,info.param.gene.post,''};
val = evalin('base','param.gene.post','[]') ;
col2 = {'pop_psot','popup',{'no','yes'},val+1,info.param.gene.post,{0,1},'param.gene.post'};
form{length(form)+1} = {colj,col1,col2};


hout=zuicreeform('Edition','ed_postrait_param','zuiedit_param_postrait_param_fct','',form);

