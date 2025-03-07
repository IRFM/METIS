function zuiedit_param_postrait_param_calc

[hfig,h] = zuiformhandle('ed_postrait_calc');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;
popupmode =      {'    off mode                    ', ...
                  '    prescribed mode             ', ...
	          '    calculated mode             ', ...
	          '    complex mode                '} ;
valeurmode = {0,1,2,3} ;
% pour mhd.stabilite
popupmode_stab = {'    void mode                    ', ...
                  '    input mode                   ', ...
	               '    calculated mode              ', ...
	               '    complex mode                 ', ...
				      '    post-processing mode        '} ;
valeurmode_stab = {0,1,2,3,4} ;

% pour la mhd
popupmode_m = {' no calculation        ', ...
               ' calculated mode       ', ...
 	            ' complex mode          '} ;
valeurmode_m = {0,1,2} ;

% pour la neo
popupmode_neo = {' input mode             ', ...
                 ' calculated mode        ', ...
                 ' complex mode           '} ;
valeurmode_neo= {1,2,3} ;
 

% pour les glacons
popupmode_g = {' no pellet         ', ...
               ' pellet            ', ...
 	            ' complex mode      '} ;
valeurmode_g = {0,1,2} ;


fonction = evalin('base','param.fonction') ;

% post
% ----
form={};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};
col1 = {'libel','text@full','Post-treatment',[],''};
form{length(form)+1} = {col1};
col1 = {'libel','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};
form{length(form)+1} = {sepa};

col1 = {'post','text','post',10,info.param.fonction.post,''} ;
fct = getfield(fonction,'post') ; if isempty(fct) fct=' '; end
col2 = {'fct_post','text',fct,12,info.param.fonction.post,''} ;
col3 = {'chng_post','radio','modify',0,' modify ...'} ;
col4 = {'para_post','radio','parameters',0,'edit function parameters'} ;
val = evalin('base','data.mode.post') ;
ind = fctval(val)+1;
col5 = {'mode_post','popup',popupmode_m,min(ind,3),'mode',valeurmode_m,''} ;
if ind~=4
	col6 = {'edit_post','radio','edit  mode',0,'edit compute mode',[],'','Enable','off'} ;
else
	col6 = {'edit_post','radio','edit mode',0,'edit compute mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;


hout=zuicreeform('Edition','ed_postrait_calc','zuiedit_param_postrait_calc_fct','',form);
