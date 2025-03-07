function zuiedit_param_modext

[hfig,h] = zuiformhandle('ed_graphique_calc');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;
popupmode =      {'    void mode                    ', ...
                  '    input mode                   ', ...
	               '    calculated mode              ', ...
	               '    complex mode                 '} ;
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
%--------------------------------------------------------------------
% le formulaire

nb1   = 10;  % longueur titre module
nb2   = 19;  % longueur nom de la fonction

form = {};
sepa ={'separation_comm','frame','',3,''};                                      form{length(form)+1} = {sepa};
col1 = {'libel_loadfile','text@full','Graphics',[],''};
form{length(form)+1} = {col1};
col1 = {'libel_loadfile','text@full','calculation mode',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% plot
% ----
col1 = {'plot','text','plot',10,info.param.fonction.plot,''} ;
fct = getfield(fonction,'plot') ; if isempty(fct) fct=' '; end
col2 = {'fct_plot','text',fct,12,info.param.fonction.plot,''} ;
col3 = {'chng_plot','radio','modify',0,' modify ...'} ;
col4 = {'para_plot','radio','parameters',0,'parameter function'} ;
val = evalin('base','data.mode.plot') ;
ind = fctval(val) + 1 ;
col5 = {'mode_plot','popup',popupmode,ind,'mode',valeurmode,''} ;
if ind~=4
	col6 = {'edit_plot','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
else
	col6 = {'edit_plot','radio','edition mode',0,'edition mode',[],''} ;
end
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',1,''};
col7 = {'void','frame','void',3,''};
col8 = {'void','frame','void',3,''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8} ;


hout=zuicreeform('Edition','ed_graphique_calc','zuiedit_graphique_calc_fct','',form);

