function zuiedit_transport_coeff_calc

[hfig,h] = zuiformhandle('ed_transport_coeff_calc');
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
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'libel_loadfile','text@full','Transport coefficient model',[],''};
form{length(form)+1} = {col1};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

%  coefficients de transport
% --------------------------
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'transport','text@full','transport',1,'transport equation coefficient',''};
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4} ;

% coefa
% -----
col1 = {'coefa','text','coefa',10,info.param.fonction.coefa,''} ;
fct = getfield(fonction,'coefa') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefa','text',fct,12,info.param.fonction.coefa,''} ;
col3 = {'chng_coefa','radio','modify',0,' modify ...'} ;
col4 = {'para_coefa','radio','parameters',0,'parameter function'} ;
%val = evalin('base','data.mode.coefa','[]') ;
%if ~isempty(val)
%	ind = val(1)+1 ;
%	col5 = {'mode_coefa','popup',popupmode,ind,'mode',valeurmode,'data.mode.coefa'} ;
%else 
%	col5 = {'mode_coefa','popup',popupmode,1,'mode',valeurmode,''} ;
%end
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefa','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% coefb
% -----
col1 = {'coefb','text','coefb',10,info.param.fonction.coefb,''} ;
fct = getfield(fonction,'coefb') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefb','text',fct,12,info.param.fonction.coefb,''} ;
col3 = {'chng_coefb','radio','modify',0,' modify ...'} ;
col4 = {'para_coefb','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefb','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% coefc
% -----
col1 = {'coefc','text','coefc',10,info.param.fonction.coefc,''} ;
fct = getfield(fonction,'coefc') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefc','text',fct,12,info.param.fonction.coefc,''} ;
col3 = {'chng_coefc','radio','modify',0,' modify ...'} ;
col4 = {'para_coefc','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefc','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% coefd
% -----
col1 = {'coefd','text','coefd',10,info.param.fonction.coefd,''} ;
fct = getfield(fonction,'coefd') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefd','text',fct,12,info.param.fonction.coefd,''} ;
col3 = {'chng_coefd','radio','modify',0,' modify ...'} ;
col4 = {'para_coefd','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefd','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% coefe
% -----
col1 = {'coefe','text','coefe',10,info.param.fonction.coefe,''} ;
fct = getfield(fonction,'coefe') ; if isempty(fct) fct=' '; end
col2 = {'fct_coefe','text',fct,12,info.param.fonction.coefe,''} ;
col3 = {'chng_coefe','radio','modify',0,' modify ...'} ;
col4 = {'para_coefe','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coefe','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% coeff
% -----
col1 = {'coeff','text','coeff',10,info.param.fonction.coeff,''} ;
fct = getfield(fonction,'coeff') ; if isempty(fct) fct=' '; end
col2 = {'fct_coeff','text',fct,12,info.param.fonction.coeff,''} ;
col3 = {'chng_coeff','radio','modify',0,' modify ...'} ;
col4 = {'para_coeff','radio','parameters',0,'parameter function'} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
%col6 = {'edit_coeff','radio','edition mode',0,'edition mode',[],'','Enable','off'} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4} ;

hout=zuicreeform('Edition','ed_transport_coeff_calc','zuiedit_transport_coeff_calc_fct','',form);

