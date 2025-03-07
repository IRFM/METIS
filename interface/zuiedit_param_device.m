function zuiedit_param_device

[hfig,h] = zuiformhandle('ed_param_device');
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

fonction = evalin('base','param.fonction') ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'libel_loadfile','text@full','Device',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% machine
% ----
col1 = {'machine','text','device',10,info.param.fonction.machine,''} ;
fct = getfield(fonction,'machine') ; if isempty(fct) fct=' '; end
col2 = {'fct_machine','text',fct,12,info.param.fonction.machine,''} ;
col3 = {'chng_machine','radio','modify',0,' changer ...'} ;
col4 = {'para_machine','radio','parameters',0,'edit function parameters '} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

hout=zuicreeform('Edition','ed_param_device','zuiedit_param_device_fct','',form);

