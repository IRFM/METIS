function zuiedit_neo

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_neo');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Neoclassical',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'neo_calcmode','radio' ,'calculation mode',0,''};
form{length(form)+1} = {col1};

col1 = {'radio_jboot','radio','jboot',0,info.data.mode.jboot,[],''} ;
form{length(form)+1} = {col1};

col1 = {'radio_eta'  ,'radio','eta'  ,0,info.data.mode.eta,[],''} ;
form{length(form)+1} = {col1};

col1 = {'radio_qei'   ,'radio','qei'   ,0,info.data.mode.qei,[],''} ;
form{length(form)+1} = {col1};

col1 = {'radio_qneo'  ,'radio','qneo'  ,0,info.data.mode.qneo,[],''} ;
form{length(form)+1} = {col1};


%--------------------------------------------------------------------
% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Close',0,''};

hout=zuicreeform('Edition','ed_neo','zuiedit_neo_fct','',form,comm);
