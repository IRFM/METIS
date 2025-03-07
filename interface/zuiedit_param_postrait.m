function zuiedit_param_postrait

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_postrait');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Post-treatment',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'postrait_param','radio' ,'calculation parameter' ,0,''};
form{length(form)+1} = {col1};

col1 = {'postrait_calc','radio' ,'calculation mode',0,''};
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

hout=zuicreeform('Edition','ed_postrait','zuiedit_param_postrait_fct','',form,comm);
