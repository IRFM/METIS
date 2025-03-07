function zuiedit_graphique

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_graphique');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Graphics',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'graphique_calc','radio' ,'calculation mode',0,''};
form{length(form)+1} = {col1};

col1 = {'graphique_plot','radio' ,'plot' ,0,''};
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

hout=zuicreeform('Edition','ed_graphique','zuiedit_graphique_fct','',form,comm);
