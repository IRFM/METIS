function hout=zuiedit_edit_pellet

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('edit_refval_pellet') ;
if ishandle(hout)
        zuiformvisible(hout) ;
        return
end

% infos pour les tooltips
info = zinfo ;

% formulaire avec la sous structure from
form={};

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text@full','Reference values for pellet',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};



% data.cons.glacon
% ----------------
nb = evalin('base','param.nombre.glacon') ;
for i=1:nb
        if i==1
                col1 = {'pellet','text','pellet',10,info.data.cons.glacon,''} ;
        else
                col1 = {'pellet','text',' '  ,10,'',''} ;
        end
        st   = sprintf('channel %s',num2str(i)) ;
        col2 = {'num' ,'text',st,10,'number ?',''} ;
        st   = strcat('edit_module_glacon',num2str(i)) ;
        col3 = {st,'radio','Edit',0,'Edit the initial value',[],''} ;
        st   = strcat('prec_glacon',num2str(i)) ;
        if i>1
                col4 = {st,'radio','previous channel',0,'Put the previous value',[],''} ;
        else
                col4 = {st,'text',' ',0,'',[],''} ;
        end
        st   = strcat('import_glacon',num2str(i)) ;
        col5 = {st,'radio','load',0,'load the initial value',[],''} ;
        form{length(form)+1} = {col1,col2,col3,col4,col5} ;
end

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

% Formulaire
% ----------
hout=zuicreeform('Edition','edit_refval_pellet','zuiedit_source_pellet_fct','',form,comm) ;

