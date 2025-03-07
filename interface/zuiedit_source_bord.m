function hout=zuiedit_edit_pellet

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('edit_refval_bord') ;
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

col1 = {'libel','text@full','Reference values for edge',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.cons.c
% -----------
nb = evalin('base','param.gene.nbg') ;
for i=1:nb
        if i==1
                col1 = {'c','text','c',10,info.data.cons.c,''} ;
        else
                col1 = {'c','text',' '  ,10,'',''} ;
        end
        st   = sprintf('channel %s',num2str(i)) ;
        col2 = {'num' ,'text',st,10,'Channel ?',''} ;
        st   = strcat('edit_module_c',num2str(i)) ;
        col3 = {st,'radio','Edit',0,'Edit the initial value',[],''} ;
        st   = strcat('prec_c',num2str(i)) ;
        if i>1
                col4 = {st,'radio','previous channel',0,'Put the previous value',[],''} ;
        else
                %col4 = {st,'radio','canal précédent',0,'Recopie le canal précédent',[],'','enable','off'} ;
                col4 = {st,'text',' ',0,'',[],''} ;
        end
        st   = strcat('import_c',num2str(i)) ;
        col5 = {st,'radio','load',0,'load initila value',[],''} ;
        form{length(form)+1} = {col1,col2,col3,col4,col5} ;
end

% data.cons.pomp
% --------------
col1 = {'pomp','text','pomp',10,info.data.cons.pomp,''} ;
col2 = {'txt' ,'text',' ',0,'',[],''} ;
col3 = {'edit_module_pomp','radio','Edit',0,'Edit the initial value',[],''} ;
col4 = {'prec_pomp','text',' ',0,'',[],''} ;
col5 = {'import_pomp','radio','load',0,'Put the previous value',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;


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
hout=zuicreeform('Edition','edit_refval_bord','zuiedit_source_bord_fct','',form,comm) ;

