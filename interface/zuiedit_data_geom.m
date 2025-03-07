% ZUIEDIT_DATA_GEOM   formulaire d'�ition des consignes de g�m�rie
%-------------------------------------------------------------------------------
% fichier zuiedit_data_geom.m ->  
%		    zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	formulaire d'�ition des consignes de g�m�rie
%
% syntaxe  :
%	hout=zuiedit_data_geom;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function hout=zuiedit_data_geom

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('data_geom') ;
if ishandle(hout)
        zuiformvisible(hout) ;
	return
end

% infos pour les tooltips
info = zinfo ;

% formulaire avec la sous structure from
form={};
liste_tag = {} ;

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text@full','separatrix references',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.geo.mode
% -------------
popupmode = {'           symetrical             ', ...
             '           asymetrical           ', ...
             '    (R,Z) separatrix   ', ...
             '    free boundary   '} ;
valeurmode = {0,1,2,3} ;
val = evalin('base','data.geo.mode') ;
ind = val(1)+1 ;
col1 = {'pop_mode','popup@center',popupmode,ind,info.data.geo.mode,valeurmode,''} ;
form{length(form)+1} = {col1} ;

% Separation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.geo.r0
% -----------
col1 = {'r0'       ,'text' ,'r0'      ,10,info.data.geo.r0,''} ;
col2 = {'edit_r0'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_r0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.z0
% -----------
col1 = {'z0'       ,'text' ,'z0'      ,10,info.data.geo.z0,''} ;
col2 = {'edit_z0'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_z0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.a
% ----------
col1 = {'a'       ,'text' ,'a'       ,10,info.data.geo.a,''} ;
col2 = {'edit_a'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_a','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.e1
% -----------
col1 = {'e1'       ,'text' ,'e1'      ,10,info.data.geo.e1,''} ;
col2 = {'edit_e1'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_e1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.trh1
%--------------
col1 = {'trh1'       ,'text' ,'trh1'    ,10,info.data.geo.trh1,''} ;
col2 = {'edit_trh1'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_trh1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.trb1
% -------------
col1 = {'trb1'       ,'text' ,'trb1'    ,10,info.data.geo.trb1,''} ;
col2 = {'edit_trb1'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_trb1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.ind1
% -------------
col1 = {'ind1'       ,'text' ,'ind1'    ,10,info.data.geo.ind1,''} ;
col2 = {'edit_ind1'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_ind1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.geo.b0
% -----------
col1 = {'b0'       ,'text' ,'b0'      ,10,info.data.geo.b0,''} ;
col2 = {'edit_b0'  ,'radio','edit'  ,0,'edit the value',[],''} ;
col3 = {'import_b0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% Separation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% import(R,Z)
% -----------
col1 = {'import_rz','radio@center','load (R,Z)',0,'load the separatrix (R,Z)',[],''} ;
form{length(form)+1} = {col1} ;

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

% Bouton Quit
% ----------
%comm{1}={'btn_quit','radio@center','Fermer',0,'Pour fermer la fen�re'};

% Formulaire
% ----------
hout=zuicreeform('Edition','data_geom','zuiedit_data_geom_fct','',form,'') ;

setappdata(hout,'liste_tag',liste_tag) ;
