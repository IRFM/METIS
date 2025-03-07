% ZUIEDIT_DATA_CONS_LIMITES formulaire d'édition des consignes aux limites
%-------------------------------------------------------------------------
% fichier zuiedit_data_cons_limites.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	formulaire d'édition des consignes aux limites
%
% syntaxe  :
%	hout = zuiedit_data_cons_limites;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 30/08/2001.
% 
% liste des modifications : 
%	* 30/08/2001 -> ajout de hform en sortie (J-F Artaud)
%
%--------------------------------------------------------------
function hform=zuiedit_data_cons_limites

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_cons_limites') ;
if ishandle(hform)
        zuiformvisible(hform) ;
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

col1 = {'libel','text@full','boundary consigns',[],''};
form{length(form)+1} = {col1};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.cons.ip
% ------------
col1 = {'ip'       ,'text' ,'ip'    ,10,info.data.cons.ip,''} ;
col2 = {'edit_ip'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_ip','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.vloop
% --------------
col1 = {'vloop'       ,'text' ,'vloop'    ,10,info.data.cons.vloop,''} ;
col2 = {'edit_vloop'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_vloop','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.flux
% --------------
col1 = {'flux'       ,'text' ,'flux'    ,10,info.data.cons.flux,''} ;
col2 = {'edit_flux'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_flux','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.ne1
% -------------
col1 = {'ne1'        ,'text' ,'ne1'    ,10,info.data.cons.ne1,''} ;
col2 = {'edit_ne1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_ne1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.ge1
% -------------
col1 = {'ge1'        ,'text' ,'ge1'    ,10,info.data.cons.ge1,''} ;
col2 = {'edit_ge1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_ge1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.te1
% -------------
col1 = {'te1'        ,'text' ,'te1'    ,10,info.data.cons.te1,''} ;
col2 = {'edit_te1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_te1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.qe1
% -------------
col1 = {'qe1'        ,'text' ,'qe1'    ,10,info.data.cons.qe1,''} ;
col2 = {'edit_qe1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_qe1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.pe1
% -------------
col1 = {'pe1'        ,'text' ,'pe1'    ,10,info.data.cons.pe1,''} ;
col2 = {'edit_pe1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_pe1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.ti1
% -------------
col1 = {'ti1'        ,'text' ,'ti1'    ,10,info.data.cons.ti1,''} ;
col2 = {'edit_ti1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_ti1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.qi1
% -------------
col1 = {'qi1'        ,'text' ,'qi1'    ,10,info.data.cons.qi1,''} ;
col2 = {'edit_qi1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_qi1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.pion1
% -------------
col1 = {'pion1'        ,'text' ,'pion1'    ,10,info.data.cons.pion1,''} ;
col2 = {'edit_pion1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_pion1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.rot1
% -------------
col1 = {'rot1'        ,'text' ,'rot1'    ,10,info.data.cons.rot1,''} ;
col2 = {'edit_rot1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_rot1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.frot1
% -------------
col1 = {'frot1'        ,'text' ,'frot1'    ,10,info.data.cons.frot1,''} ;
col2 = {'edit_frot1'  ,'radio','edit',0,'Plot the consign',[],''} ;
col3 = {'import_frot1','radio','load',0,'load a consign',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.fe1
% -------------
col1 = {'fe1'        ,'text' ,'fe1'    ,10,info.data.cons.fe1,''} ;
col2 = {'edit_fe1'  ,'radio','edit',0,'Plot the consign',[],'','enable','off'} ;
col3 = {'import_fe1','radio','load',0,'load a consign',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.ffe1
% -------------
col1 = {'ffe1'        ,'text' ,'ffe1'    ,10,info.data.cons.ffe1,''} ;
col2 = {'edit_ffe1'  ,'radio','edit',0,'Plot the consign',[],'','enable','off'} ;
col3 = {'import_ffe1','radio','load',0,'load a consign',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.fi1
% -------------
col1 = {'fi1'        ,'text' ,'fi1'    ,10,info.data.cons.fi1,''} ;
col2 = {'edit_fi1'  ,'radio','edit',0,'Plot the consign',[],'','enable','off'} ;
col3 = {'import_fi1','radio','load',0,'load a consign',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3} ;


% data.cons.ffi1
% -------------
col1 = {'ffi1'        ,'text' ,'ffi1'    ,10,info.data.cons.ffi1,''} ;
col2 = {'edit_ffi1'  ,'radio','edit',0,'Plot the consign',[],'','enable','off'} ;
col3 = {'import_ffi1','radio','load',0,'load a consign',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3} ;

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
hform=zuicreeform('Edition','data_cons_limites','zuiedit_data_cons_limites_fct','',form,comm) ;


