% ZUIPLOTCOUCHE	fonction de création de l'interface de l'outil traçant les couches
% ----------------------------------------------------------------------------
%
% fichier zuiplotcouche.m 
%			zuicreeformFCI	: creation du formulaire
%			zuiplotcouche_fct	: gestion des callbacks des uicontrols
%					  du formulaire
%
% fonction Matlab 5 :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.8, du 17/01/2002.
% 
% liste des modifications : 
%
%--------------------------------------------------------------

function zuiplotcouche

% root du programme et initilisation du path
root = getappdata(0,'root') ;
if isempty(root)
	addpath /usr/drfc/cgc/matlab5/zineb/v1.8
        zineb_path;
        root = getappdata(0,'root');
	addpath /usr/drfc/cgc/matlab5/zineb/v1.8/acces/jet/interface
end

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('jetdata');
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% Formulaire
% ----------
form = {};

% titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};
col1 = {'titre','text@full',' Accés aux données JET',[],''};
form{length(form)+1} = {col1};
sepa = {'separation_comm','frame',' ',3,''};
form{length(form)+1} = {sepa};

% Numéro de choc, lecture
% -----------------------
col1 = {'numchoc','text'  ,'Numéro choc',15,''} ;
col2 = {'edit_numchoc','edit',' ',15,'numéro du choc','',''} ;
col3 = {'radio_lecture','radio@right','Lecture',0,'gestion de la lecture des données'} ;
form{length(form)+1} = {col1,col2,col3};

% Etats des lectures
% ------------------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'titre','text@full',' Etat des lectures',[],''};
form{length(form)+1} = {col1};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'push_temp'  ,'push','  Temp   ',0,'',[],'','backgroundcolor','red'} ;
col2 = {'push_far'   ,'push',' Far/MSE ',0,'',[],'','backgroundcolor','red'} ;
col3 = {'push_prof'  ,'push','  Prof   ',0,'',[],'','backgroundcolor','red'} ;
col4 = {'push_transp','push',' Transp  ',0,'',[],'','backgroundcolor','red'} ;
form{length(form)+1} = {col1,col2,col3,col4};

col1 = {'push_idn'   ,'push','  IDN    ',0,'',[],'','backgroundcolor','red'} ;
col2 = {'push_efit'  ,'push','  EFIT   ',0,'',[],'','backgroundcolor','red'} ;
col3 = {'push_qmse'  ,'push','  qmse   ',0,'',[],'','backgroundcolor','red'} ;
col4 = {'push_qpol'  ,'push','  qpol   ',0,'',[],'','backgroundcolor','red'} ;
form{length(form)+1} = {col1,col2,col3,col4};

sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% FIT
% ---
colj = {'void','jump','void',[],''} ;

col1 = {'push_tesh'     ,'push' ,'TeSH'    ,0,'profil de Te calculé à partir des données du superhétérodyne (SH)',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_tesh','radio','FIT'     ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_tesh','radio','Retouche',0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_teth'     ,'push' ,'TeTH'    ,0,'profil de Te calculé à partir des données du LIDAR Thomson (TH)',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_teth','radio','FIT'     ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_teth','radio','Retouche',0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_teshth'     ,'push' ,'TeSH & TH',0,'profil de Te calculé à partir des données du superhétérodyne (SH) et du LIDAR Thomson (TH)',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_teshth','radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_teshth','radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_ti'         ,'push' ,'Ti'       ,0,'profil de Ti',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_ti'    ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_ti'    ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_rot'         ,'push' ,'rot'      ,0,'rotation',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_rot'    ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_rot'    ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_ne'         ,'push' ,'Ne'       ,0,'profil de Ne',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_ne'    ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_ne'    ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_zeff'       ,'push' ,'Zeff'     ,0,'profil de Zeff',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_zeff'  ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_zeff'  ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_car'         ,'push' ,'C'       ,0,'profil de Carbonne',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_car'    ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_car'    ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_nik'         ,'push' ,'Ni'       ,0,'profil de Nickel',[],'','backgroundcolor','red'} ;
col2 = {'radio_fit_nik'    ,'radio','FIT'      ,0,'',[],'','enable','off'} ;
col3 = {'radio_ret_nik'    ,'radio','Retouche' ,0,'',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

% calcul
% ------
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

%
col1 = {'push_dsmf'       ,'push' ,' DSMF' ,0,'état du choix de la DSMF',[],'','backgroundcolor','red'} ;
col2 = {'radio_choix_dsmf','radio','Choix',0,'choix de la dernière surface magnétique fermée (DSMF) en fonction de la valeur relative du flux au bord',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,colj};

%
col1 = {'push_nl'       ,'push' ,' Nl'            ,0,'état de la reconstruction des mesures de l''interférométrie (Nl)',[],'','backgroundcolor','red'} ;
col2 = {'radio_rec_nl'  ,'radio','Reconstruction' ,0,'reconstruction des mesures de l''interférométrie (Nl)',[],'','enable','off'} ;
col3 = {'radio_choix_nl','radio','Choix des voies',0,'sélection des cordes utilisées pour la renormalisation du profil de densité',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,col3};

%
col1 = {'push_compo' ,'push' ,' COMPO' ,0,'état du calcul de la composition du plasma',[],'','backgroundcolor','red'} ;
col2 = {'radio_compo','radio','Composition',0,'réglage de la composition du plasma',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,colj};

%
col1 = {'push_coher'     ,'push' ,' Cohérence' ,0,'état du calcul du bilan d''énergie',[],'','backgroundcolor','red'} ;
col2 = {'radio_cal_coher','radio','Calcul cohérence',0,'calcul du bilan d''énergie et visualisation',[],'','enable','off'} ;
form{length(form)+1} = {col1,col2,colj};

sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% Quit
% ----
comm{1}={'btn_quit','radio@center','Quitter',0,'Pour fermer la fenêtre'};


hout=zuicreeform(' ','jetdata','zuijetdata_fct','',form,comm);

