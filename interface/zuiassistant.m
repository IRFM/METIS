% ZUIASSISTANT  formulaire des 'assistants' appele depuis le formulaire principal
%---------------------------------------------------------------
% fichier zuiassistant.m ->  
%					             
%
% fonction Matlab 5 :
%	Cette fonction interface les assistant de Conos.
%
% syntaxe  :
%	zuiassistant
%  
% enetress 
%
% sorties 
%
% fonction ecrite par J-F Artaud, poste: 62-15
% version 2.1 du 17/06/2003.
%
% liste des modifications : 
%
%--------------------------------------------------------------

function zuiassistant

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('assistant');
if ishandle(hform)
        zuiformvisible(hform);
	return
end

% 1eres lignes : Nom du fichier de travail
% --------------------------------------
form = {};

col1 = {'titre','text@full','Assistants',50,'',[]};
form{1} = {col1} ;

% alignement jump
% ---------------
colj = {'void','jump','void',[],''} ;
form{length(form)+1} = {colj,colj,colj,colj,colj};

% 4eme ligne
% ----------
col1 = {'radio_jivaro'           ,'radio' ,'Data reduction' ,0,'Divided by N the time step number'};
col2 = {'radio_basetempsch'      ,'radio' ,'time resampling' ,0,'modify the time slace of CRONOS (as "extension")'};
col3 = {'radio_rayon'            ,'radio' ,'x resampling' ,0,'modify the radial coordinate x of the profiles'};
col4 = {'radio_exportascii'      ,'radio' ,'Save data ASCII' ,0,'Save the whole CRONOS datas in an ASCII file [ASCII + tar + gzip]'};
col5 = {'radio_importascii'      ,'radio' ,'load data ASCII   ' ,0,'load an ASCII file  (ASCII + tar + gzip)'};
form{length(form)+1} = {col1,col2,col3,col4,col5};


% 5eme ligne
% ----------
col1 = {'radio_etalhts'           ,'radio' ,'Eta LH TS' ,0,'adjust the LH efficiency to reproduce the flux consumption for TS'};
col2 = {'radio_rematune'          ,'radio' ,'tuning of Rema' ,0,'helper to adjust Rema parameters'};
col3 = {'radio_neident'           ,'radio' ,'nl ident' ,0,'profile identification using interferrometry based on maximum entropy algorithm'};
col4 = {'radio_eceident'          ,'radio' ,'Ece @ rema' ,0,'calculation of ECE radiative temperature with the help of REMA code'};
col5 = {'radio_0d'                ,'radio' ,'Metis' ,0,'call of Metis (fast integrated simulator of Cronos)'};
form{length(form)+1} = {col1,col2,col3,col4,col5};
% 6eme ligne
% ----------
col1 = {'radio_fcitune'           ,'radio' ,'tuning of Pion    ' ,0,'Pion simulation (minority effect, ...)'};
col2 = {'radio_kine   '           ,'radio' ,'kinezero simulation' ,0,'stability studies using kinezero at a given time'};
col3 = {'radio_nemaxent'           ,'radio' ,'ne ident' ,0,'profile identification using interferrometry and reflectometry based on maximum entropy algorithm'};
col4 = {'radio_luketune'           ,'radio' ,'tuning of LUKE    ' ,0,'LUKE simulation (n// effects, ...)'};
col5 = {'radio_neotune'           ,'radio' ,'tuning of NCLASS    ' ,0,'allow to run NCLASS from CRONOS (modification of Ti and Te allowed'};
form{length(form)+1} = {col1,col2,col3,col4,col5};


% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 11eme ligne : Bouton Quit
% ------------------------
comm{1}={'btn_quit','radio@left','Close',0,'to close the window'};
comm{2}={'aide','radio@right','Help',0,'Help'};

hout=zuicreeform('Cronos Assistants ','assistant','zuiassistant_fct','',form,comm);
