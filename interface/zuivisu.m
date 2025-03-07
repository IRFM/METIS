% ZUIVISU  formulaire de 'visulalisation' appele depuis le formulaire principal
%---------------------------------------------------------------
% fichier zuivisu.m ->  
%		zuicreeform 		   : creation du formulaire
% 		zuidirect_fonction	: fonction de gestion des callbacks des uicontrols
%					              du formulaire
%
% fonction Matlab 5 :
%	Cette fonction definie la feuille de visualisation sous 
%	la forme d'un formulaire. Le formulaire est constitue de lignes.
%	Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
%	qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
%	alignees verticalement.
%
% syntaxe  :
%	zuivisu
%  
% enetress 
%
% sorties 
%
% fonction ecrite par C. Passeron, poste 61 19
% version  3.0 , du 14/10/2005.
%
% liste des modifications : 
% * 11/12/2002 : interface en anglais
% * 14/10/2005 : interface en anglais fin

%--------------------------------------------------------------

function zuivisu

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('visu');
if ishandle(hform)
        zuiformvisible(hform);
	return
end

% 1eres lignes : Nom du fichier de travail
% --------------------------------------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'titre','text@full',' Visualisation ',50,'',[]};
form{length(form)+1} = {col1} ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

% alignement jump
% ---------------
colj = {'void','jump','void',[],''} ;
form{length(form)+1} = {colj,colj,colj,colj,colj};

% 1ere ligne
% ----------
col1 = {'radio_zdataplot'   ,'radio' ,'Main visualisation window' ,0,'display all data'};
col2 = {'fast_zdataplot'   ,'radio' ,'Fast visualisation window ' ,0,'display all data (fast window)'};
form{length(form)+1} = {colj,col1,colj,col2,colj};

%sepa ={'separation_comm','frame','',3,''};
%form{length(form)+1} = {sepa};
form{length(form)+1} = {colj,colj,colj,colj,colj};

% 2eme ligne
% ----------
col1 = {'radio_plotverif'   ,'radio' ,'plotverif   ' ,0,'plot temporal data related to the current(JET, TS)'};
col2 = {'radio_compare'     ,'radio' ,'compare     ' ,0,'compare the profiles of two different Cronos runs with TS database'};
col3 = {'radio_polarplot'   ,'radio' ,'polarplot   ' ,0,'Plot Faradays angle (TS and JET'};
col4 = {'radio_mseplot'     ,'radio' ,'mseplot     ' ,0,'Plot MSE JET (profile + time)'};
col5 = {'radio_coefverif'   ,'radio' ,'coefverif   ' ,0,'transport coefficients plot'};
form{length(form)+1} = {col1,col2,col3,col4,col5};

% 3eme ligne
% ----------
col1 = {'radio_neocohere'   ,'radio' ,'neocohere   ' ,0,'neoclassical plot (1)'};
col2 = {'radio_neoverif'    ,'radio' ,'neoverif    ' ,0,'neoclassical plot (2)'};
col3 = {'radio_etatverif'   ,'radio' ,'etaverif   ' ,0,'compare different resistivity models'};
col4 = {'radio_comparetemp' ,'radio' ,'comparetemp ' ,0,'viewgraph of scalar plasma parameters'};
col5 = {'radio_qplot'   ,'radio' ,'qplot   ' ,0,'plot various profiles (q, Te, j, ...)'};
form{length(form)+1} = {col1,col2,col3,col4,col5};

% 4eme ligne
% ----------
col1 = {'radio_power'       ,'radio' ,'Puissance   ' ,0,'Injected power for all heat sources'};
col2 = {'radio_stabmhd'     ,'radio' ,'Stab. MHD   ' ,0,'MHD plot (from Mishka)'};
col3 = {'radio_modemhd'     ,'radio' ,'Mode MHD    ' ,0,'MHD mode (from Mishka)'};
col4 = {'radio_couche'      ,'radio' ,'ICRF Harmonics  ' ,0,'ICRH layer (harmonics : 1-4, minority species : H D He3)'};
col5 = {'radio_cohere'      ,'radio' ,'Coherence TS' ,0,'data coherence for TS'};
form{length(form)+1} = {col1,col2,col3,col4,col5};

% 5eme ligne
% ----------
col1 = {'radio_wdiff'       ,'radio' ,'Data diff   ' ,0,'detects differences between the actual and the reference datasets'};
col2 = {'radio_neutron'     ,'radio' ,'Neutrons    ' ,0,'neutron flux evolution and comparison with experimental data'};
col3 = {'radio_bootstrap'   ,'radio' ,'Bootstrap   ' ,0,'Sauter law compared to CRONOS'};
col4 = {'radio_scenario'    ,'radio' ,'Scenario    ' ,0,'Display scenario'};
col5 = {'radio_postece'     ,'radio' ,'Ece         ' ,0,'Ece measured and computed'};
form{length(form)+1} = {col1,col2,col3,col4,col5};

% 6eme ligne
% ----------
col1 = {'radio_ecrh'        ,'radio' ,'ECRH source        ' ,0,'sources ECRH'};
col2 = {'radio_flux'        ,'radio' ,'Edge flux   ' ,0,'Poloidal edge flux vs time'};
col3 = {'radio_qali'         ,'radio' ,'(qa,li)' ,0,'(qa,li) stability plot'};
col4 = {'radio_itb'          ,'radio' ,'ITB ?' ,0,'Check likeliness of ITB using Tresset criterium'};
col5 = {'radio_sepa'        ,'radio' ,'T-S LCMS' ,0,'Comparison between several positions of identified  TS separatrix'};
form{length(form)+1} = {col1,col2,col3,col4,col5};

% 7eme ligne
% ----------
col1 = {'radio_catalogue'        ,'radio' ,'Summary        ' ,0,'Simulation summary'};
col2 = {'radio_li       '        ,'radio' ,'li             ' ,0,'Comparison li with TS database'};
col3 = {'radio_nustar   '        ,'radio' ,'nu* and rho*   ' ,0,'plot nu* and rho* for ion and electron'};
col4 = {'radio_st       '        ,'radio' ,'Sawteeth (TS)  ' ,0,'Check ST inversion radius'};
col5 = {'radio_movie'        ,'radio' ,'Movie' ,0,'make a movie with 4D plasma reprentation'};
form{length(form)+1} = {col1,col2,col3,col4,col5};


% 8eme ligne
% ----------
col1 = {'radio_confinement'        ,'radio' ,'Confinement   ' ,0,'confinement summary plot including H factor, betaN, qmin and non inductive current'};
col2 = {'radio_bootvalid'        ,'radio' ,'BootValid        ' ,0,'validity of perturbative calculation for bootstrap current'};
col3 = {'radio_ticxs'        ,'radio' ,'TS Ticxs' ,0,'plot Ti from cronos and Ticxs for Tore Supra'};
col4 = {'r4'        ,'radio' ,'' ,0,''};
col5 = {'r5'        ,'radio' ,'' ,0,''};
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
comm{1}={'btn_quit','radio@left','Close',0,'close the window'};
comm{2}={'btn_jeux1','radio@center','delete reference data',0,'delete the reference data '};
comm{3}={'aide','radio@right','Help',0,'Help'};

hout=zuicreeform('Visualisation Zineb','visu','zuivisu_fct','',form,comm);
