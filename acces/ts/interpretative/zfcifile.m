% ZPIONMEX profil de depot fci calcule en utilisant Pion et / ou absor
%--------------------------------------------------------------------
% fichier zfci.m ->  mexfile pionmex2 et programmeabsor
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule le depot pour FCI en FWEH ou Minoritaire.
% valbale TS et JET
% 
% syntaxe  :
%  
%      [sortie,memoire] = zfci(parametre,proto,cons,geo,equi,injection, ...
%                          prof,neo,impur,phys,composition,gene)
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.fci)
%      proto           =    prototype de la structure pour les sources, valeurs a zeros (proto = zsourceproto;)
%      cons            =    consigne de puissance par coupleur (data.cons.fci)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      injection       =    consigne d'injection de gaz (data.cons.c)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      phys             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes ( param.compo)
%      gene            =    parametres generaux (param.gne)
%      memoire         =    structure des dernieres valeurs calculees
% 
% sortie :
% 
%     sortie           =  structure de type source remplie par le module (sortie === proto)
%     memoire          =  datak.memoire.fci (valeur de reference pour le dernier calcul complet, 
%                         pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
%
% fonction ecrite par V. Basiuk , poste 61_26
% version 3.0, du 24/03/2005.
% 
% 
% liste des modifications : 
%   * 26/03/2002 -> securite psupra
%   * 17/04/2002 -> prise en compte de l'absorption directe sur les electrons (necessaire en Helium 3 minoritaire)
%   * 24/04/2002 -> correction structure output et ajout en sortie du taux de conversion de mode
%   * 10/06/2002 -> Profil de densite venant de CRONOS ou homothetique a ne
%   * 11/06/2002 -> Sauvegarde du contexte pion dans des variables MatLab
%   * 18/11/2002 -> bug flux de neutron totale (sortie sfusd), passage en m-3
%   * 11/02/2003 -> blindage absor pour le canal ionique au bord
%   * 17/03/2003 -> prise en compte de la fraction perdue dans le ripple pour TS
%   * 13/04/2004 -> cas ou le premier minoritaire n'est pas l'hydrogene
%   * 04/02/2005 -> modification de time(2)
%   * 09/03/2005 -> modification de la loi pour l'efficacite du FWCD
%   * 24/03/2005 -> suppression du parametre efficacite
%  
%
%--------------------------------------------------------------
%
function [sortie,memoire] = zfcifile(parametre,proto,cons,geo,equi,injection, ...
                             prof,neo,impur,phys,composition,gene,memoire)
                             
%
% mode initialisation 
%
langue                  = getappdata(0,'langue_cronos');
machine                 = getappdata(0,'tokamak_name');
% fonction auto declarante                             
if nargin <=1 
   if isempty(machine)
     try
       machine = evalin('base','param.from.machine');
     catch
       machine = 'TS';
     end
   end
   if nargin == 0
     if strcmp(machine,'JET')
       nbfci = 4;
     elseif strcmp(machine,'TS')
       nbfci = 3;
     else
       nbfci = 1;
     end
   else
     nbfci=parametre;       
   end
%
% description de l'antenne
%	

%
% description d'une antenne
%  	  
      nbstrap         = 2;
	  diststrap       = 0.25;
	  widthstrap      = 0.11;
	  distplas        = 0.05;
	  distwall        = 0.04;
        
        if strcmp(machine,'TS')
	  nbstrap         = 2;
	  diststrap       = 0.19;
	  widthstrap      = 0.1;
	  distplas        = 0.05;
	  distwall        = 0.16;
        end
        if strcmp(machine,'JET')
	  nbstrap         = 2;
	  diststrap       = 0.19;
	  widthstrap      = 0.1;
	  distplas        = 0.05;
	  distwall        = 0.16;
        end
        if strcmp(machine,'ITER')
	  nbstrap         = 2;
	  diststrap       = 0.19;
	  widthstrap      = 0.1;
	  distplas        = 0.05;
	  distwall        = 0.16;
        end
        if strcmp(machine,'EAST')
	  nbstrap         = 2;
	  diststrap       = 0.25;
	  widthstrap      = 0.11;
	  distplas        = 0.05;
	  distwall        = 0.04;
        end
       if strcmp(machine,'SST1')
	  nbstrap         = 2;
	  diststrap       = 0.25;
	  widthstrap      = 0.11;
	  distplas        = 0.05;
	  distwall        = 0.04;
        end
       if strcmp(machine,'FAST')
	  nbstrap         = 2;
	  diststrap       = 0.05;
	  widthstrap      = 0.056;
	  distplas        = 0.07;
	  distwall        = 0.02;
        end
        if isempty(machine)
	  nbstrap         = 2;
	  diststrap       = 0.25;
	  widthstrap      = 0.11;
	  distplas        = 0.05;
	  distwall        = 0.04;
        end

%
%
%	
	if nbfci <= 3 
	  [valeur.mode{1:nbfci}] = deal('FWEH    ');                          % mode de fonctionnement EH ou HMIN
	  valeur.frequence       = 48.*ones(1,nbfci);                     % frequence en MHz
  	  valeur.rant            = 0.02.*ones(1,nbfci);                   % rayon des antennes (m)
	  valeur.decrois         = 2e-2.*ones(1,nbfci);                   % longueur de decroissance de la densite dans la SOL (m)
	  valeur.Npas            = 100;                                   % nombre de maille
	  valeur.Nsol            = 30;                                    % nombre de pas utilise pour le calcul de l'attenuation de l'antenne
	  valeur.int             = 1;                                    % nombre de modes (maximum 81, de 0 �80)
	  valeur.auto            = 'auto';                                % choix automatique du code (pion ou absor) ou manuel (impos�en entr�)
	  valeur.code            = ones(1,nbfci);                         % lancement de pion (1) ou absorb (0) en mode manuel
	  valeur.machine         = 'TS';                                  % nature du TOKAMAK
	  valeur.profil          = 1;                                   % profil ionique CRONOS (1) ou homothetique ne (0)
	  valeur.norm            = 1;                                   % si Pabs > Pcons, 0 -> rien, 1 -> Pabs(t-dt), 2 -> normalisation Pcons 
	  valeur.rip             = 'Yes';                                   % correction ripple en cas de minoritaire 1H 
	  valeur.sup             = 'No';                                   % correction ripple en cas de minoritaire 1H 
          valeur.save             = 'No';	                  % sauvegarde des fichiers contextes
%
% description d'une antenne
%          
	  valeur.nbstrap         = nbstrap;
	  valeur.diststrap       = diststrap;
	  valeur.widthstrap      = widthstrap;
	  valeur.distplas        = distplas;
	  valeur.distwall        = distwall;
	  
	  type.mode              = 'list';                                % type liste de choix
	  type.frequence         = 'float';                               % type reel
	  type.rant              = 'float';                               % type reel
	  type.decrois           = 'float';                               % type reel
	  type.Npas              = 'integer';
	  type.Nsol              = 'integer';
	  type.int               = 'integer';
	  type.auto              = 'list';
	  type.code              = 'integer';
	  type.machine           = 'list';
	  type.profil            = 'integer';                           % type liste de choix
	  type.norm              = 'integer';                           % type liste de choix
	  type.rip             = 'string';
	  type.sup             = 'string';
	  type.save             = 'string';
%
% description d'une antenne
%          
	  type.nbstrap         = 'list';
	  type.diststrap       = 'string';
	  type.widthstrap      = 'string';
	  type.distplas        = 'string';
	  type.distwall        = 'string';

	  borne.mode             = {'HMIN_H  ','HMIN_D  ','HMIN_He ','HMIN_He3','HMIN_T  ','HMIN_2T ','FWEH    ','HARM_2H ','FWCD    ','FWCCD   '};                % valeurs possible
	  borne.frequence        = [10,100];                              % valeurs possible 
	  borne.rant             = [0,0.1];                               % valeurs possibles 
	  borne.decrois          = [1e-3,1e-1];                           % valeurs possibles
	  borne.Npas             = [50 500];                              % valeurs possibles
	  borne.Nsol             = [10 50];                               % valeurs possibles
	  borne.int              = [1 10];                               % valeurs possibles
	  borne.auto             = {'auto','manuel'};                               % valeurs possibles
	  borne.code             = {0,1};                                 % valeurs possibles
          borne.machine          = {'TS','JET','ITER','EAST','FAST'};                   % valeurs possibles
	  borne.profil           = {0,1}   ;                            % valeurs possible 
	  borne.norm             = {0,1,2}   ;                            % valeurs possible 
	  borne.rip              = {'Yes','No'};
	  borne.sup              = {'Yes','No'};
	  borne.save              = {'Yes','No'};
%
% description d'une antenne
%          
	  borne.nbstrap         = {2};
	  borne.diststrap       = [0.01 0.5];
	  borne.widthstrap      = [0.01 0.5];
	  borne.distplas        = [0.01 0.5];
	  borne.distwall        = [0.01 0.5];

	  defaut.mode            = 'FWEH    ';                            % valeurs par defaut
	  defaut.frequence       = NaN;                                   % valeurs par defaut
	  defaut.rant            = 0.02;                                  % valeur par defaut
	  defaut.decrois         = 2e-2;                                  % valeur par defaut
	  defaut.Npas            = 100;                                   % valeur par defaut
	  defaut.Nsol            = 30;                                    % valeur par defaut 
	  defaut.int             = 1;                                    % valeur par defaut
          defaut.auto            = 'auto';                                % valeur par defaut
	  defaut.code            = 0;                                     % valeurs possibles
          defaut.machine         = 'TS';                                  % valeur par defaut
	  defaut.profil          = 1;                                   % valeurs par defaut
	  defaut.norm            = 1;                                   % valeurs par defaut
	  defaut.rip             = 'No';
	  defaut.sup             = 'No';
	  defaut.save            = 'No';
%
% description d'une antenne
%          
	  defaut.nbstrap         = 2;
	  defaut.diststrap       = 0.1;
	  defaut.widthstrap      = 0.06;
	  defaut.distplas        = 0.05;
	  defaut.distwall        = 0.04;

	    info.mode              = 'FCI : FWEH, HMIN ou FWCD';        % informations
	    info.frequence         = 'wave frequency ';
	    info.rant              = '[absor code] antenna distance from the last close surface';
	    info.decrois           = '[absor code] Ln in the SOL';
	    info.Npas              = '[absor code] number of radial point in major radius';                   
	    info.Nsol              = '[absor code] number of radial point in the SOL';                                    
	    info.int               = '[absor code] step on toroidal mode (Nmin:int:Nmax)';                   
	    info.auto              = 'Internal choice for the codes(pion ou absor) or manual';                   
	    info.code              = '0 -> Absor, 1 -> Pion';                   
            info.machine           = 'TOKAMAK name (spectrum antenna)';
	    info.profil            = '[pion] Ti profile from CRONOS (1) ou Ti proportinal to ne (0)';
	    info.norm              = '[pion] if Pabs > Pcons, 0 -> do nothing, 1 -> used Pabs(t-dt), 2 -> normalization on Pcons';
	    info.rip                = 'TS ripple correction for injected power';
	    info.sup                = 'Suprathermal deduced from pion taking into account or not ';
	    info.save                = 'save context for test';
%
% description d'une antenne
%          
	  info.nbstrap         = 'number of strap (only 2 for the moment)';
	  info.diststrap       = 'distance between two straps';
	  info.widthstrap      = 'width of a strap';
	  info.distplas        = 'distance from the plasma';
	  info.distwall        = 'distance from the wall';

	  interface.ts           = '';                            % nom de la fonction d'interfacage avec les donnees TS
	  interface.jet          = '';                            % nom de la fonction d'interfacage avec les donnees Jet

	  sortie.valeur          = valeur;
	  sortie.type            = type;
	  sortie.borne           = borne;
	  sortie.defaut          = defaut;
	  sortie.info            = info;
	  sortie.interface       = interface;
	  sortie.description     = 'ICRH power deposition (PION/ABSORB)';   % description (une ligne) de la fonction
	  sortie.help            = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	  sortie.gui             = '';                            % nom de l'interface graphique specifique si elle existe
	  sortie.controle        = '';                            % nom de la fonction de controle des valeurs si elle existe
	end
	if nbfci == 4
	  [valeur.mode{1:nbfci}] = deal('HMIN_H  ');                        % mode de fonctionnement EH ou HMIN
	  valeur.frequence       = 51.*ones(1,nbfci);                   % frequence en MHz
  	  valeur.rant            = 0.02.*ones(1,nbfci);                 % rayon des antennes (m)
	  valeur.decrois         = 2e-2.*ones(1,nbfci);                 % longueur de decroissance de la densite dans la SOL (m)
	  valeur.Npas            = 100;                                 % nombre de maille
	  valeur.Nsol            = 30;                                  % nombre de pas utilise pour le calcul de l'attenuation de l'antenne
	  valeur.int             = 81;                                  % nombre de modes (maximum 81, de 0 �80)
	  valeur.auto            = 'auto';                                % choix automatique du code (pion ou absor) ou manuel (impos�en entr�)
	  valeur.code            = ones(1,nbfci);                                     % lancement de pion (1) ou absorb (0)
	  valeur.machine         = 'JET';                               % nature du TOKAMAK
	  valeur.profil          = 1;                                   % profil ionique CRONOS (1) ou homothetique ne (0)
	  valeur.norm            = 1;                                   % si Pabs > Pcons, 0 -> rien, 1 -> Pabs(t-dt), 2 -> normalisation Pcons 
	  valeur.sup             = 'No';                                   % correction ripple en cas de minoritaire 1H 
          valeur.save             = 'No';	                  % sauvegarde des fichiers contextes
%
% description d'une antenne
%          
	  valeur.nbstrap         = nbstrap;
	  valeur.diststrap       = diststrap;
	  valeur.widthstrap      = widthstrap;
	  valeur.distplas        = distplas;
	  valeur.distwall        = distwall;

	  type.mode              = 'list';                              % type liste de choix
	  type.frequence         = 'float';                             % type reel
	  type.rant              = 'float';                               % type reel
	  type.decrois           = 'float';                               % type reel
	  type.Npas              = 'integer';
	  type.Nsol              = 'integer';
	  type.int               = 'integer';
	  type.auto              = 'list';
	  type.code              = 'integer';
	  type.machine           = 'list';
	  type.profil            = 'integer';                           % type liste de choix
	  type.norm              = 'integer';                           % type liste de choix
  	  type.sup               = 'string';
	  type.save               = 'string';
%
% description d'une antenne
%          
	  type.nbstrap         = 'list';
	  type.diststrap       = 'string';
	  type.widthstrap      = 'string';
	  type.distplas        = 'string';
	  type.distwall        = 'string';

	  borne.mode             = {'HMIN_H  ','HMIN_D  ','HMIN_He ','HMIN_He3','HMIN_T  ','FWEH    ','HARM_2H ','FWCD    ','FWCCD   '};                % valeurs possible
	  borne.frequence        = [1,1000];                            % valeurs possible 
	  borne.rant             = [0,0.1];                               % valeurs possibles 
	  borne.decrois          = [1e-3,1e-1];                           % valeurs possibles
	  borne.Npas             = [50 500];                              % valeurs possibles
	  borne.Nsol             = [10 50];                               % valeurs possibles
	  borne.int              = [20 81];                               % valeurs possibles
	  borne.auto             = {'auto','manuel'};                               % valeurs possibles
	  borne.code             = {0,1};                                 % valeurs possibles
          borne.machine          = {'TS','JET','ITER','EAST','FAST'};                   % valeurs possibles
	  borne.profil           = {0,1}   ;                            % valeurs possible 
	  borne.norm             = {0,1,2}   ;                            % valeurs possible 
	  borne.sup              = {'Yes','No'};
	  borne.save             = {'Yes','No'};
%
% description d'une antenne
%          
	  borne.nbstrap         = {2};
	  borne.diststrap       = [0.01 0.5];
	  borne.widthstrap      = [0.01 0.5];
	  borne.distplas        = [0.01 0.5];
	  borne.distwall        = [0.01 0.5];

	  defaut.mode            = 'HMIN_H  ';                                % valeurs par defaut
	  defaut.frequence       = NaN;                                 % valeurs par defaut
	  defaut.rant            = 0.02;                                  % valeur par defaut
	  defaut.decrois         = 2e-2;                                  % valeur par defaut
	  defaut.Npas            = 100;                                   % valeur par defaut
	  defaut.Nsol            = 30;                                    % valeur par defaut 
	  defaut.int             = 81;                                    % valeur par defaut
          defaut.auto            = 'auto';                                % valeur par defaut
	  defaut.code            = 1;                                     % valeurs possibles
          defaut.machine         = 'JET';                                 % valeur par defaut
	  defaut.profil          = 1;                                   % valeurs par defaut
	  defaut.norm            = 1;                                   % valeurs par defaut
	defaut.sup          = 'No';
	defaut.save          = 'No';
%
% description d'une antenne
%          
	  defaut.nbstrap         = 2;
	  defaut.diststrap       = 0.1;
	  defaut.widthstrap      = 0.06;
	  defaut.distplas        = 0.05;
	  defaut.distwall        = 0.04;
 	    
	    info.mode              = 'FCI : FWEH, HMIN ou FWCD';        % informations
	    info.frequence         = 'wave frequency ';
	    info.rant              = '[absor code] antenna distance from the last close surface';
	    info.decrois           = '[absor code] Ln in the SOL';
	    info.Npas              = '[absor code] number of radial point in major radius';                   
	    info.Nsol              = '[absor code] number of radial point in the SOL';                                    
	    info.int               = '[absor code] number of toroidal mode';                   
	    info.auto              = 'Internal choice for the codes(pion ou absor) or manual';                   
	    info.code              = '0 -> Absor, 1 -> Pion';                   
       info.machine           = 'TOKAMAK name (spectrum antenna)';
	    info.profil            = '[pion] Ti profile from CRONOS (1) ou Ti proportinal to ne (0)';
	    info.norm              = '[pion] if Pabs > Pcons, 0 -> do nothing, 1 -> used Pabs(t-dt), 2 -> normalization on Pcons';
	  info.sup                = 'if No, supra deduced from pion = supra * 0';
	  info.save                = 'save context for test';
%
% description d'une antenne
%          
	  info.nbstrap         = 'number of strap (only 2 for the moment)';
	  info.diststrap       = 'distance between two straps';
	  info.widthstrap      = 'width of a strap';
	  info.distplas        = 'distance from the plasma';
	  info.distwall        = 'distance from the wall';
	  
	  interface.ts           = '';                   % nom de la fonction d'interfacage avec les donnees TS
	  interface.jet          = '';                   % nom de la fonction d'interfacage avec les donnees Jet

	  sortie.valeur          = valeur;
	  sortie.type            = type;
	  sortie.borne           = borne;
	  sortie.defaut          = defaut;
	  sortie.info            = info;
	  sortie.interface       = interface;

	  sortie.description     = 'ICRH power deposition (pion or absor)';   % description (une ligne) de la fonction

	  sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	  sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	  sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	end
	if nbfci == 1
	  [valeur.mode{1:nbfci}] = deal('HMIN_He3');                        % mode de fonctionnement EH ou HMIN
	  valeur.frequence       = 50.*ones(1,nbfci);                   % frequence en MHz
  	  valeur.rant            = 0.16.*ones(1,nbfci);                 % rayon des antennes (m)
	  valeur.decrois         = 2e-2.*ones(1,nbfci);                 % longueur de decroissance de la densite dans la SOL (m)
	  valeur.Npas            = 100;                                 % nombre de maille
	  valeur.Nsol            = 30;                                  % nombre de pas utilise pour le calcul de l'attenuation de l'antenne
	  valeur.int             = 1;                                  % nombre de modes (maximum 81, de 0 �80)
	  valeur.auto            = 'auto';                                % choix automatique du code (pion ou absor) ou manuel (impos�en entr�)
	  valeur.code            = ones(1,nbfci);                                     % lancement de pion (1) ou absorb (0)
	  valeur.machine         = 'ITER';                               % nature du TOKAMAK
	  valeur.profil          = 1;                                   % profil ionique CRONOS (1) ou homothetique ne (0)
	  valeur.norm            = 1;                                   % si Pabs > Pcons, 0 -> rien, 1 -> Pabs(t-dt), 2 -> normalisation Pcons 
	  valeur.sup             = 'No';                                   % correction ripple en cas de minoritaire 1H 
          valeur.save            = 'No';	                  % sauvegarde des fichiers contextes
%
% description d'une antenne
%          
	  valeur.nbstrap         = nbstrap;
	  valeur.diststrap       = diststrap;
	  valeur.widthstrap      = widthstrap;
	  valeur.distplas        = distplas;
	  valeur.distwall        = distwall;

	  type.mode              = 'list';                              % type liste de choix
	  type.frequence         = 'float';                             % type reel
	  type.rant              = 'float';                               % type reel
	  type.decrois           = 'float';                               % type reel
	  type.Npas              = 'integer';
	  type.Nsol              = 'integer';
	  type.int               = 'integer';
	  type.auto              = 'list';
	  type.code              = 'integer';
	  type.machine           = 'list';
	  type.profil            = 'integer';                           % type liste de choix
	  type.norm              = 'integer';                           % type liste de choix
 	  type.sup               = 'string';
 	  type.save               = 'string';
%
% description d'une antenne
%          
	  type.nbstrap         = 'list';
	  type.diststrap       = 'string';
	  type.widthstrap      = 'string';
	  type.distplas        = 'string';
	  type.distwall        = 'string';

	  borne.mode             = {'HMIN_H  ','HMIN_D  ','HMIN_He ','HMIN_He3','HMIN_T  ','HMIN_2T ','FWEH    ','HARM_2H ','FWCD    ','FWCCD   '};                % valeurs possible
	  borne.frequence        = [10,100];                            % valeurs possible 
	  borne.rant             = [0,0.3];                               % valeurs possibles 
	  borne.decrois          = [1e-3,1e-1];                           % valeurs possibles
	  borne.Npas             = [50 500];                              % valeurs possibles
	  borne.Nsol             = [10 50];                               % valeurs possibles
	  borne.int              = [1 10];                               % valeurs possibles
	  borne.auto             = {'auto','manuel'};                               % valeurs possibles
	  borne.code             = {0,1};                                 % valeurs possibles
         borne.machine           = {'TS','JET','ITER','EAST','FAST'};                   % valeurs possibles
	  borne.profil           = {0,1}   ;                            % valeurs possible 
	  borne.norm             = {0,1,2}   ;                            % valeurs possible 
	  borne.save             = {'Yes','No'};
	  borne.sup              = {'Yes','No'};
%
% description d'une antenne
%          
	  borne.nbstrap         = {2};
	  borne.diststrap       = [0.01 0.5];
	  borne.widthstrap      = [0.01 0.5];
	  borne.distplas        = [0.01 0.5];
	  borne.distwall        = [0.01 0.5];

	  defaut.mode            = 'HMIN_He3';                                % valeurs par defaut
	  defaut.frequence       = 50;                                 % valeurs par defaut
	  defaut.rant            = 0.02;                                  % valeur par defaut
	  defaut.decrois         = 2e-2;                                  % valeur par defaut
	  defaut.Npas            = 100;                                   % valeur par defaut
	  defaut.Nsol            = 30;                                    % valeur par defaut 
	  defaut.int             = 81;                                    % valeur par defaut
          defaut.auto            = 'auto';                                % valeur par defaut
	  defaut.code            = 1;                                     % valeurs possibles
          defaut.machine         = 'JET';                                 % valeur par defaut
	  defaut.profil          = 1;                                   % valeurs par defaut
	  defaut.norm            = 1;                                   % valeurs par defaut
	defaut.save          = 'No';
	defaut.sup          = 'No';
%
% description d'une antenne
%          
	  defaut.nbstrap         = 2;
	  defaut.diststrap       = 0.1;
	  defaut.widthstrap      = 0.06;
	  defaut.distplas        = 0.05;
	  defaut.distwall        = 0.04;

	    info.mode              = 'FCI : FWEH, HMIN ou FWCD';        % informations
	    info.frequence         = 'wave frequency ';
	    info.rant              = '[absor code] antenna distance from the last close surface';
	    info.decrois           = '[absor code] Ln in the SOL';
	    info.Npas              = '[absor code] number of radial point in major radius';                   
	    info.Nsol              = '[absor code] number of radial point in the SOL';                                    
	    info.int               = '[absor code] step of toroidal mode (Nmin:int:Nmax)';                   
	    info.auto              = 'Internal choice for the codes(pion ou absor) or manual';                   
	    info.code              = '0 -> Absor, 1 -> Pion';                   
            info.machine           = 'TOKAMAK name (spectrum antenna)';
	    info.profil            = '[pion] Ti profile from CRONOS (1) ou Ti proportinal to ne (0)';
	    info.norm              = '[pion] if Pabs > Pcons, 0 -> do nothing, 1 -> used Pabs(t-dt), 2 -> normalization on Pcons';
	    info.sup               = 'if No, supra deduced from pion = supra * 0';
	    info.save              = 'save context for test';
%
% description d'une antenne
%          
	  info.nbstrap         = 'number of strap (only 2 for the moment)';
	  info.diststrap       = 'distance between two straps';
	  info.widthstrap      = 'width of a strap';
	  info.distplas        = 'distance from the plasma';
	  info.distwall        = 'distance from the wall';


	  interface.ts           = '';                   % nom de la fonction d'interfacage avec les donnees TS
	  interface.jet          = '';                   % nom de la fonction d'interfacage avec les donnees Jet

	  sortie.valeur          = valeur;
	  sortie.type            = type;
	  sortie.borne           = borne;
	  sortie.defaut          = defaut;
	  sortie.info            = info;
	  sortie.interface       = interface;

	  sortie.description     = 'ICRH power deposition (pion or absor)';   % description (une ligne) de la fonction

	  sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	  sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	  sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	end
	
	
	return
end
[chemin,void]           = fileparts(gene.rapsauve);


% compatibilite entre version
if ~iscell(parametre.mode)
   modecell = {};
   for k =1:size(parametre.mode,1)
       modecell{end+1} = parametre.mode(k,:);
   end
   parametre.mode = modecell;
end
if ~isfield(parametre,'sup')

parametre.sup = 'No';

end

%
% determination du scenario de chaque antenne
%
nbfci = length(parametre.frequence);
disp (['First minority species inside CRONOS : A =',int2str(composition.a(2)),' Z =',int2str(composition.z(2))])
newcomposition      = composition;
if composition.a(2) ~= 1 & composition.z(2) ~= 1 & strcmp(parametre.auto,'auto')
  newcomposition.a(2) = composition.a(3);
  newcomposition.a(3) = composition.a(2);
  newcomposition.z(2) = composition.z(3);
  newcomposition.z(3) = composition.z(2);
end
if strcmp(parametre.auto,'auto')
  for k=1:nbfci
    freq = parametre.frequence(k);
    pos  = geo.r0 + geo.a + parametre.rant(k);
    if ~isnan(freq)
      [scenar(k),scenstr(k,:)]=scenarFCI_pccronos(geo,phys,freq,pos,newcomposition);
    end
  end
else
  scenar = parametre.code;
end

% initialisation de la sortie
sortie = proto;

if sum(scenar) > 0
%
% lancement de pion
%	
  indant = find(scenar == 1);     
   if strcmp(parametre.machine,'TS')
    if strcmp(parametre.rip,'Yes')
      [Ilost,Emean,frlost] = ripple_fci1t(cons,geo,prof,impur,phys,gene,newcomposition,parametre,equi);
      if frlost > 0.25
           %    keyboard         
      end
      Pfci                 = abs(cons) .* (1-frlost);
      phase                = angle(cons);
      Pfci(isnan(Pfci))   = 0;
      for k=1:3
        cons(:,k) = Pfci(:,k).*exp(i.*phase(:,k));
      end
      if isfield('memoire.data','ripple')
          memoire.data.ripple.frlost(end+1,:) = frlost;
          memoire.data.ripple.Ilost(end+1,:) = Ilost;
          memoire.data.ripple.Emean(end+1,:) = Emean;
          memoire.data.ripple.Prip(end+1) = abs(cons);
      else
          memoire.data.ripple.frlost = frlost;
          memoire.data.ripple.Ilost = Ilost;
          memoire.data.ripple.Emean = Emean;
          memoire.data.ripple.Prip = abs(cons);
      end
  
    end
  end
  % debut du calcul
  % petit vecteur utils
  nant   = length(indant);
  va     = ones(nant,1);
  ve     = ones(1,gene.nbrho);

  %
  % donnees concernant le temps d'analyse
  %
  times2    = gene.t;
  vpr2      = equi.vpr;
  psi2      = prof.psi;
  rhomax2   = equi.rhomax;
  dpsidx2   = prof.psid1;
  fdia2     = equi.F;
  ip2       = - 1 ./ ( 2 * pi * phys.mu0 .* rhomax2(end)) .* vpr2(end) .* equi.grho2r2(end) .*  dpsidx2(end);
  shift0    = equi.d(:,1);
  ellip    = equi.e(:,end);
  triang    = (equi.trl(:,end)+equi.trh(:,end))/2;
  q2        = equi.q;
  s2        = sqrt((psi2 - psi2(1))./(psi2(end) - psi2(1)));
  vol2      = rhomax2 .* cumtrapz(gene.x,vpr2,2);
  warning off
  rout2      = interp1(double(equi.rhoRZ),max(squeeze(double(equi.R))')',gene.x .* equi.rhomax,'spline');
  rin2       = interp1(double(equi.rhoRZ),min(squeeze(double(equi.R))')',gene.x .* equi.rhomax,'spline');
  warning on
  ne2        = prof.ne;
  te2        = prof.te;
  ni2        = prof.ni;
  ti2        = prof.ti;
  indmino    = 2;
  indmajo    = 1;
  indmino2   = 3;
  if strcmp(parametre.auto,'manuel')
    icode = find(parametre.code==1);
    if isempty(icode)
       disp('fweh')
       return
    end
    icode = icode(1);
    if strcmp(parametre.mode{icode},'HMIN_He3')
       indmino = find(composition.a == 3 & composition.z == 2);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
    if strcmp(parametre.mode{icode},'HMIN_2T ')
       indmino = find(composition.a == 3 & composition.z == 1);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
    if strcmp(parametre.mode{icode},'HMIN_He')
       indmino = find(composition.a == 4 & composition.z == 2);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
    if strcmp(parametre.mode{icode},'HMIN_H')
       indmino = find(composition.a == 1 & composition.z == 1);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
    if strcmp(parametre.mode{icode},'HMIN_D')
       indmino = find(composition.a == 2 & composition.z == 1);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
    if strcmp(parametre.mode{icode},'HMIN_T')
       indmino = find(composition.a == 3 & composition.z == 1);
       if isempty(indmino)
           error('composition not compatible with fci')
           return
       end
       indmajo = 1;
       indmino2 = 6-indmino-indmajo;
    end
  end
  njc2       = squeeze(impur.impur(1,1,[indmino indmajo indmino2 4 5]))';

  den1       = squeeze(impur.impur(:,:,indmino));
  den2       = squeeze(impur.impur(:,:,indmajo));
  den3       = squeeze(impur.impur(:,:,indmino2));
  den4       = squeeze(impur.impur(:,:,4));
  den5       = squeeze(impur.impur(:,:,5));

  lambda     = 39.1 - 0.5 * log(ne2) + log(te2 ./ 1e3);
  taus2      = max(6.27e8 .* 2 .* (te2 .^ 1.5) ./ (lambda .* ne2 ./1e6));


  if isfield(memoire.data,'input')
	  pion_data = memoire.data.input;
  else
	  pion_data = [];
  end
  if isfield(memoire.data,'output')
	  premierefois = 0;
  else
	  premierefois = 1;
  end

  % si pion_data est vide on cree la structure
  if isempty(pion_data)
  %
  % demarrage de pion avec une fonction maxwellienne
  %
    pion_data.nlanc   = 0;
  else
  %
  % demarrage de pion avec la derniere fonction de distribution calcule
  %
    pion_data.nlanc   = 1;
  end
  pion_data.dt	    = gene.dt;
  pion_data.times     = times2;
  pion_data.vpr	      = vpr2;
  pion_data.psi	    = psi2;
  pion_data.rhomax    = rhomax2;
  pion_data.dpsidx    = dpsidx2;
  pion_data.fdia      = fdia2;
  pion_data.ip	       = ip2;
  pion_data.shift0    = shift0;
  pion_data.q	       = q2;
  pion_data.ellip      = ellip;
  pion_data.triang     = triang;
  pion_data.s	       = s2;
  pion_data.vol	    = vol2;
  pion_data.rout      = rout2;
  pion_data.rin	    = rin2;
  pion_data.ne	       = ne2;
  pion_data.te	       = te2;
  pion_data.ni	       = ni2;
  pion_data.ti	       = ti2;
  pion_data.njc	    = njc2;
  pion_data.den1	    = den1;
  pion_data.den2	    = den2;
  pion_data.den3	    = den3;
  pion_data.den4	    = den4;
  pion_data.den5	    = den5;
  pion_data.cons      = cons;
  pion_data.taus      = taus2;
  %
  % si plusieurs creneaux FCI, demarrage de pion avec une fonction maxwellienne
  %
  if pion_data.dt >= 2*taus2
    pion_data.nlanc   = 0;
  end
  %
  % test sur la puissance
  %
  if all(abs(cons(indant))< 3e5) | max(prof.ti) < 500

      disp('injected power too low for PION')

  else
      disp('PION is starting')

% mis a jour de la structure memoire
    memoire.t          = gene.t;
    memoire.data.input = pion_data;
%  numero propre pour les nom de fichier
    pid = getappdata(0,'pid');
    if isempty(pid)
	 [pid,pid_sess,user]=getidprocess;
	 setappdata(0,'pid',pid);
    end

    % formatage du nombre
    nums = sprintf('%5.5d',pid);
    if length(nums) >5
       nums = nums((end-4):end);
    end

    % creation des donnees pour pionmex
    % creation des variables
    % parametre reechantillonage:
    param.ondelette        = 0;
    param.defaut.temps     = 0;
    param.defaut.espace    = 0;
    param.defaut.inf       = 0;
    param.plus             = 1;
    param.energie          = 1;
    %
    % changement en 21 points radiaux
    %
    rhopion   = linspace(0,1,21);
    Num       = str2num(nums);
    ntact     = 2;
    Nrho      = 21;
    Nspec     = gene.nbg;
    Nant      = nant;
    %
    % variables d'entree de pion
    %
    time      = pion_data.times;
    curr      = pion_data.ip./1e6;
    shift0    = pion_data.shift0;
    qax       = min(pion_data.q);
    ellip     = pion_data.ellip;
    triang    = pion_data.triang;
    PICRH     = abs(pion_data.cons(indant))./1e6;
    FREQ      = parametre.frequence(indant);
    FREQ(isnan(FREQ)) = 50*ones(size(FREQ(isnan(FREQ))));
    Aph       = mean(angle(pion_data.cons(indant)),1);
    rho       = rhopion;
    x         = gene.x;
    s         = zsample(pion_data.s',x',rho,param)';
    Vol       = zsample(pion_data.vol',x',rho,param)';
    Rout      = zsample(pion_data.rout',x',rho,param)';
    Rin       = zsample(pion_data.rin',x',rho,param)';
    Fdia      = zsample(pion_data.fdia',x',rho,param)';
    Bout      = Fdia ./ Rout;
    Bin       = Fdia ./ Rin;
    if strcmp(parametre.auto,'auto')
      Acharg    = newcomposition.z([2,1,3:length(composition.z)]);
      Amass     = newcomposition.a([2,1,3:length(composition.z)]);
    else
      Acharg    = composition.z([indmino,indmajo,indmino2,4:length(composition.z)]);
      Amass     = composition.a([indmino,indmajo,indmino2,4:length(composition.z)]);
    end
    Tempe     = zsample(pion_data.te'/1e3+eps,x',rho,param)';
    Tempi     = zsample(pion_data.ti'/1e3+eps,x',rho,param)';
    Dense     = zsample(pion_data.ne'+1e17,x',rho,param)';
    nion1     = zsample(pion_data.den1'+1e17,x',rho,param)';
    nion2     = zsample(pion_data.den2'+1e17,x',rho,param)';
    nion3     = zsample(pion_data.den3'+1e17,x',rho,param)';
    nion4     = zsample(pion_data.den4'+1e17,x',rho,param)';
    nion5     = zsample(pion_data.den5'+1e17,x',rho,param)';
    Densi0    = pion_data.njc+1e17;
    Densi0    = Densi0(:,[1,2,3:length(composition.z)]);
    %
    % Il faut deux temps dans pion, le calcul ne s'effectuant que sur 1 temps
    %
    Densi0(2,:)   = Densi0(1,:)*0.999;
    nion1(2,:)    = nion1(1,:)*0.999;
    nion2(2,:)    = nion2(1,:)*0.999;
    nion3(2,:)    = nion3(1,:)*0.999;
    nion4(2,:)    = nion4(1,:)*0.999;
    nion5(2,:)    = nion5(1,:)*0.999;
    Dense(2,:)    = Dense(1,:)*0.999;
    Tempi(2,:)    = Tempi(1,:)*0.999;
    Tempe(2,:)    = Tempe(1,:)*0.999;
    Bin(2,:)      = Bin(1,:);
    Bout(2,:)     = Bout(1,:);
    Rin(2,:)      = Rin(1,:);
    Rout(2,:)     = Rout(1,:);
    Fdia(2,:)     = Fdia(1,:);
    Vol(2,:)      = Vol(1,:);
    s(2,:)        = s(1,:);
    PICRH(2,:)    = PICRH(1,:);
    FREQ(2,:)     = FREQ(1,:);
    curr(2)       = curr(1);
    qax(2)        = qax(1);
    ellip(2)      = ellip(1);
    triang(2)     = triang(1);
    shift0(2)     = shift0(1);
    time(2)       = time(1)+gene.dt;
    %
    % dimensionnement a 1000 de toutes les variables temps
    % pour le mexfile
    %
    Densi0(1000,5)   = 0;
    Dense(1000,Nrho) = 0;
    nion1(1000,Nrho) = 0;
    nion2(1000,Nrho) = 0;
    nion3(1000,Nrho) = 0;
    nion4(1000,Nrho) = 0;
    nion5(1000,Nrho) = 0;
    Tempi(1000,Nrho) = 0;
    Tempe(1000,Nrho) = 0;
    Bin(1000,Nrho)   = 0;
    Bout(1000,Nrho)  = 0;
    Rin(1000,Nrho)   = 0;
    Rout(1000,Nrho)  = 0;
    Fdia(1000,Nrho)  = 0;
    Vol(1000,Nrho)   = 0;
    s(1000,Nrho)     = 0;
    PICRH(1000,3)    = 0;
    FREQ(1000,3)     = 0;
    curr(1000)       = 0;
    qax(1000)        = 0;
    ellip(1000)      = 0;
    triang(1000)     = 0;
    shift0(1000)     = 0;
    time(1000)       = 0;
    %
    % on rajoute des antennes fictives de puissance nulle
    %
    if Nant < 4
      PICRH(1000,4)  = 0;
      FREQ(1000,4)   = 0;
      Aph(4)         = pi;
    end
    %
    %
    %
    zvar(1)          = Num;
    zvar(2)          = ntact;
    zvar(3)          = Nrho;
    zvar(4)          = Nspec;
    zvar(5)          = Nant;
    zvar(6)          = pion_data.nlanc;
    zvar(7)          = 0;
    if isfield(parametre,'profil')
     zvar(8)         = parametre.profil;
    else
     zvar(8)         = 1;
    end
    if strcmp(parametre.machine,'ITER')
     zvar(9)         = 1;
%
% probleme puissance depose sur les ions
% on considere pion toujours premier temps en attendant mieux
%
%     zvar(6) = 0;

    else
     zvar(9)         = 0;
    end
%
% description d'une antenne (2 straps, antenne identique)
%          
    zvar(10) = parametre.diststrap;       
    zvar(11) = parametre.widthstrap;
    zvar(12) = parametre.distplas;
    zvar(13) = parametre.distwall;
    
%
% lancement de pionmex en mexfile
%

if zvar(6) == 0
   disp('first time for PION')

  ipicall = 1;
  pioncall;

%
% passage en m-3 du flux de neutron total deduit de pion
%
  sfusd = sfusd * 1e6;
elseif premierefois == 0


    disp('next PION (use the last distribution function calculate)')
 
  flo1  		 = memoire.data.output.flo1;
  flo2  		 = memoire.data.output.flo2;
  flo3  		 = memoire.data.output.flo3;
  flo4  		 = memoire.data.output.flo4;
  flo5  		 = memoire.data.output.flo5;
  flo6  		 = memoire.data.output.flo6;
  sig1  		 = memoire.data.output.sig1;
  sig2  		 = memoire.data.output.sig2;
  sig3  		 = memoire.data.output.sig3;
  sig4  		 = memoire.data.output.sig4;
  sig5  		 = memoire.data.output.sig5;
  sig6  		 = memoire.data.output.sig6;
  va1			 = memoire.data.output.va1;
  va2			 = memoire.data.output.va2;
  va3			 = memoire.data.output.va3;
  va4			 = memoire.data.output.va4;
  va5			 = memoire.data.output.va5;
  va6			 = memoire.data.output.va6;
  fno1  		 = memoire.data.output.fno1;
  fno2  		 = memoire.data.output.fno2;
  fno3  		 = memoire.data.output.fno3;
  fno4  		 = memoire.data.output.fno4;
  fno5  		 = memoire.data.output.fno5;
  fno6  		 = memoire.data.output.fno6;
  nvm			 = memoire.data.output.nvm;
  fc 			 = memoire.data.output.fc;
  fl 			 = memoire.data.output.fl;
  fr 			 = memoire.data.output.fr;
  tzna       = memoire.data.output.tzna;
  sd 			 = memoire.data.output.sd;
  tzf			 = memoire.data.output.tzf;
  nhres 		 = memoire.data.output.nhres;
  nharmi		 = memoire.data.output.nharmi;
  is 			 = memoire.data.output.is;
  nres  		 = memoire.data.output.nres;

  ipicall = 2;
  pioncall;
  
%
% passage en m-3 du flux de neutron total deduit de pion
%
  sfusd = sfusd * 1e6;
else
    if strcmp(langue,'francais')
       disp('PION idem first time (no memory structure)')
    else
      disp('lancement pion idem premier temps, pas de fichier output')
    end
  zvar(6) = 0;

  
  ipicall = 1;
  pioncall;

  
%
% passage en m-3 du flux de neutron total deduit de pion
%
  sfusd = sfusd * 1e6;
end



	 output.piglob       = piglob(2,:);
	 output.spoyfl       = spoyfl(:,2);
	 output.spdrf1       = spdrf1(:,2);
	 output.spdrf2       = spdrf2(:,2);
	 output.spdrf3       = spdrf3(:,2);
	 output.spdrfe       = spdrfe(:,2);
	 output.spdic        = spdic(:,2);
	 output.spdec        = spdec(:,2);
	 output.sw1          = sw1(:,2);
	 output.sw2          = sw2(:,2);
	 output.sw3          = sw3(:,2);
	 output.swz1         = swz1(:,2);
	 output.swz2         = swz2(:,2);
	 output.swz3         = swz3(:,2);
	 output.sfusd        = sfusd(:,2);
	 output.sfd1         = sfd1(:,2);
	 output.sfd2         = sfd2(:,2);
	 output.swdf         = swdf(:,2);
	 output.swdzf        = swdzf(:,2);
	 output.flo1         = flo1;
	 output.flo2         = flo2;
	 output.flo3         = flo3;
	 output.flo4         = flo4;
	 output.flo5         = flo5;
	 output.flo6         = flo6;
	 output.sig1         = sig1;
	 output.sig2         = sig2;
	 output.sig3         = sig3;
	 output.sig4         = sig4;
	 output.sig5         = sig5;
	 output.sig6         = sig6;
	 output.va1          = va1;
	 output.va2          = va2;
	 output.va3          = va3;
	 output.va4          = va4;
	 output.va5          = va5;
	 output.va6          = va6;
	 output.fno1         = fno1;
	 output.fno2         = fno2;
	 output.fno3         = fno3;
	 output.fno4         = fno4;
	 output.fno5         = fno5;
	 output.fno6         = fno6;
	 output.nvm          = nvm;
	 output.fc           = fc;
	 output.fl           = fl;
	 output.fr           = fr;
	 output.tzna         = tzna;
	 output.sd           = sd;
	 output.tzf          = tzf;
	 output.nhres        = nhres;
	 output.nharmi       = nharmi;
	 output.is           = is;
	 output.nres         = nres;

    sortie.ion          = interp1(rhopion,spdic(:,2)' .* (spdic(:,2) > 0)' .* 1e6,x);
%
% contribution sur les electrons : spdec (venant du minoritaire) + spdrfe (absorption directe)
%
    contri_el           = spdec(:,2) + spdrfe(:,2);
    sortie.el           = interp1(rhopion,contri_el' .* (contri_el > 0)' .* 1e6,x);
%
%  corretcion bosse au bord en attente d'etude du probleme
%
    sortie.el           = corpion(sortie.el,gene.x);


%
% protection si puissance absorbee trop forte dans pion
% on prend le dernier calcul renormalisee sur la puissance injectee
%
   Ptot                = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
   normpion            = Ptot/sum(abs(cons(indant)));
   convm               = output.spoyfl(end);
   if strcmp(parametre.save,'Yes')

       [chemin,void] = fileparts(gene.rapsauve);
       save(fullfile(chemin,sprintf('last_pion@%s',int2str(fix(gene.t*1000)))));

   end

   if normpion < 0.2 | normpion > 10
     ecart = 1;
     [chemin,void] = fileparts(gene.rapsauve);
     save(fullfile(chemin,'prob_pion'));
     disp(['problem with pion, new run as first time and prob_pion generated'])
     zvar(6) = 0;

     ipicall = 1;
     pioncall;
     
%
% passage en m-3 du flux de neutron total deduit de pion
%
     sfusd = sfusd * 1e6;
     output.piglob       = piglob(2,:);
     output.spoyfl       = spoyfl(:,2);
	 output.spdrf1       = spdrf1(:,2);
	 output.spdrf2       = spdrf2(:,2);
	 output.spdrf3       = spdrf3(:,2);
	 output.spdrfe       = spdrfe(:,2);
	 output.spdic        = spdic(:,2);
	 output.spdec        = spdec(:,2);
	 output.sw1          = sw1(:,2);
	 output.sw2          = sw2(:,2);
	 output.sw3          = sw3(:,2);
	 output.swz1         = swz1(:,2);
	 output.swz2         = swz2(:,2);
	 output.swz3         = swz3(:,2);
	 output.sfusd        = sfusd(:,2);
	 output.sfd1         = sfd1(:,2);
	 output.sfd2         = sfd2(:,2);
	 output.swdf         = swdf(:,2);
	 output.swdzf        = swdzf(:,2);
	 output.flo1         = flo1;
	 output.flo2         = flo2;
	 output.flo3         = flo3;
	 output.flo4         = flo4;
	 output.flo5         = flo5;
	 output.flo6         = flo6;
	 output.sig1         = sig1;
	 output.sig2         = sig2;
	 output.sig3         = sig3;
	 output.sig4         = sig4;
	 output.sig5         = sig5;
	 output.sig6         = sig6;
	 output.va1          = va1;
	 output.va2          = va2;
	 output.va3          = va3;
	 output.va4          = va4;
	 output.va5          = va5;
	 output.va6          = va6;
	 output.fno1         = fno1;
	 output.fno2         = fno2;
	 output.fno3         = fno3;
	 output.fno4         = fno4;
	 output.fno5         = fno5;
	 output.fno6         = fno6;
	 output.nvm          = nvm;
	 output.fc           = fc;
	 output.fl           = fl;
	 output.fr           = fr;
	 output.tzna         = tzna;
	 output.sd           = sd;
	 output.tzf          = tzf;
	 output.nhres        = nhres;
	 output.nharmi       = nharmi;
	 output.is           = is;
	 output.nres         = nres;

      sortie.ion          = interp1(rhopion,spdic(:,2)' .* (spdic(:,2) > 0)' .* 1e6,x);
%
% contribution sur les electrons : spdec (venant du minoritaire) + spdrfe (absorption directe)
%
      contri_el           = spdec(:,2) + spdrfe(:,2);
      sortie.el           = interp1(rhopion,contri_el' .* (contri_el > 0)' .* 1e6,x);
%
%  corretcion bosse au bord en attente d'etude du probleme
%
      sortie.el           = corpion(sortie.el,gene.x);
%
% protection si puissance absorbee trop forte dans pion
% on prend le dernier calcul renormalisee sur la puissance injectee
%
      Ptot                = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
      normpion            = Ptot/sum(abs(cons(indant)));
      convm               = output.spoyfl(end);
   else
     ecart = 0;
   end


   if ecart == 1 & parametre.norm > 0
       if strcmp(langue,'francais')
          disp(['probleme PION, Pabs/Pcons = ',num2str(normpion,3)])
       else
          disp(['problem PION, Pabs/Pcons = ',num2str(normpion,3)])

       end
		if parametre.norm == 1
	 %
	 % prise du temps precedent et normalisation avec la puissance sauf pour le premier lancement de pion
	 %
   	  if premierefois == 0
      	 fact   = sum(abs(cons(indant))) / sum(abs(memoire.data.output.cons(indant)));
      	 if isfield(memoire.data.output,'spdic')
      		sortie.ion          = fact*interp1(rhopion,memoire.data.output.spdic' .* (memoire.data.output.spdic > 0)' .* 1e6,x);
	 %
	 % contribution sur les electrons : spdec (venant du minoritaire) + spdrfe (absorption directe)
	 %
      		contri_el           = memoire.data.output.spdec + memoire.data.output.spdrfe;
      		sortie.el           = fact*interp1(rhopion,contri_el' .* (contri_el > 0)' .* 1e6,x);
      		sortie.neutron.dd   = fact*interp1(rhopion,memoire.data.output.sfusd' .* (memoire.data.output.sfusd > 0)',x);
      		psupratot           = fact*interp1(rhopion,memoire.data.output.swdf'  .* (memoire.data.output.swdf  > 0)' .* 1e6,x);
      		psuprapara          = fact*interp1(rhopion,memoire.data.output.swdzf' .* (memoire.data.output.swdzf > 0)' .* 1e6,x);
	   		sortie.err          = (1-memoire.data.output.spoyfl(end))+ i * fact;
      	 else
	 %
	 % attente d'un pion correct, normalisation forcee sur la valeur de consigne
	 %
	   		fact                = normpion;
      		sortie.ion          = normpion*sortie.ion;
      		sortie.el           = normpion*sortie.el;
      		sortie.neutron.dd   = normpion*interp1(rhopion,output.sfusd' .* (output.sfusd > 0)',x);
      		psupratot           = normpion*interp1(rhopion,output.swdf'  .* (output.swdf  > 0)' .* 1e6,x);
      		psuprapara          = normpion*interp1(rhopion,output.swdzf' .* (output.swdzf > 0)' .* 1e6,x) ;
	   		sortie.err          = (1-output.spoyfl(end)) + i * fact;

      	 end
		  else
	 fact                = normpion;
      	 sortie.ion          = normpion*sortie.ion;
      	 sortie.el           = normpion*sortie.el;
      	 sortie.neutron.dd   = normpion*interp1(rhopion,output.sfusd' .* (output.sfusd > 0)',x);
      	 psupratot           = normpion*interp1(rhopion,output.swdf'  .* (output.swdf  > 0)' .* 1e6,x);
      	 psuprapara          = normpion*interp1(rhopion,output.swdzf' .* (output.swdzf > 0)' .* 1e6,x) ;
	   	 sortie.err          = (1-output.spoyfl(end)) + i * fact;

		  end
		else
	 %
	 % normalisation sur la valeur de consigne
	 %
	 fact                = normpion;
      	 sortie.ion          = normpion*sortie.ion;
      	 sortie.el           = normpion*sortie.el;
      	 sortie.neutron.dd   = normpion*interp1(rhopion,output.sfusd' .* (output.sfusd > 0)',x);
      	 psupratot           = normpion*interp1(rhopion,output.swdf'  .* (output.swdf  > 0)' .* 1e6,x);
      	 psuprapara          = normpion*interp1(rhopion,output.swdzf' .* (output.swdzf > 0)' .* 1e6,x) ;
	 sortie.err          = (1-output.spoyfl(end)) + i * fact;


		end
	 else
	 %
	 % pas de probleme sur la puissance absorbe, mise en memoire de pion
	 %
		if ecart == 0
   	  disp(['PION ok, Pabs/Pcons = ',num2str(normpion,3)])
		else
    if strcmp(langue,'francais')
   	  disp(['PION nok, Pabs/Pcons = ',num2str(normpion,3),' Pabs conserve par choix'])
    else
   	  disp(['PION nok, Pabs/Pcons = ',num2str(normpion,3),' Pabs=Pcons imposed'])
    end

		end
		memoire.data.output     = output;
		memoire.data.output.cons = cons;
		sortie.neutron.dd   = interp1(rhopion,sfusd(:,2)' .* (sfusd(:,2) > 0)',x);
		psupratot           = interp1(rhopion,swdf(:,2)'  .* (swdf(:,2)  > 0)' .* 1e6,x);
		psuprapara          = interp1(rhopion,swdzf(:,2)' .* (swdzf(:,2) > 0)' .* 1e6,x) ;
		sortie.err          = (1-spoyfl(end,2)) + i * normpion;

	 end

	 sortie.psupra       = psupratot - psuprapara;
%
% blindage pression supra negative
%
    sortie.psupra(sortie.psupra<=0) = psupratot(sortie.psupra<=0)*2/3;
    sortie.psupra(sortie.psupra<=0) = 0;
    sortie.paniso       = psuprapara - (1/2) .* sortie.psupra;
    Psup                = zintvol(sortie.psupra,gene.x,equi.vpr,equi.rhomax)/1e6;
    if Psup > 0.15
      %keyboard
    end
    %
    % probleme psupra (test 12 septembre 2001, V. Basiuk)
    %
    if strcmp(parametre.sup,'No')
      sortie.psupra     = 0*sortie.psupra;
      sortie.paniso     = 0*sortie.paniso;
    end
    % commutateurs de mode


  end

end

%
% lancement d'absor
%
if sum(scenar) < nbfci

  indant = find(scenar == 0);
  % compatibilite
  parametre.freq = parametre.frequence(indant);

  % test sur la puissance
  if all(abs(cons(indant))< 3e5)

     disp('injected power too low for ABSORB ')


  else


    % debut du calcul
    % petit vecteur utils
    nant   = length(indant);
    va     = ones(nant,1);
    ve     = ones(1,gene.nbrho);

    % positions antennes
    rant   = (geo.r0+geo.a+parametre.rant(indant));    % rayon des antennes (m)

    % appel de absorb-TS
    rho    = equi.a ./ equi.a(end);
    nion1  = squeeze(impur.impur(1,:,1));
    nion2  = squeeze(impur.impur(1,:,2));
    rm2    = geo.r0 + equi.d;
    for k = 1: nant
	 % traitement des donnees
	 %
	 % machine = 1=TS,2=JET,3=ITER,...
	 % Rant = position de l'antenne en m
	 % phi0 = phase en radian
	 % Puisantenne = puissance de l'antenne en MW
	 % freq = frequence en MHz
	 % phys = structure des constantes physiques
	 % compo = structure de la composition plasma
	 % rho = rho geometrique normalise
	 % Te = temperature electronique en keV
	 % Tion1 = temperature ionique epsece majoritaire en keV
	 % Tion2 = temperature ionique epsece minoritaire en keV
	 % ne    = densite electronique en m-3
	 % nion1 = densite ionique epsece majoritaire en m-3
	 % nion2 = densite ionique epsece minoritaire en m-3
	 % rm2   = centre des surfaces magnetiques (R0+d0(rho))
	 % a     = petit rayon (m)
	 % R0    = grand rayon (m)
	 % decrois = decroissance de la densite dans la SOL (m)
	 % bcentre = B0 (T)
	 %
	 % sortie
	 % PDP1 = profil de puissance absorbe par les electrons
	 % PDP2 = profil de puissance absorbe par les ions majoritaires
	 % PDP3 = profil de puissance absorbe par les ions minoritaires
	 % Pd1  = profil de puissance volumique absorbe par les electrons (MW/m3)
	 % Pd2  = profil de puissance volumique absorbe par les ions majoritaires (MW/m3)
	 % Pd3  = profil de puissance volumique absorbe par les ions minoritaires (MW/m3)
	 % P    = puissance integree absorbee (MW) (1-> electron, 2 -> majoritaire, 3 -> minoritaire)

	 kant   = indant(k);
	 if strcmp(parametre.machine,'TS')
	   machine = 1;
	 end
	 if strcmp(parametre.machine,'JET')
	   machine = 2;
	 end
	 if strcmp(parametre.machine,'ITER')
	   machine = 3;
	 end
	 option.Npas = parametre.Npas;
	 option.Nsol = parametre.Nsol;
	 option.int  = parametre.int;
         phi = angle(cons(kant));
	 [PDP1,PDP2,PDP3,Pd1,Pd2,Pd3,P]  = ...
	   zdepotas(machine,rant(kant),phi,abs(cons(kant))./1e6, ...
                      parametre.freq(kant),phys,newcomposition,rho,prof.te./1e3, ...
                      prof.ti./1e3,prof.ti./1e3,prof.ne,nion1,nion2, ...
                      rm2,geo.a,geo.r0,parametre.decrois(kant),geo.b0,option);

	sortie.ion        = sortie.ion + Pd2 .* 1e6 + Pd3 .*1e6;
%
% blindage depot ionique au bord
%
   if sortie.ion(end) > 3e6
      sortie.ion(end) = 1e-5;
      sortie.ion(end-1) = 2e-5;
   end
	sortie.el         = sortie.el  + Pd1 .* 1e6;
    end
   Ptot    = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
   normab  = sum(abs(cons(indant)))/Ptot;
    if strcmp(parametre.save,'Yes')

       [chemin,void] = fileparts(gene.rapsauve);
       save(fullfile(chemin,sprintf('last_absor@%s',int2str(fix(gene.t*1000)))));

    end

   if normab < 0.98 | normab > 1.02
     ecart = 1;
    if strcmp(langue,'francais')
     disp(['probleme Pabs d''absorb ecart=',num2str(normab)])
     else
     disp(['problem Pabs d''ABSORB difference=',num2str(normab)])
     end

     sortie.el   = sortie.el .* normab;
     sortie.ion  = sortie.ion .* normab;    
   else
     ecart = 0;
   end

    % commutateurs de mode
    f_fwcd = strcmp(deblank(upper(parametre.mode(indant))),'FWCD')'*ve  -  strcmp(deblank(upper(parametre.mode(indant))),'FWCCD')'*ve;
    % repartition par antenne
    p = (va*sortie.el).* ((abs(cons(indant)')./sum(abs(cons(indant))))*ve); 

    % generation de courant 
    pp      = [-0.00025162783491,0.01039466555003,-0.00111295969404];
    % etaj = etaI *2 * pi (passage de la densite de puissance et de courant au courant et a la puissance)
    gamma   = va * (polyval(pp,prof.te ./ 1e3) .* 6 ./ (5 + prof.zeff)./ prof.ne .* 1e20  .* 2 .* pi);
    sortie.j =sortie.j + sum( gamma .* p .* f_fwcd,1);
  
  end
	
end
% integralle volumique 
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax   
function s=zintvol(e,x,vpr,rhomax)   

  s = rhomax.*trapz(x,vpr .* e,2);
