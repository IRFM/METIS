% ZINIT initialisation des structures de donnees de zineb a 0
%--------------------------------------------------------------
% fichier zinit.m ->  init
%
%
% fonction Matlab 5 :
%
% Initialise les structures de donnees du programme zineb a 0
% (sauf le temps). Les donnes sont sauvegardees dans un fichier
% en vue de l'appel de zineb
%
% syntaxe  :
%  
% [cr,data,param]  = zinit(mode,temps,nbrho,file,tdeb,tfin, ...
%                         {nbfci,nbfce,nbidn,nbhyb,nbglacon,nbg});
%
%
% entrees :
%
%     mode  =   mode d'initialisation (chaine):
%                    'TS'  ->   prepare les structure de donnees
%                               pour Tore-Supra.
%
%                    'JET' ->   prepare les structure de donnees
%                               pour Jet.
%
%                    ''    ->   prepare les structure de donnees
%                               pour une machine quelconque.
%                               
%     temps   = base temps de la simulation 
%               (cette base de temps peut etre modifier ensuite en fonction  des evenements)
%             
%     nbrho   = nombre de points radiaux
%
%     file    =  nom du fichier de sauvegarde des resultats 
%               (avec le chemin complet d'acces)
%
%     tdeb    =  temps de debut pour les calculs , si vide premier temps
%
%     tfin    =  temps de fin pour les calculs , si vide dernier temps
%
%
% si mode == '' :
% 
%     nbfci      =  nombre de coupleurs fci
%     nbfce      =  nombre de coupleurs fce
%     nbidn      =  nombre d'injecteurs idn
%     nbhyb      =  nombre de coupleurs hybride
%     nbglacon   =  nombre de ligne d'injection de glacons
%     nbg        =  nombre de gaz et d'impuretee
% 
% sorties :
% 
%     cr         =  compte rendu d'execution
% 
%     data       =  structure de donnees contenant les variables 
%                   dependant du temps.
%                  
%     param      = structure de donnees contenant les variables 
%                  independantes du temps. 
%
% fonction ecrite par J-F Artaud , poste 62-15
% version 3.0, du 21/03/2005.
% 
% 
% liste des modifications : 
%
%   * 08/06/2001 -> ajout de la structure data.source.externe
%   * 19/06/2001 -> remplissage uniquement des coefficients des equations utilisees (ajout du parametre gene.nbeq_mode)
%   * 25/06/2001 -> ajout de param.source_bord (controle du calcul du recyclage)
%   * 12/07/2001 -> ajout de data.equi.raxe 
%   * 24/07/2001 -> ajout de consignes (pseudo consigne) pour l'edition des donnees
%   * 24/07/2001 -> ajout de consignes pour les asservissements (ti0 et  ti1)
%   * 28/08/2001 -> ajout des champ option et createur dans param.from
%   * 14/09/2001 -> ajout du champ param.gene.filetype
%   * 14/09/2001 -> mise en place de la fonction zinebversion
%   * 20/09/2001 -> ajout de data.equi.shear
%   * 21/09/2001 -> ajout de la paroi dans la struture param.from
%   -------------------- derniere mise a jour de la liste des variable le 27/09/2001 ---------------
% 
%   * 28/09/2001 -> ajout de la date d'execution param.gene.date_exec
%   * 03/10/2001 -> ajout du flag param.profile.memoire
%   * 09/10/2001 -> ajout du flag param.gene.guido
%   * 10/10/2001 -> ajout de param.gene.creux
%   * 23/10/2001 -> ajout de source.xxx.q et source.xxx.err
%   * 05/11/2001 -> suppression de data.gene.wmag et data.gene.dwmagdt
%   * 05/11/2001 -> suppression de data.gene.sne (doublon avec nadd)
%   * 14/11/2001 -> changement definition bord.fluxXXbord
%   * 14/11/2001 -> valeur par defaut de lambda = 5/2 (et qneo non pris en compte)
%   * 13/12/2001 -> ajout de zinit_perso
%   -------------------- derniere mise a jour de la liste des variable le 21/12/2001 ---------------
%   * 21/12/2001 -> ajout de data.neo.gammae 
%   * 21/12/2001 -> ajout de data.prof.alpha 
%   * 21/01/2002 -> modification de valeurs possible de param.gene.adiabatic
%   * 22/01/2002 -> ajout des derivee temporelle 3 pts
%   * 30/01/2002 -> ajout de param.gene.force
%   * 30/01/2002 -> ajout de param.gene.nonneg
%   * 06/02/2002 -> modification des valeur par defaut de force,nonneg,fast , neomulti.ve et neomulti.vi
%   * 14/02/2002 -> ajout de param.gene.corrae
%   * 14/02/2002 -> ajout de param.gene.coefplat
%   * 14/02/2002 -> ajout de param.gene.coefbord
%   * 01/03/2002 -> valeur par defaut de psiequi =0
%   * 01/03/2002 -> ajout option gene.force = 2
%   -------------------- derniere mise a jour de la liste des variable le 11/03/2002 ---------------
%   * 14/03/2002 -> ajout des donnees pour la mhd : 
%                      * strcuture mhd
%                      * mode stabilite
%                      * fonction stabilite
%                      * consigne fonction stabilite (parametre de la fonction)
%                      * consigne stabilite (nombre de mode toroidaux)
%  * 14/03/2002 -> changement des conssignes d'asservissement (minimale)
%  * 18/03/2002 -> la structure mhd a ete mise dans equi
%  * 21/03/2002 -> ajout de la variable erreurmode dans porfile
%  * 21/03/2002 -> ajout des variables varlog1, varlog2, varlog3 dans profile
%  * 17/04/2002 -> changement de dtmin 1e-6 -> 1e-5
%  * 27/09/2002 -> changement des modules par defaut pour fci et fce
%  * 17/10/2002 -> ajout de la variable ifce
%  * 11/12/2002 -> ajout de la variable param.gene.ti_invar
%   -------------------- nouvelles variables v2.2------------------------------
%  * 01/09/2003 -> changement de definition des grandeurs de la rotation plasma
%  * 01/09/2003 -> ajout de la source ripple thermique data.source.rip
%  * 01/09/2003 -> modification de la strcuture neo
%  * 03/09/2003 -> ajout des donnees pour le module externe "ripple"
%  * 12/09/2003 -> changement de signification des fonctions mhd et de leur enchainement
%  * 12/09/2003 -> ajout de equi.mhd.deltap et equi.mhd.m_deltap
%  * 15/09/2003 -> ajout de la matrice et des etat de charge de chaque atome, passe du module d'impuretes au module neoclassique
%  * 18/09/2003 -> correction bug mode.rot defini 2 fois
%  * 22/09/2003 -> retrait de source.xx.q, ajout de source.xx.wb et source.xx.qb
%  * 22/09/2003 -> ajout de param.gene.factti et param.gene.nbforce
%  * 22/09/2003 -> ajout de data.neo.fail
%  * 10/11/2003 -> ajout du parametrage machine dependant
%  * 10/11/2003 -> ajout de la gestion de la fonction de post traitement (fonction,parametre et mode)
%  * 10/03/2004 -> ajout de l'efficacite LH dependant du temps pour zhybsimple
%  * 31/03/2004 -> ajout d'une donnees dans equi pour le calcul de jmoy 
%  * 10/05/2004 -> ajout des donnes transp (correction d'un ancien bug) 
%  * 21/06/2004 -> ajout de la donnee neo.g
%  * 21/07/2004 -> ajout de la consigne de spectre continu pour l'hybride sur TS data.cons.hybspec
% ----------------------- version 3.0 de cronos ---------------------------------------------------
%  * 14/12/2004 -> ajout de la structure simu (pour simulink) dans param
%  * 10/01/2005 -> ajout du module externe cyclo 
%  * 10/01/2005 -> supression de cyclo de la structure impur
%  * 11/01/2005 -> ajout de nombre.impur
%  * 13/01/2005 -> ajout de noe.force1, neo.force2 et neo.force3
%  * 14/01/2005 -> ajout de impur.conv, impur.fail et impur.neofail
%  * 24/01/2005 -> ajout de memoire.neo
%  * 31/01/2005 -> ajout de impur.pradsol
%  * 31/01/2005 -> ajout de impur.alpha
%  * 09/02/2005 -> ajout dans equi de df2RZ, dprRZ et frmode pour le reclacul rapide de l'equilibre
%  * 09/02/2005 -> ajout param.gene.ifus
%  * 21/03/2005 -> ajout de param.gene.computer & rapsauvemem
%  * 19/04/2006 -> ajout de data.cons.asser.qmin, data.gene.qmin & data.gene.rhoqmin
%---------------------------------------------------------------------------------------------------
%
function [cr,data,param]=zinit(mode,temps,nbrho,file,tdeb,tfin,nbfci,nbfce,nbidn,nbhyb,nbglacon,nbg)

% cr si sortie ok
cr =0;

% valeurs par defaut
nbcell         = 100;           % par defaut il est prevu 100 parametres par fonction externe
%fastonoff      = 0;             % rapide non pris en compte

% test des entrees
if nargin <6
	disp(' Nombre d''arguments incorrect');
	cr = -1 ;
	return
end

% commutation selon le mode
if isempty(mode)
	if nargin < 12
		disp(' Nombre d''arguments incorrect');
		cr = -2 ;
		return
	end
	if isempty(nbfci)
		nbfci=0;
	end
	if isempty(nbfce)
		nbfce=0;
	end
	if isempty(nbidn)
		nbidn=0;
	end
	if isempty(nbhyb)
		nbhyb=0;
	end
	if isempty(nbglacon)
		nbglacon=0;
	end
	if isempty(nbg)
		nbg=0;
	end
elseif strcmp(upper(mode),'TS')
	nbfci    = 3;
	nbfce    = 1;
	nbhyb    = 2;
	nbidn    = 1;
	nbglacon = 3;
        nbg      = 5;
elseif strcmp(upper(mode),'JET')
	disp(' Le mode Jet n''est pas encore implante');
	cr = 1 ;
	return
else
	disp(' mode inconnu');
	cr = -3 ;
	return
end

% securite sur le vecteur de temps
indbad = find(diff(temps) <=0);
if ~isempty(indbad)
	disp('negtive time step will be corrected');
	while(~isempty(indbad))
		temps(indbad) =[];
		indbad = find(diff(temps) <=0);
	end	
end

% nombre de temps
nbt =length(temps);

% parametre de la grille RZ (pour helena par defaut, redimensionner au moment de la connexion des modules)
nbrhorz     = 101;
nbthetarz   = 35;
nbsepa      = 250;
nbmode      = 128;
% pour l'equilibre a frontiere libre
nbprobes       = 1024;         % number of probes
nbgaps         = 101;          % number of gap measurements
nbcoils        = 101;          % number of independant coils currents
nbcoilelements = 101;          % number of coil elements (each coil can be compound of many elements)


% parametre du nombre de bobines poloidales (par defaut 11, redimensionner au moment de la connexion des modules)
nbpfcoils   = 11;

% temps de debut et de fin
if isempty(tdeb)
	tdeb=min(temps);
	k=1;
	kmin=1;
elseif tdeb < min(temps)
	disp('temps de debut de calcul < min(temps)');
	cr = -10 ;
	return
elseif tdeb > max(temps)
	disp('temps de debut de calcul > max(temps)');
	cr = -11 ;
	return
else
	k=min(find(temps>=tdeb));
	kmin=k;
end
	
if isempty(tfin)
	tfin=max(temps);
	kmax=nbt-1;
elseif tfin < min(temps)
	disp('temps de fin de calcul < min(temps)');
	cr = -12 ;
	return
elseif tfin > max(temps)
	disp('temps de fin de calcul > max(temps)');
	cr = -13 ;
	return
else
	kmax=max(find(temps<=tfin));
end

if kmin > kmax
	disp('tdeb > tfin')
	cr = -14;
	return
end

% theta :
% theta =linspace(0,2*pi,nbtheta);
% theta=theta(1:nbtheta);

% debut de la creation des structures
%<DEBUT_PARAM>  (balise ne pas enlever ) 


% 1 - parametres
% a - generaux   (gene)
gene.nbrho = nbrho;              % number of radial points
gene.nbt   = nbt;                % number of time points
gene.nbeq  = 7;                  % number of equations in the solver
gene.tdeb  = tdeb;               % beginning time of calculation
gene.tfin  = tfin;               % ending time of calculation
gene.t     = tdeb;               % current time
gene.dt    = 0;                  % current time variation between two indices
gene.tnp    = Inf .* ones(1,nbglacon);               % time for next pellet injection in frequency mode
gene.kmin  = kmin;               % beginning time index of calculation
gene.kmax  = kmax;               % ending time index of calculation
gene.k     = k;                  % current time index
gene.dk    = 1;                  % number of time steps calculated or added during last loop
gene.nbfci = nbfci;              % number of ICRF launchers
gene.nbfce = nbfce;              % number of ECRF launchers
gene.nbhyb = nbhyb;              % number of LH launchers
gene.nbidn = nbidn;              % number of NBI launchers
gene.nbglacon = nbglacon;        % number of pellet launchers

gene.x  = linspace(0,1,nbrho);   % normalised radial position vector : rho = rhomax(t) * x;
gene.dx = mean(diff(gene.x));    % normalised radial step
gene.nbg   = nbg;                % number of gaz in plasma composition
gene.nbrhorz   = nbrhorz;        % number of radial (rho) points in the RZ grid of the magnetic surfaces
gene.nbthetarz = nbthetarz;      % number of angular (theta) points in the RZ grid of the magnetic surfaces
gene.nbmoderz  = nbmode;         % number of coefficient in FOurrier transform of the LCMS
gene.nbsepa    = nbsepa;         % number points to define the separatrix
gene.nbpfcoils = nbpfcoils;         % number points to define the separatrix

%
% le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
% lorsque le tokamak est regarde depuis le haut
%
gene.signe.ip  = 1;              % sign of plasma current (positive for trigonometric direction when looking at the tokamak from above)
gene.signe.b0  = 1;              % sign of toroidal field (positive for trigonometric direction when looking at the tokamak from above)

%
% lambda =  permet de prendre en compte le flux de chleur convectif associe a la diffusion des particules :
%   Qe = lambda * Ge * Te ou Qi = lambda * Gi * Ti.
%   La composante neoclassique de ce flux est incluse dans les coefficient de transport neoclassique.
%   Si la composante anormale est incluse dans les coefficient de transport anormaux, alors lambda =0.
%   Sinon lambda = 3/2.
%   Si les multiplicateurs des vitesses de convections neoclassiques sont mis a zeros et que l'on veut 
%   traiter un probleme purement neoclassique, alors lambda = 5/2
%
gene.lambda   = 5/2;           % lambda factor in the contribution from particule flux to the heat fluxes {0,3/.2,5/2}
gene.modecoef = 1;             % calculation mode of the coefficients in the transport equations  0 -> all, 1 -> convective + diagonal
gene.self  = 0;                % fully self consistent mode if = 1, 0 -> coefficient + neoclassical self consistent , -1 -> like + sources at the end of each internal time step
gene.fast  = 1;                % 0 -> standard and rigorous mode , 1-> no neoclassical sources, no self consistent neoclassical coeff. 2 -> optimised mode for current diffusion
gene.source_bord = 0;          % recycling calculation control  : 1 = recycling calculated at every internal time step, 0 = as the other sources
gene.ti_invar = 0;             % with data.mode.pion = 1, 1 = Ti  preserved & 0 = Pion preserved (in equation pion = ni * ti)
gene.guido   = 1;               % if = 1, psi from current diffusion is forced to be monotonic (no negative current in the centre)
gene.qdds    = 0;              % q value for sawteeth trigger; if = 0 , no effect. if safety factor drop under qdds value, the safety factor is clamp to qdds value (this is a emulation of mean effect of sawteeth).
gene.xidds   = 0;              % value of diffusivity in the region where the safety factor is clamp for sawteeth emulation (m^2/s).
gene.force   = 1;              % 0-> does nothing, 1-> forces convergence on non-linearities if slow convergence (does not concern psi), 2-> forces convergence on non-linearities in any case (does not concern psi)
gene.nonneg   = 1;              % if =1, prevents Pe, Pion and Ne from being negative during convergence
gene.corrae   = 1;               % if = 1 imposes ae so that dni/dx has same sign as dne/dx
gene.coefplat = 5;               % if > 0 imposes Ke, Ki, Krot and D to be constant at the centre on coefplat points (allows to treat correctly the inverse problem)
gene.coefbord = 1;               % if =1 regularisation of Ke, Ki, Krot and D at the edge
gene.factti   = 30;              % link between Ti and toroidal rotation (use in zneoclass if rot undefined)
gene.nbforce  = 1;              % maximale number of iteration in radial electic field calculation : 1 = no convergence, set at 100  to used convergence mode

% attention cette variable doit etre mise a jour si changement de version :
gene.nbeq_mode = 4;            % optimisation of number of equations in the solver : 0 -> standard; 1  -> current diffusion only; {2,3,4} -> Psi, Pe, Pion, Ne; 5 -> Psi, Pe, Pion, Ne and Rotation
%
%
% le mode de fonctionnenemnt du solver :
% gene.slef = 1  ->  les equations sont resolues de maniere auto consistante pour l'ensemble des grandeurs (sources,
%                    equilibre, neoclassique , coefficient, bord, impurete et rayonnement ...)
%             0  ->  les equations sont resolues de maniere auto consistante uniquement pour les coefficients de transports
%                    et les grandeurs neoclassiques. Les sources ne sont calculees que pour les temps de la base temps donnee
%                    en entree.
%             -1 ->  comme 0 sauf que les sources sont evaluees a chaque fin de calcul des sous temps de split
%
% gene.fast  = 0 -> mode standart correct, l'equilibre neoclassique est self consistant
%              1 -> pas de sources neoclassique, ni de coefficient neo self consistante. ces garndeurs sont calculer comme les autres sources 
%                   (on commet une petite erreur proportionnelle au pas de temps de la base temps du fichier)
%              2 -> mode optimiser pour la diffusion du courant. L'equilibre neoclassique est appele qu'une seule fois en debut
%                   de boucle de convergence. Ce mode doit etre utiliser uniquement en mode diffusion de courant. Les autres champs
%                   doivent etre donnes. (on fait une petite erreur sur E.B en entree)
%
%
%
% gene.nmax = nombre de boucle maximum pour prendre en compte les non linearite pour un pas de temps elementaire
% gene.cn   = coefficient de melange implicite/ explicite du solver. Ce coefficient intervient sur la stabilite du calcul.
%             un premier appel est automatiquent fait en mode explicite (cn =1) pour avoir un premiere estimation des valeurs 
%             des coefficients (et eventuellement des sources, ...) au temps suivant avant d'appeler en mode implicite.
%             la valeur de cn par defaut est 0.5 (shema de C-N). Pour certain probleme a forte non linearite la valeur
%             peut etre diminuee. Toute valeur superieure a 0.5 est instable numeriquement. Si le probleme doit etre vu
%             comme une suite d'equilibre, on peut mettre cn = 0 (shema implicite pure).
%             si cn = -1, le solveur utilise une methode de type "exponential integrator"
%
% gene.amorti = coefficient d'amortissement des oscillations des coefficients de transport dans la boucle de convergence du solveur
%               pour les non linearite : 
%                  Dnew  = amorti * Dcalcule + (1-amorti) * Dold
%               si gene.amorti = 0  utilise la regle Dnew = 2/3 .* (Dcalcule + Dold) - (1/3) .* sign(Dcalcule + Dold) .* real(sqrt(Dcalcule * Dold));
%
% gene.psiequi = 0 -> le psi diffusion n'est pas modifie apres le calcul de l'equilibre
%                1 -> recopie la variable psi de l'equilibre dans la variable psi de la diffusion du courant
%                2 -> recopie la variable psi de l'equilibre dans la variable psi de la diffusion du courant 
%                     si il y a une grande difference en le Jmoy equilibre et diffusion sqrt(djmoy)
%                ce mode sert a rendre plus stable le couplage diffusion-equilibre 
%                ce mode peut avoir des effet sur le resultat de la diffusion du courant
%
%
gene.psiequi   = 0;              % 0,1,2 -> no copy / systematic copy / copy if different {2}
gene.adiabatic = 1;              % first evaluation of the fields at t+dt : 0 => explicit method, 1 => adiabatic calculation from equilibrium , 2 => elaborated estimation
gene.delta_adia = 0.05;          % maximum of relative adiabatic variation
gene.nmax  = 15;                 % maximum convergence loops number in solver for non-linearities (0 = no convergence)
gene.nequi_ini = 50;             % maximum convergence loops number for initial equilibrium (0 = no convergence)
gene.nequi = 50;                 % maximum convergence loops number for equilibrium (0 = no convergence)
gene.dpsi_ini = 1e-2;            % tolerance on initial psi at first time step (condition for exit of loop)
gene.djmoy    = 1e-2;            % tolerance on jmoy in solver (condition for exit of loop)
gene.creux    = 0.3;             % tunes current density in the centre of jmoy profile is too hollow (j0 = creux * mean(jmoy(jmoy >0)) [0.3] [0.01 1]
gene.cn    =  -1;               % coefficient f in solver (implicit - explicit); if cn = -1, the solver uses a metod type exponential integrator (more stable, order 2 in time)
gene.amorti =  0.5;              % coefficient for damping oscillations of transport coefficients in solver convergence loop (0.5)
gene.mjmoy  = 0.7;               % coefficient for damping oscillations of jmoy in equilibrium convergence loop {0.5/0.7}
gene.verbose =1;                 % if 1, writes information on the state of the calculation
gene.evx_inter =0;               % if 1, inserts MHD events in time base

gene.critere.ne   = 1e-3;        % convergence criterion (relative deviation) for ne
gene.critere.pe   = 1e-3;        % convergence criterion (relative deviation) for pour Pe
gene.critere.pion = 1e-3;        % convergence criterion (relative deviation) for pour Pi
gene.critere.psi  = 1e-6;        % convergence criterion (relative deviation) for psi
gene.critere.fluce    = 1e-3;    % convergence criterion (relative deviation) for fluce
gene.critere.flucion  = 1e-3;    % convergence criterion (relative deviation) for flucion
gene.critere.rot  = 1e-3;        % convergence criterion (relative deviation) for rotation

gene.file         = file;        % name of result file   (with full path)
gene.origine      = file;        % name of source file (with full path)
gene.nbsauve      = 0;           % number of time steps between two save of global result file (0 = saves the global result file only at the end of calculation)
gene.rapsauve     = '';          % name of temporary save folder + temporary save rootname, set to empty do deactivate the temporary save
gene.rebuilt      = 0;           % if 1, automatic reconstruction of result file if an error occurs during the calculation
gene.post         = 0;           % if 1, automatic run of post-treatments after the calculation

[zver,zdate]        = zinebversion;
gene.version_zineb  = zver;      % CRONOS version
gene.date_zineb     = zdate;     % date of installation of version
gene.date_exec      = NaN;       % date of simulation run
gene.computer       = '';        % name of the computer used to run CRONOS
gene.memrapsauve    = '';        % initial name of temporary save folder is case of redirection

gene.filetype       = 'source';  % file type (source ou result)

% donnees pour la gestion des fichiers cronos dans la base de donnees
gene.cronos_db.simulation_id       = '';                                 % unique identifier of the simulation file, changed at each file save.(Null = not save)
gene.cronos_db.direct_ancester_id  = '';                                 % previous identifier of the simulation file. (Null = no ancester)
gene.cronos_db.logical_ancester_id = '';                                 % for result simulation file, identifier of the source simulation file.(Null = no ancester)
gene.cronos_db.status              = 0;                                  % set to 1 when the simulation file add a reccord in database
gene.cronos_db.fail                = 0;                                  % error flag, set to non zero if a error occurred during the communication with database
gene.cronos_db.text                = 0;                                  % text of the log of the communication with database
gene.cronos_db.sim                 = struct([]);                         % sim java object contening data to reccord in the database.
gene.cronos_db.sql                 = '';                                 % text of the SQL query
gene.cronos_db.source2result       = 0;                                  % set to 1 if the function save a result after a cronos computation



% pour la generation automatique de code (redondance)
nombre.fci = nbfci;              % number of ICRF launchers
nombre.fce = nbfce;              % number of ECRF launchers
nombre.hyb = nbhyb;              % number of LH launchers
nombre.idn = nbidn;              % number of NBI launchers
nombre.impur = nbg;              % nombre d'impuretes
nombre.glacon = nbglacon;        % number of pellet launchers

% composition du gaz
%      * composition du gaz (compo)
compo.z     =  NaN .* ones(1,nbg);  % ion charge
compo.a     =  NaN .* ones(1,nbg);  % atomic mass (A.M.U.)


% nom des fonction associees au mode de calcul
% (pde)
fonction.impur      = 'zinebcompo';       % function for calculation of plasma composition + prad  (zeff, ae, nj ...) and impurity transport
fonction.equi       = 'zequi_helena';     % function for calculation of plasma equilibrium + analytical stabillity limit
fonction.neo        = 'zneo';             % function for calculation neoclassical quantities
fonction.rip        = 'zripple_therm';     % function for calculation thermal ripple sources  and other quantities

% (mhd)
fonction.mhd.dds    = 'zddscrash';        % first function for calculation of mhd reconnexion or crash, associate to evx.dds = 1, can trigger an secondary event evx.elm =1
fonction.mhd.elm    = '';                 % second function for calculation of  mhd reconnexion or crash, associate to evx.elm = 1
fonction.mhd.limite = 'zlims1';           % function for calculation of simple stability limits (sawtooth and elm threshold)
fonction.mhd.stab   = '';                 % function for calculation of stability of low toroidal number MHD modes


% (sources)
fonction.fci        = 'zfcifile';             % function for calculation of ICRF sources
fonction.fce        = 'zremafile';            % function for calculation of ECRF sources
fonction.hyb        = 'zhybsimple';       % function for calculation of LH sources
fonction.idn        = 'zsinbad2temps';          % function for calculation of NBI sources
fonction.n0         = 'zneutres';         % function for calculation of edge neutrals sources
fonction.bord       = 'zrecycle';         % function for calculation of wall and gas puff sources
fonction.glacon     = 'zglaquelc';        % function for calculation of pellet sources
fonction.fus        = 'zfusion';          % function for calculation of fusion power and alpha particules
fonction.cyclo      = 'zcytran77';         % function for calculation of cyclotronic  radiative losses


% coefficients des equation de transport
fonction.coefa      = 'zbgbs';            % function for calculation of coefficients in the transport equations "a"
fonction.coefb      = '';                 % function for calculation of coefficients in the transport equations "b"
fonction.coefc      = '';                 % function for calculation of coefficients in the transport equations "c"
fonction.coefd      = '';                 % function for calculation of coefficients in the transport equations "d"
fonction.coefe      = '';                 % function for calculation of coefficients in the transport equations "e"
fonction.coeff      = '';                 % function for calculation of coefficients in the transport equations "f"

% autres fonctions
fonction.plot      =  'zplot';       % function for plotting variables during interactive run
fonction.asser     =  '';            % function defining feedback controls
fonction.machine   =  '';            % shadow funtion use to declare device dependent parameters
fonction.post   =  '';               % post processing function name (use defautlt value if isempty, as in previus version)

% memoire pour le declenchement des calculs de sources lorsque le calcul est long
% structure pous les memoires
ms.t               = NaN;       % time of last data storage
ms.data            = [];        % storage structure
% pour impur
msi                = ms;
% etat de charge des impuretes
msi.etatcharge     = NaN .* ones(nbrho,101 .* nbg);     % state charge of impurities, internal communication impur -> neo

% impuretes
memoire.impur      = msi;        % storage structure for plasma composition (zeff, ae, nj ...) function

% bord
memoire.bord       = ms;        % storage structure for wall and gaz puff function

% (sources)
memoire.fci        = ms;        % storage structure for ICRF sources function
memoire.fce        = ms;        % storage structure for ECRF sources function
memoire.hyb        = ms;        % storage structure for LH sources function
memoire.idn        = ms;        % storage structure for NBI sources function
memoire.n0         = ms;        % storage structure for neutrals sources function
memoire.fus        = ms;        % storage structure for fusion power and alpha particules function
memoire.cyclo      = ms;        %  storage structure for cyclotronic radiative losses function
memoire.neo        = ms;	%  storage structure for bootstrap calculation (NCLASS)

% stabilite mhd
memoire.stab       = ms;        % storage structure for MHD stability function

%equilibre 
memoire.equi       = ms;        % storage structure for equilibrium function
memoire.neo        = ms;        % storage structure for neoclassical function


% mode batch pour les calculs de sources lorsque le calcul est long
% 0 -> pas de batch
% 1 -> mode batch

% impuretee
batch.impur         = 0;        % batch mode for calculation of plasma composition (zeff, ae, nj ...)

% (sources)
batch.fci        = 0;        % batch mode for calculation of ICRF sources
batch.fce        = 0;        % batch mode for calculation of ECRF sources
batch.hyb        = 0;        % batch mode for calculation of LH sources
batch.idn        = 0;        % batch mode for calculation of NBI sources
batch.n0         = 0;        % batch mode for calculation of neutrals sources
batch.fus        = 0;        % batch mode for calculation of fusion power + alpha particules sources
batch.cyclo      = 0;        %  batch mode for calculation of cyclotronic radiative losses function



% parametre des fonctions externes (cons)
%  il sont donnees dans une structure,
%  la fonction doit retourner une structure modele
%  si elle recoit comme premier argument la chaine 
%  'init', ainsi qu'une tructure defaut, une structure
%  definissant les valeurs possibles et une structure 
%  d'information : 
%  
%    [modele, defaut,possible,info]=zbidon('init');
%    
%  toutes les structures ont au moins les champs de modele. 
%  les champs de la structure possible sont compose de
%  2 sous champs :
%      .type    -> type du parametre
%      .val     -> valeurs possible :
%                    [a,b] = intervalle
%                    {a,b,c, ...} = valeurs discrete  
%  
%  la structure info a en plus les champs :
%      .nom  -> nom de la focntion
%      .label -> label utiliser dans les interface graphiques
%      .gene  -> description courte de la fonction
%      .gui   -> nom du module specifique (si besoin) de l'interface graphique
%      

% structure generique :
scg = [];


% (pde)
cons.impur      = scg;                 % parameters of plasma composition function (zeff, ae, nj ...)
cons.equi       = scg;                 % parameters of plasma equilibrium function
cons.neo        = scg;                 % parameters of neoclassical quantities function
cons.rip        = scg;                 % parameters of thermal ripple quantities function

% (mhd)
cons.mhd.dds    = scg;                 % parameters of sawteeth function
cons.mhd.elm    = scg;                 % parameters of elm function
cons.mhd.limite = scg;                 % parameters of simple elm and sawtooth stability function
cons.mhd.stab   = scg;                 % parameters of MHD stability (low toroidal mode numbers) function

% (sources)
cons.fci        = scg;                 % parameters of ICRH sources function
cons.fce        = scg;                 % parameters of ECRH sources function
cons.hyb        = scg;                 % parameters of LH sources function
cons.idn        = scg;                 % parameters of NBI sources function
cons.n0         = scg;                 % parameters of neutrals sources function
cons.bord       = scg;                 % parameters of wall and gaz puff function
cons.glacon     = scg;                 % parameters of pellet sources function
cons.fus        = scg;                 % parameters of fusion power and alpha particules function
cons.cyclo      = scg;                 % parameters of cyclotronic radiative losses function
cons.fast       = scg;                 % parameters of fast particules function


% coefficient des fonctions des coef de transport
cons.coefa        =  scg;                 % parameters of transport coefficient function "a"
cons.coefb        =  scg;                 % parameters of transport coefficient function "b"
cons.coefc        =  scg;                 % parameters of transport coefficient function "c"
cons.coefd        =  scg;                 % parameters of transport coefficient function "d"
cons.coefe        =  scg;                 % parameters of transport coefficient function "e"
cons.coeff        =  scg;                 % parameters of transport coefficient function "f"

%coefficient multiplicateur des variables neoclassiques (utilise dans les equations de transport)
cons.neomulti.ee       = 1;    % neoclassical electron thermal diffusion coefficient multiplyer
cons.neomulti.ve       = 0;    % neoclassical electron thermal convection speed multiplyer
cons.neomulti.ii       = 1;    % neoclassical ion thermal diffusion coefficient multiplyer
cons.neomulti.vi       = 0;    % neoclassical ion thermal convection speed multiplyer
cons.neomulti.nn       = 1;    % neoclassical particule (el.) diffusion coefficient multiplyer
cons.neomulti.vn       = 1;    % neoclassical particule (el.) convection speed multiplyer
cons.neomulti.rot      = 1;    % neoclassical ion velocity diffusion coefficient multiplyer
cons.neomulti.rotv     = 1;    % neoclassical ion velocity convection speed multiplyer

% autres
cons.plot         =  scg;       % parameters of plot function during interactive calculations
cons.asser        =  scg;       % parameters of feedback controls function
cons.machine      =  scg;       % parameters device dependent.
cons.post         =  scg;       % parameters of posprocessing function

% constante physique (phys)
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)

% control du plot
plot.onoff        =  0;                         % 1 if plot activated
plot.intervalle   =  1;                         % plots once every 'intervalle' time steps, uses mode.plot if 0
plot.run          =  0;                         % 1 if program is running, 0 if paused
plot.fin          =  0;                         % 1 if 'stop' button is activated
plot.pause        =  0;                         % pause between 2 calculations in s (if =0 uses drawnow)
plot.figure1      = [];                         % handle of figure 1 (control)
plot.figure2      = [];                         % handle of figure 2 (time)
plot.figure3      = [];                         % handle of figure 3 (profile)
plot.h.run        = [];                         % handle of 'stop' button
plot.h.pause      = [];                         % handle of 'pause' button
plot.h.keyboard   = [];                         % handle of 'keyboard' button
plot.h.step       = [];                         % handle of 'step by step' button in paused mode
plot.h.fin        = [];                         % handle of 'end' button
plot.h.info       = [];                         % handle of information zone
plot.order        = [];                         % default colormap

% control du decoupage en temps (split)
%
% split.equi :  
%    si la decoupe temporelle est activee (split.onoff =1), split.equi permet de controler l'appel a l'equilibre :
%     split.equi = 0  -> l'equilibre est appele a chaque sous pas de temps
%     split.equi = dt -> l'equilibre est appele une fois a la fin de chaque pas de temps et toutes les dt (s) . La derniere valeur calculer est
%                        utilisee pour les sous pas de temps intermediaires.
% 
%
split.onoff       = 1;              % time-splitting authorised if = 1
split.mode        = 1;              % time-splitting mode : 1 -> automatic on convergence criterion, 0 -> imposed
split.dtmax       = 1e-2;           % maximum time variation between two steps (s)
split.dtmin       = 1e-5;           % minimum time variation between two steps (s)
split.nb          = 3;              % number of points in the interval if imposed
split.equi        = 1e-2;           % time between two equilibrium calcuations (s)

% parametre d'echantillonage (cf. zsample.m)
%    1 - pour un signal simple :
signal.ondelette        = 1;
signal.defaut.temps     = NaN;
signal.defaut.espace    = 0;
signal.defaut.inf       = [];
signal.plus             = 0;

%    2 - pour un groupe de signaux :
groupe.ondelette         = 0;
groupe.energie           = 0.01;
groupe.defaut.temps     = NaN;
groupe.defaut.espace    = 0;
groupe.defaut.inf       = [];
groupe.plus             = 0;

% origine des donnees (from)
from.machine         = mode;        % name of server producing data
from.shot.num        = [];          % shot number
from.shot.date       = [];          % date + time of shot
from.shot.info       = {''};        % information on shot
from.creation.date   = [];          % date + time of creation of CRONOS source file
from.creation.user   = 'nobody';    % name of the user who created the CRONOS source file
from.creation.info   = {''};        % information on CRONOS source file generation (automatic)
from.creation.com    = '';          % information written by user during CRONOS source file generation
from.source.desc     = {};          % description of data sources (names of data sources)
from.source.info     = {};          % information about data used in Cronos (format variable with machine)
from.sample.signal   = signal;      % sampling parameter for time dependent data =s(time,1)
from.sample.groupe   = groupe;      % sampling parameter for time and space dependent data =g(time,space)
from.createur        = '';          % name of function used to create source file
from.option          = [];          % options of source file creation function
from.paroi.R         = [];          % R vector defining vaccum chamber shape
from.paroi.Z         = [];          % Z vector defining vaccum chamber shape

% structrures pour les asservissements  (asser)
asser.data =[];
asser.delais=[];

% intervalle de temps (pour la modification de variable de maniere groupees)
% cette structure est modifiee lors de l'edition des donnees
intervalle.temps  =[min(temps),max(temps)];    % time base limits
intervalle.calcul =[tdeb,tfin];                % beginning and ending calculation times

% structure utiliser par l'interface graphique
edit.currentfile ='';              % name of file loaded in workspace


% structure pour le "profiler"
profile.onoff      = 0;           % if 1, triggers matlab "profiler" function
profile.rapport    = 0;           % if 1, plots profiler report
profile.data       = [];          % data structure of "profiler" function
profile.memoire    = 0;           % if 1, saves every time step the "memoire" structure (temporary save must be activated)
profile.erreurmode = 0;           % if 1, saves the full CRONOS workspace if an error occurs
profile.varlog1    = '';          % name of variable 1 saved at every internal time step (temporary save must be activated)
profile.varlog2    = '';          % name of variable 2 saved at every internal time step (temporary save must be activated)
profile.varlog3    = '';          % name of variable 3 saved at every internal time step (temporary save must be activated)


% structure pour simulink
simu.onoff        = 0;            % on/off flag for simulink in Cronos : if = 1, use simulink in Cronos simulation.
simu.mdlname      = '';           % name of the simulink model used in the Cronos simulation
simu.mdltxt       = '';           % source text of the simulink model used in the Cronos simulation
simu.mdlinput     = [];           % stucture in which are saved workspace parameters used by the simulink model
simu.mdloutput    = [];           % stucture in which are saved workspace output used by the simulink model
simu.simparam     = [];           % simulink parameter
simu.mdlorigine   = '';            % original path to model file (.mdl)
simu.mdlloc       = '';           % actual path to model file (.mdl)

% regroupement
param.gene=gene;
clear gene

param.compo =compo;
clear compo

param.fonction = fonction;
clear fonction

param.memoire = memoire;
clear memoire

param.batch = batch;
clear batch

param.cons = cons;
clear cons

param.phys = phys;
clear phys

param.plot = plot;
clear plot

param.split = split;
clear split

param.from = from;
clear from

param.asser = asser;
clear asser

param.nombre = nombre;
clear nombre

param.intervalle = intervalle;
clear inter

param.edit = edit;
clear edit

param.profile = profile;
clear profile

param.simu = simu;
clear simu


%<FIN_PARAM>  (balise ne pas enlever ) 



% matrices de 0 generiques
vt    = NaN .* ones(nbt,1);
mt    = NaN .* ones(nbt,nbrho);

% 2 - data
% dans ce qui suit, les champs marque d'une "*" sont obligatoires pour le premier temps
% et, les champs marque d'une "**" sont obligatoires pour tous les temps

%<DEBUT_DATA>  (balise ne pas enlever ) 

% a - surface magnetiques et flux (magn)

%      * geometrie (geo)
% valeur de geo.mode :
%        0 = plasma symetrique, donnee de a, e1, tr1 (avec tr1 = trh1 = trb1)
%        1 = plasma asym�rique, donnee de a, e1, trh1, trb1
%        2 = plasma asym�rique, donnee de (R,Z) de la derniere surface magnetique
%
geo.mode  = vt;                   % description mode of last closed flux surface : 0 -> symmetric, 1 -> asymmetric, 2 -> RZ, 3 -> free boundary **
geo.r0    = vt;                   % major radius (m)  **
geo.z0    = vt;                   % Z shift of last closed flux surface (m)  **
geo.a     = vt;                   % minor radius (m)  **
geo.e1    = vt;                   % ellipticity of last closed flux surface  **
geo.trh1  = vt;                   % upper triangularity of last closed flux surface  **
geo.trb1  = vt;                   % lower triangularity of last closed flux surface **
geo.ind1  = vt;                   % indentation of last closed flux surface **
geo.b0    = vt;                   % toroidal field on axis (r0) without plasma (T) **
geo.R     = single(NaN .* ones(nbt,nbsepa));       % R vector defining last closed flux surface (m)
geo.Z     = single(NaN .* ones(nbt,nbsepa));       % Z defining last closed flux surface (not recentered) (m)

%      * equilibre (equi)
equi.rhomax    = vt;              % toroidal flux coordinate of last closed flux surface (m)
equi.drhomaxdt = vt;              % time derivative of rhomax (m*s^-1)
equi.phi       = mt;              % toroidal flux;
equi.phid1     = mt;              % first radial derivative of toroidal flux;
equi.phid2     = mt;              % second radial derivative of toroidal flux;
equi.dphidt    = mt;              % time derivative of toroidal flux;
equi.psi       = mt;              % poloidal flux (calculated by equilibrium for checking);
equi.rhog      = mt;              % normalised geometrical radius (r/a)
equi.d         = mt;              % Shafranov shift (m)
equi.e         = mt;              % Ellipticity
equi.a         = mt;              % minor radius of the flux surface (m)
equi.raxe      = mt;              % major radius of the centre of the flux surfaces; raxe(end) = r0_equilibre (m)
equi.zaxe      = mt;              % Z coordonate of the centre of the flux surfaces; (m)
equi.trh       = mt;              % upper triangularity
equi.trl       = mt;              % lower triangularity
equi.indh      = mt;              % upper idendation
equi.indl      = mt;              % lower idendation
equi.grho2     = mt;              % <|grad(rho)|^2>
equi.grho      = mt;              % <|grad(rho)|>
equi.vpr       = mt;              % V', <RJ>/4/pi^2 in Houlberg
equi.spr       = mt;              % S'
equi.dvprdt    = mt;              % time derivative of vpr
equi.dsprdt    = mt;              % time derivative of spr
equi.F         = mt;              % diamagnetic function F
equi.grho2r2   = mt;              % <|grad(rho)|^2/R^2>
equi.c2c       = mt;              % V'*<|gradient(rho)|^2/R^2> compute to find the correct jmoy at the edge and ip
equi.grhor     = mt;              % <|grad(rho)|/R>
equi.ri        = mt;              % <1/R>
equi.b2        = mt;              % <B^2>
equi.b2i       = mt;              % <1/B^2>
equi.psi0      = vt;              % poloidal flux at plasma centre
equi.q         = mt;              % q calculated by equilibrium for checking
equi.shear     = mt;              % magnetic shear associated to equi.q
equi.rmoy      = mt;              % <R> = <sqrt(g)>
equi.ftrap     = mt;              % trapped electrons fraction
equi.jmoy      = mt;              % current density recalculated from equilibrium
equi.ptot      = mt;              % pressure recalculated from equilibrium
equi.ip        = vt;              % plasma current from equilibrium (A)
equi.li        = vt;              % internal inductance
equi.betap     = vt;              % poloidal normalised pressure
equi.betat     = vt;              % toroidal normalised pressure
equi.betan     = vt;              % current normalised pressure
equi.q0        = vt;              % safety factor central value
equi.volume    = mt;              % plasma volume (m^3)
equi.surface   = mt;              % plasma surface (m^2)
equi.r2        = mt;              % <R^2> (m^2) 
equi.r2i       = mt;              % <1/R2> (m^-2)
equi.r2tau2    = mt;              % <R^2/tau^2> (m^-6)
equi.grho2b2   = mt;              % < |grad(rho)|^2/B^2 > (m^2.T^-2)
equi.r3tau3    = mt;              % <R^3/tau^3> (m^-3)
equi.r3tau     = mt;              % <R^3/tau> (m)
equi.bnorme    = vt;              % normalisation used in equilibrium
equi.conv      = vt;              % number of loops for convergence on jmoy ( <0 means no convergence)
equi.oscil     = vt;              % 1 if convergence exited on oscillation
equi.fail      = vt;              % 0 if convergence is ok
equi.errcur    = vt;              % final relative error on current
equi.errit     = vt;              % final relative error on equilibriu determination
equi.amix      = vt;              % final value of convergence mixing parameter



% la grille des surface magnetique
equi.R         = single(NaN .* ones(nbt,nbrhorz,nbthetarz));       % R of magnetic flux surfaces (m)
equi.Z         = single(NaN .* ones(nbt,nbrhorz,nbthetarz));       % Z of magnetic flux surfaces (m)
equi.BR        = single(NaN .* ones(nbt,nbrhorz,nbthetarz));       % Br of magnetic flux surfaces (T)
equi.BZ        = single(NaN .* ones(nbt,nbrhorz,nbthetarz));       % Bz of magnetic flux surfaces (T)
equi.BPHI      = single(NaN .* ones(nbt,nbrhorz,nbthetarz));       % Bphi of magnetic flux surfaces (T)
equi.rhoRZ     = single(NaN .* ones(nbt,nbrhorz));                 % toroidal flux coordinate of magnetic flux surfaces (m)
equi.psiRZ     = single(NaN .* ones(nbt,nbrhorz));                 % poloidal flux of magnetic flux surfaces (m)
equi.df2RZ     = single(NaN .* ones(nbt,nbrhorz));                 % FF' of magnetic flux surfaces
equi.dprRZ     = single(NaN .* ones(nbt,nbrhorz));                 % P' of magnetic flux surfaces
equi.frmode    = single(NaN .* ones(nbt,nbmode));                    % Fourrier coefficient of LCMS transform.
%     * stabilite mhd (mhd)
% les criteres sont chosi pour > 0 => instable
equi.mhd.gamma1       = vt + i .* vt;  % growth rate for mhd mode 1
equi.mhd.gamma2       = vt + i .* vt;  % growth rate for mhd mode 2
equi.mhd.gamma3       = vt + i .* vt;  % growth rate for mhd mode 3
equi.mhd.vper1        = mt;            % for mhd mode 1, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)
equi.mhd.vper2        = mt;            % for mhd mode 2, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)
equi.mhd.vper3        = mt;            % for mhd mode 3, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)
equi.mhd.ballooning   = mt;            % stability criterion of ballooning modes ((1- FM) .* (FM ~= 0))
equi.mhd.mercier      = mt;            % Mercier stability criterion ( -DR)
equi.mhd.ideal        = mt;            % Ideal mhd stability criterion ( 0.25 - DI)
equi.mhd.fail         = vt ;           % error return code : if = 0, equilibrium calculation for mishka is ok
equi.mhd.conv1        = vt ;           % for mhd mode 1, number of convergence loops (ok if < 21, verify result if = 21)
equi.mhd.conv2        = vt ;           % for mhd mode 2, number of convergence loops (ok if < 21, verify result if = 21)
equi.mhd.conv3        = vt ;           % for mhd mode 3, number of convergence loops (ok if < 21, verify result if = 21)
equi.mhd.erreur1      = vt ;           % for mhd mode 1, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)
equi.mhd.erreur2      = vt ;           % for mhd mode 2, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)
equi.mhd.erreur3      = vt ;           % for mhd mode 3, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)
%  q(k) = m(k)/n(k)
equi.mhd.deltap       = mt;            % vector of delta' values for rationnal equi.q  ( >0 instable)
equi.mhd.m_deltap     = mt;            % vector of m associate to delta', with equi.q(k) = equi.mhd.m_deltap(k) / n(k)


% donnees specifique aux equilibres a frontiere libre (on reserve la memoire pour 101 bobines et 1024 mesures magnetiques et 101 gaps)
equi.free.pfcur      =  NaN .* ones(nbt,nbcoils);           % poloidal coils currents (A)
equi.free.vext       =  NaN .* ones(nbt,nbcoils);           % poloidal coils applied voltage (V)
equi.free.probe.br   =  NaN .* ones(nbt,nbprobes);           % magnetic radial field measured at probe point (T)
equi.free.probe.bz   =  NaN .* ones(nbt,nbprobes);           % magnetic vertical field measured at probe point (T)
equi.free.probe.psi  =  NaN .* ones(nbt,nbprobes);           % magnetic vertical field measured at probe point (T)
equi.free.probe.gaps =  NaN .* ones(nbt,nbgaps);             % vector of control gaps values (m)
equi.free.strength.bmax     = NaN .* ones(nbt,nbcoilelements);             % maximal value of B on the coil (T)
equi.free.strength.fr       = NaN .* ones(nbt,nbcoilelements);             % total radial force on coil (N)
equi.free.strength.fz       = NaN .* ones(nbt,nbcoilelements);             % total vertical force on coil (N)
equi.free.strength.flux     = NaN .* ones(nbt,nbcoilelements);             % poloidal flux at the center of the coil (Wb, cronos definition) (N)

% b - donnees scalaires
% consignes   generales  = conditions au limites pour le pde (cons)
% une seule des deux consignes est prise en compte en fonction de
% la variable mode.cons.xx
cons.ip         = vt;             % plasma current (A) * {used to initialise the starting equilibrium}
cons.vloop      = vt;             % loop voltage at plasma surface (V)
cons.flux       = vt;             % edge poloidal flux (Wb)

cons.ne1        = vt;             % edge electron density (m^-3)
cons.ge1        = vt;             % edge particule flux (m^-2*s^-1)

cons.te1        = vt;             % edge electron temperature (eV)
cons.qe1        = vt;             % edge electron heat flux (W*m^-2*s^-1)
cons.pe1        = vt;             % edge electron pressure (Pa)

cons.ti1        = vt;             % edge ion temperature (eV)
cons.qi1        = vt;             % edge ion heat flux (W*m^-2*s^-1)
cons.pion1      = vt;             % edge ion pressure (Pa)

cons.ffe1       = vt;             % edge electron "turbulence" flux
cons.fe1        = vt;             % amplitude of edge electron "turbulence"

cons.ffi1       = vt;             % edge ion "turbulence" flux
cons.fi1        = vt;             % amplitude of edge ion "turbulence"

cons.frot1       = vt;            % edge rotation flux
cons.rot1        = vt;            % edge rotation

% conditition aux limites pour les impuretees
cons.c         = NaN .* ones(nbt,nbg);  % edge gas injection (atoms/s)
% pour le module de base :
% concentration de gaz majoritaire      -> 1   (D ou He)
% concentration de minoriatire 1        -> 2   (H,T)
% concentration de minoriatire 2        -> 3   (H,He)
% concentration d'impurete 1            -> 4   C
% concentration d'impurete 2            -> 5   O

% consigne de pompage (vitesse en atome /s)
cons.pomp    = vt;     % pumping speed (atoms/s)

% consignes de puissance pour les sources
%comlexe -> abs(cons) =puissance en W et angle(cons) = phase de l'antenne en radian
% sauf pour idn ou real(cons) = puissance et imag(cons) = derivee temporelle de la puissance
%% MS 25.08.08: addition of tests number>0 for initalization of each power
%  if nbfci>0
%    cons.fci     = (NaN + i .* NaN) .* ones(nbt,nbfci);   % ICRH power
%  else
%    cons.fci     = (NaN + i .* NaN) .* ones(nbt,1);   % ICRH power
%  end
%  if nbfce>0
%    cons.fce     = (NaN + i .* NaN) .* ones(nbt,nbfce);   % ECRH power
%  else
%    cons.fce     = (NaN + i .* NaN) .* ones(nbt,1);   % ECRH power
%  end
%  if nbhyb>0
%    cons.hyb     = (NaN + i .* NaN) .* ones(nbt,nbhyb);   % LH power
%  else
%    cons.hyb     = (NaN + i .* NaN) .* ones(nbt,1);   % LH power
%  end
%  if nbidn>0
%    cons.idn     = (NaN + i .* NaN) .* ones(nbt,nbidn);   % NBI power
%  else
%    cons.idn     = (NaN + i .* NaN) .* ones(nbt,1);   % NBI power
%  end
cons.fci     = (NaN + i .* NaN) .* ones(nbt,nbfci);   % ICRH power
cons.fce     = (NaN + i .* NaN) .* ones(nbt,nbfce);   % ECRH power
cons.hyb     = (NaN + i .* NaN) .* ones(nbt,nbhyb);   % LH power
cons.idn     = (NaN + i .* NaN) .* ones(nbt,nbidn);   % NBI power
cons.ext     = (NaN + i .* NaN) .* ones(nbt,1);       % external pseudo-power

% consigne d'injection de glacons
cons.glacon = NaN .* ones(nbt,nbglacon); % pellet injection

% consigne pour le zeffm quand il est donnee a la place du zeff
cons.zeffm   = vt;   % reference value for average zeff when zeffm is used (no profile given)
cons.nhnd    = (NaN + i .* NaN) .* ones(nbt,1);   % reference value for nH/nD (real part) and nT/nD (imaginary part)

% consigne pour le calcul de la stabilite mhd
cons.stab    = vt;   % number of toroidal modes in mhd stability calculations

% consignes utilisee par les asservissements
% cette partie est personnalisable avec init_perso et M a J code dynamique
cons.asser.nl0        = vt;                    % reference for central line density
cons.asser.ne0        = vt;                    % reference for central density
cons.asser.ne1        = vt;                    % reference for edge density
cons.asser.nemoy      = vt;                    % reference for averaged density
cons.asser.nbar       = vt;                    % refernce for central line avaraged density (nbar)
cons.asser.beta       = vt;                    % reference for beta
cons.asser.fgr        = vt;                    % greenwald fraction reference
cons.asser.te0        = vt;                    % reference for central electron temperature
cons.asser.te1        = vt;                    % reference for edge electron temperature
cons.asser.ti0        = vt;                    % reference for central ion temperature
cons.asser.ti1        = vt;                    % reference for edge ion temperature
cons.asser.li         = vt;                    % reference for internal inductance
cons.asser.q0         = vt;                    % reference for q0
cons.asser.qmin       = vt;                    % reference for qmin
cons.asser.c          = NaN .* ones(nbt,nbg);  % reference for gas injection (precise which gas)
cons.asser.vloop      = vt;                    % reference for vloop
cons.asser.ip         = vt;                    % reference for ip
cons.asser.pfcur      = NaN .* ones(nbt,nbpfcoils);   % references for poloidal coil current

%
% evolution temporelle de eta (pour modification de n// de l'onde lors d'un choc)
%
cons.asser.etalh = 5e18*ones(nbt,1);        % reference of LH efficiency depending on time (for zhybsimple)

%      * donnees generiques du plasma (gene)
gene.temps      = temps;          % time array **
gene.dt         = vt;             % internal time step used for the calculation of the time derivatives (s)
gene.conv       = vt;             % number of convergence loops done by the solver (<0 if no convergence)
gene.nbsplit    = vt;             % number of splits during a time step for convergence
gene.flops      = vt;             % number of floating operations in matlab
gene.cputime    = vt;             % total CPU time (s) + i * user CPU time (s)
gene.memory     = vt;             % real memory used (Mo) + i * virtual memory (Mbytes)
gene.datation   = vt;             % elapsed time since beginning of run (s)

gene.surface    = vt;             % surface of plasma cross section
gene.volume     = vt;             % plasma volume

gene.zeffm      = vt;             % line averaged zeff
gene.zeffne     = vt;             % averaged zeff, weighted by density

gene.netot      = vt;             % number of electrons in plasma
gene.nemoy      = vt;             % volume averaged electron density
gene.dnetotdt   = vt;             % time derivative of the number of electrons in plasma
gene.bilan_e    = vt;             % particule balance (electrons/s)
gene.piqne      = vt;             % electron density peaking ne0/nemoy
gene.nbar       = vt;             % nl0/2/rhomax (central line averaged density)

gene.temoy      = vt;             % volume averaged electron temperature
gene.piqte      = vt;             % electron temperature peaking te0/temoy

gene.timoy      = vt;             % volume averaged ion temperature
gene.piqti      = vt;             % ion temperature peaking ti0/timoy

gene.qmin       = vt;             % value of the minimum of the q profile
gene.rhoqmin    = vt;             % position of the minimum of the q profile

gene.we         = vt;             % thermal electrons energy
gene.wion       = vt;             % thermal ions energy
gene.wdia       = vt;             % diamagnetic energy (thermal+supra)
gene.wth        = vt;             % total thermal energy
gene.wbp        = vt;             % poloidal energy
%gene.wmag       = vt;             % energie magnetique du plama

gene.dwedt      = vt;             % time derivative of thermal electrons energy
gene.dwiondt    = vt;             % time derivative of thermal ions energy
gene.dwdiadt    = vt;             % time derivative of diamagnetic energy (thermal+supra)
gene.dwthdt     = vt;             % time derivative of total thermal energy
gene.dwbpdt     = vt;             % time derivative of poloidal energy
%gene.dwmagdt    = vt;             % derivee temporelle de l'energie magnetique du plama

gene.pel        = vt;             % power balance on electrons (additional +alpha +ohm -ei -brem -rad -en0 -synch )
gene.pion       = vt;             % power balance on ions (additional +alpha +ohm -ei -brem -rad -en0 -synch )
gene.padde      = vt;             % total heating power on electrons (additional +alpha +ohm)
gene.paddion    = vt;             % total heating power on ions (additional +alpha +ohm)
gene.paddtot    = vt;             % total heating power (additional +alpha +ohm)
gene.paddfci    = vt;             % ICRH power deposited on thermals
gene.paddhyb    = vt;             % LH power deposited on thermals
gene.paddfce    = vt;             % ECRH power deposited on thermals
gene.paddidn    = vt;             % NBI power deposited on thermals
gene.paddohm    = vt;             % ohmic heating power
gene.paddfus    = vt;             % alpha heating power
gene.paddext    = vt;             % external heat power
gene.paddrip    = vt;             % ripple heat loss
gene.pbrem      = vt;             % bremsstrahlung power loss
gene.prad       = vt;             % radiative power loss (without brem/sync)
gene.pcyclo     = vt;             % synchrotron power loss
gene.prip       = vt;             % thermal ripple power losses
gene.pn0        = vt;             % charge exchange with cold neutrals power loss
gene.nadd       = vt;             % number of electrons injected in plasma
gene.ploss      = vt;             % plasma power loss (total heating power - Pbrem, Psynch and dwdia/dt)
gene.qei        = vt;             % electron-ion exchange power due to collisions (equipartition)
gene.qneo       = vt;             % electron-ion exchange power due to neoclassical effects

%gene.sne        = vt;            % source totale d'electron
gene.sne_bord   = vt;             % electron source due to edge neutrals
gene.sne_idn    = vt;             % electron source due to NBI
gene.snions     = NaN*ones(nbt,nbg);  % source of the various ion species
gene.evolution  = NaN*ones(nbt,nbg);  % evolution of the 0d content for each specie (used for calculation of zeff0d)
gene.zeff0d     = vt;             % average zeff calculated from 0d data
gene.neutralite = vt;             % electroneutrality
gene.nions      = NaN*ones(nbt,nbg); % ion species content

gene.taue       = vt;             % electron energy confinement time
gene.tauion     = vt;             % ion energy confinement time
gene.taune      = vt;             % particules confinement time
gene.tau        = vt;             % total confinement time
gene.tauth      = vt;             % thermal energy confinement time
gene.tauj       = vt;             % current diffusion time
gene.tauei      = vt;             % equipartition time

gene.betap      = vt;             % beta poloidal
gene.beta       = vt;             % beta total
gene.betastar   = vt;             % beta *  (linked to fusion power)
gene.betadia    = vt;             % beta shafranof (diamagnetic)

gene.ip         = vt;             % plasma current (analytical formula) (A)
gene.ipepar     = vt;             % plasma current (from E/eta +Jni) (A)
gene.ipohm      = vt;             % ohmic current (E/eta) (A)
gene.ipjmoy     = vt;             % plasma current (from jmoy) (A)
gene.icd        = vt;             % additionnal driven current (A)
gene.ihyb       = vt;             % LH driven current (A)
gene.ifce       = vt;             % ECRF driven current (A)
gene.ifci       = vt;             % ICRF driven current (A)
gene.iidn       = vt;             % NBI driven current (A)
gene.ifus       = vt;             % Alpha (fusion) induce current (A)
gene.iboot      = vt;             % bootstrap current (A)
gene.ini        = vt;             % non-inductive current (add. + boot.) (A)
gene.li         = vt;             % internal inductance (H) {formula with bpol}
gene.liip       = vt;             % internal inductance (H) {formula with bpol and ip}
gene.licr       = vt;             % internal inductance (H) {formula with q}
gene.vloop      = vt;             % loop voltage on fixed loop  (V)
gene.vsurf      = vt;             % loop voltage on plasma surface (V)
gene.vres       = vt;             % average resistive loop voltage U = Pohm/Ip (V)

%      * mode de calcul des donnees
%      
%  pour les coef et les sources  :
%     
%      0  = coef ou source a 0
%      1  = coef ou source lu dans la structure d'entree
%      2  = coef ou source calcule avec le module designe par fonction.x
%      3  = report des valeurs obtenues au temps precedent (supporte uniquement par fci,idn,hyb,fce,n0 et bord; pour les autres = 1)
%
%  pour les evenements (glacon ,dds et elm)
%  
%      0  = pas  de prise en compte des evenements
%      1  = prise en compte des venements
%      
%  pour les pde  sauf ne :
%
%      0  = mis a 0
%      1  = lu dans la structure d'entree
%      2  = calcule a l'aide de l'equation differentielle
%
%  pour  pde  sur ne  :
%     
%      0  = mis a 0
%      1  = lu dans la structure d'entree
%      2  = calcule a l'aide de l'equation differentielle
%      3  = forme donnees en entree, densite totale controlee
%      4  = forme donnees en entree, densite totale controlee + piquage loi en n_gr/nbar
%      5  = forme donnees en entree, densite totale controlee + piquage loi en de H. Wiesen
%      
%  pour l'equilibre :
%      0 et 1 = lu dans la structure d'entree
%      2      = calculer
%      3      = extrapoler pour le temps k+1 du temps k     
%      
% (pde)
mode.impur      = vt;                 % caculation mode of plasma composition (zeff, ae, nj ...)
mode.psi        = vt;                 % calculation mode of poloidal flux (current diffusion)
mode.nel        = vt;                 % calculation mode of electron density
mode.pe         = vt;                 % calculation mode of electron pressure (temperature)
mode.pion       = vt;                 % calculation mode of ion pressure (temperature)
mode.equi       = vt;                 % calculation mode of magnetic equilibrium
mode.neo        = vt;                 % calculation mode neoclassical quantities
mode.rip        = zeros(size(vt));   % calculation mode thermal ripple quantities (default value =0)
mode.fluce      = vt;                 % calculation mode of electron turbulence amplitude
mode.flucion    = vt;                 % calculation mode of ion turbulence amplitude
mode.rot        = vt;                 % calculation mode of plasma rotation

% mode des consigne pour les conditions aux limites
% mode = 0 -> valeur (ne,Te et Ti) ou Ip donne
% mode = 1 -> flux (des equations sur ne, Pe et Pion)  ou Vloop donne
% mode = 2 -> valeur (Pe ou Pion) ou Psi (W)  au bord donnee 
mode.cons.psi        = vt;            % boundary condition for current diffusion : Ip (0), Vloop (1), edge psi (2), edge psi computed by a free boundary equilibrium 
mode.cons.ne         = vt;            % boundary condition for electron density equation : edge density (0), edge particule flux (1)
mode.cons.pe         = vt;            % boundary condition for electron temperature equation : edge temperature (0), edge heat flux (1), edge pressure (2)
mode.cons.pion       = vt;            % boundary condition for ion temperature equation : edge temperature (0), edge heat flux (1), edge pressure (2)
mode.cons.fluce      = vt;            % boundary condition for electron turbulence equation : edge turbulence (0), edge turbulence flux (1)
mode.cons.flucion    = vt;            % boundary condition for ion turbulence equation : edge turbulence (0), edge turbulence flux (1)
mode.cons.rot        = vt;            % boundary condition for plasma rotation equation : edge rotation (0), edge rotation flux (1)

% source des consignes pour le bord
mode.consbord.ne     = vt;            % edge electron density given as input (0) or by edge module (1)
mode.consbord.te     = vt;            % edge electron temperature given as input (0) or by edge module (1)
mode.consbord.ti     = vt;            % edge ion temperature given as input (0) or by edge module (1)
mode.consbord.fluxge = vt;            % edge electron particule flux given as input (0) or by edge module (1)
mode.consbord.fluxqe = vt;            % edge electron heat flux given as input (0) or by edge module (1)
mode.consbord.fluxqi = vt;            % edge ion heat flux given as input (0) or by edge module (1)

% (mhd)
mode.mhd.dds    = vt;                 % calculation mode of sawteeth
mode.mhd.elm    = vt;                 % calculation mode of elms
mode.mhd.limite = vt;                 % calculation mode of stability limits
mode.mhd.stab   = vt;                 % calculation mode of mhd modes stability

% (sources)
mode.fci        = vt;                 % calculation mode of H/CD sources due to ICRH
mode.fce        = vt;                 % calculation mode of H/CD sources due to ECRH
mode.hyb        = vt;                 % calculation mode of H/CD sources due to LH
mode.idn        = vt;                 % calculation mode of H/CD sources due to NBI
mode.n0         = vt;                 % calculation mode of edge neutrals source
mode.bord       = vt;                 % calculation mode of edge module
mode.glacon     = vt;                 % calculation mode of pellet deposition
mode.fus        = vt;                 % calculation mode of fusion + alpha power
mode.ohm        = vt;                 % calculation mode of ohmic power
mode.qneo       = vt;                 % calculation mode of equipartition due to ion pressure gradient
mode.qei        = vt;                 % calculation mode of ion - electron equipartition
mode.prad       = vt;                 % calculation mode of radiative power losses
mode.brem       = vt;                 % calculation mode of bremsstrahlung power lossess
mode.cyclo      = vt;                 % calculation mode of synchrotron power losses

% coefficient de transport
mode.ee        =  vt;                 % calculation mode coef in equation Pe in front of Pe (Te)
mode.ei        =  vt;                 % calculation mode coef in equation Pe in front of Pi
mode.en        =  vt;                 % calculation mode coef in equation Pe in front of ne
mode.ej        =  vt;                 % calculation mode coef in equation Pe in front of j
mode.ve        =  vt;                 % calculation mode coef in equation Pe, convection speed
mode.ep        =  vt;                 % calculation mode coef in equation Pe in front of Pe (Pe)

mode.ie        =  vt;                 % calculation mode coef in equation Pi in front of Pe
mode.ii        =  vt;                 % calculation mode coef in equation Pi in front of Pi (Ti)
mode.in        =  vt;                 % calculation mode coef in equation Pi in front of ne
mode.ij        =  vt;                 % calculation mode coef in equation Pi in front of j
mode.vi        =  vt;                 % calculation mode coef in equation Pi, convection speed
mode.ip        =  vt;                 % calculation mode coef in equation Pi in front of Pi (Pi)

mode.ne        =  vt;                 % calculation mode coef in equation Ne in front of Pe
mode.ni        =  vt;                 % calculation mode coef in equation Ne in front of Pi
mode.nn        =  vt;                 % calculation mode coef in equation Ne in front of ne
mode.nj        =  vt;                 % calculation mode coef in equation Ne in front of j
mode.vn        =  vt;                 % calculation mode coef in equation Ne, convection speed

mode.fefe      =  vt;                 % calculation mode coef in equation fluce in front of fluce
mode.fev       =  vt;                 % calculation mode coef in equation fluce, convection speed

mode.fifi      =  vt;                 % calculation mode coef in equation flucion in front of flucion
mode.fiv       =  vt;                 % calculation mode coef in equation flucion, convection speed

mode.rotc      =  vt;                 % calculation mode coef in equation rotation plasma in front of rot
mode.rotv      =  vt;                 % calculation mode coef in equation rotation plasma, convection speed

%mode de variables neo
mode.eta       = vt;    % calculation mode of resistivity
mode.jboot     = vt;    % calculation mode of bootstrap current

% autres
mode.plot          = vt;    % parametre of plot function during interactive run
mode.zeff          = vt;    % zeff calculated or given as input
mode.ae            = vt;    % ae calculated or given as input
mode.asser         = vt;    % feedback mode (0,1) = reference values given as input / 2 = reference values calculated using feedback module
mode.premiertemps  = vt;    % parameter for the initialisation of some modules
mode.post          = vt;    % post processing mode (effect depend of the value of param.gene.post)


% mode.cons.zeffm    = vt;    % zeffm est donne en consigne externe (0) ou par le modele 0d (1)
% mode.zeffm         = vt;    % zeffm est donne en consigne externe (0) ou par le modele 0d (1)

% c - profils  (prof) 
% des equations 
prof.psi        = mt;             % poloidal flux profile (Wb) *
prof.ne         = mt;             % electron density profile (m^-3) *
prof.pe         = mt;             % electron pressure profile (Pa) *
prof.pion       = mt;             % ion pressure profile (Pa) *
prof.ae         = mt;             % ni/ne profile
prof.fluce      = mt;             % electron turbulence amplitude profile (su)
prof.flucion    = mt;             % ion turbulence amplitude profile  (su)
prof.rot        = mt;             % toroidal moment rotation profile  of the plasma : sum(m_k n_k <R Vphi_k>) (kg m^-1 s^-1)

% les derivees spatiales premieres (en x pas en rho)
prof.psid1      = mt;             % 1st radial derivative of poloidal flux
prof.ned1       = mt;             % 1st radial derivative of electron density
prof.ped1       = mt;             % 1st radial derivative of electron pressure
prof.piond1     = mt;             % 1st radial derivative of ion pressure
prof.aed1       = mt;             % 1st radial derivative of ni/ne
prof.fluced1    = mt;             % 1st radial derivative of electron turbulence amplitude profile
prof.fluciond1  = mt;             % 1st radial derivative of ion turbulence amplitude profile
prof.rotd1      = mt;             % 1st radial derivative of toroidal moment rotation

% les derivees spatiales seconde  (en x pas en rho)
prof.psid2      = mt;             % 2nd radial derivative of poloidal flux
prof.ned2       = mt;             % 2nd radial derivative of electron density
prof.ped2       = mt;             % 2nd radial derivative of electron pressure
prof.piond2     = mt;             % 2nd radial derivative of ion pressure
prof.aed2       = mt;             % 2nd radial derivative ni/ne
prof.fluced2    = mt;             % 2nd radial derivative of electron turbulence amplitude profile
prof.fluciond2  = mt;             % 2nd radial derivative of ion turbulence amplitude profile
prof.rotd2      = mt;             % 2nd radial derivative of toroidal moment rotation

% les profils de Te et Ti (car souvent utilis)
prof.te         = mt;             % electron temperature profile (eV)
prof.ti         = mt;             % ion temperature profile (eV)
prof.ni         = mt;             % total ion density profile (m-3)
prof.gte        = mt;             % gradient of electron temperature (eV/m)
prof.gti        = mt;             % gradient of ion temperature  (eV/m)
prof.gne        = mt;             % gradient of electron density  (m^-2)
prof.gni        = mt;             % gradient of ion density  (m^-2)
prof.alpha      = mt;             % profile of alpha = - q^2 R grad(Beta)

% les flux
prof.flux.qe    = mt;             % electron heat flux profile (W * m^-2)
prof.flux.qi    = mt;             % ion heat flux profile    (W * m^-2)
prof.flux.ge    = mt;             % electron flux profile (m^-2 * s^-1)
prof.flux.gi    = mt;             % ion flux profile (m^-2 * s^-1)
prof.flux.fluce = mt;             % electron turbulence amplitude flux profile
prof.flux.fluci = mt;             % ion turbulence amplitude flux profile
prof.flux.rot   = mt;             % toroidal moment rotation flux profile

% les coefficients deduits
% effectif -> flux/gradient
% anormal -> corrige des coefficient neoclassiques
prof.flux.keeff = mt;             % effective electron heat diffusion coefficient
prof.flux.kieff = mt;             % effective ion heat diffusion coefficient
prof.flux.deeff = mt;             % effective electron diffusion coefficient
prof.flux.dieff = mt;             % effective ion diffusion coefficient
prof.flux.kean  = mt;             % effective anomalous electron heat diffusion coefficient
prof.flux.kian  = mt;             % effective anomalous ion heat diffusion coefficient
prof.flux.dean  = mt;             % effective anomalous electron diffusion coefficient
prof.flux.dian  = mt;             % effective anomalous ion diffusion coefficient

% flux sortant d'une surface magnetique
prof.flux.sortant.ge = mt;  % electron flux coming through magnetic surface (electrons/s)
prof.flux.sortant.gi = mt;  % ion flux coming through magnetic surface (ions/s) [ supppose totalement ionis]
prof.flux.sortant.qe = mt;  % electron heat flux coming through magnetic surface (W)
prof.flux.sortant.qi = mt;  % ion heat flux coming through magnetic surface (W)

% les derivees temporelles
prof.dpsidt     = mt;             % time derivative of the poloidal flux
prof.dnedt      = mt;             % time derivative of the electron density
prof.dpedt      = mt;             % time derivative of the electron pressure
prof.dpiondt    = mt;             % time derivative of the ion pressure
prof.daedt      = mt;             % time derivative of ni/ne
prof.dflucedt   = mt;             % time derivative of electron turbulence amplitude
prof.dfluciondt = mt;             % time derivative of ion turbulence amplitude
prof.drotdt     = mt;             % time derivative of toroidal moment rotation

% les derivees temporelles 3 pts
prof.dpsidt3p     = mt;             % time derivative of the poloidal flux (3 points)
prof.dnedt3p      = mt;             % time derivative of the electron density (3 points)
prof.dpedt3p      = mt;             % time derivative of the electron pressure(3 points)
prof.dpiondt3p    = mt;             % time derivative of the ion pressure (3 points)
prof.daedt3p      = mt;             % time derivative of ni/ne (3 points)
prof.dflucedt3p   = mt;             % time derivative of electron turbulence amplitude (su) (3 points)
prof.dfluciondt3p = mt;             % time derivative of ion turbulence amplitude (3 points)
prof.drotdt3p     = mt;             % time derivative of toroidal moment rotation (3 points)

% les profils de pressions supra et totale
% prof.psupra     = source.totale.psupra;             % composante suprathermique de la pression totale
prof.ptot       = mt;             % total pressure (including suprathermals) *

% les longueurs de gradient les plus utilisee
prof.shear      = mt;             % magnetic shear
prof.lte        = mt;             % electron temperature gradient length (Te / grad(Te))
prof.lti        = mt;             % ion temperature gradient length (Ti / grad(Ti))
prof.lne        = mt;             % electron density gradient length (ne / grad(ne))
prof.lni        = mt;             % ion density gradient length (ni / grad(ni))
prof.rhostarte  = mt;             % "rhostar" criterion for electron ITBs

% electromagnetisme
prof.dpsidt_phi = mt;             % time derivative of psi (poloidal flux) at constant phi (toroidal flux)
prof.jphi       = mt;             % effective toroidal current density (A/m^2)
prof.jpol       = mt;             % effective poloidal current density (A/m^2)
prof.jeff       = mt;             % effective current density (A/m^2)
prof.jmoy       = mt;             % current density (magnetic surface averaged)(A/m^2) *
prof.q          = mt;             % safety factor
prof.epar       = mt;             % parellel electric field E// = <E.B>/B0   (V/m)
%prof.eparlis    = mt;            % profil de E// = <E.B>/B0   lisser (lissage glissant sur tauj/10 s,V/m)
prof.bpolm      = mt;             % poloidal field (magnetic surface averaged) <|Bp|> (T)
prof.bpol       = mt;             % poloidal field (magnetic surface averaged) sqrt(<Bp.Bp>)  (T) -> relie a l'energie
prof.bphi       = mt;             % toroidal field (magnetic surface averaged) <B.grad(phi)> (T)
prof.ej         = mt;             % ohmic power

% profils divers
prof.zeff       = mt;             % effective charge
prof.xdur       = mt;             % hard-X ray profile (from diagnostic, on Tore Supra)
prof.vtor_exp   = mt;             % toroidal angular speed (rad/s)
prof.vpol_exp   = mt;             % poloidal angular speed (rad/s)
prof.mhd_cd     = mt;             % factor applied to current drive efficiency when MHD occurs (Jcd = mhd_cd * Jcd0)

% d -sources (source)
% les flux de neutrons de diverses origines
neutron.dt        = mt;    % neutron source from reaction  D + T -> He4 + n
neutron.tt        = mt;    % neutron source from reaction  T + T -> He4 + 2n
neutron.dd        = mt;    % neutron source from reaction  D + D -> He3 + ne
neutron.the3      = mt;    % neutron source from reaction  T + He3 -> He4 + p + n

% les flux de protons de diverses origines
proton.dd         = mt;    % proton source from reaction D + D -> T + p
proton.dhe3       = mt;    % proton source from reaction D + He3 -> He4 + p


% 1 - structure generique pour les sources de type chauffage
sg.el             = mt;    % power density coupled to electrons
sg.ion            = mt;    % power density coupled to ions
sg.ne             = mt;    % particle source
sg.j              = mt;    % parallel current density source ( j = <J.B>/Bo)
sg.w              = mt;    % rotation moment source <R Fk.ephi> due to heat source, toroidal component ,  unity = kg m^-1 s^-2
sg.wb             = mt;    % rotation moment source <F1k.B> due to heat source, parallel B component for NClass module ,  unity = T N m^-3
sg.q              = mt;    % extrenal heat force source <F2k.B> due to heat source, parallel B component for NClass module  ,  unity =  T N m^-3
sg.fluce          = mt;    % electron turbulence amplitude source
sg.flucion        = mt;    % ion turbulence amplitude source
sg.psupra         = mt;    % total suprathermal pressure (perpendicular) (Pa)
sg.paniso         = mt;    % anisotropic component of suprathermal pressure = p// - 1/2 * pperp (Pa)
sg.neutron        = neutron;
sg.proton         = proton;
sg.err            = vt;    % return code of source module
sg.synergie       = mt;    % synergie factor for LHCD & ECCD

% remarque :
%  Psupra  =  densite d'energie supra thermique perpendiculaire (= pperp)
%             si la pression suprathermique est isotrope, la densite totale d'energie
%             supra thermique est : Wsupra = 3/2 Psupra,
%             sinon Wsupra = 3/2 Psupra + Paniso.
%  Paniso  =  anisotropie de la pression suprathermique dans la direction // a B
%             Paniso = Wsupra - 3/2 * Psupra.
%             soit encore en terme de pression // :
%             Paniso = W// - 1/2 * Psupra = p// - 1/2 * pperp.

% 2 - les differentes 

% les chauffages

source.fci        = sg;    % ICRF sources
source.fce        = sg;    % ECRF sources
source.idn        = sg;    % NBI sources
source.hyb        = sg;    % LH sources
source.fus        = sg;    % sources due to alpha particles
source.rip        = sg;    % sources due to ripple losses of  thermal particules
source.ext        = sg;    % other external sources

% sources impuretees 
source.n0         =  sg;                             % sources due to plasma interaction with neutrals coming from the edge
source.n.n0       =  NaN .* ones(nbt,nbrho,nbg);     % sources of matter per specy due to cold neutrals (atoms)
source.n.fus      =  NaN .* ones(nbt,nbrho,nbg);     % sources of matter per specy due to fusion reaction (atoms)
source.n.idn      =  NaN .* ones(nbt,nbrho,nbg);     % sources of matter per specy due to NBI (atoms)
source.n.rip      =  NaN .* ones(nbt,nbrho,nbg);     % sources of matter per specy due to ripple (atoms)
source.n.n0th     =  mt;                             % thermal neutrals density coming from the edge (equivalent electron)

% sources rayonnements
source.prad       = mt;    % radiated power source/sink
source.brem       = mt;    % bremsstrahlung power source/sink
source.cyclo      = mt;    % cyclotron power source/sink

% sources collisions
source.ohm        = mt;    % ohmic power deposition
source.qei        = mt;    % electron-ion heat exchange (equipartition)
source.qneo       = mt;    % neoclassical electron-ion heat exchange
source.jboot      = mt;    % bootstrap current density (A/m^2)

% turbulence
source.fluce      = sg;
source.flucion    = sg;

% sommation des termes sources
source.totale     = sg;    % sum of all sources

% les coefficents (coef) 
% coefficient de transport
coef.eta        = mt;                 % resistivity (ohm * m)

coef.ee        =  mt;                 % coefficient of Pe (Te) in Pe equation  (m^-1 * s^-1)
coef.ei        =  mt;                 % coefficient of Pi in Pe equation  (m^-1 * s^-1)
coef.en        =  mt;                 % coefficient of ne in Pe equation  (m^-1 * s^-1)
coef.ej        =  mt;                 % coefficient of j in Pe equation   (?)
coef.ve        =  mt;                 % convective speed in Pe equation   (m*s^-1)
coef.ep        =  mt;                 % coefficient of Pe (Pe) in Pe equation (m^-1 * s^-1)

coef.ie        =  mt;                 % coefficient of Pe in Pi equation  (m^-1 * s^-1)
coef.ii        =  mt;                 % coefficient of Pi (Ti) in Pi equation (m^-1 * s^-1)
coef.in        =  mt;                 % coefficient of ne in Pi equation  (m^-1 * s^-1)
coef.ij        =  mt;                 % coefficient of j in Pi equation     (?)
coef.vi        =  mt;                 % convective speed in Pi equation (m*s^-1)
coef.ip        =  mt;                 % coefficient of Pi (Pi) in Pi equation (m^-1 * s^-1)

coef.ne        =  mt;                 % coefficient of Pe in ne equation   (m^2 * s^-1)
coef.ni        =  mt;                 % coefficient of Pi in ne equation   (m^2 * s^-1)
coef.nn        =  mt;                 % coefficient of ne in ne equation   (m^2 * s^-1)
coef.nj        =  mt;                 % coefficient of j in ne equation     (?)
coef.vn        =  mt;                 % convective speed in ne equation (m*s^-1)

coef.fev       =  mt;                 % convective speed in fluce equation
coef.fefe      =  mt;                 % diffusion coefficient in fluce equation

coef.fiv       =  mt;                 % convective speed in flucion equation
coef.fifi      =  mt;                 % diffusion coefficient in flucion equation

coef.rotv      =  mt;                 % convective speed in toroidal rotation equation
coef.rot       =  mt;                 % diffusion coefficient in toroidal rotation equation

% les grandeur neoclassique (neo)
neo.eta           = mt;    % neoclassical resistivity (ohm * m)
neo.jboot         = mt;    % bootstrap current density (A/m^2)
neo.coef.ee       = mt;    % neoclassical electron heat diffusivity  (m/s)
neo.coef.ve       = mt;    % neoclassical electron heat convection
neo.coef.ii       = mt;    % neoclassical ion heat diffusivity
neo.coef.vi       = mt;    % neoclassical ion heat convection
neo.coef.nn       = mt;    % neoclassical particle diffusivity
neo.coef.vn       = mt;    % neoclassical particle convection
neo.coef.rotv     = mt;    % neoclassical rotation convection (m s^-1)
neo.coef.rot      = mt;    % neoclassical rotation diffusivity (m^-1 s^-1)
neo.vtor          = NaN .* ones(nbt,nbrho,nbg);  % neoclassical toroidal rotation speed (m/s), at R = Rmax of each flux surface, for each ion specy
neo.vtheta        = NaN .* ones(nbt,nbrho,nbg);  % neoclassical poloidal rotation speed (m/s), at R = Rmax of each flux surface, for each ion specy
neo.mach          = NaN .* ones(nbt,nbrho,nbg);  % neoclassical mach number, at R = Rmax of each flux surface, for each ion specy
neo.er            = mt;    % neoclassical radial electric field (V/m) = Er / gradient(rho).
neo.gammae        = mt;    % neoclassical rotation rotation shearing (s^-1), at R = Rmax of each flux surface.
neo.g             = mt;    % this parameter measures the magnitude of the // friction force relative to the // pressure gradient; abs(G)<< 1 to have valid NClass results.
neo.flux.ne       = mt;    % neoclassical electron flux
neo.flux.nion     = mt;    % neoclassical ion flux (sum over species)
neo.flux.qe       = mt;    % neoclassical electron heat flux
neo.flux.qion     = mt;    % neoclassical ion heat flux (sum over species)
neo.qei           = mt;    % electron-ion heat exchange (equipartition)
neo.qneo          = mt;    % electron-ion heat exchange (neoclassical)
neo.w_ion         = NaN .* ones(nbt,nbrho,nbg); % 2*pi * plasma solid rotation frequency in toroidal direction , for each ion specy
neo.w_e           = mt;    % 2*pi * plasma electron solid rotation frequency in toroidal direction
neo.utheta_i      = NaN .* ones(nbt,nbrho,nbg); % fluid velocity, theta component :<V_k . theta> / <B . theta> , for each ion specy (cf NCLASS)
neo.utheta_e      = mt;    % electron fluid velocity, theta component :<V_e . theta> / <B . theta>  (cf NCLASS)
%  le moment d'ordre 1 est la quantite de mouvement (v), 2 -> la chaleur (v^3), 3- > en v^5
%  ref : W. A. Houlberg, Phys. PLasmas 4 (9), september 1997, p 3230- (eq 4 et suivantes)
%  l'usage normale : fexizpr1 et fexizpr2 doivent etre fournis pour IDN et fexizpr3 = 0;
%  fexizpr1 = datak.source.totale.wb
%  fexizpr2 = datak.source.totale.q
neo.force1        = mt;    % first moment of external parallel force on species i,z (T*j/m**3)
neo.force2        = mt;    % second moment of external parallel force on species i,z (T*j/m**3)
neo.force3        = mt;    % third moment of external parallel force on species i,z (T*j/m**3)
neo.fail          = vt;    % convergence failure flag ( force ball ou Nclass)

% les donnees d'un run transp
transp.jnbtrx   = mt;    % neutral beam current
transp.jbstrx   = mt;    % bootstrap current
transp.jtottrx  = mt;    % total current
transp.titrx    = mt;    % ion temperature
transp.tetrx    = mt;    % electron temperature
transp.netrx    = mt;    % electron density
transp.nitrx    = mt;    % ion density
transp.qohtrx   = mt;    % ohmic heating source
transp.qrfetrx  = mt;    % electron heating source due to RF
transp.qrfitrx  = mt;    % ion heating source due to RF
transp.qnbetrx  = mt;    % electron heating source due to NBI
transp.qnbitrx  = mt;    % ion heating source due to NBI
transp.qbthtrx  = mt;    %
transp.zefftrx  = mt;    % Zeff profile
transp.rotortrx = mt;    % toroidal velocity
transp.chietrx  = mt;    % transport coefficient (ion-electron)
transp.chiitrx  = mt;    % transport coefficient (ion-ion)
transp.qietrx   = mt;    %  electron-ion heat exchange (equipartition)
transp.qtrx     = mt;    % q profile
transp.ttransp  = vt;    % slice time of transp run
% les impuretees
impur.prad       = mt;    % radiated power
impur.brem       = mt;    % bremsstrahlung
%impur.cyclo      = mt;    % cyclotron radiation
impur.zeff       = mt;    % Zeff
impur.impur      = NaN .* ones(nbt,nbrho,nbg);  % impurities density
impur.ae         = mt;    % ni/ne
impur.conv       = vt;    % number of convergence loops done with the impurities transport code and the neoclassical module to obtain consitentcy.
impur.fail       = vt;    % convergence failure flag (impurities or consitency with neoclassical code)
impur.neofail    = vt;    %convergence failure flag (Nclass)
impur.pradsol    = vt;    % radiated power in the SOL (W)
impur.alpha      = mt;    % electroneutrality normalisation a the ouptut of impurities transport module

% le bord (structure remplie par le module de bord)
% attention cette structure ne doit pas avoir de sous structure
% sinon il faut modifier zbord.m !
bord.temp         = vt;    % wall temperature (K)
bord.nb           = vt;    % number of atoms in wall (electron equivalent)
bord.pref         = vt;    % power loss due to cooling (W)
bord.pchauffe     = vt;    % heating power (W)
bord.reflex       = vt;    % recycling coefficient
bord.fluxpompe    = vt;    % pumped particle flux (electrons * s^-1)
bord.fluxplasma   = vt;    % particle flux coming from plasma (electrons *  s^-1)
bord.fluxgazout   = NaN*ones(nbt,nbg);    % particle flux coming from plasma, for each specy (atoms  * s^-1)
bord.fluxgaz      = NaN*ones(nbt,nbg);    % fuelling, for each specy (atoms * s^-1)
bord.fluxmur_c    = vt;    % hot neutral flux recycled towards the plasma (electrons * s^-1)
bord.fluxmur_f    = vt;    % cold neutral flux from the wall towards the plasma (fuelling +degasing + sputtering, electrons * s^-1)
bord.fluxgazin_c  = NaN*ones(nbt,nbg);    % hot neutral flux recycled towards the plasma, for each specy (atoms * s^-1)
bord.fluxgazin_f  = NaN*ones(nbt,nbg);    % cold neutral flux from the wall towards the plasma, for each specy (fuelling +degasing + sputtering, atoms * s^-1)
bord.nebord       = vt;    % edge electron density (depending on limiter parameters) (m^-3)
bord.nibord       = NaN*ones(nbt,nbg);    % edge ion density, for each specy (depending on limiter parameters) (m^-3)
bord.tebord       = vt;    % edge electron temperature (depending on limiter parameters) (eV)
bord.tibord       = vt;    % edge ion temperature (depending on limiter parameters) (eV)
bord.fluxgebord   = vt;    % edge particle flux (electrons * s^-1)
bord.fluxqebord   = vt;    % edge electron heat flux (W)
bord.fluxqibord   = vt;    % edge ion heat flux(W)

% description detaillee de la SOL
bord.x            = vt * linspace(1,1.5,51);      % extended normalized coordonnate param.gene.x in the SOL
bord.ne           = NaN .* bord.x;    % electron density profile in the SOL (m^-3)
bord.ni           = NaN .* bord.x;    % ion density profile in the SOL (m^-3)
bord.te           = NaN .* bord.x;    % electron temperature profile in the SOL (ev)
bord.ti           = NaN .* bord.x;    % ion temperature profile in the SOL (ev)
bord.impur        = NaN .* ones(param.gene.nbt,size(bord.x,2),param.gene.nbg);    % ions and impurities density  in the SOL (m^-3)

% les evennents (evx)
evx.dds      = vt;      % sawtooth event triggering
evx.elm      = vt;      % elm event triggering


% les donnees experimentale pour le plot
exp.ne0      = vt;     % measured central electron density (m^-3)
exp.nea      = vt;     % measured edge electron density (m^-3)
exp.te0      = vt;     % measured central electron temperature(eV)
exp.tea      = vt;     % measured edge electron temperature (eV)
exp.ni0      = vt;     % measured central ion density (m^-3)
exp.nia      = vt;     % measured edge ion density (m^-3)
exp.ti0      = vt;     % measured central ion temperature (eV)
exp.tia      = vt;     % measured edge ion temperature (eV)
exp.j0       = vt;     % measured central current density (A*m^-2)
exp.ip       = vt;     % measured plasma current (A)
exp.li       = vt;     % measured internal inductance
exp.vloop    = vt;     % measured loop voltage (V)
exp.q0       = vt;     % measured central safety factor
exp.qa       = vt;     % measured edge safety factor
exp.betadia  = vt;     % measured diamagnetic beta
exp.hybspec.npar = NaN .* ones(nbt,2,nbhyb);      % npar of continuous LH spectrum
exp.hybspec.pow  = NaN .* ones(nbt,2,nbhyb);      % power of continuous LH spectrum


% regroupement
data.geo = geo;
clear geo

data.equi = equi;
clear equi

% data.mhd = mhd;
% clear mhd

data.cons = cons;
clear cons

data.gene = gene;
clear gene

data.mode = mode;
clear mode

data.prof = prof;
clear prof

data.source = source;
clear source

data.coef = coef;
clear coef

data.neo = neo;
clear neo

data.transp = transp;
clear transp

data.impur = impur;
clear impur

data.bord = bord;
clear bord

data.evx = evx;
clear evx

data.exp = exp;
clear exp

%<FIN_DATA>  (balise ne pas enlever ) 

% structure vide pour le fichier
post =[];


% appel de la fonction de personalisation de la structure des donnees
[data,param] = zinit_perso(data,param);

if ~isempty(file) & (nargout < 3)
	try
	   % compactage des donnees
           data=zreduit(param,data,'compact');
           % sauvegarde
 		if  verLessThan('matlab','7.0')
					savedb(file,'param','data','post','-V6');
		else
					savedb(file,'param','data','post');
		end
           %save(file,'param','data','post');
           % compression du fichier
           zgzip(file,'compress');
           %save(file,'param','data');
           data=zreduit(param,data,'uncompact');
	catch
	   cr =-100
	   disp('error during the saving')
	   return
	end
end
