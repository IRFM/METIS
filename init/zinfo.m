% ZINFO fonction creant la structure info
%----------------------------------------------------------------
% fichier zinfo.m ->  zinfo
%
%
% fonction Matlab 5 :
%
% Cree la structure info de zineb
%
% syntaxe  :
% 
%   info=zinfo;
%      -> retourne l'aide complete (celle des modules externe aussi) 
%
%   info=zinfo('noconnexion');
%      -> retourne l'aide incomplete (ne lit pas celle des modules externes) 
% 
% sortie :
% 
%     info    =  structure info
% 
%
% fonction gfenere par zcreeinfo
%
%--------------------------------------------------------------
%
function info=zinfo(noconnexion)

if nargin < 1
	noconnexion ='';
end



                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% 1 - parametres
% a - generaux   (gene)
gene.nbrho =' number of radial points';
gene.nbt   =' number of time points';
gene.nbeq  =' number of equations in the solver';
gene.tdeb  =' beginning time of calculation';
gene.tfin  =' ending time of calculation';
gene.t     =' current time';
gene.dt    =' current time variation between two indices';
gene.tnp    =' time for next pellet injection in frequency mode';
gene.kmin  =' beginning time index of calculation';
gene.kmax  =' ending time index of calculation';
gene.k     =' current time index';
gene.dk    =' number of time steps calculated or added during last loop';
gene.nbfci =' number of ICRF launchers';
gene.nbfce =' number of ECRF launchers';
gene.nbhyb =' number of LH launchers';
gene.nbidn =' number of NBI launchers';
gene.nbglacon =' number of pellet launchers';
                                                                                                                                                                                                                     
gene.x  =' normalised radial position vector : rho = rhomax(t) * x;';
gene.dx =' normalised radial step';
gene.nbg   =' number of gaz in plasma composition';
gene.nbrhorz   =' number of radial (rho) points in the RZ grid of the magnetic surfaces';
gene.nbthetarz =' number of angular (theta) points in the RZ grid of the magnetic surfaces';
gene.nbmoderz  =' number of coefficient in FOurrier transform of the LCMS';
gene.nbsepa    =' number points to define the separatrix';
gene.nbpfcoils =' number points to define the separatrix';
                                                                                                                                                                                                                     
%
% le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
% lorsque le tokamak est regarde depuis le haut
%
gene.signe.ip  =' sign of plasma current (positive for trigonometric direction when looking at the tokamak from above)';
gene.signe.b0  =' sign of toroidal field (positive for trigonometric direction when looking at the tokamak from above)';
                                                                                                                                                                                                                     
%
% lambda =  permet de prendre en compte le flux de chleur convectif associe a la diffusion des particules :
%   Qe = lambda * Ge * Te ou Qi = lambda * Gi * Ti.
%   La composante neoclassique de ce flux est incluse dans les coefficient de transport neoclassique.
%   Si la composante anormale est incluse dans les coefficient de transport anormaux, alors lambda =0.
%   Sinon lambda = 3/2.
%   Si les multiplicateurs des vitesses de convections neoclassiques sont mis a zeros et que l''on veut
%   traiter un probleme purement neoclassique, alors lambda = 5/2
%
gene.lambda   =' lambda factor in the contribution from particule flux to the heat fluxes {0,3/.2,5/2}';
gene.modecoef =' calculation mode of the coefficients in the transport equations  0 -> all, 1 -> convective + diagonal';
gene.self  =' fully self consistent mode if = 1, 0 -> coefficient + neoclassical self consistent , -1 -> like + sources at the end of each internal time step';
gene.fast  =' 0 -> standard and rigorous mode , 1-> no neoclassical sources, no self consistent neoclassical coeff. 2 -> optimised mode for current diffusion';
gene.source_bord =' recycling calculation control  : 1 = recycling calculated at every internal time step, 0 = as the other sources';
gene.ti_invar =' with data.mode.pion = 1, 1 = Ti  preserved & 0 = Pion preserved (in equation pion = ni * ti)';
gene.guido   =' if = 1, psi from current diffusion is forced to be monotonic (no negative current in the centre)';
gene.qdds    =' q value for sawteeth trigger; if = 0 , no effect. if safety factor drop under qdds value, the safety factor is clamp to qdds value (this is a emulation of mean effect of sawteeth).';
gene.xidds   =' value of diffusivity in the region where the safety factor is clamp for sawteeth emulation (m^2/s).';
gene.force   =' 0-> does nothing, 1-> forces convergence on non-linearities if slow convergence (does not concern psi), 2-> forces convergence on non-linearities in any case (does not concern psi)';
gene.nonneg   =' if =1, prevents Pe, Pion and Ne from being negative during convergence';
gene.corrae   =' if = 1 imposes ae so that dni/dx has same sign as dne/dx';
gene.coefplat =' if > 0 imposes Ke, Ki, Krot and D to be constant at the centre on coefplat points (allows to treat correctly the inverse problem)';
gene.coefbord =' if =1 regularisation of Ke, Ki, Krot and D at the edge';
gene.factti   =' link between Ti and toroidal rotation (use in zneoclass if rot undefined)';
gene.nbforce  =' maximale number of iteration in radial electic field calculation : 1 = no convergence, set at 100  to used convergence mode';
                                                                                                                                                                                                                     
% attention cette variable doit etre mise a jour si changement de version :
gene.nbeq_mode =' optimisation of number of equations in the solver : 0 -> standard; 1  -> current diffusion only; {2,3,4} -> Psi, Pe, Pion, Ne; 5 -> Psi, Pe, Pion, Ne and Rotation';
%
%
% le mode de fonctionnenemnt du solver :
% gene.slef = 1  ->  les equations sont resolues de maniere auto consistante pour l''ensemble des grandeurs (sources,
%                    equilibre, neoclassique , coefficient, bord, impurete et rayonnement ...)
%             0  ->  les equations sont resolues de maniere auto consistante uniquement pour les coefficients de transports
%                    et les grandeurs neoclassiques. Les sources ne sont calculees que pour les temps de la base temps donnee
%                    en entree.
%             -1 ->  comme 0 sauf que les sources sont evaluees a chaque fin de calcul des sous temps de split
%
% gene.fast  = 0 -> mode standart correct, l''equilibre neoclassique est self consistant
%              1 -> pas de sources neoclassique, ni de coefficient neo self consistante. ces garndeurs sont calculer comme les autres sources
%                   (on commet une petite erreur proportionnelle au pas de temps de la base temps du fichier)
%              2 -> mode optimiser pour la diffusion du courant. L''equilibre neoclassique est appele qu''une seule fois en debut
%                   de boucle de convergence. Ce mode doit etre utiliser uniquement en mode diffusion de courant. Les autres champs
%                   doivent etre donnes. (on fait une petite erreur sur E.B en entree)
%
%
%
% gene.nmax = nombre de boucle maximum pour prendre en compte les non linearite pour un pas de temps elementaire
% gene.cn   = coefficient de melange implicite/ explicite du solver. Ce coefficient intervient sur la stabilite du calcul.
%             un premier appel est automatiquent fait en mode explicite (cn =1) pour avoir un premiere estimation des valeurs
%             des coefficients (et eventuellement des sources, ...) au temps suivant avant d''appeler en mode implicite.
%             la valeur de cn par defaut est 0.5 (shema de C-N). Pour certain probleme a forte non linearite la valeur
%             peut etre diminuee. Toute valeur superieure a 0.5 est instable numeriquement. Si le probleme doit etre vu
%             comme une suite d''equilibre, on peut mettre cn = 0 (shema implicite pure).
%
% gene.amorti = coefficient d''amortissement des oscillations des coefficients de transport dans la boucle de convergence du solveur
%               pour les non linearite :
%                  Dnew  = amorti * Dcalcule + (1-amorti) * Dold
%
% gene.psiequi = 0 -> le psi diffusion n''est pas modifie apres le calcul de l''equilibre
%                1 -> recopie la variable psi de l''equilibre dans la variable psi de la diffusion du courant
%                2 -> recopie la variable psi de l''equilibre dans la variable psi de la diffusion du courant
%                     si il y a une grande difference en le Jmoy equilibre et diffusion sqrt(djmoy)
%                ce mode sert a rendre plus stable le couplage diffusion-equilibre
%                ce mode peut avoir des effet sur le resultat de la diffusion du courant
%
%
gene.psiequi   =' 0,1,2 -> no copy / systematic copy / copy if different {2}';
gene.adiabatic =' first evaluation of the fields at t+dt : 0 => explicit method, 1 => adiabatic calculation from equilibrium , 2 => elaborated estimation';
gene.delta_adia =' maximum of relative adiabatic variation';
gene.nmax  =' maximum convergence loops number in solver for non-linearities (0 = no convergence)';
gene.nequi_ini =' maximum convergence loops number for initial equilibrium (0 = no convergence)';
gene.nequi =' maximum convergence loops number for equilibrium (0 = no convergence)';
gene.dpsi_ini =' tolerance on initial psi at first time step (condition for exit of loop)';
gene.djmoy    =' tolerance on jmoy in solver (condition for exit of loop)';
gene.creux    =' tunes current density in the centre of jmoy profile is too hollow (j0 = creux * mean(jmoy(jmoy >0)) [0.3] [0.01 1]';
gene.cn    =' coefficient f in solver (implicit - explicit)';
gene.amorti =' coefficient for damping oscillations of transport coefficients in solver convergence loop (0.5)';
gene.mjmoy  =' coefficient for damping oscillations of jmoy in equilibrium convergence loop {0.5/0.7}';
gene.verbose =' if 1, writes information on the state of the calculation';
gene.evx_inter =' if 1, inserts MHD events in time base';
                                                                                                                                                                                                                     
gene.critere.ne   =' convergence criterion (relative deviation) for ne';
gene.critere.pe   =' convergence criterion (relative deviation) for pour Pe';
gene.critere.pion =' convergence criterion (relative deviation) for pour Pi';
gene.critere.psi  =' convergence criterion (relative deviation) for psi';
gene.critere.fluce    =' convergence criterion (relative deviation) for fluce';
gene.critere.flucion  =' convergence criterion (relative deviation) for flucion';
gene.critere.rot  =' convergence criterion (relative deviation) for rotation';
                                                                                                                                                                                                                     
gene.file         =' name of result file   (with full path)';
gene.origine      =' name of source file (with full path)';
gene.nbsauve      =' number of time steps between two save of global result file (0 = saves the global result file only at the end of calculation)';
gene.rapsauve     =' name of temporary save folder + temporary save rootname, set to empty do deactivate the temporary save';
gene.rebuilt      =' if 1, automatic reconstruction of result file if an error occurs during the calculation';
gene.post         =' if 1, automatic run of post-treatments after the calculation';
                                                                                                                                                                                                                     
[zver,zdate]        = zinebversion;
gene.version_zineb  =' CRONOS version';
gene.date_zineb     =' date of installation of version';
gene.date_exec      =' date of simulation run';
gene.computer       =' name of the computer used to run CRONOS';
gene.memrapsauve    =' initial name of temporary save folder is case of redirection';
                                                                                                                                                                                                                     
gene.filetype       =' file type (source ou result)';
                                                                                                                                                                                                                     
% donnees pour la gestion des fichiers cronos dans la base de donnees
gene.cronos_db.simulation_id       =' unique identifier of the simulation file, changed at each file save.(Null = not save)';
gene.cronos_db.direct_ancester_id  =' previous identifier of the simulation file. (Null = no ancester)';
gene.cronos_db.logical_ancester_id =' for result simulation file, identifier of the source simulation file.(Null = no ancester)';
gene.cronos_db.status              =' set to 1 when the simulation file add a reccord in database';
gene.cronos_db.fail                =' error flag, set to non zero if a error occurred during the communication with database';
gene.cronos_db.text                =' text of the log of the communication with database';
gene.cronos_db.sim                 =' sim java object contening data to reccord in the database.';
gene.cronos_db.sql                 =' text of the SQL query';
gene.cronos_db.source2result       =' set to 1 if the function save a result after a cronos computation';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% pour la generation automatique de code (redondance)
nombre.fci =' number of ICRF launchers';
nombre.fce =' number of ECRF launchers';
nombre.hyb =' number of LH launchers';
nombre.idn =' number of NBI launchers';
nombre.impur =' nombre d''impuretes';
nombre.glacon =' number of pellet launchers';
                                                                                                                                                                                                                     
% composition du gaz
%      * composition du gaz (compo)
compo.z     =' ion charge';
compo.a     =' atomic mass (A.M.U.)';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% nom des fonction associees au mode de calcul
% (pde)
fonction.impur      =' function for calculation of plasma composition + prad  (zeff, ae, nj ...) and impurity transport';
fonction.equi       =' function for calculation of plasma equilibrium + analytical stabillity limit';
fonction.neo        =' function for calculation neoclassical quantities';
fonction.rip        =' function for calculation thermal ripple sources  and other quantities';
                                                                                                                                                                                                                     
% (mhd)
fonction.mhd.dds    =' first function for calculation of mhd reconnexion or crash, associate to evx.dds = 1, can trigger an secondary event evx.elm =1';
fonction.mhd.elm    =' second function for calculation of  mhd reconnexion or crash, associate to evx.elm = 1';
fonction.mhd.limite =' function for calculation of simple stability limits (sawtooth and elm threshold)';
fonction.mhd.stab   =' function for calculation of stability of low toroidal number MHD modes';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% (sources)
fonction.fci        =' function for calculation of ICRF sources';
fonction.fce        =' function for calculation of ECRF sources';
fonction.hyb        =' function for calculation of LH sources';
fonction.idn        =' function for calculation of NBI sources';
fonction.n0         =' function for calculation of edge neutrals sources';
fonction.bord       =' function for calculation of wall and gas puff sources';
fonction.glacon     =' function for calculation of pellet sources';
fonction.fus        =' function for calculation of fusion power and alpha particules';
fonction.cyclo      =' function for calculation of cyclotronic  radiative losses';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% coefficients des equation de transport
fonction.coefa      =' function for calculation of coefficients in the transport equations "a"';
fonction.coefb      =' function for calculation of coefficients in the transport equations "b"';
fonction.coefc      =' function for calculation of coefficients in the transport equations "c"';
fonction.coefd      =' function for calculation of coefficients in the transport equations "d"';
fonction.coefe      =' function for calculation of coefficients in the transport equations "e"';
fonction.coeff      =' function for calculation of coefficients in the transport equations "f"';
                                                                                                                                                                                                                     
% autres fonctions
fonction.plot      =' function for plotting variables during interactive run';
fonction.asser     =' function defining feedback controls';
fonction.machine   =' shadow funtion use to declare device dependent parameters';
fonction.post   =' post processing function name (use defautlt value if isempty, as in previus version)';
                                                                                                                                                                                                                     
% memoire pour le declenchement des calculs de sources lorsque le calcul est long
% structure pous les memoires
ms.t               =' time of last data storage';
ms.data            =' storage structure';
% pour impur
msi                = ms;
% etat de charge des impuretes
msi.etatcharge     =' state charge of impurities, internal communication impur -> neo';
                                                                                                                                                                                                                     
% impuretes
memoire.impur      =' storage structure for plasma composition (zeff, ae, nj ...) function';
                                                                                                                                                                                                                     
% bord
memoire.bord       =' storage structure for wall and gaz puff function';
                                                                                                                                                                                                                     
% (sources)
memoire.fci        =' storage structure for ICRF sources function';
memoire.fce        =' storage structure for ECRF sources function';
memoire.hyb        =' storage structure for LH sources function';
memoire.idn        =' storage structure for NBI sources function';
memoire.n0         =' storage structure for neutrals sources function';
memoire.fus        =' storage structure for fusion power and alpha particules function';
memoire.cyclo      ='  storage structure for cyclotronic radiative losses function';
memoire.neo        = ms;	%  storage structure for bootstrap calculation (NCLASS)
                                                                                                                                                                                                                     
% stabilite mhd
memoire.stab       =' storage structure for MHD stability function';
                                                                                                                                                                                                                     
%equilibre
memoire.equi       =' storage structure for equilibrium function';
memoire.neo        =' storage structure for neoclassical function';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% mode batch pour les calculs de sources lorsque le calcul est long
% 0 -> pas de batch
% 1 -> mode batch
                                                                                                                                                                                                                     
% impuretee
batch.impur         =' batch mode for calculation of plasma composition (zeff, ae, nj ...)';
                                                                                                                                                                                                                     
% (sources)
batch.fci        =' batch mode for calculation of ICRF sources';
batch.fce        =' batch mode for calculation of ECRF sources';
batch.hyb        =' batch mode for calculation of LH sources';
batch.idn        =' batch mode for calculation of NBI sources';
batch.n0         =' batch mode for calculation of neutrals sources';
batch.fus        =' batch mode for calculation of fusion power + alpha particules sources';
batch.cyclo      ='  batch mode for calculation of cyclotronic radiative losses function';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% parametre des fonctions externes (cons)
%  il sont donnees dans une structure,
%  la fonction doit retourner une structure modele
%  si elle recoit comme premier argument la chaine
%  ''init'', ainsi qu''une tructure defaut, une structure
%  definissant les valeurs possibles et une structure
%  d''information :
%
%    [modele, defaut,possible,info]=zbidon(''init'');
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
%      .gui   -> nom du module specifique (si besoin) de l''interface graphique
%
                                                                                                                                                                                                                     
% structure generique :
scg = [];
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% (pde)
cons.impur      =' parameters of plasma composition function (zeff, ae, nj ...)';
cons.equi       =' parameters of plasma equilibrium function';
cons.neo        =' parameters of neoclassical quantities function';
cons.rip        =' parameters of thermal ripple quantities function';
                                                                                                                                                                                                                     
% (mhd)
cons.mhd.dds    =' parameters of sawteeth function';
cons.mhd.elm    =' parameters of elm function';
cons.mhd.limite =' parameters of simple elm and sawtooth stability function';
cons.mhd.stab   =' parameters of MHD stability (low toroidal mode numbers) function';
                                                                                                                                                                                                                     
% (sources)
cons.fci        =' parameters of ICRH sources function';
cons.fce        =' parameters of ECRH sources function';
cons.hyb        =' parameters of LH sources function';
cons.idn        =' parameters of NBI sources function';
cons.n0         =' parameters of neutrals sources function';
cons.bord       =' parameters of wall and gaz puff function';
cons.glacon     =' parameters of pellet sources function';
cons.fus        =' parameters of fusion power and alpha particules function';
cons.cyclo      =' parameters of cyclotronic radiative losses function';
cons.fast       =' parameters of fast particules function';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% coefficient des fonctions des coef de transport
cons.coefa        =' parameters of transport coefficient function "a"';
cons.coefb        =' parameters of transport coefficient function "b"';
cons.coefc        =' parameters of transport coefficient function "c"';
cons.coefd        =' parameters of transport coefficient function "d"';
cons.coefe        =' parameters of transport coefficient function "e"';
cons.coeff        =' parameters of transport coefficient function "f"';
                                                                                                                                                                                                                     
%coefficient multiplicateur des variables neoclassiques (utilise dans les equations de transport)
cons.neomulti.ee       =' neoclassical electron thermal diffusion coefficient multiplyer';
cons.neomulti.ve       =' neoclassical electron thermal convection speed multiplyer';
cons.neomulti.ii       =' neoclassical ion thermal diffusion coefficient multiplyer';
cons.neomulti.vi       =' neoclassical ion thermal convection speed multiplyer';
cons.neomulti.nn       =' neoclassical particule (el.) diffusion coefficient multiplyer';
cons.neomulti.vn       =' neoclassical particule (el.) convection speed multiplyer';
cons.neomulti.rot      =' neoclassical ion velocity diffusion coefficient multiplyer';
cons.neomulti.rotv     =' neoclassical ion velocity convection speed multiplyer';
                                                                                                                                                                                                                     
% autres
cons.plot         =' parameters of plot function during interactive calculations';
cons.asser        =' parameters of feedback controls function';
cons.machine      =' parameters device dependent.';
cons.post         =' parameters of posprocessing function';
                                                                                                                                                                                                                     
% constante physique (phys)
phys.c           =' speed of light in vacuum (m/s)  (definition)';
phys.h           =' Planck constant (J*s) (+/- 0.0000052e-34)';
phys.e           =' electron charge (C)   (+/- 0.000000063e-19)';
phys.mu0         =' permeablity of vacuum (H/m) (definition)';
phys.epsi0       =' permitivity of vacuum (F/m)  (definition)';
phys.g           =' gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)';
phys.k           =' Boltzmann constant (J/K)  (+/- 0.0000024e-23)';
phys.alpha       =' fine structure constant (+/- 0.000000027e-3 )';
phys.me          =' electron mass (at rest) (kg) (+/- 0.00000079e-31)';
phys.mp          =' proton mass (at rest) (kg)';
phys.ua          =' Atomic mass unit (kg) (+/- 1.00000013e-27)';
phys.avo         =' Avogadro number (mol^-1) (+/- 0.00000047e23)';
phys.sigma       =' Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)';
                                                                                                                                                                                                                     
% control du plot
plot.onoff        =' 1 if plot activated';
plot.intervalle   =' plots once every ''intervalle'' time steps, uses mode.plot if 0';
plot.run          =' 1 if program is running, 0 if paused';
plot.fin          =' 1 if ''stop'' button is activated';
plot.pause        =' pause between 2 calculations in s (if =0 uses drawnow)';
plot.figure1      =' handle of figure 1 (control)';
plot.figure2      =' handle of figure 2 (time)';
plot.figure3      =' handle of figure 3 (profile)';
plot.h.run        =' handle of ''stop'' button';
plot.h.pause      =' handle of ''pause'' button';
plot.h.keyboard   =' handle of ''keyboard'' button';
plot.h.step       =' handle of ''step by step'' button in paused mode';
plot.h.fin        =' handle of ''end'' button';
plot.h.info       =' handle of information zone';
plot.order        =' default colormap';
                                                                                                                                                                                                                     
% control du decoupage en temps (split)
%
% split.equi :
%    si la decoupe temporelle est activee (split.onoff =1), split.equi permet de controler l''appel a l''equilibre :
%     split.equi = 0  -> l''equilibre est appele a chaque sous pas de temps
%     split.equi = dt -> l''equilibre est appele une fois a la fin de chaque pas de temps et toutes les dt (s) . La derniere valeur calculer est
%                        utilisee pour les sous pas de temps intermediaires.
%
%
split.onoff       =' time-splitting authorised if = 1';
split.mode        =' time-splitting mode : 1 -> automatic on convergence criterion, 0 -> imposed';
split.dtmax       =' maximum time variation between two steps (s)';
split.dtmin       =' minimum time variation between two steps (s)';
split.nb          =' number of points in the interval if imposed';
split.equi        =' time between two equilibrium calcuations (s)';
                                                                                                                                                                                                                     
% parametre d''echantillonage (cf. zsample.m)
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
from.machine         =' name of server producing data';
from.shot.num        =' shot number';
from.shot.date       =' date + time of shot';
from.shot.info       =' information on shot';
from.creation.date   =' date + time of creation of CRONOS source file';
from.creation.user   =' name of the user who created the CRONOS source file';
from.creation.info   =' information on CRONOS source file generation (automatic)';
from.creation.com    =' information written by user during CRONOS source file generation';
from.source.desc     =' description of data sources (names of data sources)';
from.source.info     =' information about data used in Cronos (format variable with machine)';
from.sample.signal   =' sampling parameter for time dependent data =s(time,1)';
from.sample.groupe   =' sampling parameter for time and space dependent data =g(time,space)';
from.createur        =' name of function used to create source file';
from.option          =' options of source file creation function';
from.paroi.R         =' R vector defining vaccum chamber shape';
from.paroi.Z         =' Z vector defining vaccum chamber shape';
                                                                                                                                                                                                                     
% structrures pour les asservissements  (asser)
asser.data =[];
asser.delais=[];
                                                                                                                                                                                                                     
% intervalle de temps (pour la modification de variable de maniere groupees)
% cette structure est modifiee lors de l''edition des donnees
intervalle.temps  =' time base limits';
intervalle.calcul =' beginning and ending calculation times';
                                                                                                                                                                                                                     
% structure utiliser par l''interface graphique
edit.currentfile =' name of file loaded in workspace';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% structure pour le "profiler"
profile.onoff      =' if 1, triggers matlab "profiler" function';
profile.rapport    =' if 1, plots profiler report';
profile.data       =' data structure of "profiler" function';
profile.memoire    =' if 1, saves every time step the "memoire" structure (temporary save must be activated)';
profile.erreurmode =' if 1, saves the full CRONOS workspace if an error occurs';
profile.varlog1    =' name of variable 1 saved at every internal time step (temporary save must be activated)';
profile.varlog2    =' name of variable 2 saved at every internal time step (temporary save must be activated)';
profile.varlog3    =' name of variable 3 saved at every internal time step (temporary save must be activated)';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% structure pour simulink
simu.onoff        =' on/off flag for simulink in Cronos : if = 1, use simulink in Cronos simulation.';
simu.mdlname      =' name of the simulink model used in the Cronos simulation';
simu.mdltxt       =' source text of the simulink model used in the Cronos simulation';
simu.mdlinput     =' stucture in which are saved workspace parameters used by the simulink model';
simu.mdloutput    =' stucture in which are saved workspace output used by the simulink model';
simu.simparam     =' simulink parameter';
simu.mdlorigine   =' original path to model file (.mdl)';
simu.mdlloc       =' actual path to model file (.mdl)';
                                                                                                                                                                                                                     
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
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     

                                                                                                                                                                                                                     
% a - surface magnetiques et flux (magn)
                                                                                                                                                                                                                     
%      * geometrie (geo)
% valeur de geo.mode :
%        0 = plasma symetrique, donnee de a, e1, tr1 (avec tr1 = trh1 = trb1)
%        1 = plasma asym�rique, donnee de a, e1, trh1, trb1
%        2 = plasma asym�rique, donnee de (R,Z) de la derniere surface magnetique
%
geo.mode  =' description mode of last closed flux surface : 0 -> symmetric, 1 -> asymmetric, 2 -> RZ, 3 -> free boundary **';
geo.r0    =' major radius (m)  **';
geo.z0    =' Z shift of last closed flux surface (m)  **';
geo.a     =' minor radius (m)  **';
geo.e1    =' ellipticity of last closed flux surface  **';
geo.trh1  =' upper triangularity of last closed flux surface  **';
geo.trb1  =' lower triangularity of last closed flux surface **';
geo.ind1  =' indentation of last closed flux surface **';
geo.b0    =' toroidal field on axis (r0) without plasma (T) **';
geo.R     =' R vector defining last closed flux surface (m)';
geo.Z     =' Z defining last closed flux surface (not recentered) (m)';
                                                                                                                                                                                                                     
%      * equilibre (equi)
equi.rhomax    =' toroidal flux coordinate of last closed flux surface (m)';
equi.drhomaxdt =' time derivative of rhomax (m*s^-1)';
equi.phi       =' toroidal flux;';
equi.phid1     =' first radial derivative of toroidal flux;';
equi.phid2     =' second radial derivative of toroidal flux;';
equi.dphidt    =' time derivative of toroidal flux;';
equi.psi       =' poloidal flux (calculated by equilibrium for checking);';
equi.rhog      =' normalised geometrical radius (r/a)';
equi.d         =' Shafranov shift (m)';
equi.e         =' Ellipticity';
equi.a         =' toroidal flux coordinate (m)';
equi.raxe      =' major radius of the centre of the flux surfaces; raxe(end) = r0_equilibre (m)';
equi.zaxe      =' Z coordonate of the centre of the flux surfaces; (m)';
equi.trh       =' upper triangularity';
equi.trl       =' lower triangularity';
equi.indh      =' upper idendation';
equi.indl      =' lower idendation';
equi.grho2     =' <|grad(rho)|^2>';
equi.grho      =' <|grad(rho)|>';
equi.vpr       =' V'', <RJ>/4/pi^2 in Houlberg';
equi.spr       =' S''';
equi.dvprdt    =' time derivative of vpr';
equi.dsprdt    =' time derivative of spr';
equi.F         =' diamagnetic function F';
equi.grho2r2   =' <|grad(rho)|^2/R^2>';
equi.c2c       =' V''*<|gradient(rho)|^2/R^2> compute to find the correct jmoy at the edge and ip';
equi.grhor     =' <|grad(rho)|/R>';
equi.ri        =' <1/R>';
equi.b2        =' <B^2>';
equi.b2i       =' <1/B^2>';
equi.psi0      =' poloidal flux at plasma centre';
equi.q         =' q calculated by equilibrium for checking';
equi.shear     =' magnetic shear associated to equi.q';
equi.rmoy      =' <R> = <sqrt(g)>';
equi.ftrap     =' trapped electrons fraction';
equi.jmoy      =' current density recalculated from equilibrium';
equi.ptot      =' pressure recalculated from equilibrium';
equi.ip        =' plasma current from equilibrium (A)';
equi.li        =' internal inductance';
equi.betap     =' poloidal normalised pressure';
equi.betat     =' toroidal normalised pressure';
equi.betan     =' current normalised pressure';
equi.q0        =' safety factor central value';
equi.volume    =' plasma volume (m^3)';
equi.surface   =' plasma surface (m^2)';
equi.r2        =' <R^2> (m^2)';
equi.r2i       =' <1/R2> (m^-2)';
equi.r2tau2    =' <R^2/tau^2> (m^-6)';
equi.grho2b2   =' < |grad(rho)|^2/B^2 > (m^2.T^-2)';
equi.r3tau3    =' <R^3/tau^3> (m^-3)';
equi.r3tau     =' <R^3/tau> (m)';
equi.bnorme    =' normalisation used in equilibrium';
equi.conv      =' number of loops for convergence on jmoy ( <0 means no convergence)';
equi.oscil     =' 1 if convergence exited on oscillation';
equi.fail      =' 0 if convergence is ok';
equi.errcur    =' final relative error on current';
equi.errit     =' final relative error on equilibriu determination';
equi.amix      =' final value of convergence mixing parameter';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% la grille des surface magnetique
equi.R         =' R of magnetic flux surfaces (m)';
equi.Z         =' Z of magnetic flux surfaces (m)';
equi.BR        =' Br of magnetic flux surfaces (T)';
equi.BZ        =' Bz of magnetic flux surfaces (T)';
equi.BPHI      =' Bphi of magnetic flux surfaces (T)';
equi.rhoRZ     =' toroidal flux coordinate of magnetic flux surfaces (m)';
equi.psiRZ     =' poloidal flux of magnetic flux surfaces (m)';
equi.df2RZ     =' FF'' of magnetic flux surfaces';
equi.dprRZ     =' P'' of magnetic flux surfaces';
equi.frmode    =' Fourrier coefficient of LCMS transform.';
%     * stabilite mhd (mhd)
% les criteres sont chosi pour > 0 => instable
equi.mhd.gamma1       =' growth rate for mhd mode 1';
equi.mhd.gamma2       =' growth rate for mhd mode 2';
equi.mhd.gamma3       =' growth rate for mhd mode 3';
equi.mhd.vper1        =' for mhd mode 1, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)';
equi.mhd.vper2        =' for mhd mode 2, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)';
equi.mhd.vper3        =' for mhd mode 3, square root of the sum over the poloidal modes of the square of the flux surface perpendicular speed (m/s)';
equi.mhd.ballooning   =' stability criterion of ballooning modes ((1- FM) .* (FM ~= 0))';
equi.mhd.mercier      =' Mercier stability criterion ( -DR)';
equi.mhd.ideal        =' Ideal mhd stability criterion ( 0.25 - DI)';
equi.mhd.fail         =' error return code : if = 0, equilibrium calculation for mishka is ok';
equi.mhd.conv1        =' for mhd mode 1, number of convergence loops (ok if < 21, verify result if = 21)';
equi.mhd.conv2        =' for mhd mode 2, number of convergence loops (ok if < 21, verify result if = 21)';
equi.mhd.conv3        =' for mhd mode 3, number of convergence loops (ok if < 21, verify result if = 21)';
equi.mhd.erreur1      =' for mhd mode 1, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)';
equi.mhd.erreur2      =' for mhd mode 2, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)';
equi.mhd.erreur3      =' for mhd mode 3, error  ( ok if < 1e-6; if < 1e-5 ask a specialist)';
%  q(k) = m(k)/n(k)
equi.mhd.deltap       =' vector of delta'' values for rationnal equi.q  ( >0 instable)';
equi.mhd.m_deltap     =' vector of m associate to delta'', with equi.q(k) = equi.mhd.m_deltap(k) / n(k)';
                                                                                                                                                                                                                     
% b - donnees scalaires
% consignes   generales  = conditions au limites pour le pde (cons)
% une seule des deux consignes est prise en compte en fonction de
% la variable mode.cons.xx
cons.ip         =' plasma current (A) * {used to initialise the starting equilibrium}';
cons.vloop      =' loop voltage at plasma surface (V)';
cons.flux       =' edge poloidal flux (Wb)';
                                                                                                                                                                                                                     
cons.ne1        =' edge electron density (m^-3)';
cons.ge1        =' edge particule flux (m^-2*s^-1)';
                                                                                                                                                                                                                     
cons.te1        =' edge electron temperature (eV)';
cons.qe1        =' edge electron heat flux (W*m^-2*s^-1)';
cons.pe1        =' edge electron pressure (Pa)';
                                                                                                                                                                                                                     
cons.ti1        =' edge ion temperature (eV)';
cons.qi1        =' edge ion heat flux (W*m^-2*s^-1)';
cons.pion1      =' edge ion pressure (Pa)';
                                                                                                                                                                                                                     
cons.ffe1       =' edge electron "turbulence" flux';
cons.fe1        =' amplitude of edge electron "turbulence"';
                                                                                                                                                                                                                     
cons.ffi1       =' edge ion "turbulence" flux';
cons.fi1        =' amplitude of edge ion "turbulence"';
                                                                                                                                                                                                                     
cons.frot1       =' edge rotation flux';
cons.rot1        =' edge rotation';
                                                                                                                                                                                                                     
% conditition aux limites pour les impuretees
cons.c         =' edge gas injection (atoms/s)';
% pour le module de base :
% concentration de gaz majoritaire      -> 1   (D ou He)
% concentration de minoriatire 1        -> 2   (H,T)
% concentration de minoriatire 2        -> 3   (H,He)
% concentration d''impurete 1            -> 4   C
% concentration d''impurete 2            -> 5   O
                                                                                                                                                                                                                     
% consigne de pompage (vitesse en atome /s)
cons.pomp    =' pumping speed (atoms/s)';
                                                                                                                                                                                                                     
% consignes de puissance pour les sources
%comlexe -> abs(cons) =puissance en W et angle(cons) = phase de l''antenne en radian
% sauf pour idn ou real(cons) = puissance et imag(cons) = derivee temporelle de la puissance
%% MS 25.08.08: addition of tests number>0 for initalization of each power
%  if nbfci>0
%    cons.fci     =' ICRH power';
%  else
%    cons.fci     =' ICRH power';
%  end
%  if nbfce>0
%    cons.fce     =' ECRH power';
%  else
%    cons.fce     =' ECRH power';
%  end
%  if nbhyb>0
%    cons.hyb     =' LH power';
%  else
%    cons.hyb     =' LH power';
%  end
%  if nbidn>0
%    cons.idn     =' NBI power';
%  else
%    cons.idn     =' NBI power';
%  end
cons.fci     =' ICRH power';
cons.fce     =' ECRH power';
cons.hyb     =' LH power';
cons.idn     =' NBI power';
cons.ext     =' external pseudo-power';
                                                                                                                                                                                                                     
% consigne d''injection de glacons
cons.glacon =' pellet injection';
                                                                                                                                                                                                                     
% consigne pour le zeffm quand il est donnee a la place du zeff
cons.zeffm   =' reference value for average zeff when zeffm is used (no profile given)';
cons.nhnd    =' reference value for nH/nD (real part) and nT/nD (imaginary part)';
                                                                                                                                                                                                                     
% consigne pour le calcul de la stabilite mhd
cons.stab    =' number of toroidal modes in mhd stability calculations';
                                                                                                                                                                                                                     
% consignes utilisee par les asservissements
% cette partie est personnalisable avec init_perso et M a J code dynamique
cons.asser.nl0        =' reference for central line density';
cons.asser.ne0        =' reference for central density';
cons.asser.ne1        =' reference for edge density';
cons.asser.nemoy      =' reference for averaged density';
cons.asser.nbar       =' refernce for central line avaraged density (nbar)';
cons.asser.beta       =' reference for beta';
cons.asser.fgr        =' greenwald fraction reference';
cons.asser.te0        =' reference for central electron temperature';
cons.asser.te1        =' reference for edge electron temperature';
cons.asser.ti0        =' reference for central ion temperature';
cons.asser.ti1        =' reference for edge ion temperature';
cons.asser.li         =' reference for internal inductance';
cons.asser.q0         =' reference for q0';
cons.asser.qmin       =' reference for qmin';
cons.asser.c          =' reference for gas injection (precise which gas)';
cons.asser.vloop      =' reference for vloop';
cons.asser.ip         =' reference for ip';
cons.asser.pfcur      =' references for poloidal coil current';
                                                                                                                                                                                                                     
%
% evolution temporelle de eta (pour modification de n// de l''onde lors d''un choc)
%
cons.asser.etalh =' reference of LH efficiency depending on time (for zhybsimple)';
                                                                                                                                                                                                                     
%      * donnees generiques du plasma (gene)
gene.temps      =' time array **';
gene.dt         =' internal time step used for the calculation of the time derivatives (s)';
gene.conv       =' number of convergence loops done by the solver (<0 if no convergence)';
gene.nbsplit    =' number of splits during a time step for convergence';
gene.flops      =' number of floating operations in matlab';
gene.cputime    =' total CPU time (s) + i * user CPU time (s)';
gene.memory     =' real memory used (Mo) + i * virtual memory (Mbytes)';
gene.datation   =' elapsed time since beginning of run (s)';
                                                                                                                                                                                                                     
gene.surface    =' surface of plasma cross section';
gene.volume     =' plasma volume';
                                                                                                                                                                                                                     
gene.zeffm      =' line averaged zeff';
gene.zeffne     =' averaged zeff, weighted by density';
                                                                                                                                                                                                                     
gene.netot      =' number of electrons in plasma';
gene.nemoy      =' volume averaged electron density';
gene.dnetotdt   =' time derivative of the number of electrons in plasma';
gene.bilan_e    =' particule balance (electrons/s)';
gene.piqne      =' electron density peaking ne0/nemoy';
gene.nbar       =' nl0/2/rhomax (central line averaged density)';
                                                                                                                                                                                                                     
gene.temoy      =' volume averaged electron temperature';
gene.piqte      =' electron temperature peaking te0/temoy';
                                                                                                                                                                                                                     
gene.timoy      =' volume averaged ion temperature';
gene.piqti      =' ion temperature peaking ti0/timoy';
                                                                                                                                                                                                                     
gene.qmin       =' value of the minimum of the q profile';
gene.rhoqmin    =' position of the minimum of the q profile';
                                                                                                                                                                                                                     
gene.we         =' thermal electrons energy';
gene.wion       =' thermal ions energy';
gene.wdia       =' diamagnetic energy (thermal+supra)';
gene.wth        =' total thermal energy';
gene.wbp        =' poloidal energy';
%gene.wmag       =' energie magnetique du plama';
                                                                                                                                                                                                                     
gene.dwedt      =' time derivative of thermal electrons energy';
gene.dwiondt    =' time derivative of thermal ions energy';
gene.dwdiadt    =' time derivative of diamagnetic energy (thermal+supra)';
gene.dwthdt     =' time derivative of total thermal energy';
gene.dwbpdt     =' time derivative of poloidal energy';
%gene.dwmagdt    =' derivee temporelle de l''energie magnetique du plama';
                                                                                                                                                                                                                     
gene.pel        =' power balance on electrons (additional +alpha +ohm -ei -brem -rad -en0 -synch )';
gene.pion       =' power balance on ions (additional +alpha +ohm -ei -brem -rad -en0 -synch )';
gene.padde      =' total heating power on electrons (additional +alpha +ohm)';
gene.paddion    =' total heating power on ions (additional +alpha +ohm)';
gene.paddtot    =' total heating power (additional +alpha +ohm)';
gene.paddfci    =' ICRH power deposited on thermals';
gene.paddhyb    =' LH power deposited on thermals';
gene.paddfce    =' ECRH power deposited on thermals';
gene.paddidn    =' NBI power deposited on thermals';
gene.paddohm    =' ohmic heating power';
gene.paddfus    =' alpha heating power';
gene.paddext    =' external heat power';
gene.paddrip    =' ripple heat loss';
gene.pbrem      =' bremsstrahlung power loss';
gene.prad       =' radiative power loss (without brem/sync)';
gene.pcyclo     =' synchrotron power loss';
gene.prip       =' thermal ripple power losses';
gene.pn0        =' charge exchange with cold neutrals power loss';
gene.nadd       =' number of electrons injected in plasma';
gene.ploss      =' plasma power loss (total heating power - Pbrem, Psynch and dwdia/dt)';
gene.qei        =' electron-ion exchange power due to collisions (equipartition)';
gene.qneo       =' electron-ion exchange power due to neoclassical effects';
                                                                                                                                                                                                                     
%gene.sne        =' source totale d''electron';
gene.sne_bord   =' electron source due to edge neutrals';
gene.sne_idn    =' electron source due to NBI';
gene.snions     =' source of the various ion species';
gene.evolution  =' evolution of the 0d content for each specie (used for calculation of zeff0d)';
gene.zeff0d     =' average zeff calculated from 0d data';
gene.neutralite =' electroneutrality';
gene.nions      =' ion species content';
                                                                                                                                                                                                                     
gene.taue       =' electron energy confinement time';
gene.tauion     =' ion energy confinement time';
gene.taune      =' particules confinement time';
gene.tau        =' total confinement time';
gene.tauth      =' thermal energy confinement time';
gene.tauj       =' current diffusion time';
gene.tauei      =' equipartition time';
                                                                                                                                                                                                                     
gene.betap      =' beta poloidal';
gene.beta       =' beta total';
gene.betastar   =' beta *  (linked to fusion power)';
gene.betadia    =' beta shafranof (diamagnetic)';
                                                                                                                                                                                                                     
gene.ip         =' plasma current (analytical formula) (A)';
gene.ipepar     =' plasma current (from E/eta +Jni) (A)';
gene.ipohm      =' ohmic current (E/eta) (A)';
gene.ipjmoy     =' plasma current (from jmoy) (A)';
gene.icd        =' additionnal driven current (A)';
gene.ihyb       =' LH driven current (A)';
gene.ifce       =' ECRF driven current (A)';
gene.ifci       =' ICRF driven current (A)';
gene.iidn       =' NBI driven current (A)';
gene.ifus       =' Alpha (fusion) induce current (A)';
gene.iboot      =' bootstrap current (A)';
gene.ini        =' non-inductive current (add. + boot.) (A)';
gene.li         =' internal inductance (H) {formula with bpol}';
gene.liip       =' internal inductance (H) {formula with bpol and ip}';
gene.licr       =' internal inductance (H) {formula with q}';
gene.vloop      =' loop voltage on fixed loop  (V)';
gene.vsurf      =' loop voltage on plasma surface (V)';
gene.vres       =' average resistive loop voltage U = Pohm/Ip (V)';
                                                                                                                                                                                                                     
%      * mode de calcul des donnees
%
%  pour les coef et les sources  :
%
%      0  = coef ou source a 0
%      1  = coef ou source lu dans la structure d''entree
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
%      1  = lu dans la structure d''entree
%      2  = calcule a l''aide de l''equation differentielle
%
%  pour  pde  sur ne  :
%
%      0  = mis a 0
%      1  = lu dans la structure d''entree
%      2  = calcule a l''aide de l''equation differentielle
%      3  = forme donnees en entree, densite totale controlee
%      4  = forme donnees en entree, densite totale controlee + piquage loi en n_gr/nbar
%      5  = forme donnees en entree, densite totale controlee + piquage loi en de H. Wiesen
%
%  pour l''equilibre :
%      0 et 1 = lu dans la structure d''entree
%      2      = calculer
%      3      = extrapoler pour le temps k+1 du temps k
%
% (pde)
mode.impur      =' caculation mode of plasma composition (zeff, ae, nj ...)';
mode.psi        =' calculation mode of poloidal flux (current diffusion)';
mode.nel        =' calculation mode of electron density';
mode.pe         =' calculation mode of electron pressure (temperature)';
mode.pion       =' calculation mode of ion pressure (temperature)';
mode.equi       =' calculation mode of magnetic equilibrium';
mode.neo        =' calculation mode neoclassical quantities';
mode.rip        =' calculation mode thermal ripple quantities (default value =0)';
mode.fluce      =' calculation mode of electron turbulence amplitude';
mode.flucion    =' calculation mode of ion turbulence amplitude';
mode.rot        =' calculation mode of plasma rotation';
                                                                                                                                                                                                                     
% mode des consigne pour les conditions aux limites
% mode = 0 -> valeur (ne,Te et Ti) ou Ip donne
% mode = 1 -> flux (des equations sur ne, Pe et Pion)  ou Vloop donne
% mode = 2 -> valeur (Pe ou Pion) ou Psi (W)  au bord donnee
mode.cons.psi        =' boundary condition for current diffusion : Ip (0), Vloop (1), edge psi (2), edge psi computed by a free boundary equilibrium';
mode.cons.ne         =' boundary condition for electron density equation : edge density (0), edge particule flux (1)';
mode.cons.pe         =' boundary condition for electron temperature equation : edge temperature (0), edge heat flux (1), edge pressure (2)';
mode.cons.pion       =' boundary condition for ion temperature equation : edge temperature (0), edge heat flux (1), edge pressure (2)';
mode.cons.fluce      =' boundary condition for electron turbulence equation : edge turbulence (0), edge turbulence flux (1)';
mode.cons.flucion    =' boundary condition for ion turbulence equation : edge turbulence (0), edge turbulence flux (1)';
mode.cons.rot        =' boundary condition for plasma rotation equation : edge rotation (0), edge rotation flux (1)';
                                                                                                                                                                                                                     
% source des consignes pour le bord
mode.consbord.ne     =' edge electron density given as input (0) or by edge module (1)';
mode.consbord.te     =' edge electron temperature given as input (0) or by edge module (1)';
mode.consbord.ti     =' edge ion temperature given as input (0) or by edge module (1)';
mode.consbord.fluxge =' edge electron particule flux given as input (0) or by edge module (1)';
mode.consbord.fluxqe =' edge electron heat flux given as input (0) or by edge module (1)';
mode.consbord.fluxqi =' edge ion heat flux given as input (0) or by edge module (1)';
                                                                                                                                                                                                                     
% (mhd)
mode.mhd.dds    =' calculation mode of sawteeth';
mode.mhd.elm    =' calculation mode of elms';
mode.mhd.limite =' calculation mode of stability limits';
mode.mhd.stab   =' calculation mode of mhd modes stability';
                                                                                                                                                                                                                     
% (sources)
mode.fci        =' calculation mode of H/CD sources due to ICRH';
mode.fce        =' calculation mode of H/CD sources due to ECRH';
mode.hyb        =' calculation mode of H/CD sources due to LH';
mode.idn        =' calculation mode of H/CD sources due to NBI';
mode.n0         =' calculation mode of edge neutrals source';
mode.bord       =' calculation mode of edge module';
mode.glacon     =' calculation mode of pellet deposition';
mode.fus        =' calculation mode of fusion + alpha power';
mode.ohm        =' calculation mode of ohmic power';
mode.qneo       =' calculation mode of equipartition due to ion pressure gradient';
mode.qei        =' calculation mode of ion - electron equipartition';
mode.prad       =' calculation mode of radiative power losses';
mode.brem       =' calculation mode of bremsstrahlung power lossess';
mode.cyclo      =' calculation mode of synchrotron power losses';
                                                                                                                                                                                                                     
% coefficient de transport
mode.ee        =' calculation mode coef in equation Pe in front of Pe (Te)';
mode.ei        =' calculation mode coef in equation Pe in front of Pi';
mode.en        =' calculation mode coef in equation Pe in front of ne';
mode.ej        =' calculation mode coef in equation Pe in front of j';
mode.ve        =' calculation mode coef in equation Pe, convection speed';
mode.ep        =' calculation mode coef in equation Pe in front of Pe (Pe)';
                                                                                                                                                                                                                     
mode.ie        =' calculation mode coef in equation Pi in front of Pe';
mode.ii        =' calculation mode coef in equation Pi in front of Pi (Ti)';
mode.in        =' calculation mode coef in equation Pi in front of ne';
mode.ij        =' calculation mode coef in equation Pi in front of j';
mode.vi        =' calculation mode coef in equation Pi, convection speed';
mode.ip        =' calculation mode coef in equation Pi in front of Pi (Pi)';
                                                                                                                                                                                                                     
mode.ne        =' calculation mode coef in equation Ne in front of Pe';
mode.ni        =' calculation mode coef in equation Ne in front of Pi';
mode.nn        =' calculation mode coef in equation Ne in front of ne';
mode.nj        =' calculation mode coef in equation Ne in front of j';
mode.vn        =' calculation mode coef in equation Ne, convection speed';
                                                                                                                                                                                                                     
mode.fefe      =' calculation mode coef in equation fluce in front of fluce';
mode.fev       =' calculation mode coef in equation fluce, convection speed';
                                                                                                                                                                                                                     
mode.fifi      =' calculation mode coef in equation flucion in front of flucion';
mode.fiv       =' calculation mode coef in equation flucion, convection speed';
                                                                                                                                                                                                                     
mode.rotc      =' calculation mode coef in equation rotation plasma in front of rot';
mode.rotv      =' calculation mode coef in equation rotation plasma, convection speed';
                                                                                                                                                                                                                     
%mode de variables neo
mode.eta       =' calculation mode of resistivity';
mode.jboot     =' calculation mode of bootstrap current';
                                                                                                                                                                                                                     
% autres
mode.plot          =' parametre of plot function during interactive run';
mode.zeff          =' zeff calculated or given as input';
mode.ae            =' ae calculated or given as input';
mode.asser         =' feedback mode (0,1) = reference values given as input / 2 = reference values calculated using feedback module';
mode.premiertemps  =' parameter for the initialisation of some modules';
mode.post          =' post processing mode (effect depend of the value of param.gene.post)';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% mode.cons.zeffm    =' zeffm est donne en consigne externe (0) ou par le modele 0d (1)';
% mode.zeffm         =' zeffm est donne en consigne externe (0) ou par le modele 0d (1)';
                                                                                                                                                                                                                     
% c - profils  (prof)
% des equations
prof.psi        =' poloidal flux profile (Wb) *';
prof.ne         =' electron density profile (m^-3) *';
prof.pe         =' electron pressure profile (Pa) *';
prof.pion       =' ion pressure profile (Pa) *';
prof.ae         =' ni/ne profile';
prof.fluce      =' electron turbulence amplitude profile (su)';
prof.flucion    =' ion turbulence amplitude profile  (su)';
prof.rot        =' toroidal moment rotation profile  of the plasma : sum(m_k n_k <R Vphi_k>) (kg m^-1 s^-1)';
                                                                                                                                                                                                                     
% les derivees spatiales premieres (en x pas en rho)
prof.psid1      =' 1st radial derivative of poloidal flux';
prof.ned1       =' 1st radial derivative of electron density';
prof.ped1       =' 1st radial derivative of electron pressure';
prof.piond1     =' 1st radial derivative of ion pressure';
prof.aed1       =' 1st radial derivative of ni/ne';
prof.fluced1    =' 1st radial derivative of electron turbulence amplitude profile';
prof.fluciond1  =' 1st radial derivative of ion turbulence amplitude profile';
prof.rotd1      =' 1st radial derivative of toroidal moment rotation';
                                                                                                                                                                                                                     
% les derivees spatiales seconde  (en x pas en rho)
prof.psid2      =' 2nd radial derivative of poloidal flux';
prof.ned2       =' 2nd radial derivative of electron density';
prof.ped2       =' 2nd radial derivative of electron pressure';
prof.piond2     =' 2nd radial derivative of ion pressure';
prof.aed2       =' 2nd radial derivative ni/ne';
prof.fluced2    =' 2nd radial derivative of electron turbulence amplitude profile';
prof.fluciond2  =' 2nd radial derivative of ion turbulence amplitude profile';
prof.rotd2      =' 2nd radial derivative of toroidal moment rotation';
                                                                                                                                                                                                                     
% les profils de Te et Ti (car souvent utilis)
prof.te         =' electron temperature profile (eV)';
prof.ti         =' ion temperature profile (eV)';
prof.ni         =' total ion density profile (m-3)';
prof.gte        =' gradient of electron temperature (eV/m)';
prof.gti        =' gradient of ion temperature  (eV/m)';
prof.gne        =' gradient of electron density  (m^-2)';
prof.gni        =' gradient of ion density  (m^-2)';
prof.alpha      =' profile of alpha = - q^2 R grad(Beta)';
                                                                                                                                                                                                                     
% les flux
prof.flux.qe    =' electron heat flux profile (W * m^-2)';
prof.flux.qi    =' ion heat flux profile    (W * m^-2)';
prof.flux.ge    =' electron flux profile (m^-2 * s^-1)';
prof.flux.gi    =' ion flux profile (m^-2 * s^-1)';
prof.flux.fluce =' electron turbulence amplitude flux profile';
prof.flux.fluci =' ion turbulence amplitude flux profile';
prof.flux.rot   =' toroidal moment rotation flux profile';
                                                                                                                                                                                                                     
% les coefficients deduits
% effectif -> flux/gradient
% anormal -> corrige des coefficient neoclassiques
prof.flux.keeff =' effective electron heat diffusion coefficient';
prof.flux.kieff =' effective ion heat diffusion coefficient';
prof.flux.deeff =' effective electron diffusion coefficient';
prof.flux.dieff =' effective ion diffusion coefficient';
prof.flux.kean  =' effective anomalous electron heat diffusion coefficient';
prof.flux.kian  =' effective anomalous ion heat diffusion coefficient';
prof.flux.dean  =' effective anomalous electron diffusion coefficient';
prof.flux.dian  =' effective anomalous ion diffusion coefficient';
                                                                                                                                                                                                                     
% flux sortant d''une surface magnetique
prof.flux.sortant.ge =' electron flux coming through magnetic surface (electrons/s)';
prof.flux.sortant.gi =' ion flux coming through magnetic surface (ions/s) [ supppose totalement ionis]';
prof.flux.sortant.qe =' electron heat flux coming through magnetic surface (W)';
prof.flux.sortant.qi =' ion heat flux coming through magnetic surface (W)';
                                                                                                                                                                                                                     
% les derivees temporelles
prof.dpsidt     =' time derivative of the poloidal flux';
prof.dnedt      =' time derivative of the electron density';
prof.dpedt      =' time derivative of the electron pressure';
prof.dpiondt    =' time derivative of the ion pressure';
prof.daedt      =' time derivative of ni/ne';
prof.dflucedt   =' time derivative of electron turbulence amplitude';
prof.dfluciondt =' time derivative of ion turbulence amplitude';
prof.drotdt     =' time derivative of toroidal moment rotation';
                                                                                                                                                                                                                     
% les derivees temporelles 3 pts
prof.dpsidt3p     =' time derivative of the poloidal flux (3 points)';
prof.dnedt3p      =' time derivative of the electron density (3 points)';
prof.dpedt3p      =' time derivative of the electron pressure(3 points)';
prof.dpiondt3p    =' time derivative of the ion pressure (3 points)';
prof.daedt3p      =' time derivative of ni/ne (3 points)';
prof.dflucedt3p   =' time derivative of electron turbulence amplitude (su) (3 points)';
prof.dfluciondt3p =' time derivative of ion turbulence amplitude (3 points)';
prof.drotdt3p     =' time derivative of toroidal moment rotation (3 points)';
                                                                                                                                                                                                                     
% les profils de pressions supra et totale
% prof.psupra     =' composante suprathermique de la pression totale';
prof.ptot       =' total pressure (including suprathermals) *';
                                                                                                                                                                                                                     
% les longueurs de gradient les plus utilisee
prof.shear      =' magnetic shear';
prof.lte        =' electron temperature gradient length (Te / grad(Te))';
prof.lti        =' ion temperature gradient length (Ti / grad(Ti))';
prof.lne        =' electron density gradient length (ne / grad(ne))';
prof.lni        =' ion density gradient length (ni / grad(ni))';
prof.rhostarte  =' "rhostar" criterion for electron ITBs';
                                                                                                                                                                                                                     
% electromagnetisme
prof.dpsidt_phi =' time derivative of psi (poloidal flux) at constant phi (toroidal flux)';
prof.jphi       =' effective toroidal current density (A/m^2)';
prof.jpol       =' effective poloidal current density (A/m^2)';
prof.jeff       =' effective current density (A/m^2)';
prof.jmoy       =' current density (magnetic surface averaged)(A/m^2) *';
prof.q          =' safety factor';
prof.epar       =' parellel electric field E// = <E.B>/B0   (V/m)';
%prof.eparlis    =' profil de E// = <E.B>/B0   lisser (lissage glissant sur tauj/10 s,V/m)';
prof.bpolm      =' poloidal field (magnetic surface averaged) <|Bp|> (T)';
prof.bpol       =' poloidal field (magnetic surface averaged) sqrt(<Bp.Bp>)  (T) -> relie a l''energie';
prof.bphi       =' toroidal field (magnetic surface averaged) <B.grad(phi)> (T)';
prof.ej         =' ohmic power';
                                                                                                                                                                                                                     
% profils divers
prof.zeff       =' effective charge';
prof.xdur       =' hard-X ray profile (from diagnostic, on Tore Supra)';
prof.vtor_exp   =' toroidal angular speed (rad/s)';
prof.vpol_exp   =' poloidal angular speed (rad/s)';
prof.mhd_cd     =' factor applied to current drive efficiency when MHD occurs (Jcd = mhd_cd * Jcd0)';
                                                                                                                                                                                                                     
% d -sources (source)
% les flux de neutrons de diverses origines
neutron.dt        =' neutron source from reaction  D + T -> He4 + n';
neutron.tt        =' neutron source from reaction  T + T -> He4 + 2n';
neutron.dd        =' neutron source from reaction  D + D -> He3 + ne';
neutron.the3      =' neutron source from reaction  T + He3 -> He4 + p + n';
                                                                                                                                                                                                                     
% les flux de protons de diverses origines
proton.dd         =' proton source from reaction D + D -> T + p';
proton.dhe3       =' proton source from reaction D + He3 -> He4 + p';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% 1 - structure generique pour les sources de type chauffage
sg.el             =' power density coupled to electrons';
sg.ion            =' power density coupled to ions';
sg.ne             =' particle source';
sg.j              =' parallel current density source ( j = <J.B>/Bo)';
sg.w              =' rotation moment source <R Fk.ephi> due to heat source, toroidal component ,  unity = kg m^-1 s^-2';
sg.wb             =' rotation moment source <F1k.B> due to heat source, parallel B component for NClass module ,  unity = T N m^-3';
sg.q              =' extrenal heat force source <F2k.B> due to heat source, parallel B component for NClass module  ,  unity =  T N m^-3';
sg.fluce          =' electron turbulence amplitude source';
sg.flucion        =' ion turbulence amplitude source';
sg.psupra         =' total suprathermal pressure (perpendicular) (Pa)';
sg.paniso         =' anisotropic component of suprathermal pressure = p// - 1/2 * pperp (Pa)';
sg.neutron        = neutron;
sg.proton         = proton;
sg.err            =' return code of source module';
sg.synergie       =' synergie factor for LHCD & ECCD';
                                                                                                                                                                                                                     
% remarque :
%  Psupra  =  densite d''energie supra thermique perpendiculaire (= pperp)
%             si la pression suprathermique est isotrope, la densite totale d''energie
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
source.n.n0       =' sources of matter per specy due to cold neutrals (atoms)';
source.n.fus      =' sources of matter per specy due to fusion reaction (atoms)';
source.n.idn      =' sources of matter per specy due to NBI (atoms)';
source.n.rip      =' sources of matter per specy due to ripple (atoms)';
source.n.n0th     =' thermal neutrals density coming from the edge (equivalent electron)';
                                                                                                                                                                                                                     
% sources rayonnements
source.prad       =' radiated power source/sink';
source.brem       =' bremsstrahlung power source/sink';
source.cyclo      =' cyclotron power source/sink';
                                                                                                                                                                                                                     
% sources collisions
source.ohm        =' ohmic power deposition';
source.qei        =' electron-ion heat exchange (equipartition)';
source.qneo       =' neoclassical electron-ion heat exchange';
source.jboot      =' bootstrap current density (A/m^2)';
                                                                                                                                                                                                                     
% turbulence
source.fluce      = sg;
source.flucion    = sg;
                                                                                                                                                                                                                     
% sommation des termes sources
source.totale     = sg;    % sum of all sources
                                                                                                                                                                                                                     
% les coefficents (coef)
% coefficient de transport
coef.eta        =' resistivity (ohm * m)';
                                                                                                                                                                                                                     
coef.ee        =' coefficient of Pe (Te) in Pe equation  (m^-1 * s^-1)';
coef.ei        =' coefficient of Pi in Pe equation  (m^-1 * s^-1)';
coef.en        =' coefficient of ne in Pe equation  (m^-1 * s^-1)';
coef.ej        =' coefficient of j in Pe equation   (?)';
coef.ve        =' convective speed in Pe equation   (m*s^-1)';
coef.ep        =' coefficient of Pe (Pe) in Pe equation (m^-1 * s^-1)';
                                                                                                                                                                                                                     
coef.ie        =' coefficient of Pe in Pi equation  (m^-1 * s^-1)';
coef.ii        =' coefficient of Pi (Ti) in Pi equation (m^-1 * s^-1)';
coef.in        =' coefficient of ne in Pi equation  (m^-1 * s^-1)';
coef.ij        =' coefficient of j in Pi equation     (?)';
coef.vi        =' convective speed in Pi equation (m*s^-1)';
coef.ip        =' coefficient of Pi (Pi) in Pi equation (m^-1 * s^-1)';
                                                                                                                                                                                                                     
coef.ne        =' coefficient of Pe in ne equation   (m^2 * s^-1)';
coef.ni        =' coefficient of Pi in ne equation   (m^2 * s^-1)';
coef.nn        =' coefficient of ne in ne equation   (m^2 * s^-1)';
coef.nj        =' coefficient of j in ne equation     (?)';
coef.vn        =' convective speed in ne equation (m*s^-1)';
                                                                                                                                                                                                                     
coef.fev       =' convective speed in fluce equation';
coef.fefe      =' diffusion coefficient in fluce equation';
                                                                                                                                                                                                                     
coef.fiv       =' convective speed in flucion equation';
coef.fifi      =' diffusion coefficient in flucion equation';
                                                                                                                                                                                                                     
coef.rotv      =' convective speed in toroidal rotation equation';
coef.rot       =' diffusion coefficient in toroidal rotation equation';
                                                                                                                                                                                                                     
% les grandeur neoclassique (neo)
neo.eta           =' neoclassical resistivity (ohm * m)';
neo.jboot         =' bootstrap current density (A/m^2)';
neo.coef.ee       =' neoclassical electron heat diffusivity  (m/s)';
neo.coef.ve       =' neoclassical electron heat convection';
neo.coef.ii       =' neoclassical ion heat diffusivity';
neo.coef.vi       =' neoclassical ion heat convection';
neo.coef.nn       =' neoclassical particle diffusivity';
neo.coef.vn       =' neoclassical particle convection';
neo.coef.rotv     =' neoclassical rotation convection (m s^-1)';
neo.coef.rot      =' neoclassical rotation diffusivity (m^-1 s^-1)';
neo.vtor          =' neoclassical toroidal rotation speed (m/s), at R = Rmax of each flux surface, for each ion specy';
neo.vtheta        =' neoclassical poloidal rotation speed (m/s), at R = Rmax of each flux surface, for each ion specy';
neo.mach          =' neoclassical mach number, at R = Rmax of each flux surface, for each ion specy';
neo.er            =' neoclassical radial electric field (V/m) = Er / gradient(rho).';
neo.gammae        =' neoclassical rotation rotation shearing (s^-1), at R = Rmax of each flux surface.';
neo.g             =' this parameter measures the magnitude of the // friction force relative to the // pressure gradient; abs(G)<< 1 to have valid NClass results.';
neo.flux.ne       =' neoclassical electron flux';
neo.flux.nion     =' neoclassical ion flux (sum over species)';
neo.flux.qe       =' neoclassical electron heat flux';
neo.flux.qion     =' neoclassical ion heat flux (sum over species)';
neo.qei           =' electron-ion heat exchange (equipartition)';
neo.qneo          =' electron-ion heat exchange (neoclassical)';
neo.w_ion         =' 2*pi * plasma solid rotation frequency in toroidal direction , for each ion specy';
neo.w_e           =' 2*pi * plasma electron solid rotation frequency in toroidal direction';
neo.utheta_i      =' fluid velocity, theta component :<V_k . theta> / <B . theta> , for each ion specy (cf NCLASS)';
neo.utheta_e      =' electron fluid velocity, theta component :<V_e . theta> / <B . theta>  (cf NCLASS)';
%  le moment d''ordre 1 est la quantite de mouvement (v), 2 -> la chaleur (v^3), 3- > en v^5
%  ref : W. A. Houlberg, Phys. PLasmas 4 (9), september 1997, p 3230- (eq 4 et suivantes)
%  l''usage normale : fexizpr1 et fexizpr2 doivent etre fournis pour IDN et fexizpr3 = 0;
%  fexizpr1 = datak.source.totale.wb
%  fexizpr2 = datak.source.totale.q
neo.force1        =' first moment of external parallel force on species i,z (T*j/m**3)';
neo.force2        =' second moment of external parallel force on species i,z (T*j/m**3)';
neo.force3        =' third moment of external parallel force on species i,z (T*j/m**3)';
neo.fail          =' convergence failure flag ( force ball ou Nclass)';
                                                                                                                                                                                                                     
% les donnees d''un run transp
transp.jnbtrx   =' neutral beam current';
transp.jbstrx   =' bootstrap current';
transp.jtottrx  =' total current';
transp.titrx    =' ion temperature';
transp.tetrx    =' electron temperature';
transp.netrx    =' electron density';
transp.nitrx    =' ion density';
transp.qohtrx   =' ohmic heating source';
transp.qrfetrx  =' electron heating source due to RF';
transp.qrfitrx  =' ion heating source due to RF';
transp.qnbetrx  =' electron heating source due to NBI';
transp.qnbitrx  =' ion heating source due to NBI';
transp.qbthtrx  ='';
transp.zefftrx  =' Zeff profile';
transp.rotortrx =' toroidal velocity';
transp.chietrx  =' transport coefficient (ion-electron)';
transp.chiitrx  =' transport coefficient (ion-ion)';
transp.qietrx   ='  electron-ion heat exchange (equipartition)';
transp.qtrx     =' q profile';
transp.ttransp  =' slice time of transp run';
% les impuretees
impur.prad       =' radiated power';
impur.brem       =' bremsstrahlung';
%impur.cyclo      =' cyclotron radiation';
impur.zeff       =' Zeff';
impur.impur      =' impurities density';
impur.ae         =' ni/ne';
impur.conv       =' number of convergence loops done with the impurities transport code and the neoclassical module to obtain consitentcy.';
impur.fail       =' convergence failure flag (impurities or consitency with neoclassical code)';
impur.neofail    ='convergence failure flag (Nclass)';
impur.pradsol    =' radiated power in the SOL (W)';
impur.alpha      =' electroneutrality normalisation a the ouptut of impurities transport module';
                                                                                                                                                                                                                     
% le bord (structure remplie par le module de bord)
% attention cette structure ne doit pas avoir de sous structure
% sinon il faut modifier zbord.m !
bord.temp         =' wall temperature (K)';
bord.nb           =' number of atoms in wall (electron equivalent)';
bord.pref         =' power loss due to cooling (W)';
bord.pchauffe     =' heating power (W)';
bord.reflex       =' recycling coefficient';
bord.fluxpompe    =' pumped particle flux (electrons * s^-1)';
bord.fluxplasma   =' particle flux coming from plasma (electrons *  s^-1)';
bord.fluxgazout   =' particle flux coming from plasma, for each specy (atoms  * s^-1)';
bord.fluxgaz      =' fuelling, for each specy (atoms * s^-1)';
bord.fluxmur_c    =' hot neutral flux recycled towards the plasma (electrons * s^-1)';
bord.fluxmur_f    =' cold neutral flux from the wall towards the plasma (fuelling +degasing + sputtering, electrons * s^-1)';
bord.fluxgazin_c  =' hot neutral flux recycled towards the plasma, for each specy (atoms * s^-1)';
bord.fluxgazin_f  =' cold neutral flux from the wall towards the plasma, for each specy (fuelling +degasing + sputtering, atoms * s^-1)';
bord.nebord       =' edge electron density (depending on limiter parameters) (m^-3)';
bord.nibord       =' edge ion density, for each specy (depending on limiter parameters) (m^-3)';
bord.tebord       =' edge electron temperature (depending on limiter parameters) (eV)';
bord.tibord       =' edge ion temperature (depending on limiter parameters) (eV)';
bord.fluxgebord   =' edge particle flux (electrons * s^-1)';
bord.fluxqebord   =' edge electron heat flux (W)';
bord.fluxqibord   =' edge ion heat flux(W)';
                                                                                                                                                                                                                     
% description detaillee de la SOL
bord.x            =' extended normalized coordonnate param.gene.x in the SOL';
bord.ne           =' electron density profile in the SOL (m^-3)';
bord.ni           =' ion density profile in the SOL (m^-3)';
bord.te           =' electron temperature profile in the SOL (ev)';
bord.ti           =' ion temperature profile in the SOL (ev)';
bord.impur        =' ions and impurities density  in the SOL (m^-3)';
                                                                                                                                                                                                                     
% les evennents (evx)
evx.dds      =' sawtooth event triggering';
evx.elm      =' elm event triggering';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
% les donnees experimentale pour le plot
exp.ne0      =' measured central electron density (m^-3)';
exp.nea      =' measured edge electron density (m^-3)';
exp.te0      =' measured central electron temperature(eV)';
exp.tea      =' measured edge electron temperature (eV)';
exp.ni0      =' measured central ion density (m^-3)';
exp.nia      =' measured edge ion density (m^-3)';
exp.ti0      =' measured central ion temperature (eV)';
exp.tia      =' measured edge ion temperature (eV)';
exp.j0       =' measured central current density (A*m^-2)';
exp.ip       =' measured plasma current (A)';
exp.li       =' measured internal inductance';
exp.vloop    =' measured loop voltage (V)';
exp.q0       =' measured central safety factor';
exp.qa       =' measured edge safety factor';
exp.betadia  =' measured diamagnetic beta';
exp.hybspec.npar =' npar of continuous LH spectrum';
exp.hybspec.pow  =' power of continuous LH spectrum';
                                                                                                                                                                                                                     
                                                                                                                                                                                                                     
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
                                                                                                                                                                                                                     
param.fonction.impur = [ param.fonction.impur,' {zinebcompo = Calcul simple de alphae, du rayonnement et de la densite d''impurtes (+zeff)}'];
if ~isempty('zinebcompo') & isempty(noconnexion)
    try
       par = feval('zinebcompo');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.impur = par.info;
   end
end
param.fonction.equi = [ param.fonction.equi,' {zequi_helena = Equilibrium computation}'];
if ~isempty('zequi_helena') & isempty(noconnexion)
    try
       par = feval('zequi_helena');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.equi = par.info;
   end
end
param.fonction.neo = [ param.fonction.neo,' {zneo = Neoclassical calcul using NCLASS}'];
if ~isempty('zneo') & isempty(noconnexion)
    try
       par = feval('zneo');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.neo = par.info;
   end
end
param.fonction.rip = [ param.fonction.rip,' {zripple_therm = Calcul des flux ripple thermique + sources associes}'];
if ~isempty('zripple_therm') & isempty(noconnexion)
    try
       par = feval('zripple_therm');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.rip = par.info;
   end
end
param.fonction.fci = [ param.fonction.fci,' {zfcifile = ICRH power deposition (PION/ABSORB)}'];
if ~isempty('zfcifile') & isempty(noconnexion)
    try
       par = feval('zfcifile');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.fci = par.info;
   end
end
param.fonction.fce = [ param.fonction.fce,' {zremafile = ECRH source term, using REMA}'];
if ~isempty('zremafile') & isempty(noconnexion)
    try
       par = feval('zremafile');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.fce = par.info;
   end
end
param.fonction.hyb = [ param.fonction.hyb,' {zhybsimple = Calcul simple du depot de puisance Hybride}'];
if ~isempty('zhybsimple') & isempty(noconnexion)
    try
       par = feval('zhybsimple');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.hyb = par.info;
   end
end
param.fonction.idn = [ param.fonction.idn,' {zsinbad2temps = calcul des sources dues a l''IDN avec le code Simbad}'];
if ~isempty('zsinbad2temps') & isempty(noconnexion)
    try
       par = feval('zsinbad2temps');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.idn = par.info;
   end
end
param.fonction.n0 = [ param.fonction.n0,' {zneutres = Calcul de l''equilibre}'];
if ~isempty('zneutres') & isempty(noconnexion)
    try
       par = feval('zneutres');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.n0 = par.info;
   end
end
param.fonction.bord = [ param.fonction.bord,' {zrecycle = Calcul tres simple du bord et de l''etat du mur}'];
if ~isempty('zrecycle') & isempty(noconnexion)
    try
       par = feval('zrecycle');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.bord = par.info;
   end
end
param.fonction.glacon = [ param.fonction.glacon,' {zglaquelc = Calcul du depot de matiere de l''injection de glacon utilisant GLAQUELC (Parks/Fois)}'];
if ~isempty('zglaquelc') & isempty(noconnexion)
    try
       par = feval('zglaquelc');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.glacon = par.info;
   end
end
param.fonction.fus = [ param.fonction.fus,' {zfusion = Calcul simple des produits de fuison et du chauffage du a alpha}'];
if ~isempty('zfusion') & isempty(noconnexion)
    try
       par = feval('zfusion');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.fus = par.info;
   end
end
param.fonction.cyclo = [ param.fonction.cyclo,' {zcytran77 = Compute the cycltronic radiation with CYTRAN code}'];
if ~isempty('zcytran77') & isempty(noconnexion)
    try
       par = feval('zcytran77');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.cyclo = par.info;
   end
end
param.fonction.coefa = [ param.fonction.coefa,' {zbgbs = Modele de coefficients de transport Bohm/Gyro Bohm }'];
if ~isempty('zbgbs') & isempty(noconnexion)
    try
       par = feval('zbgbs');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.coefa = par.info;
   end
end
param.fonction.plot = [ param.fonction.plot,' {zplot = Plot interratif pour Cronos-Zineb}'];
if ~isempty('zplot') & isempty(noconnexion)
    try
       par = feval('zplot');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.plot = par.info;
   end
end
param.fonction.mhd.dds = [ param.fonction.mhd.dds,' {zddscrash = Calcul du crash des DDS }'];
if ~isempty('zddscrash') & isempty(noconnexion)
    try
       par = feval('zddscrash');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.mhd.dds = par.info;
   end
end
param.fonction.mhd.limite = [ param.fonction.mhd.limite,' {zlims1 = limite de stabilite des DDS (modele du shear critique F.Porcelli)}'];
if ~isempty('zlims1') & isempty(noconnexion)
    try
       par = feval('zlims1');
    catch
        par.info='';
    end
   if ~isempty(par.info)
      param.cons.mhd.limite = par.info;
   end
end
info.param =param;
clear param;
info.data  =data;
clear data;
