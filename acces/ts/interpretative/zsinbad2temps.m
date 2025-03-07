% ZSINBAD calcul des sources dues a l'IDN avec le code Simbad
%--------------------------------------------------------------
% fichier zsinbad2temps.m ->  zsinbad
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule les sources dues a l'IDN avec le code Sinbad.
% 
% syntaxe  :
%  
%      [sortie,matiere,memoire] = zsinbad2temps(param,proto,cons,geo,equi,injection, ...
%                                     prof,neo,impur,phy,composition,gene,memoire);
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
%      phy             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes ( param.compo)
%      gene            =    parametres generaux (param.gne)
%      memoire         =    structure des dernieres valeurs calculees
% 
% sortie :
% 
%     sortie           =  structure de type source remplie par le module (sortie === proto)
%     matiere          =  source de matiere par espece
%     memoire          =  datak.memoire.idn (valeur de reference pour le dernier calcul complet, 
%                         pas utiliser dans cette fonction, reserve pour d'autres modules)
%
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.0 du 05/02/2003.
% 
% 
% liste des modifications : 
% 19/11/2002 -> ajout en sortie des flux de neutrons
% 09/12/2002 -> interface en anglais
% 14/13/2003 -> debuggage du fokkerplanck (sortie des fonctions s2d, s2t, ....)
% 15/13/2003 -> introduction ASDEX
% 05/02/2003 -> nouvelle version fokker
%
%--------------------------------------------------------------
%
function [sortie,particules,memoire] = zsinbad2temps(param,proto,cons,geo,equi,injection, ...
                                          prof,neo,impur,phys,composition,gene,memoire)
                             
                                                         
% mode initialisation 
% fonction auto declarante
langue                  = getappdata(0,'langue_cronos');
if nargin <=1
	if nargin ==0
		nbidn=1;
	else
		nbidn=param;
	end


	valeur.align            = zeros(1,nbidn);         %  position des injecteurs pour Jet: 0 -> standard, 1 -> upshift
	valeur.energie          = 80e3 .* ones(1,nbidn);  %  energie du faiceau en eV
	valeur.fraction1        = 0.8 .* ones(1,nbidn);   %  fraction a l'energie entiere
	valeur.fraction2        = 0.15 .* ones(1,nbidn);  %  fraction a l'energie 1/2
	valeur.fraction3        = 0.05 .* ones(1,nbidn);  %  fraction a l'energie 1/3
	valeur.charge           = 1 .* ones(1,nbidn);     %  nombre de charge de l'espece injecte
	valeur.masse            = 2 .* ones(1,nbidn);     %  masee de l'espece injecte (ua)
	valeur.type             = 3;                      %  type d'interraction faisceau plasma : 3 -> D-D
	valeur.geom             = 2;                      %  geometrie utilisee [0 -> R0,a,d,e 1-> Volume CRONOS, Surface R0, a 2-> Volume et surface CRONOS] {2}
	valeur.backcur          = 'LIN-LIU';                      %  geometrie utilisee [0 -> R0,a,d,e 1-> Volume CRONOS, Surface R0, a 2-> Volume et surface CRONOS] {2}
	valeur.angle            = 2.268;                  %  angle d'injectiojn variable pour ITER
	valeur.debut            = 'Before';               %  debut de simulation pendant IdN
	valeur.fokker           = 'old';               %  debut de simulation pendant IdN
	valeur.depot           = 'full';               %  debut de simulation pendant IdN
	valeur.fileversion      = 0;    % flag=0 to run with mexfiles, =1 for I/O via files
        valeur.spot_use         = 'No'; % flag to choose if SINBAD used with SPOT or not
	valeur.spot_rlong       = 21;   % size of output profiles in SPOT
	valeur.spot_nproc       = 1;    % number of processors used in SPOT
	valeur.spot_ncreated    = 10;   % number of MC particles created at each SPOT time step
	valeur.spot_verbose     = 'No'; % do you want SPOT to be talkative?
	valeur.spot_init        = 'Zero';
	valeur.spot_mpi         = 'No';
	valeur.spot_save        = 'No';
	valeur.spot_anomalous   = 'No';
	valeur.spot_ano_dcoef   = 0.01;
	
	if exist('machine_nom')
	  valeur.machine_nom          = machine_nom; %  nom de la machine_nom (utiliser pour chosir la description geometrique des injecteurs)
 	else
	  valeur.machine_nom          = 'JET'; %  nom de la machine_nom (utiliser pour chosir la description geometrique des injecteurs)
	end
	valeur.geometrie        = 1;
	valeur.save             = 'No';	% sauvegarde des fichiers contextes
	valeur.orbit             = 'No';% sauvegarde des fichiers contextes
	valeur.zinterpolation    = '';	% interpolation de sinbad pour asservissment

	type.align              = 'integer'; % type du parametre
	type.energie            = 'float';
	type.fraction1          = 'float';
	type.fraction2          = 'float';
	type.fraction3          = 'float';
	type.charge             = 'integer';
	type.masse              = 'integer';
	type.type               = 'integer';
	type.geom               = 'integer';
	type.backcur              = 'string';
	type.angle              = 'float';
	type.debut              = 'string';
	type.fokker              = 'string';
	type.depot              = 'string';
	type.fileversion        = 'integer';
	type.spot_use           = 'string';
	type.spot_rlong         = 'real';    % size of output profiles in SPOT
	type.spot_nproc         = 'integer'; % nb of processors used in SPOT
	type.spot_ncreated      = 'real';    % nb of MC particles created at each SPOT time step
	type.spot_verbose       = 'string';  % do you want SPOT to be talkative?		
	type.spot_init          = 'string'; 
	type.spot_mpi           = 'string'; 
	type.spot_save          = 'string';
	type.spot_anomalous     = 'string';
	type.spot_ano_dcoef     = 'real';
	type.machine_nom        = 'string';
	type.geometrie          = 'integer';
	type.save               = 'string';
	type.orbit              = 'string';
	type.zinterpolation     = 'string';

	borne.align             = {0,1};     % bornes ou ensemble de valeur du parametre
	borne.energie           = [1,2e6];
	borne.fraction1         = [0,1];
	borne.fraction2         = [0,1];
	borne.fraction3         = [0,1];
	borne.charge            = {1,2};
	borne.masse             = {1,2,3,4};
	borne.type              = {1,2,3,4,5};
	borne.geom              = {0,1,2};
	borne.backcur             = {'old','LIN-LIU'};
 	borne.angle             = [2.268 3.35];
	borne.debut             = {'Before','During'};
	borne.fokker             = {'old','new'};
	borne.depot             = {'full','half'};
	borne.fileversion       = {0,1}; % flag=0 to run with mexfiles, =1 for I/O via files
        borne.spot_use          = {'Yes','No'}; % choose if SINBAD used vith/wo SPOT
	borne.spot_rlong        = [10 101];    % size of output profiles in SPOT
	borne.spot_nproc        = [1,100];     % nb of processors used in SPOT
	borne.spot_ncreated     = [1,1000];  % nb of MC particles created at each SPOT time step
	borne.spot_verbose      = {'Yes','No'};% do you want SPOT to be talkative?		
	borne.spot_init         = {'Steady-state','Zero'};
	borne.spot_mpi          = {'Yes','No'};% is MPI installed on your machine?
	borne.spot_save         = {'Yes','No'};
	borne.spot_anomalous    = {'Yes','No'};
	borne.spot_ano_dcoef    = [0,10];
	borne.machine_nom       = {'JET','TS','ITER','ASDEX','DIIID','DEMO'};
	borne.geometrie         = {1,2};
	borne.save              = {'Yes','No'};
	borne.orbit             = {'Yes','No'};
	borne.zinterpolation    = {'','zsinbad2temps','zinterpsinbad'};

	defaut.align            = 0;                      % valeur par defautdu parametre
	defaut.energie          = 80;
	defaut.fraction1        = 0.8;
	defaut.fraction2        = 0.15;
	defaut.fraction3        = 0.05;
	defaut.charge           = 1;
	defaut.masse            = 2;
	defaut.type             = 3;
	defaut.geom             = 2;
	defaut.backcur            = 'LIN-LIU';
	defaut.angle            = 2.268;
	defaut.debut            = 'Before';
	defaut.fokker            = 'old';
	defaut.depot            = 'full';
	defaut.fileversion      = 0;
	defaut.spot_use         = 'No';
	defaut.spot_rlong       = 21;
	defaut.spot_nproc       = 1;
	defaut.spot_ncreated    = 10;
	defaut.spot_verbose     = 'Yes';
	defaut.spot_init        = 'Zero';
	defaut.spot_mpi         = 'No';
	defaut.spot_save        = 'No';
	defaut.spot_anomalous   = 'No';
	defaut.spot_ano_dcoef   = 0.01;
	defaut.machine_nom      = 'JET';
	defaut.geometrie        = 1;
	defaut.save             = 'No';
	defaut.orbit            = 'No';
	defaut.zinterpolation   = '';


	info.align              = 'upshift of the PINIS: 0 -> standard, 1 -> upshift'; % information sur le parametre
	info.energie            = 'beam eneergy (eV)';
	info.fraction1          = 'fraction of plain energy ions';
	info.fraction2          = 'fraction of half energy ions';
	info.fraction3          = 'fraction of third energy ions';
	info.charge             = 'Z of injected ion';
	info.masse              = 'A of injected ion';
	info.type               = 'beam plasma interaction : 1 D+T => 4HE+N, 2 D+D => T+P, 3 D+D => 3HE+N, 4 T+T => 4HE+2N, 5 D+3HE => 4HE+P';
	info.geom               = 'geometry [0 -> simple one (R0,a,d,e) 1-> Volume deduced from CRONOS, Sipmle area (R0, a) 2-> Geometry of CRONOS {2}';
	info.backcur               = 'back current law, Lin-Liu new formula for general tokamak equilibria';
	info.angle               = 'for ITER geometric parameters, T. Oikawa, tilting angle';
	info.debut              = 'NBI is present before the simulation or not';
	info.fokker              = 'old version, 11 radial steps';
	info.depot              = 'half, symetrical deposition in Z';
	info.machine_nom        = 'Tokamak name ';
	info.geometrie          = 'for the geometry of the injector, only for ITER (1-> standard, 2 -> allow change of tiltingangle)';
	info.fileversion        = 'flag=0 to run with mexfiles, flag=1 for I/O via files';
	info.spot_use           = 'YES to run SINBAD with SPOT, NO to run without SPOT';
	info.spot_rlong         = 'Size of radial output array'; 
	info.spot_nproc         = 'Number of processors used to run SPOT (if MPI is available)';
	info.spot_ncreated      = 'Number of MC particles generated at each SPOT time step';
	info.spot_verbose       = 'Flag to get a talkative SPOT or not';
	info.spot_init          = 'Status of fast ion particle distribution function for the first call to the module';
	info.spot_mpi           = 'YES if MPI installed on your machine, NO in the other case';
	info.spot_save          = 'Save context of SPOT for tests (WARNING: it fills the disk very quickly....)';
	info.spot_anomalous     = 'YES to account for the fast ion anomalous transport in SPOT';
	info.spot_ano_dcoef     = 'Diffusion coefficient for fast ion anomalous transport in SPOT';
	info.machine_nom        = 'Tokamak name (for the geometry of the injector)';
	info.save               = 'save context for test';
	info.orbit              = 'smooth the current on magnetic axis due to orbit width effect';
	info.zinterpolation     = 'interpolation de zsinbad2temps pour asservissement';
	
	interface.ts = '';   % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';  % nom de la fonction d'interfacage avec les donnees Jet

	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'calcul des sources dues a l''IDN avec le code Simbad';   % description (une ligne) de la fonction

	sortie.help = '';     % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';      % nom de l'interface graphique specifique si elle existe
	sortie.controle = ''; % nom de la fonction de controle des valeurs si elle existe

	return
end

%% If we want SINBAD to be used with SPOT, we have to use the
%% SINBAD version with I/O via files (i.e. fileversion = 1)
if strcmp(param.spot_use,'Yes')
  param.fileversion = 1;
end

if isfield(param,'machine_nom')
           machine_nom = param.machine_nom;
else
          machine_nom = [];
end
if ~isfield(param,'zinterpolation')
  param.zinterpolation='';
end
if ~isfield(param,'geometrie')
	  param.geometrie = 1;
end

% debut du calcul
sortie     = proto;
particules = zeros(1,gene.nbrho,gene.nbg);

% petits vecteurs utils
va = ones(gene.nbhyb,1);
ve = ones(1,gene.nbrho);


% donnees le dernier temps
times   = gene.t;
dt      = gene.dt;
% la geometrie est fixe sauf rhomax :
R0      = geo.r0;
Z       = geo.z0;
a       = equi.a(end);

% patch provisioir sinbad 
e0      = min(1.8,equi.e(1));

e1      = equi.e(end);
d0      = equi.d(1);
tri     = 0.5 .* (equi.trl(end) + equi.trh(end)); % traingularite en m


% courant plasma
Ip      = equi.ip ./ 1e6;  % A -> MA
% fluxn -> pas utiliser
fluxn   = 1e13;  
% calcul des valeur avec evolution (zevols) a partir de la donee et de sa derivee
ne      = prof.ne./1e19;                                         % densite electronique (1019 m-3)
ae      = prof.ae;                                               % rapport ne/ni
ni      = ae .* ne;                                              % densite ionique totale
rot     = prof.rot;
%
% attention pas de NaN dans rot (pour les abrutis dixit Mireille
% V. Basiuk, 6 septembre 2005, apres une pause cafe
%
if isnan(rot)
  rot = zeros(size(rot));
end
% rotation toroidal echange de chargedu plasma (m/s)
Te      = prof.te ./ 1e3 ;                                       % temperature electronique en keV
Te      = zbornes(Te,13.6e-3,100,13.6e-3);
Ti      = prof.ti ./ 1e3 ;                                       % temperature ionique en keV
Ti      = zbornes(Ti,13.6e-3,100,13.6e-3);
% facteur d'evolution des especes ionique dans le plasma (meme rapport de composition et zeff)
ni      = squeeze(impur.impur);
n1      = squeeze(impur.impur(1,:,1)) ./ 1e19;    % densite ionique premiere espece
n2      = squeeze(impur.impur(1,:,2)) ./ 1e19;    % densite ionique seconde espece
n3      = squeeze(impur.impur(1,:,3)) ./ 1e19;    % densite ionique troisieme espece
n4      = squeeze(impur.impur(1,:,4)) ./ 1e19;    % densite ionique quatrieme espece
Zeff    = impur.zeff;                             % profil de la charge effective
ftrap   = equi.ftrap;
% calcul des consigne avec evolution
p0      = real(cons);
dp0dt   = imag(cons);
% temps de thermalisation
lambda    = 39.1 - 0.5 * log(prof.ne) + log(prof.te ./ 1e3);
Ad        = mean(prof.ne) .* phys.e^4 .* mean(lambda) / 2 / pi / phys.epsi0^2 / (phys.mp*param.masse(1))^2;
tause     = 3*sqrt(2*pi) ./ sqrt(phys.me) ./ (phys.mp*param.masse(1)) ./ Ad .* mean(prof.te*phys.e)^(3/2);
zmoy      = (sum(mean(ni).*(compo.z).^2./compo.a)./mean(prof.ne)).^(2/3);
Ecrit     = 14.8 * param.masse(1) * mean(prof.te) * zmoy;
taus      = tause/3*log(1+(mean(param.energie)./Ecrit)^(3/2));
%taus      = max(6.27e8 .* 2 .* (prof.te .^ 1.5) ./ (lambda .* prof.ne ./1e6));

% rayon externe
rext = interp1(double(equi.rhoRZ),max(double(squeeze(equi.R))')',gene.x .* equi.rhomax,'pchip');
% recupere les donees globale de pion
%sinbad_data = getappdata(0,'sinbad_data');
if isfield(memoire.data,'input')
	sinbad_data = memoire.data.input;
else
	sinbad_data = [];
end
%
%   application de taus pour le dt
%
if dt > taus
 
  dt        = taus/5;
  ncoup     = 1; 
  
else

  ncoup     = 1;
end 

% si sinbad_data est vide on cree la structure
if isempty(sinbad_data)
   % creation de la base temps fictive sur ntemps
   ntemps                    = 2;
   sinbad_data.times         = [times times+dt];
   sinbad_data.dt(1:2)       = dt;
   sinbad_data.R0            = ones(ntemps,1) * R0;
   sinbad_data.Z             = ones(ntemps,1) * Z;
   sinbad_data.a             = ones(ntemps,1) * a;
   sinbad_data.e0            = ones(ntemps,1) * e0;
   sinbad_data.e1            = ones(ntemps,1) * e1;
   sinbad_data.d0            = ones(ntemps,1) * d0;
   sinbad_data.tri           = ones(ntemps,1) * tri;
   sinbad_data.Ip            = ones(ntemps,1) * Ip;
   sinbad_data.fluxn         = ones(ntemps,1) * fluxn;
   sinbad_data.ne            = ones(ntemps,1) * ne;
   sinbad_data.ni            = ones(ntemps,1) * ni;
   sinbad_data.Te            = ones(ntemps,1) * Te;
   sinbad_data.rot           = ones(ntemps,1) * rot;
   sinbad_data.Ti            = ones(ntemps,1) * Ti;
   sinbad_data.n1            = ones(ntemps,1) * n1;
   sinbad_data.n2            = ones(ntemps,1) * n2;
   sinbad_data.n3            = ones(ntemps,1) * n3;
   sinbad_data.n4            = ones(ntemps,1) * n4;
   sinbad_data.Zeff          = ones(ntemps,1) * Zeff;
   sinbad_data.ftrap         = ones(ntemps,1) * ftrap;
%   sinbad_data.p0            = ones(ntemps,1) * p0 + (sinbad_data.times - max(sinbad_data.times)) * dp0dt;
%   sinbad_data.p0(sinbad_data.p0 < 0) = 0;
   if ~isfield(param,'Before')
	  param.debut = 'Before';
   end
   if strcmp(lower(param.debut),lower('Before'))
	  sinbad_data.p0            = [zeros(ntemps-1,1)*p0;p0];
   elseif strcmp(lower(param.debut),lower('During'))
	  sinbad_data.p0            = [ones(ntemps-1,1)*p0;p0];
   end

   sinbad_data.taus          = ones(ntemps,1) * taus;
   sinbad_data.rext          = ones(ntemps,1) * rext;

else
   sinbad_data.times(end+1)  = sinbad_data.times(end)+dt;
   sinbad_data.dt(end+1)     = dt;
   sinbad_data.R0(end+1)     = R0;
   sinbad_data.Z(end+1)      = Z;
   sinbad_data.a(end+1)      = a;
   sinbad_data.e0(end+1)     = e0;
   sinbad_data.e1(end+1)     = e1;
   sinbad_data.d0(end+1)     = d0;
   sinbad_data.tri(end+1)    = tri;
   sinbad_data.Ip(end+1)     = Ip;
   sinbad_data.fluxn(end+1)  = fluxn;
   sinbad_data.ne(end+1,:)   = ne;
   sinbad_data.ni(end+1,:)   = ni;
   sinbad_data.Te(end+1,:)   = Te;
   sinbad_data.rot(end+1,:)  = rot;
   sinbad_data.Ti(end+1,:)   = Ti;
   sinbad_data.n1(end+1,:)   = n1;
   sinbad_data.n2(end+1,:)   = n2;
   sinbad_data.n3(end+1,:)   = n3;
   sinbad_data.n4(end+1,:)   = n4;
   sinbad_data.Zeff(end+1,:) = Zeff;
   sinbad_data.ftrap(end+1,:) = ftrap;
   sinbad_data.p0(end+1,:)   = p0;
   sinbad_data.taus(end+1)   = taus;
   sinbad_data.rext(end+1,:) = rext;
end
   % et on retire les vieux temps
nmax  = 2;
if length(sinbad_data.times) > nmax
   sinbad_data.times = sinbad_data.times(2:end);
   sinbad_data.dt    = sinbad_data.dt(2:end);
   sinbad_data.R0    = sinbad_data.R0(2:end);
   sinbad_data.Z     = sinbad_data.Z(2:end);
   sinbad_data.a     = sinbad_data.a(2:end);
   sinbad_data.e0    = sinbad_data.e0(2:end);
   sinbad_data.e1    = sinbad_data.e1(2:end);
   sinbad_data.d0    = sinbad_data.d0(2:end);
   sinbad_data.tri   = sinbad_data.tri(2:end);
   sinbad_data.Ip    = sinbad_data.Ip(2:end);
   sinbad_data.fluxn = sinbad_data.fluxn(2:end);
   sinbad_data.ne    = sinbad_data.ne(2:end,:);
   sinbad_data.ni    = sinbad_data.ni(2:end,:);
   sinbad_data.Te    = sinbad_data.Te(2:end,:);
   sinbad_data.rot   = sinbad_data.rot(2:end,:);
   sinbad_data.Ti    = sinbad_data.Ti(2:end,:);
   sinbad_data.n1    = sinbad_data.n1(2:end,:);
   sinbad_data.n2    = sinbad_data.n2(2:end,:);
   sinbad_data.n3    = sinbad_data.n3(2:end,:);
   sinbad_data.n4    = sinbad_data.n4(2:end,:);
   sinbad_data.Zeff  = sinbad_data.Zeff(2:end,:);
   sinbad_data.ftrap  = sinbad_data.ftrap(2:end,:);
   sinbad_data.p0    = sinbad_data.p0(2:end,:);
   sinbad_data.taus  = sinbad_data.taus(2:end);
   sinbad_data.rext  = sinbad_data.rext(2:end,:);
end
%
% champ magnetique total fonction de x et theta
%
tempsB     = times(end);
tcron      = gene.t;
xcron      = gene.x;
indcron    = iround(tcron,tempsB);
Rcron      = double(squeeze(equi.R(indcron,:,:)));
Zcron      = double(squeeze(equi.Z(indcron,:,:)));
Rcentre    = mean(Rcron,2)*ones(1,size(Rcron,2));
Zcentre    = mean(Zcron,2)*ones(1,size(Rcron,2));
thetacron  = angle((Rcron-Rcentre)+i*(Zcron-Zcentre));
thetacron  = unwrap(thetacron')';
BRcron     = double(squeeze(equi.BR(indcron,:,:)));
BZcron     = double(squeeze(equi.BZ(indcron,:,:)));
Bphicron   = double(squeeze(equi.BPHI(indcron,:,:)));
Bcron      = sqrt(BRcron.^2+BZcron.^2+Bphicron.^2);
BPcron     = sqrt(BRcron.^2+BZcron.^2);
Boutcron   = Bcron(:,1);

if sum(real(cons)) < 3e5
  sinbad_data.option.init = 0;
elseif ~isfield(sinbad_data,'option')
  if ~isfield(memoire.data,'input')
    if strcmp(lower(param.debut),lower('During'))
      sinbad_data.option.init = 2;
    elseif strcmp(lower(param.debut),lower('Before'))
      sinbad_data.option.init = 0;
    end

  end
end


% memorisation des donees globale de pion
memoire.t           = gene.t;
memoire.data.input  = sinbad_data;
%setappdata(0,'sinbad_data',sinbad_data);


% test presence de puissance
if sum(real(cons)) < 3e5
     disp('injected Power too low for SINBAD')
     % pas de puissance
     return
end

% creation des donnees pour sinbad
% creation des variables
% parametre reechantillonage:
par.ondelette        = 0;
par.defaut.temps     = 0;
par.defaut.espace    = 0;
par.defaut.inf       = 0;
par.plus             = 1;
par.energie          = 1;
%
% base temps pour sinbad
%
ntact               = 2;
%
% choix d'intervalle temporelle <= a 100 ms
% blindage du nombre d'intervalle <= 1000 pour SINDAB
%

temps   = [sinbad_data.times(1) sinbad_data.times(2)];

vt      = ones(size(temps));
  
% la geometrie est fixe sauf rhomax :
R0      = sinbad_data.R0';
Z       = sinbad_data.Z';
a       = sinbad_data.a';
e0      = sinbad_data.e0';
e1      = sinbad_data.e1';
d0      = sinbad_data.d0';
tri     = sinbad_data.tri';
% courant plasma
Ip      = sinbad_data.Ip';
% composition du plasma
mass    = composition.a(1:4);
charg   = composition.z(1:4);
% fluxn -> pas utiliser
fluxn   = sinbad_data.fluxn';
% coordonnees (x)
rho     = gene.x;
rhoin     = gene.x;
% calcul des valeur avec evolution (zevols) a partir de la donee et de sa derivee
ne      = sinbad_data.ne';
ni      = sinbad_data.ni';
Te      = sinbad_data.Te';
rot     = sinbad_data.rot';
Ti      = sinbad_data.Ti';

% facteur d'evolution des especes ionique dans le plasma (meme rapport de composition et zeff)
n1      = sinbad_data.n1';
n2      = sinbad_data.n2';
n3      = sinbad_data.n3';
n4      = sinbad_data.n4';
Zeff    = sinbad_data.Zeff';
ftrap    = sinbad_data.ftrap';
% rayon externe des surfaces de flux
rext    = sinbad_data.rext';


% calcul des consigne avec evolution
vpini   = 1:size(sinbad_data.p0,2);
pcons   = sinbad_data.p0;

% selon la machine_nom
if strcmp(upper(param.machine_nom),'JET')
	% premier groupe
	align1      =  param.align(1:8);      % alignement
	en1         =  param.energie(1:8)./1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z1          =  param.charge(1:8);        % nombre de charge des especes injectees
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	% deuxieme groupe
	align2      =  param.align(9:16);      % alignement
	en2         =  param.energie(9:16)./1e3;    % energie du faisceau
	fr2         =  cat(2,param.fraction1(9:16)',param.fraction2(9:16)',param.fraction3(9:16)'); % fraction de 1, 1/2 1/3 de l'energie
	A2          =  param.masse(9:16);        % masse des especes injectee (ua)
	Z2          =  param.charge(9:16);        % nombre de charge des especes injectees
	pin2        =  (pcons(:,9:16) ./ 1e6)' ;     % consigne de puisance en MW
	% aurtres parametres
	type        =  param.type;             % type interraction (3= DD)
	neu         =  1;                      % calcul des neutrons
   option.iter = 0;
elseif strcmp(upper(param.machine_nom),'TS')
    error ('La machine_nom TS n''est pas deja implantee dans sinbad')
elseif strcmp(upper(param.machine_nom),'ITER') & param.geometrie == 1
	align1      =  param.align(1:8);      % alignement
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	en1         =  param.energie(1:8)./1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z1          =  param.charge(1:8);        % nombre de charge des especes injectees
	align2      =  param.align(9:16);      % alignement
	A2          =  param.masse(9:16);        % masse des especes injectee (ua)
	Z2          =  param.charge(9:16);        % nombre de charge des especes injectees
	en2         =  param.energie(9:16)./1e3;    % energie du faisceau
	fr2         =  cat(2,param.fraction1(9:16)',param.fraction2(9:16)',param.fraction3(9:16)'); % fraction de 1, 1/2 1/3 de l'energie
	pin2        =  (pcons(:,9:16) ./ 1e6)' ;     % consigne de puisance en MW
	type        =  param.type;             % type interraction (1= DD)
	neu         =  1;                      % calcul des neutrons
   option.iter = 1;
elseif strcmp(upper(param.machine_nom),'ASDEX')
	% premier injecteur
	align1      =  param.align(1:8);      % alignement
	en1         =  param.energie(1:8)./1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1);        % masse des especes injectee (ua)
	Z1          =  param.charge(1);        % nombre de charge des especes injectees
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	% deuxieme injecteur
	align2      =  param.align(9:16);      % alignement
	en2         =  param.energie(9:16)./1e3;    % energie du faisceau
	fr2         =  cat(2,param.fraction1(9:16)',param.fraction2(9:16)',param.fraction3(9:16)'); % fraction de 1, 1/2 1/3 de l'energie
	A2          =  param.masse(9);        % masse des especes injectee (ua)
	Z2          =  param.charge(9);        % nombre de charge des especes injectees
	pin2        =  (pcons(:,9:16) ./ 1e6)' ;     % consigne de puisance en MW
	% aurtres parametres
	type        =  param.type;             % type interraction (3= DD)
	neu         =  1;                      % calcul des neutrons
   option.iter = 2;
elseif strcmp(upper(param.machine_nom),'DIIID')
	% premier injecteur
	align1      =  zeros(1,8);      % alignement
	en1         =  param.energie(1:8)/1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z1          =  param.charge(1:8);        % nombre de charge des especes injectees
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	% deuxieme injecteur
	align2      =  zeros(1,8);      % alignement
	en2         =  zeros(1,8)/1e3;    % energie du faisceau
	fr2         =  0*fr1; % fraction de 1, 1/2 1/3 de l'energie
	A2          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z2          =  param.charge(1:8);        % nombre de charge des especes injectees
	pin2        =  0*pin1;     % consigne de puisance en MW
	% aurtres parametres
	type        =  param.type;             % type interraction (3= DD)
	neu         =  1;                      % calcul des neutrons
     	option.iter = 3;
elseif strcmp(upper(param.machine_nom),'DEMO')
	align1      =  param.align(1:8);      % alignement
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	en1         =  param.energie(1:8)./1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z1          =  param.charge(1:8);        % nombre de charge des especes injectees
	align2      =  param.align(9:16);      % alignement
	A2          =  param.masse(9:16);        % masse des especes injectee (ua)
	Z2          =  param.charge(9:16);        % nombre de charge des especes injectees
	en2         =  param.energie(9:16)./1e3;    % energie du faisceau
	fr2         =  cat(2,param.fraction1(9:16)',param.fraction2(9:16)',param.fraction3(9:16)'); % fraction de 1, 1/2 1/3 de l'energie
	pin2        =  (pcons(:,9:16) ./ 1e6)' ;     % consigne de puisance en MW
	type        =  param.type;             % type interraction (1= DD)
	neu         =  1;                      % calcul des neutrons
    	option.iter = 4;
elseif strcmp(upper(param.machine_nom),'ITER') & param.geometrie == 2
	align1      =  param.align(1:8);      % alignement
	pin1        =  (pcons(:,1:8) ./ 1e6)' ;     % consigne de puisance en MW
	en1         =  param.energie(1:8)./1e3;    % energie du faisceau
	fr1         =  cat(2,param.fraction1(1:8)',param.fraction2(1:8)',param.fraction3(1:8)'); % fraction de 1, 1/2 1/3 de l'energie
	A1          =  param.masse(1:8);        % masse des especes injectee (ua)
	Z1          =  param.charge(1:8);        % nombre de charge des especes injectees
	align2      =  param.align(9:16);      % alignement
	A2          =  param.masse(9:16);        % masse des especes injectee (ua)
	Z2          =  param.charge(9:16);        % nombre de charge des especes injectees
	en2         =  param.energie(9:16)./1e3;    % energie du faisceau
	fr2         =  cat(2,param.fraction1(9:16)',param.fraction2(9:16)',param.fraction3(9:16)'); % fraction de 1, 1/2 1/3 de l'energie
	pin2        =  (pcons(:,9:16) ./ 1e6)' ;     % consigne de puisance en MW
	type        =  param.type;             % type interraction (1= DD)
	neu         =  1;                      % calcul des neutrons
    	option.iter = 5;
else
    error('Not avialable tokamak')
end


% -- consignes injection de neutres --
% pin1  : puissance de l'injecteur numero 1 de JET (8 sources)
%         ou source de neutres TS (MW)(base temps), dim (source,temps)
% -- parametre injection de neutres --
% en1   : energie pleine injectee (8 pour JET, 1 pour TS)
% fr1   : fraction des energies pleines, 1/2 et 1/3 injectee en %
%         (8*3 pour JET, 1*3 pour TS)
% A1    : masse du gaz injecte
% Z1    : charge du gaz injecte
% align1: seulement pour JET, alignement standard ou upshifted des faisceaux
%
% type  : 3 -> reactions D-D
% neu    : 0 -> pas de flux de neutron, 1-> flux de neutrons

% sauvegarde des donnees de sinbad pour debbugage
cs       = 0;
% appel de sinbad
%
% donnees concernant l'equilibre sur le dernier temps
%
equis.R      	= double(squeeze(equi.R(indcron,:,:)));
equis.Z      	= double(squeeze(equi.Z(indcron,:,:)));
equis.vpr    	= squeeze(equi.vpr(indcron,:,:));
equis.rhoRZ  	= double(equi.rhoRZ(indcron,:));
equis.a      	= equi.a(indcron,:);
equis.rhomax 	= equi.rhomax(indcron);
geos.z0      	= geo.z0(indcron);
genes.x      	= gene.x;
[chemin,void] 	= fileparts(gene.rapsauve);
genes.chemin 	= chemin;
genes.orbit 	= param.orbit;

option.geo   	= param.geom;
if ~isfield(sinbad_data,'option')
  sinbad_data.option.init=0;
end
option.init  	= sinbad_data.option.init;
if option.iter == 4
  option.init = 0;
end
if option.init == 0
  dis	= [];
  dist	= [];
  xv	= [];
  ev	= [];
  mat	= [];
  s2d_1	= [];
  s2d_2	= [];
  s2d_3	= [];
  s2d_4	= [];
  s2d_5	= [];
  s2d_6	= [];
  s2d_7	= [];
  s2d_8	= [];
  s2d_9	= [];
  s2d_10= [];
  s2t	= [];
  divers= [];
%
% mis en forme des donnees
%
  plas.temps 	= temps;
  plas.R0	= R0;
  plas.a 	= a;
  plas.Z	= Z;
  plas.e0 	= e0;
  plas.e1	= e1;
  plas.d0	= d0;
  plas.tri	= tri;
  plas.Ip	= Ip;
  plas.mass	= mass;
  plas.charg	= charg;
  plas.fluxn	= fluxn;
  plas.rho 	= rho;
  plas.Ti 	= Ti;
  plas.Te 	= Te;
  plas.ne 	= ne;
  plas.n1	= n1;
  plas.n2	= n2;
  plas.n3	= n3;
  plas.n4	= n4;
  plas.Zeff 	= Zeff;
  plas.ftrap 	= ftrap;
  plas.ftlaw 	= param.backcur;
  fais.pin1	= pin1;
  fais.en1	= en1;
  fais.fr1	= fr1;
  fais.A1	= A1;
  fais.Z1	= Z1;
  fais.align1	= align1;
  fais.pin2	= pin2;
  fais.en2	= en2;
  fais.fr2	= fr2;
  fais.A2	= A2;
  fais.Z2	= Z2;
  fais.align2 	= align2;
  fais.type	= type;
  fais.neu	= neu;
  fais.rext	= rext;
  fais.cs	= cs;
  fais.angle    = param.angle;
  plas.Bcron	= Bcron;
  plas.Boutcron = Boutcron;
  plas.Bphicron	= Bphicron;
  plas.BPcron	= BPcron;
  plas.thetacron= thetacron;
  plas.xcron 	= xcron;
  plas.Rcron 	= Rcron;
  para.rot	= rot;
  para.option	= option;
  plas.equis	= equis;
  para.genes	= genes;
  plas.geo 	= geos;
  plas.rext	= rext;
  plas.fok      = param.fokker;
  plas.dep      = param.depot;
  RAD		= [0 0.1 0.15 0.2 0.25 0.35 0.5 0.65 0.75 0.9 1];
  plas.vpr      = interp1(gene.x,equi.vpr,RAD);
  plas.rhomax   = equi.rhomax;
  sort.dis 	= dis;
  sort.dist 	= dist;
  sort.xv	= xv;
  sort.ev	= ev;
  sort.mat 	= mat;
  sort.s2d_1 	= s2d_1;
  sort.s2d_2 	= s2d_2;
  sort.s2d_3 	= s2d_3;
  sort.s2d_4 	= s2d_4;
  sort.s2d_5 	= s2d_5;
  sort.s2d_6 	= s2d_6;
  sort.s2d_7 	= s2d_7;
  sort.s2d_8 	= s2d_8;
  sort.s2d_9 	= s2d_9;
  sort.s2d_10 	= s2d_10;
  sort.s2t	= s2t;
  sort.divers 	= divers;
  if option.iter == 3 | option.iter == 4
    sort.past   = taus/ncoup;
  else
    sort.past   = 0.8;
  end
  sort.dtcron   = gene.dt;
%  disp('zsinbad2temps1')
%  keyboard
  sor  		= nbi_sinbad_2tempssat(plas,fais,para,sort,equi,prof,impur,composition,geo,param,gene,memoire,cons);
  memoire.data.memoryspot = sor.memoryspot;
  memoire.data.input.option.init = 1;
  sor1		= sor.DEP1;
  sor2		= sor.DEP2;
  sor3		= sor.DEP3;
  sor4		= sor.DEP4;
  sor5		= sor.DEP5;
  sor6		= sor.DEP6;
  pion1		= sor.pion1;
  pion2		= sor.pion2;
  pion3		= sor.pion3;
  pion4		= sor.pion4;
  pelec		= sor.pelec;
  ploss		= sor.ploss;
  Jidn		= sor.Jidn;
  RAD		= sor.RAD;
  wpar		= sor.Wpar;
  wtot		= sor.W;
  wth		= sor.Wth;
  ang		= sor.ang;
  dis		= sor.dis;
  dist		= sor.dist;
  xv		= sor.xv;
  ev		= sor.ev;
  neu1		= sor.neu1;
  neu2		= sor.neu2;
  neu3		= sor.neu3;
  neu4		= sor.neu4;
  neu5		= sor.neu5;
  neu6		= sor.neu6;
  mat		= sor.mat;
  s2d_1		= sor.s2d_1;
  s2d_2		= sor.s2d_2;
  s2d_3		= sor.s2d_3;
  s2d_4		= sor.s2d_4;
  s2d_5		= sor.s2d_5;
  s2d_6		= sor.s2d_6;
  s2d_7		= sor.s2d_7;
  s2d_8		= sor.s2d_8;
  s2d_9		= sor.s2d_9;
  s2d_10	= sor.s2d_10;
  s2t		= sor.s2t;
  divers	= sor.divers;

else
  oldpui = sum(memoire.data.input.p0(:))/1e6/2;
  newpui = sum(pin1(:)+pin2(:))/2;
  if oldpui ~= newpui
%    option.init = 2;
  end

  if option.init == 2
    dis		= [];
    dist	= [];
    xv		= [];
    ev		= [];
    mat		= [];
    s2d_1	= [];
    s2d_2	= [];
    s2d_3	= [];
    s2d_4	= [];
    s2d_5	= [];
    s2d_6	= [];
    s2d_7	= [];
    s2d_8	= [];
    s2d_9	= [];
    s2d_10	= [];
    s2t		= [];
    divers	= [];

    option2      = option;
    option2.init = 0;
%
% mis en forme des donnees
%
    plas.temps 		= temps;
    plas.R0		= R0;
    plas.a 		= a;
    plas.Z		= Z;
    plas.e0 		= e0;
    plas.e1		= e1;
    plas.d0		= d0;
    plas.tri		= tri;
    plas.Ip		= Ip;
    plas.mass		= mass;
    plas.charg		= charg;
    plas.fluxn		= fluxn;
    plas.rho 		= rho;
    plas.Ti 		= Ti;
    plas.Te 		= Te;
    plas.ne 		= ne;
    plas.n1		= n1;
    plas.n2		= n2;
    plas.n3		= n3;
    plas.n4		= n4;
    plas.Zeff 		= Zeff;
    plas.ftrap 		= ftrap;
    plas.ftlaw 		= param.backcur;
    fais.pin1		= pin1;
    fais.en1		= en1;
    fais.fr1		= fr1;
    fais.A1		= A1;
    fais.Z1		= Z1;
    fais.align1		= align1;
    fais.pin2		= pin2;
    fais.en2		= en2;
    fais.fr2		= fr2;
    fais.A2		= A2;
    fais.Z2		= Z2;
    fais.align2 	= align2;
    fais.type		= type;
    fais.neu		= neu;
    fais.rext		= rext;
    fais.cs		= cs;
    fais.angle   	= param.angle;
    plas.Bcron		= Bcron;
    plas.Boutcron 	= Boutcron;
    plas.Bphicron 	= Bphicron;
    plas.BPcron		= BPcron;
    plas.thetacron	= thetacron;
    plas.xcron 		= xcron;
    plas.Rcron 		= Rcron;
    para.rot		= rot;
    para.option		= option2;
    plas.equis		= equis;
    para.genes		= genes;
    plas.geo 		= geos;
    plas.rext		= rext;
    plas.fok            = param.fokker
    plas.dep      = param.depot;
    RAD		        = [0 0.1 0.15 0.2 0.25 0.35 0.5 0.65 0.75 0.9 1];
    plas.vpr            = interp1(gene.x,equi.vpr,RAD);
    plas.rhomax         = equi.rhomax;
    sort.dis 		= dis;
    sort.dist 		= dist;
    sort.xv		= xv;
    sort.ev		= ev;
    sort.mat 		= mat;
    sort.s2d_1 		= s2d_1;
    sort.s2d_2 		= s2d_2;
    sort.s2d_3 		= s2d_3;
    sort.s2d_4 		= s2d_4;
    sort.s2d_5 		= s2d_5;
    sort.s2d_6 		= s2d_6;
    sort.s2d_7 		= s2d_7;
    sort.s2d_8 		= s2d_8;
    sort.s2d_9 		= s2d_9;
    sort.s2d_10 	= s2d_10;
    sort.s2t		= s2t;
    sort.divers 	= divers;
    if option.iter == 3
    sort.past     = taus/ncoup;
  else 
    sort.past             = 0.8;
  end
  sort.dtcron   = gene.dt;
    sor  		= nbi_sinbad_2tempssat(plas,fais,para,sort,equi,prof,impur,composition,geo,param,gene,memoire,cons);
    memoire.data.memoryspot = sor.memoryspot;
    memoire.data.input.option.init = 1;
    sor1		= sor.DEP1;
    sor2		= sor.DEP2;
    sor3		= sor.DEP3;
    sor4		= sor.DEP4;
    sor5		= sor.DEP5;
    sor6		= sor.DEP6;
    pion1		= sor.pion1;
    pion2		= sor.pion2;
    pion3		= sor.pion3;
    pion4		= sor.pion4;
    pelec		= sor.pelec;
    ploss		= sor.ploss;
    Jidn		= sor.Jidn;
    RAD			= sor.RAD;
    wpar		= sor.Wpar;
    wtot		= sor.W;
    wth			= sor.Wth;
    ang			= sor.ang;
    dis			= sor.dis;
    dist		= sor.dist;
    xv			= sor.xv;
    ev			= sor.ev;
    neu1		= sor.neu1;
    neu2		= sor.neu2;
    neu3		= sor.neu3;
    neu4		= sor.neu4;
    neu5		= sor.neu5;
    neu6		= sor.neu6;
    mat			= sor.mat;
    s2d_1		= sor.s2d_1;
    s2d_2		= sor.s2d_2;
    s2d_3		= sor.s2d_3;
    s2d_4		= sor.s2d_4;
    s2d_5		= sor.s2d_5;
    s2d_6		= sor.s2d_6;
    s2d_7		= sor.s2d_7;
    s2d_8		= sor.s2d_8;
    s2d_9		= sor.s2d_9;
    s2d_10		= sor.s2d_10;
    s2t			= sor.s2t;
    divers		= sor.divers;

  end
  if isfield(memoire.data,'output')
    dis             = memoire.data.output.dis ;
    dist            = memoire.data.output.dist ;
    xv              = memoire.data.output.xv ;
    ev              = memoire.data.output.ev ;
    mat             = memoire.data.output.mat;
    s2d_1           = memoire.data.output.s2d_1;
    s2d_2           = memoire.data.output.s2d_2;
    s2d_3           = memoire.data.output.s2d_3;
    s2d_4           = memoire.data.output.s2d_4;
    s2d_5           = memoire.data.output.s2d_5;
    s2d_6           = memoire.data.output.s2d_6;
    s2d_7           = memoire.data.output.s2d_7;
    s2d_8           = memoire.data.output.s2d_8;
    s2d_9           = memoire.data.output.s2d_9;
    s2d_10          = memoire.data.output.s2d_10;
    divers          = memoire.data.output.divers;
    s2t             = memoire.data.output.s2t;
  end
%
% mis en forme des donnees
%
  plas.temps 		= temps;
  plas.R0		= R0;
  plas.a 		= a;
  plas.Z		= Z;
  plas.e0 		= e0;
  plas.e1		= e1;
  plas.d0		= d0;
  plas.tri		= tri;
  plas.Ip		= Ip;
  plas.mass		= mass;
  plas.charg		= charg;
  plas.fluxn		= fluxn;
  plas.rho 		= rho;
  plas.Ti 		= Ti;
  plas.Te 		= Te;
  plas.ne 		= ne;
  plas.n1		= n1;
  plas.n2		= n2;
  plas.n3		= n3;
  plas.n4		= n4;
  plas.Zeff 		= Zeff;
  plas.ftrap 		= ftrap;
  plas.ftlaw 		= param.backcur;
  fais.pin1		= pin1;
  fais.en1		= en1;
  fais.fr1		= fr1;
  fais.A1		= A1;
  fais.Z1		= Z1;
  fais.align1		= align1;
  fais.pin2		= pin2;
  fais.en2		= en2;
  fais.fr2		= fr2;
  fais.A2		= A2;
  fais.Z2		= Z2;
  fais.align2 		= align2;
  fais.type		= type;
  fais.neu		= neu;
  fais.rext		= rext;
  fais.cs		= cs;
  fais.angle   	= param.angle;
  plas.Bcron		= Bcron;
  plas.Boutcron 	= Boutcron;
  plas.Bphicron 	= Bphicron;
  plas.BPcron		= BPcron;
  plas.thetacron	= thetacron;
  plas.xcron 		= xcron;
  plas.Rcron 		= Rcron;
  para.rot		= rot;
  para.option		= option;
  plas.equis		= equis;
  para.genes		= genes;
  plas.geo 		= geos;
  plas.rext		= rext;
  plas.fok              = param.fokker
  plas.dep      	= param.depot;
  RAD			= [0 0.1 0.15 0.2 0.25 0.35 0.5 0.65 0.75 0.9 1];
  plas.vpr      	= interp1(gene.x,equi.vpr,RAD);
  plas.rhomax   	= equi.rhomax;
  sort.dis 		= dis;
  sort.dist 		= dist;
  sort.xv		= xv;
  sort.ev		= ev;
  sort.mat 		= mat;
  sort.s2d_1 		= s2d_1;
  sort.s2d_2 		= s2d_2;
  sort.s2d_3 		= s2d_3;
  sort.s2d_4 		= s2d_4;
  sort.s2d_5 		= s2d_5;
  sort.s2d_6 		= s2d_6;
  sort.s2d_7 		= s2d_7;
  sort.s2d_8 		= s2d_8;
  sort.s2d_9 		= s2d_9;
  sort.s2d_10 		= s2d_10;
  sort.s2t		= s2t;
  sort.divers 		= divers;
    if option.iter == 3
    sort.past     = taus/ncoup;
  else 
    sort.past             = 0.8;
  end
  sort.dtcron   = gene.dt;
  sor  			= nbi_sinbad_2tempssat(plas,fais,para,sort,equi,prof,impur,composition,geo,param,gene,memoire,cons);
  memoire.data.memoryspot = sor.memoryspot;
  memoire.data.input.option.init = 1;
  sor1			= sor.DEP1;
  sor2			= sor.DEP2;
  sor3			= sor.DEP3;
  sor4			= sor.DEP4;
  sor5			= sor.DEP5;
  sor6   		= sor.DEP6;
  pion1			= sor.pion1;
  pion2			= sor.pion2;
  pion3			= sor.pion3;
  pion4			= sor.pion4;
  pelec			= sor.pelec;
  ploss			= sor.ploss;
  Jidn			= sor.Jidn;
  RAD			= sor.RAD;
  wpar			= sor.Wpar;
  wtot			= sor.W;
  wth			= sor.Wth;
  ang			= sor.ang;
  dis			= sor.dis;
  dist			= sor.dist;
  xv			= sor.xv;
  ev			= sor.ev;
  neu1			= sor.neu1;
  neu2			= sor.neu2;
  neu3			= sor.neu3;
  neu4			= sor.neu4;
  neu5			= sor.neu5;
  neu6			= sor.neu6;
  mat			= sor.mat;
  s2d_1			= sor.s2d_1;
  s2d_2			= sor.s2d_2;
  s2d_3			= sor.s2d_3;
  s2d_4			= sor.s2d_4;
  s2d_5			= sor.s2d_5;
  s2d_6			= sor.s2d_6;
  s2d_7			= sor.s2d_7;
  s2d_8			= sor.s2d_8;
  s2d_9			= sor.s2d_9;
  s2d_10		= sor.s2d_10;
  s2t			= sor.s2t;
  divers		= sor.divers;

end

output.sor1            =  sor1;
output.sor2            =  sor2;
output.sor3            =  sor3;
output.sor4            =  sor4 ;
output.sor5            =  sor5 ;
output.sor6            =  sor6 ;
output.pion1           = pion1 ;
output.pion2           = pion2 ;
output.pion3           = pion3;
output.pion4           = pion4;
output.pelec           = pelec;
output.ploss           = ploss;
output.Jidn            = Jidn;
output.RAD             = RAD ;
output.wpar            = wpar ;
output.wtot            = wtot;
output.wth             = wth  ;
output.ang             = ang ;
output.dis             = dis ;
output.dist            = dist ;
output.xv              = xv ;
output.neutron.thth    = neu1 ;
output.neutron.beth    = neu2 ;
output.neutron.bebe    = neu3 ;
output.neutron.itot    = neu4 ;
output.neutron.abeb    = neu5 ;
output.neutron.atot    = neu6 ;
output.ev              = ev ;
output.divers          = divers;
output.mat             = mat;
output.s2d_1           = s2d_1;
output.s2d_2           = s2d_2;
output.s2d_3           = s2d_3;
output.s2d_4           = s2d_4;
output.s2d_5           = s2d_5;
output.s2d_6           = s2d_6;
output.s2d_7           = s2d_7;
output.s2d_8           = s2d_8;
output.s2d_9           = s2d_9;
output.s2d_10          = s2d_10;
output.s2t             = s2t;

memoire.data.output     = output;

%
% -- entrees --
% temps : temps d'analyse (t1:0.05:t1+0.3) en seconde
% -- donnees geometriques
% R0    : grand rayon en m (base temps)
% a     : petit rayon en m (base temps)
% Z     : hauteur du centre du plasma (en m)
% e0    : elongation cemtrale
% e1    : elongation au bord
% d0    : decentrement de Shafranov au centre (m)
% tri   : triangularite (en m)
% Rext  : grand rayon externe de la surface de flux (m) [rho,temps]
% cs    : 0 [defaut] -> sinus de l'angle d'injection, 1 -> cosinus
% -- donnees plasma --
% Ip        : courant plasma (MA)
% mass      : masse des especes ioniques (au maximum 4 especes)
% char      : charges associees
% exemples  : especes Deuterium, Helium, Carbone, Oxygene
%
% mass   = [2 4 12 18];
% charg  = [1 2 6 8];
% Bcron     : champ magnetique totale le long d'une surface de flux [x,theta]
% Boutcron  : champ magnetique externe de la surface de flux [x,theta]
% thetacron : angle poloidal associe a Rcron [x,theta]
% xcron     : rayon normalise des differentes surfaces de flux [1,x]
% Rcron     : position en R de la surface de flux [x,theta]
% fluxn : flux de neutrons (/s)
% -- profils --
% rho   : rayon normalise cronos
% Te    : temperature electronique en keV [rho,temps]
% Ti    : temperature ionique en keV
% ne    : densite electronique (1019 m-3)
% n1    : densite ionique premiere espece
% n2    : densite ionique seconde espece
% n3    : densite ionique troisieme espece
% n4    : densite ionique quatrieme espece
% Zeff  : profil de la charge effective
% ftrap  : profil de la fonction de trappée
% -- donnees injection de neutres --
% pin1  : puissance de l'injecteur numero 1 de JET (8 sources)
%         ou source de neutres TS (MW)(base temps), dim (source,temps)
% en1   : energie pleine injectee (8 pour JET, 1 pour TS)
% fr1   : fraction des energies pleines, 1/2 et 1/3 injectee en %
%         (8*3 pour JET, 1*3 pour TS)
% A1    : masse du gaz injecte
% Z1    : charge du gaz injecte
% align1: seulement pour JET, alignement standard ou upshifted des faisceaux
%
% type  : 0 -> NO REACTION
%         1 -> D+T  ==> 4HE+N
%         2 -> D+D  ==>  T+P
%         3 -> D+D  ==> 3HE+N
%         4 -> T+T  ==> 4HE+2N
%         5 -> D+3HE==> 4HE+P.
% ne    : 0 -> pas de flux de neutron, 1-> flux de neutrons
% -- sorties -- particules / m3
% sor1  : depots des ions suprathermiques injecteur 1 JET, pleine energie
% sor2  : depots des ions suprathermiques injecteur 1 JET, energie moitie
% sor3  : depots des ions suprathermiques injecteur 1 JET, energie tiers
% sor4  : depots des ions suprathermiques injecteur 2 JET, pleine energie
% sor5  : depots des ions suprathermiques injecteur 2 JET, energie moitie
% sor6  : depots des ions suprathermiques injecteur 2 JET, energie tiers
%
%
% -- sortie du Fokker Planck (W/m3 pour la chaleur et A/m2 pour le courant)
% pion1 : depot sur l'espece ionique 1
% pion2 : depot sur l'espece ionique 2
% pion3 : depot sur l'espece ionique 3
% pion4 : depot sur l'espece ionique 4
% pelec : depot sur les electrons
% ploss : se rajoute a la somme des pionX (x=1-4) pour le depot total sur les ions
% Jidn  : courant genere
% RAD   : rayon normalise de sortie
% Wpar  : densite d'energie parallele des suprathermiques
% W     : densite d'energie totale des suprathermiques
% Wth   : densite d'energie totale thermique
% nombre de pini
%
%
% V. Basiuk, 17 octobre 2001
% modifications
% V. Basiuk, 21 janvier 2002
% problemes resolus : calcul de sinbad meme si pas de puissance sur tous les temps
% 11 octobre 2006
% version preliminaire avec une bonne gestion de 2 boites d'IdN
%
% calcul des sources zineb
% on prend l'avant dernier temps
if ndims(pion1) == 2
  pion   =  pion1(:,end)' + pion2(:,end)' +pion3(:,end)' +pion4(:,end)' + ploss(:,end)';
  pel    =  pelec(:,end)';
  %
  % pour la puissance totale depose sur ions rapides
  %
  depr1  =  sor1(:,end)*mean(en1(en1>10)) + ...
	    sor2(:,end)*mean(en1(en1>10))/2 + ...
	    sor3(:,end)*mean(en1(en1>10))/3;
  if ~all(en2==0)
    depr2  =  sor4(:,end)*mean(en2(en2>10)) + ...
	      sor5(:,end)*mean(en2(en2>10))/2 + ...
	      sor6(:,end)*mean(en2(en2>10))/3;
  else
    depr2 = zeros(11,1);
  end
  depr1(~isfinite(depr1))=0;
  depr2(~isfinite(depr2))=0;
  deprt  =  depr1+depr2;
  js     =  Jidn(:,end)';
  nes    =  Z1(1) .* (sor1(:,end)'+sor2(:,end)'+sor3(:,end)') + ...
	    Z2(1) .* (sor4(:,end)'+sor5(:,end)'+sor6(:,end)');
  
  par1   =  (sor1(:,end)'+sor2(:,end)'+sor3(:,end)');
  par2   =  (sor4(:,end)'+sor5(:,end)'+sor6(:,end)');
  
  psupra =  wtot(:,end)' .* 1e6 - wpar(:,end)' .* 1e6;
  psupra(psupra<0) = 0;
  paniso =  wpar(:,end)' .* 1e6 - (1/2) .* psupra;
  
  
  % 03/02/2005 :  correction provisoire  de la puissance de bord pour simulation ITER
  % la renormalisation n'est pas utile
  %Vpr_loc   = interp1(rho,equi.vpr,RAD);
  %Pel_loc   = equi.rhomax .* trapz(RAD,Vpr_loc .* pel,2);
  pel(end) = 0;
  %Pion_loc  = equi.rhomax .* trapz(RAD,Vpr_loc .* pion,2);
  pion(end) = 0;

  
  % reechantillonage du  profil
  sortie.el           = zbornes(interp1(RAD,pel,rho,'pchip'),0,inf,0);
  %Pel_fin             = equi.rhomax .* trapz(rho,equi.vpr .* sortie.el,2);    % ajout du 03/02/2005
  %sortie.el           = sortie.el .* Pel_loc ./ Pel_fin;    % ajout du 03/02/2005  
  sortie.ion          = zbornes(interp1(RAD,pion,rho,'pchip'),0,inf,0);
  %Pion_fin            = equi.rhomax .* trapz(rho,equi.vpr .* sortie.ion,2);    % ajout du 03/02/2005
  %sortie.ion          = sortie.el .* Pion_loc ./ Pion_fin;    % ajout du 03/02/2005  
  sortie.ne           = zbornes(interp1(RAD,nes,rho,'pchip'),0,inf,0);
  sortie.j            = zbornes(interp1(RAD,js,rho,'pchip'),0,inf,0);
  sortie.psupra       = zbornes(interp1(RAD,psupra,rho,'pchip'),0,inf,0);
  sortie.paniso       = zbornes(interp1(RAD,paniso,rho,'pchip'),-inf,inf,0);
  deprtx              = zbornes(interp1(RAD,deprt,rho,'pchip'),-inf,inf,0);

%  keyboard
  % Replace SINBAD ouput by SPOT output
  if strcmp(param.spot_use,'Yes')
    rhospot       = sor.memoryspot.inprhoout / max(sor.memoryspot.inprhoout);
    sortie.el     = zbornes(interp1(rhospot,sor.memoryspot.elpower1d,rho),0,inf,0);
    sortie.ion    = zbornes(interp1(rhospot,sor.memoryspot.ionpower1d,rho),0,inf,0);
    sortie.j      = zbornes(interp1(rhospot,sor.memoryspot.icurpara1d(:,1),rho),0,inf,0);
    sortie.psupra = zbornes(interp1(rhospot,sor.memoryspot.pressperp1d,rho),0,inf,0);
    sortie.paniso = zbornes(interp1(rhospot,sor.memoryspot.presspara1d,rho),0,inf,0)-0.5*zbornes(interp1(rhospot,sor.memoryspot.pressperp1d,rho),0,inf,0);
  end
  
  
  %sortie.neutron.thth = zbornes(interp1(RAD,neu1(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.beth = zbornes(interp1(RAD,neu2(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.bebe = zbornes(interp1(RAD,neu3(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.itot = zbornes(interp1(RAD,neu4(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.abeb = zbornes(interp1(RAD,neu5(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.atot = zbornes(interp1(RAD,neu6(:,end),rho,'linear'),-inf,inf,0);
  % puissance de chauffage IDN
  paddidn            = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
  jidn               = zintsurf(sortie.j,gene.x,equi.spr,equi.rhomax);
  depidn             = zintvol(deprtx*phys.e,gene.x,equi.vpr,equi.rhomax)*1e3;
  pmaxidn            = sum(pcons(2,:),2);
  if strcmp(param.save,'Yes')
    [chemin,void] = fileparts(gene.rapsauve);
    save(fullfile(chemin,sprintf('last_sinbad@%s',int2str(fix(gene.t*1000)))));
  end
  rap               = paddidn/pmaxidn;
  sortie.err        = rap;
  disp([' J nbi = ',num2str(jidn/1e6,3),' MA'])
  disp([' absorbed power on fast ion = ',num2str(depidn/1e6,3),' MW'])
  disp([' total absorbed power = ',num2str(paddidn/1e6,3),' MW'])

  if paddidn < pmaxidn*0.9
    disp(['rap =',num2str(rap,3)])
  end
  if (rap > 1.1 | rap < 0.9) & strcmp(param.spot_use,'No')
    
    disp(['absorbed power different from the injected one, renormalisation, rap=',num2str(rap,3)])
    memoire.data.input.option.init = 0; 
    sortie.el         = sortie.el  / rap;
    sortie.ion        = sortie.ion / rap;
    sortie.psupra     = sortie.psupra/rap;
    sortie.j     = sortie.j/rap;
    [chemin,void] = fileparts(gene.rapsauve);
    if strcmp(param.save,'Yes')
      save(fullfile(chemin,sprintf('prob_sinbad@%s',int2str(fix(gene.t*1000)))));
    end
  end
  val                   = neu3(:,2) <= 0;
  neu3(val,2)           = 0;
  
  if any(neu3(:,2) >= 0)
    neutrons            = zbornes(interp1(RAD,neu2(:,2),rho,'linear'),-inf,inf,0) + ...
	zbornes(interp1(RAD,neu3(:,2),rho,'linear'),-inf,inf,0);
  else
    neutrons            = zbornes(interp1(RAD,neu2(:,2),rho,'linear'),-inf,inf,0);
  end
  %
  % passage en m-3
  %
  sortie.neutron.dd    = neutrons*1e6;
  
  % source de matiere par espece
  par1   = zbornes(interp1(RAD,par1,rho,'linear'),0,inf,0);
  par2   = zbornes(interp1(RAD,par2,rho,'linear'),0,inf,0);
  
  % recherche de lindice de l'espece
  ind1       = find((composition.z == Z1(1)) &(composition.a == A1(1)));
  ind2       = find((composition.z == Z2(1)) &(composition.a == A2(1)));
  particules = zeros(1,gene.nbrho,gene.nbg);
  
  if ~isempty(ind1)
    particules(1,:,ind1)   = par1(:);
  end
  if ~isempty(ind2)
    particules(1,:,ind2)   = par2(:);
  end
else
  % recherche de lindice de l'espece
  if any(diff(A1)) == 0
    indesp1   = 1:8;
    indesp2   = [];
  else
    indesp1   = find(A1 == A1(1));
    indesp2   = find(A1 ~= A1(1));
  end
  if any(diff(A2)) == 0
    indesp3   = 1:8;
    indesp4   = [];
  else
    indesp3   = find(A2 == A1(1));
    indesp4   = find(A2 ~= A1(1));
  end
  pion = zeros(1,11);
  pel = zeros(1,11);
  js = zeros(1,11);
  nes = zeros(1,11);
  psupra = zeros(1,11);
  paniso = zeros(1,11);
  for kesp=1:2
    eval(['indv1=indesp',int2str(kesp),';']);
    eval(['indv2=indesp',int2str(kesp+2),';']);
    pion   =  pion + pion1(:,end,kesp)' + pion2(:,end,kesp)' +pion3(:,end,kesp)' +pion4(:,end,kesp)' + ploss(:,end,kesp)';
    pel    =  pel  + pelec(:,end,kesp)';
    js     =  js   + Jidn(:,end,kesp)';
    if ~isempty(indv1)
      nes    =  nes + Z1(indv1(1)) .* (sor1(:,end,kesp)'+sor2(:,end,kesp)'+sor3(:,end,kesp)');
    end
    if ~isempty(indv2)
      nes   = nes + Z2(indv2(1)) .* (sor4(:,end,kesp)'+sor5(:,end,kesp)'+sor6(:,end,kesp)');
    end	  
    par1(kesp,:)   =  (sor1(:,end,kesp)'+sor2(:,end,kesp)'+sor3(:,end,kesp)');
    par2(kesp,:)   =  (sor4(:,end,kesp)'+sor5(:,end,kesp)'+sor6(:,end,kesp)');
    
    psupra =  psupra + wtot(:,end,kesp)' .* 1e6 - wpar(:,end,kesp)' .* 1e6;
    paniso =  paniso + wpar(:,end,kesp)' .* 1e6 - (1/2) .* psupra;
  end
  
  % 03/02/2005 :  correction provisoire  de la puissance de bord pour simulation ITER
  pel(end) = 0;
  pion(end) = 0;
  
  % reechantillonage du  profil
  sortie.el           = zbornes(interp1(RAD,pel,rho,'pchip'),0,inf,0);
  sortie.ion          = zbornes(interp1(RAD,pion,rho,'pchip'),0,inf,0);
  sortie.ne           = zbornes(interp1(RAD,nes,rho,'pchip'),0,inf,0);
  sortie.j            = zbornes(interp1(RAD,js,rho,'pchip'),0,inf,0);
  sortie.psupra       = zbornes(interp1(RAD,psupra,rho,'pchip'),0,inf,0);
  sortie.paniso       = zbornes(interp1(RAD,paniso,rho,'pchip'),-inf,inf,0);
  %sortie.neutron.thth = zbornes(interp1(RAD,neu1(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.beth = zbornes(interp1(RAD,neu2(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.bebe = zbornes(interp1(RAD,neu3(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.itot = zbornes(interp1(RAD,neu4(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.abeb = zbornes(interp1(RAD,neu5(:,end),rho,'linear'),-inf,inf,0);
  %sortie.neutron.atot = zbornes(interp1(RAD,neu6(:,end),rho,'linear'),-inf,inf,0);
  % puissance de chauffage IDN

  % Replace SINBAD ouput by SPOT output
%  keyboard
  if strcmp(param.spot_use,'Yes')
    rhospot       = sor.memoryspot.inprhoout / max(sor.memoryspot.inprhoout);
    sortie.el     = zbornes(interp1(rhospot,sor.memoryspot.elpower1d,rho),0,inf,0);
    sortie.ion    = zbornes(interp1(rhospot,sor.memoryspot.ionpower1d,rho),0,inf,0);
    sortie.j      = zbornes(interp1(rhospot,sor.memoryspot.icurpara1d(:,1),rho),0,inf,0);
    sortie.psupra = zbornes(interp1(rhospot,sor.memoryspot.pressperp1d,rho),0,inf,0);
    sortie.paniso = zbornes(interp1(rhospot,sor.memoryspot.presspara1d,rho),0,inf,0)-0.5*zbornes(interp1(rhospot,sor.memoryspot.pressperp1d,rho),0,inf,0);
  end
  
  paddidn             = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
  pmaxidn            = sum(sum(pcons,2))/2;
  if strcmp(param.save,'Yes')
    [chemin,void] = fileparts(gene.rapsauve);
    save(fullfile(chemin,sprintf('last_sinbad@%s',int2str(fix(gene.t*1000)))));
  end
  rap               = paddidn/pmaxidn;
  sortie.err        = rap;
  
  if paddidn > pmaxidn & strcmp(param.spot_use,'No')
    
    disp(['absorbed power greater than the injected one, renormalisation, rap=',num2str(rap,3)])
    
    sortie.el         = sortie.el  / rap;
    sortie.ion        = sortie.ion / rap;
    [chemin,void] = fileparts(gene.rapsauve);
    if rap > 100
      %      save probsinbad
    end
  end
  val                   = neu3(:,2) <= 0;
  neu3(val,2)           = 0;
  
  if any(neu3(:,2) >= 0)
    neutrons            = zbornes(interp1(RAD,neu2(:,2),rho,'linear'),-inf,inf,0) + ...
	zbornes(interp1(RAD,neu3(:,2),rho,'linear'),-inf,inf,0);
  else
    neutrons            = zbornes(interp1(RAD,neu2(:,2),rho,'linear'),-inf,inf,0);
  end
  %
  % passage en m-3
  %
  sortie.neutron.dd    = neutrons*1e6;
  
  % source de matiere par espece
  for kesp=1:2
    npar1(kesp,:)   = zbornes(interp1(RAD,par1(kesp,:),rho,'linear'),0,inf,0);
    npar2(kesp,:)   = zbornes(interp1(RAD,par2(kesp,:),rho,'linear'),0,inf,0);
  end
  particules = zeros(1,gene.nbrho,gene.nbg);
  
  for kesp = 1:2
    eval(['indv1=indesp',int2str(kesp),';']);
    eval(['indv2=indesp',int2str(kesp+2),';']);
    if ~isempty(indv1)
      ind1       = find((composition.z == Z1(indv1(1))) &(composition.a == A1(indv1(1))));
    else
      ind1 = [];
    end
    if ~isempty(indv2)
      ind2       = find((composition.z == Z2(indv2(1))) &(composition.a == A2(indv2(1))));
    else
      ind2 = [];
    end
    if ~isempty(ind1)
      particules(1,:,ind1)   = particules(1,:,ind1)+npar1(kesp,:);
    end
    if ~isempty(ind2)
      particules(1,:,ind2)   = particules(1,:,ind2)+npar2(kesp,:);
    end
  end
  
end

memoire.data.sortie = sortie;
memoire.data.particules = particules;
memoire.data.consin = cons;

% integralle volumique
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax
function s=zintvol(e,x,vpr,rhomax)

  s = rhomax.*trapz(x,vpr .* e,2);


% the end

% integralle surfacique (d'une grandeur independante de theta)
%  s = integrale de surface
%  e = valeur a integree
%  x = coordonnees normalisee
%  sp = datak.equi.sp
%  rhomax = datak.equi.rhomax
function s=zintsurf(e,x,sp,rhomax)

    s = rhomax .* trapz(x,sp .* e,2);

