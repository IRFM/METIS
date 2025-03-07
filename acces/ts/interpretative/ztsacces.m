function [cr,data,param]=ztsacces(numchoc,chemin,temps,option)
% [cr,data,param]=ztsacces(numchoc,chemin,temps,option);
% acces aux donnees TS, pr???aration du fichier CRONOS
%
%
% Auteur : J.F. Artaud (46 78)
%
% dernieres modifications :
% V. Basiuk, 12 ferier 2002
% nouveau calcul du sc???ario avec choix de lancement de pion ou absor
% 19 septembre 2003
% blindage si pas itor
% 13 fevrier 2004
% correction scenariop fci si hydrogene pas premier minoritaire
% * 21/07/2004 -> ajout de la consigne de spectre continu de l'hybride
% * 21/01/2005 -> correction de la deuxieme dimension de hybspec si pas de signal GHYBSPEC (2 au lieu de 1, comme dans zinit)

if nargin <=1 
	valeur.nbrho       = 101;     % nombre de points radiaux [101]
	valeur.tdebut      = [];      % temps de debut du calcul , si vide min(temps)
	valeur.tfin        = [];      % temps de fin du calcul, si vide max(temps)
	valeur.fjli        = 1;       % profil initial de Jmoy (ancien choc date <= 1999) : 0 -> coupole , 1 -> ip et li (dpolo)
	valeur.tiop        = 0;       % valeur de Ti : 0 -> tiprof, 1 -> a*ne, 2 -> a*Te, 3 -> =Te, 4 et repli -> scaling [4] 
	valeur.facte       = 1;       % facteur multiplicatif de Te [1]
	valeur.faczeff     = 1;       % facteur multiplicatif de Zeffm [1]
	valeur.psimode     = 2;       % mode Psi : 1 -> interpretatif, 2 -> predictif
	valeur.pemode      = 1;       % mode Pe : 1 -> interpretatif, 2 -> predictif
	valeur.pionmode    = 1;       % mode Pion : 1 -> interpretatif, 2 -> predictif
	valeur.nelmode     = 1;       % mode Ne : 1 -> interpretatif, 2 -> predictif
	valeur.psilim      = 0;       % condition au limite sur Psi : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)
	valeur.lambda      = 0;       % facteur lambda dans les flux de chaleur provenant du flux de particules {0,3/.2,5/2}
	valeur.modecoef    = 1;       % mode de calcul des coefficients des equations de transport  0 -> tous, 1 -> convectif +diagonaux
	%valeur.self        = 0;       % mode de fonctionnement completement  auto consistant si = 1, 0 -> coefficient +neoclassique auto consistant 
	valeur.fast        = 0;       % 0 -> mode standart correct, 1-> pas de sources neoclassique, ni de coefficient neo self consistante, 2 -> mode optimiser pour la diffusion du courant
	valeur.source_bord = 0;       % controle du clacul de recyclage  : 1 = recyclage recalculer a chaque sous pas de temps, 0 = comme les sources
	valeur.cn          = 0.3;     % mode solveur : 0 -> implicite, 0.5 -> Cranck-Nickolson
	valeur.plotonoff   = 0;       % mode interractif avec plot si 1 
	valeur.verbose     = 0;       % 1 -> commentaires , 0 -> pas de commentaires
	valeur.rebuilt     = 1;       % 1 -> reconstruction automatique du fichier resultat en fin d'execution 
	valeur.post        = 1;       % 1 -> postprocessing en fin d'execution
	
	type.nbrho       = 'integer';    
	type.tdebut      = 'float';      
	type.tfin        = 'float';    
	type.fjli        = 'integer';       
	type.tiop        = 'integer';       
	type.facte       = 'float';       
	type.faczeff     = 'float';       
	type.psimode     = 'integer';       
	type.pemode      = 'integer';       
	type.pionmode    = 'integer';       
	type.nelmode     = 'integer';       
	type.psilim      = 'integer';       
	type.lambda      = 'integer';     
	type.modecoef    = 'integer';      
	type.self        = 'integer';      
	type.fast        = 'integer';       
	type.source_bord = 'integer'; 
	type.cn          = 'float';
	type.plotonoff   = 'integer';       
	type.verbose     = 'integer'; 
	type.rebuilt     = 'integer';       
	type.post        = 'integer'; 
   
	borne.nbrho       = {21,51,101,201,501,1001};    
	borne.tdebut      = [0,inf];      
	borne.tfin        = [0,inf];    
	borne.fjli        = {0,1};       
	borne.tiop        = {0,1,2,3,4,5};       
	borne.facte       = [0.1,10];       
	borne.faczeff     = [0.1,10];       
	borne.psimode     = {1,2};       
	borne.pemode      = {1,2};       
	borne.pionmode    = {1,2};      
	borne.nelmode     = {1,2};       
	borne.psilim      = {0,1,2};      
	borne.lambda      = {0,3/2,5/2};     
	borne.modecoef    = {0,1};      
	borne.self        = {0,1};       
	borne.fast        = {0,1,2};       
	borne.source_bord = {0,1};       
	borne.cn          = {0,0.3,1/2};
	borne.plotonoff   = {0,1};      
	borne.verbose     = {0,1};       
 	borne.rebuilt     = {0,1};      
	borne.post        = {0,1};       
  
	defaut.nbrho       = 101;    
	defaut.tdebut      = [];      
	defaut.tfin        = [];    
	defaut.fjli        = 1;       
	defaut.tiop        = 0;       
	defaut.facte       = 1;       
	defaut.faczeff     = 1;       
	defaut.psimode     = 2;       
	defaut.pemode      = 1;       
	defaut.pionmode    = 1;       
	defaut.nelmode     = 1;       
	defaut.psilim      = 0;       
	defaut.lambda      = 0;     
	defaut.modecoef    = 1;      
	defaut.self        = 0;      
	defaut.fast        = 0;       
	defaut.source_bord = 0;       
	defaut.cn          = 0.3;
 	defaut.plotonoff   = 0;       
	defaut.verbose     = 0;       
 	defaut.rebuilt     = 1;       
	defaut.post        = 1;       
  
	info.nbrho       = 'radial points [101]';
	info.tdebut      = 'initial time (if empty, then min(temps) is used';
	info.tfin        = 'final time (if empty, then max(temps) is used';
	info.fjli        = 'first current profile (old shot <= 1999, otherwise used TPROF) : 0 -> coupole, 1 -> TPROF';
	info.tiop        = 'ion temperature calculation : 0 -> tiprof, 1 -> a*ne, 2 -> a*Te, 3 -> =Te, 4 and default -> scaling [4],5 -> Tiprof and neutron ';
	info.facte       = 'multiplier of the electron temperature [1]';
	info.faczeff     = 'fmultiplier of the Zeff [1]';
	info.psimode     = 'Psi mode: 1 -> interpretatif, 2 -> predictif';
	info.pemode      = 'Pe mode : 1 -> interpretatif, 2 -> predictif';
	info.pionmode    = 'Pion mode: 1 -> interpretatif, 2 -> predictif';
	info.nelmode     = 'Ne mode: 1 -> interpretatif, 2 -> predictif';
	info.psilim      = 'boundary on Psi : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)';
	info.lambda      = 'lambda factor in heat flux from density {0,3/.2,5/2}';
	info.modecoef    = 'matrix transport coefficient  0 -> all, 1 -> convectif + diagonal';
	info.self        = 'auto coherent if  = 1, 0 -> only neoclassic auto coherent ';
	info.fast        = '0 -> correct mode, 1-> no neoclassical source, neither self consistant neoclassical coefficient, 2 -> fast mode for current diffusion only';
	info.source_bord = 'recycling mode  : 1 = recycling in the internal time step, 0 = as sources';
	info.cn          = 'solver mode : 0 -> implicite, 0.5 -> Cranck-Nickolson 0.3 -> recommended';
	info.plotonoff   = 'show plot if = 1 ';
	info.verbose     = '1-> comments , 0 -> no comments';
	info.rebuilt     = '1 -> reconstruction of the output file from intermediate file'; 
	info.post        = '1 -> postprocessing at the end (reconstruction of some diagnostics)';

	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Access Module to Tore Supra data';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end

% initialisation variable vide
contenu     = [];
remplissage = [];
consigne    = [];

% compatibilite
nbrho     = option.nbrho;
tdebut    = option.tdebut;
tfin      = option.tfin;
fjli      = option.fjli;
tiop      = option.tiop;
facte     = option.facte;
faczeff   = option.faczeff;


% cr par defaut
cr =0;
% choc charniere ciel
charniere = 28353;

% le bon profil de ni
niprofok =[];

% variables contenant les information sur les donnees
infotsdata =[];
desctsdata ={};
% gestion des entrees
if nargin <4
	cr = 1;
	disp('Nombre d''arguments incorrect !');
end

% lecture des donnees
%
% probleme dans cgcgettrait pour tcoupol (traitement inexistant)
% on evite le probleme en mettant une charniere plus petite
% V. Basiuk, 30 mars 2004
%
charniere2 = 28000;
if fix(numchoc) <= charniere2 
	[bile,sbile]   = cgcgettrait(numchoc,'tprof');
	infotsdata = zinfocat(infotsdata,sbile,'tprof');
   desctsdata{end+1} = 'trpof';
   [coup,scoup]   = cgcgettrait(numchoc,'tcoupol');
   desctsdata{end+1} = 'tcoupol';
	infotsdata = zinfocat(infotsdata,scoup,'tcoupol');
	%[boot,sboot]   = cgcgettrait(numchoc,'tbootstrap');
elseif 	fix(numchoc) > charniere2 
	[bile,sbile]   = cgcgettrait(numchoc,'tprof');
   desctsdata{end+1} = 'tprof';
	infotsdata = zinfocat(infotsdata,sbile,'tprof');
	coup.fluxp =[];
end
if isempty(bile.times)
	[bile,sbile]   = zcgcgetprof(numchoc);
	infotsdata = zinfocat(infotsdata,sbile,'prof');
   desctsdata{end+1} = 'prof';
	if isempty(bile)
		disp('Pas de donnees pour ce choc ...');
		return
		cr = -3;
	end
   if ~isfield(bile,'jmoy') & isfield(bile,'j')
     bile.jmoy = bile.j;
   elseif ~isfield(bile,'jmoy')
     bile.jmoy = [];
   end	
end

% pour les plot
zassignin('base','bile',bile);

% lecture de tiprof	
file  = sprintf('/usr/drfc/cgc/matlab5/tcron/TS/%d/Ti%d.mat',fix(numchoc),fix(numchoc));

if exist(file,'file')
  tiprof=load(file);
  stiprof.tiprof.commentaire = file;
  desctsdata{end+1} = 'tiprof';
  infotsdata = zinfocat(infotsdata,stiprof,'tiprof');
else
  disp('fichier de donnees tiprof indisponible')
  tiprof=[];
end


% si presence hybride
phyb  = [];
tphyb =[];
if any(bile.plh >0.5)
  if fix(numchoc)<20000   
      [phyb,tphyb,void,cert] = tsbase(fix(numchoc),'ghyb');
	   infotsdata = zinfocat(infotsdata,'ghyb','thyb',cert);
      desctsdata{end+1} = 'thyb';
  else
      [phyb,tphyb,void,cert] = tsbase(fix(numchoc),'gphyb');
	   infotsdata = zinfocat(infotsdata,'gphyb','dhyb',cert);
      desctsdata{end+1} = 'dhyb';
      if size(tphyb,2) == 1
         tphyb = tphyb * ones(1,size(phyb,2));
      end
  end
  if (fix(numchoc)<20000) & isempty(phyb) 
      [phyb,tphyb,void,cert] = tsbase(fix(numchoc),'gphyb');
	   infotsdata = zinfocat(infotsdata,'gphyb','dhyb',cert);
      desctsdata{end+1} = 'dhyb';
      if size(tphyb,2) == 1
         tphyb = tphyb * ones(1,size(phyb,2));
      end
  end 
end
 
% pour tous les chocs
%[pecrh1,tecrh1]= tsbase(fix(numchoc),'sicata1');
%[pecrh2,tecrh2]= tsbase(fix(numchoc),'sicata2');
%if ~isempty(pecrh1);
%	pecrh1(pecrh1<5) = 0;
%   pecrh1 = pecrh1 - min(pecrh1);
%   pecrh1 = pecrh1 ./ max(pecrh1+eps) .* 0.300;
%end 
%if ~isempty(pecrh2);
%	pecrh2(pecrh2<5) = 0;
%   pecrh2 = pecrh2 - min(pecrh2);
%   pecrh2 = pecrh2 ./ max(pecrh2+eps) .* 0.500;
%end 
 
% valeurs par default
if isempty(temps)
	temps=bile.times;
elseif size(temps,2)>1
	temps=temps';
end
if isempty(tdebut)
	tdebut =min(temps);
end
if isempty(tfin)
	tfin = max(temps);
end
if isempty(nbrho)
	nbrho=length(bile.rhofit);
end

charniere = 28353;

if fix(numchoc) <= charniere

    [gfa,tgfa,ygfa,cgfa]             = tsbase(fix(numchoc),'gfa');
	  infotsdata = zinfocat(infotsdata,'gfa','tempete',cgfa);
    [gfteta,tgfteta,ygfteta,cgfteta] = tsbase(fix(numchoc),'gfteta'); 
	  infotsdata = zinfocat(infotsdata,'gfteta','tempete',cgfteta);
     desctsdata{end+1} = 'tempete';
    [z0,times]                       = tsbase(numchoc,'sprofzpos');
    [r0,times]                       = tsbase(numchoc,'sprofrmaj');
    [a0,times]                       = tsbase(numchoc,'sprofamin');
    zt                               = interp1(times,z0,temps)'; 
    Rt                               = interp1(times,r0,temps)'; 
    at                               = interp1(times,a0,temps)'; 
    gfa                              = tsample(gfa,tgfa,temps,'fen');
    gfteta                           = tsample(gfteta,tgfa,temps,'fen');  
    teta  = linspace(0,2*pi,101);
    comp  = ones(1,size(teta,2));
    teta  = ones(length(temps),1)*teta;
    rayon = (at'*comp)+(gfa(:,1)*comp).*cos(teta-(gfteta(:,1)*comp))+...
	  (gfa(:,2)*comp).*cos(2*(teta-(gfteta(:,2)*comp)))+...
	  (gfa(:,3)*comp).*cos(3*(teta-(gfteta(:,3)*comp)))+...
	  (gfa(:,4)*comp).*cos(4*(teta-(gfteta(:,4)*comp)));
    %
    % derniere surface
    %
    Rext  = rayon.*cos(teta)+(Rt'*comp);
    Zext  = rayon.*sin(teta)+(zt'*comp);
	 %
	 % le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
	 % lorsque le tokamak est regarde depuis le haut
	 %
	 signe.ip  = 1;              % signe du courant plasma
	 signe.b0  = 1;              % signe du champ toroidal
else

    [grho,tgrho,void,cert]   = tsbase(fix(numchoc),'grho');
	 infotsdata = zinfocat(infotsdata,'grho','dpolo',cert);
    desctsdata{end+1} = 'dpolo';
    tgrho          = tgrho(:,1);
    % attente d'acces base
    zaxe           = 0;
    raxe           = 2.42;
    
    alpha          = (0:15:345) ./ 180 .* pi;
    vt             = ones(size(grho,1),1);
    rr             = raxe + grho .* cos(vt*alpha);
    zz             = zaxe + grho .* sin(vt*alpha);
    rr(:,end+1)    = rr(:,1);
    zz(:,end+1)    = zz(:,1);
    alpha(end+1)   = 2*pi;
    %
    % derniere surface
    %
    teta           = linspace(0,2*pi,201);
    Rext           = tsplinet(vt*alpha,rr,vt*teta);
    Zext           = tsplinet(vt*alpha,zz,vt*teta);
    Rext(:,end)    = Rext(:,1);
    Zext(:,end)    = Zext(:,1);
    
    % recalcul des parametres pour verification
    ve    = ones(1,size(Rext,2));
    rmin  = min(Rext,[],2);
    rmax  = max(Rext,[],2);
    ra    = 0.5 .* (rmin + rmax);   
    a     = 0.5 .* (rmax - rmin);
    zmin  = min(Zext,[],2);
    zmax  = max(Zext,[],2);
    za    = (zmin + zmax) ./ 2;
    b     = 0.5 .* (zmax - zmin);
    k     = b ./ a;
    mask1 = (Zext == (max(Zext,[],2)*ve));
    mask2 = (Zext == (min(Zext,[],2)*ve));
   
    rzmax = max(Rext .* mask1,[],2);
    rzmin = max(Rext .* mask2,[],2);
    cl    = ra - rzmin;
    cu    = ra - rzmax;
    d     = (cl+cu) ./2 ./ a;
    
    % calcul des parametres pour helena
    hr0     = ra;
    hz0     = za;
    ha      = a;
    he1     = k;
    htrl    = - asin(cl ./ a); 
    htrh    =   asin(cu ./ a); 
	 %
	 % le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
	 % lorsque le tokamak est regarde depuis le haut
	 %
	 signe.ip  = -1;              % signe du courant plasma
	 signe.b0  = -1;              % signe du champ toroidal
    
end

% calcul de R*B0                              

[tori,titor,cert] = tsbase(fix(numchoc),'stori');
[torvar,titor,cert] = tsbase(fix(numchoc),'storvar');
if ~isempty(tori) &  ~isempty(torvar)
	itor = mean(tori((titor<-10) & (titor>= -32))) + torvar - mean(torvar((titor<-30) & (titor>=-32)));
else
	[itor,titor,cert] = tsbase(fix(numchoc),'sitor');
end
infotsdata = zinfocat(infotsdata,'itor','generic',cert);
desctsdata{end+1} = 'generic';
if isempty(itor)
	[itor,titor,cert] = tsbase(fix(numchoc),'stori');
   infotsdata = zinfocat(infotsdata,'itor','dmag',cert);
   desctsdata{end+1} = 'dmag';
end
if isempty(itor)
   itor=tsmat(fix(numchoc),'EXP=T=S;GENERAL;ITOR');
   itor = itor*ones(size(temps));
   titor = temps;
   infotsdata = zinfocat(infotsdata,'itor','EXP-T-S;GENERAL;ITOR','');
   desctsdata{end+1} = 'EXP=T=S';
end
	
ind          = find(titor >=0);
itor         = itor(ind);
titor        = titor(ind);
[b,a]        = butter(11,0.1);
itor         = filtfilt(b,a,itor);
rb0          = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi;
trb0         = titor;

% consigne de courant plasma
[sipmes,tsipmes,cert] = tsbase(fix(numchoc),'sipmes');
infotsdata = zinfocat(infotsdata,'sipmes','dpolo',cert);
desctsdata{end+1} = 'dpolo';

% valeurs pour TS
nbfci    = 3;
nbfce    = 3;
nbhyb    = 2;
nbidn    = 1;
nbglacon = 3;
nbg      = 5;

% parametre d'echantillnage (cf. zsample.m)
%    1 - pour un signal simple :
signal.ondelette        = 0;  %1
signal.defaut.temps     = NaN;
signal.defaut.espace    = 0;
signal.defaut.inf       = [];
signal.plus             = 0;

%    2 - pour un groupe de signaux :
groupe.ondelette         = 0;
groupe.energie           = 1;  %0.01;
groupe.defaut.temps     = NaN;
groupe.defaut.espace    = 0;
groupe.defaut.inf       = [];
groupe.plus             = 0;


% pour tiprof 
mfl =0;

% creation du nom du fichier
if isempty(chemin)
	chemin=strcat(getenv('HOME'),'/zineb/data');
end
fichier = strcat(chemin,'/zineb',int2str(fix(numchoc)),'ts',int2str(round(rem(numchoc,1)*10)), ...
                 'de',int2str(round(tdebut*1000)),'a',int2str(round(tfin*1000)));

% initialisation de la structure
[cr,data,param]=zinit('',temps,nbrho,fichier,tdebut,tfin,nbfci,nbfce,nbidn,nbhyb,nbglacon,nbg);
if cr ~=0
	return
end
vt = ones(param.gene.nbt,1);
ve = ones(1,param.gene.nbrho);

% parametres generaux
param.gene.modecoef     = option.modecoef;   % modele diffussif/convectif - pas de terme non diagonaux
%param.gene.self         = option.self;
param.gene.lambda       = option.lambda;     
param.gene.fast         = option.fast;       
param.gene.source_bord  = option.source_bord;       
param.gene.cn           = option.cn;       
param.plot.onoff        = option.plotonoff;       
param.gene.verbose      = option.verbose;     
param.gene.rebuilt      = option.rebuilt;     
param.gene.post         = option.post;     
if any([option.nelmode,option.pemode,option.pionmode] > 1)
	param.gene.nbeq_mode    = 4;        
else
	param.gene.nbeq_mode    = 1;       
	
end
param.gene.signe = signe;

% connexion des modules externes
[cr,data,param] = zconnexion(data,param);
if cr ~=0
	return
end


% changement du modele de transport par defaut
param.fonction.coefa = 'zbgbs_TS';
infocoefa = zbgbs_TS;
param.cons.coefa = infocoefa.valeur;

% 1- preparation du changement de rho
if (fix(numchoc) <= charniere)
	vtin = ones(size(bile.times));
	vein = ones(size(bile.rhofit)); 
	rhogin =  vtin * bile.rhofit;
	shiftin = (bile.d0 * vein) .* ( 1 - rhogin .^ (bile.piqd * vein));
	delta = ((bile.rmaj * vein) + shiftin) .^ 2 - (bile.amin * vein).^2 .* rhogin.^2;
	phiin = 2 .* pi .* (bile.rmaj * vein) .* (bile.btor * vein) .* ...
	((bile.rmaj * vein) + shiftin - sqrt(delta));
	
	% coordonnee de flux dans la base temps d'origine        
	rhoin = sqrt(phiin ./ pi ./ (bile.btor * vein));   % coordonnee sqrt(Phi/pi/B0) pour les donnees bile
	rhomaxin = max(rhoin')';   % rhomax dans bile
	
	% coordonnee de flux dans la base temps finale
	rhomax = zsample(rhomaxin,bile.times,data.gene.temps,signal);  % rhomax pour zineb
	rho = rhomax * param.gene.x;  % rho sqrt(Phi/pi/b0) pour zineb
	
	% appel de l'equilibre de TS
	[rhophi_eq,rhog_eq,d,vpr,grho2,r2i,ri,grho2r2,psi_eq,phi,fdia,ptoteq,jmoy,q] =  ...
	equilTS(bile.rmaj,bile.amin,bile.ip,bile.beli,bile.wdia,-bile.btor);
	vteq = ones(size(rhog_eq,1),1);
	veeq = ones(1,size(rhog_eq,2));
	xx = rhophi_eq ./ (max(rhophi_eq')'*veeq);   % coordonnees  normalisee pour l'equilibre (correspond a rhofit pour bile)
	
	% coordonnee de reechantillonage (pour eviter les probleme d'arrondis)
	xz  = vt * param.gene.x;
	xb  = rhoin ./ (max(rhoin')'*vein);           % coordonnee de flux normalisee de bile
	rsa    = zsample(rhogin,bile.times,xb,temps,xz,groupe); % equivalent de rhofit pour cronos (r/a)
	
else
	vtin = ones(size(bile.times));
	vein = ones(size(bile.rhofit));
	rhogin =  vtin * bile.rhofit;
	% coordonnee de flux dans la base temps d'origine        
	rhoin = sqrt(bile.phi ./ pi ./ (bile.btor * vein));   % coordonnee sqrt(Phi/pi/B0) pour les donnees bile
	rhomaxin = max(rhoin')';   % rhomax dans bile
	% coordonnee de flux dans la base temps finale
	rhomax = zsample(rhomaxin,bile.times,data.gene.temps,signal);  % rhomax pour zineb
	rho = rhomax * param.gene.x;  % rho sqrt(Phi/pi/b0) pour zineb
	% coordonnee de reechantillonage (pour eviter les probleme d'arrondis)
	xz  = vt * param.gene.x;
	xb  = rhoin ./ (max(rhoin')'*vein);           % coordonnee de flux normalisee de bile
	% utilisation de l'equilibre tprof
	rsa  = zsample(rhogin,bile.times,xb,temps,xz,groupe); % coordonnees geometrique normalisee pour cronos (correspond a rhofit pour bile)
end

% ajout dans data.equi pour utilisation ulterieure :
% attention ceci ne doit pas etre utilise comme un sortie de  l'equilibre
data.equi.a = - rsa;

% profil de xdur
[xdur,sxdur]   = zgetxdur_base(fix(numchoc),temps,rsa,bile.ip);
infotsdata = zinfocat(infotsdata,sxdur,'thxrprof');
desctsdata{end+1} = 'thxrprof';

% 2 - remplissage des parametres
% les informations sur les donnees
param.from.machine = 'TS';
param.from.shot.num = numchoc;
param.from.shot.date = eval(bile.info.bile.date);
param.from.shot.info = bile.info;
param.from.creation.date =clock;
[s,whoami] = unix('whoami');
param.from.creation.user = whoami(whoami > ' ');

signal.plus =0;
groupe.plus=0;
param.from.sample.signal   = signal;      % parametre d'echantillonage  pour un signal simple =s(temps,1)
param.from.sample.groupe   = groupe;      % parametre d'echantillonage  pour un groupe =g(temps,espace)
param.from.createur = 'ztsacces';
param.from.option  = option;
% paroi pour ciel
if fix(numchoc) < 28590
	[zp,rp,void,cert] = tsbase(fix(numchoc),'sparoi');
   infotsdata = zinfocat(infotsdata,'sparoi','tempete',cert);
   desctsdata{end+1} = 'tempete';
else
	x = tsmat(fix(numchoc),'APOLO;+Parametres;Paroi');
   sdata.paroi.commentaire  = 'APOLO;+Parametres;Paroi';
   infotsdata = zinfocat(infotsdata,sdata,'apolo');
   desctsdata{end+1} = 'apolo';

	rp = x(:,1);
	zp = x(:,2);
end
param.from.paroi.R         = rp;          % vecteur R decrivant la parois de la machine 
param.from.paroi.Z         = zp;          % vecteur Z decrivant la parois de la machine 


% la composition du plasma

% on commence par mettre des 0
param.compo.z = zeros(size(param.compo.z));
param.compo.a = zeros(size(param.compo.a));

param.compo.z(1)=bile.zmain(1);
param.compo.a(1)=bile.mmain(1);
if length(bile.zmain)== 2
	param.compo.z(2)=bile.zmain(2);
	param.compo.a(2)=bile.mmain(2);
elseif length(bile.zmain)==3
	param.compo.z(2)=bile.zmain(2);
	param.compo.a(2)=bile.mmain(2);
	param.compo.z(3)=bile.zmain(3);
	param.compo.a(3)=bile.mmain(3);
end
while (any(param.compo.z == 0))
      ind = min(find(param.compo.z == 0));
      
      % 1 - minoritaire H
      if ~ any((param.compo.z == 1) & (param.compo.a == 1))
        param.compo.z(ind)=1;
	param.compo.a(ind)=1;
      % 2 - minoritaire D
      elseif ~any((param.compo.z == 1) & (param.compo.a == 2))
        param.compo.z(ind)=1;
	param.compo.a(ind)=2;
      % 3 - minoritaire He
      elseif ~any((param.compo.z == 2) & (param.compo.a == 4))
        param.compo.z(ind)=2;
	param.compo.a(ind)=4;
      % 3 - minoritaire He3
      elseif ~any((param.compo.z == 2) & (param.compo.a == 3))
        param.compo.z(ind)=2;
	param.compo.a(ind)=3;
      % 5 - repli T
      elseif ~any((param.compo.z == 1) & (param.compo.a == 3))
        param.compo.z(ind)=1;
	param.compo.a(ind)=3;
      % 6 - repli Li
      elseif ~any((param.compo.z == 3) & (param.compo.a == 6))
        param.compo.z(ind)=3;
	param.compo.a(ind)=6;
      end
end

param.compo.z(4)=bile.zimp(1);
param.compo.a(4)=bile.mimp(1);
if length(bile.zimp)>1
	param.compo.z(5)=bile.zimp(2);
	param.compo.a(5)=bile.mimp(2);
elseif bile.zimp(1)==6
	param.compo.z(5)=8;
	param.compo.a(5)=16;
else
	param.compo.z(5)=6;
	param.compo.a(5)=12;
end

if strcmp(param.fonction.impur,'zinebcompo')
	if length(bile.cmain)==2
	   if sum(bile.cmain) > 1
			param.cons.impur.cmin1    =  0.05;
			param.cons.impur.cmin2    =  0;
		else
			param.cons.impur.cmin1    =  bile.cmain(2)./ bile.cmain(1);
			param.cons.impur.cmin2    =  0;
		end
	elseif length(bile.cmain)==3
	   if sum(bile.cmain) > 1
			param.cons.impur.cmin1    =  0.1;
			param.cons.impur.cmin2    =  0.05;
		else
			param.cons.impur.cmin1    =  bile.cmain(2)./bile.cmain(1);
			param.cons.impur.cmin2    =  bile.cmain(3)./bile.cmain(1);
		end
	else
		param.cons.impur.cmin1    =  0;
		param.cons.impur.cmin2    =  0;
	end
	param.cons.impur.rimp     =  0.3;
	param.cons.impur.zeff     =  1;
end

if ~isfinite(param.cons.impur.cmin1)
	param.cons.impur.cmin1=0;
end
if ~isfinite(param.cons.impur.cmin2)
	param.cons.impur.cmin2=0;
end

% 3 - remplissage des signaux simple (pas de rho)

% le temps
data.gene.temps = temps;

% la coordonnees r/a (preevaluation, pour l'importation du profil des xdurs)
data.equi.rhog = rsa;

% la geometrie
if fix(numchoc) <= charniere
    signal.plus =1;
    data.geo.r0     = Rt';
    data.geo.z0     = zt';
    data.geo.a      = at';
    data.geo.e1     = zsample(bile.elong,bile.times,data.gene.temps,signal);
    data.geo.b0     = zsample(rb0,trb0,data.gene.temps,signal) ./ data.geo.r0;
    data.geo.trh1   = zsample(bile.triang,bile.times,data.gene.temps,signal);
    data.geo.trb1   = zsample(bile.triang,bile.times,data.gene.temps,signal);
    data.geo.ind1   = zeros(size(temps));
    data.geo.mode   = 2.*ones(size(temps));
    data.geo.R      = Rext;
    data.geo.Z      = Zext;
else
    signal.plus     = 1;
    data.geo.r0     = zsample(hr0,tgrho,data.gene.temps,signal);
    data.geo.z0     = zsample(hz0,tgrho,data.gene.temps,signal);
    data.geo.a      = zsample(ha,tgrho,data.gene.temps,signal);
    data.geo.e1     = zsample(he1,tgrho,data.gene.temps,signal);
    data.geo.b0     = zsample(rb0,trb0,data.gene.temps,signal) ./ data.geo.r0;
    signal.plus     = 0;
    data.geo.trh1   = zsample(htrh,tgrho,data.gene.temps,signal);
    data.geo.trb1   = zsample(htrl,tgrho,data.gene.temps,signal);
    data.geo.ind1   = zeros(size(temps));
    data.geo.mode   = 2.*ones(size(temps));
    groupe.plus     = 0;
    void            = 1:size(Rext,2);
    data.geo.R      = zsample(Rext,tgrho,void,data.gene.temps,void,groupe);
    data.geo.Z      = zsample(Zext,tgrho,void,data.gene.temps,void,groupe);
end
% consignes
data.cons.ip      = zsample(sipmes.*1e6,tsipmes,data.gene.temps,signal);
data.cons.vloop   = medfilt1(zsample(bile.vs,bile.times,data.gene.temps,signal),3);
data.cons.zeffm   = faczeff .* zsample(bile.zeff,bile.times,data.gene.temps,signal);

% control de zeffm
z1   = max(param.compo.z);
ind1 = find(param.compo.z < z1);
zmax = 2/3 .* max(param.compo.z(ind1)) + 1/3 .* z1;
zmin = 1.1 .* min(param.compo.z);
ind = find(data.cons.zeffm > zmax);
if ~isempty(ind)
  data.cons.zeffm(ind) = zmax .* ones(1,length(ind));
end 
ind = find(data.cons.zeffm < zmin);
if ~isempty(ind)
  data.cons.zeffm(ind) = zmin .* ones(1,length(ind));
end 

if ~isempty(tiprof)
   data.cons.nhnd  = zsample(tiprof.nH ./ tiprof.nD,tiprof.t,data.gene.temps,signal);
   
   % parametre de concentration (zinebcompo)
   if strcmp(param.fonction.impur,'zinebcompo')
        if mfl >=3
            nh       = zsample(medfilt1(tiprof.nH,mfl),tiprof.t,data.gene.temps,signal);
            nd       = zsample(medfilt1(tiprof.nD,mfl),tiprof.t,data.gene.temps,signal);
        else
            nh       = zsample(tiprof.nH,tiprof.t,data.gene.temps,signal);
            nd       = zsample(tiprof.nD,tiprof.t,data.gene.temps,signal);
        end
        ni       = zsample(bile.nifit(:,1),bile.times,data.gene.temps,signal);
	ne       = zsample(bile.nefit(:,1),bile.times,data.gene.temps,signal);
	zeff     = data.cons.zeffm;
	z        = param.compo.z;
	rimp     = param.cons.impur.rimp;
	if (param.compo.z(1) ==1)
	   if (param.compo.a(1) ==1)
	      ni = nh;
	   else
	      ni = nd;
	   end
	elseif ( any((param.compo.z == 1) & (param.compo.a == 1)) & ...
	          any((param.compo.z == 1) & (param.compo.a == 2)))
	    ni = ((ne .* zeff - nd - nh) .* (z(4) + rimp .* z(5)) - ...
	          (ne - nd -nh) .* (z(4)^2  + rimp .* z(5) ^ 2)) ./ ...
		  (z(1) ^ 2 .* (z(4) + rimp .* z(5)) - z(1) .* ...
		  (z(4)^2  + rimp .* z(5) ^ 2));  
        else
	    ni = ni ./ bile.cmain(1);  
	end    
        nhnmain  = nh ./ ni;
	cminh    = mean(nhnmain(isfinite(nhnmain)));
        ndnmain  = nd ./ ni;
	cmind    = mean(ndnmain(isfinite(ndnmain)));
	ind   = find((param.compo.z == 1) & (param.compo.a == 1));
	if ~isempty(ind) & (ind == 2)
	   if ~isempty(cminh)
	      if isfinite(cminh)
	         param.cons.impur.cmin1=max(cminh,0);
	      end 
	   end
	elseif ~isempty(ind) & (ind == 3)
	   if ~isempty(cminh)
	      if isfinite(cminh)
	         param.cons.impur.cmin2=max(cminh,0);
	      end 
	   end
	end
	ind   = find((param.compo.z == 1) & (param.compo.a == 2));
	if ~isempty(ind) & (ind == 2)
	   if ~isempty(cmind)
	      if isfinite(cmind)
	         param.cons.impur.cmin1=max(cmind,0);
	      end 
	   end
	elseif ~isempty(ind) & (ind == 3)
	   if ~isempty(cmind)
	      if isfinite(cmind)
	         param.cons.impur.cmin2=max(cmind,0);
	      end 
	   end
	  end
   end 
end


% securite zeff
ind = find(data.cons.zeffm < min(param.compo.z));
if ~isempty(ind)
	data.cons.zeffm(ind) = min(param.compo.z) .* ones(1,length(ind));
end
ind = find(data.cons.zeffm > max(param.compo.z));
if ~isempty(ind)
	data.cons.zeffm(ind) = max(param.compo.z) .* ones(1,length(ind));
end

% 4 - remplissage des donnees dependants de rho et du temps

% les profils
if (fix(numchoc) > charniere) 
	groupe.plus =0;
    % on recupere le flux au bord du plasma
    [fluxbord,tfluxbord,void,cert] = tsbase(fix(numchoc),'gfluxnoy%4');
    infotsdata = zinfocat(infotsdata,'gfluxnoy','dpolo',cert);
    desctsdata{end+1} = 'dpolo';

    signal.plus      = 0;
    fluxbord         =  zsample(fluxbord,tfluxbord,data.gene.temps,signal) ./ 2 ./ pi;
    psi              =  zsample(bile.psi,bile.times,xb,data.gene.temps,xz,groupe) ./ 2 ./ pi;
    psi              =  psi + (fluxbord - psi(:,end)) * ones(1,size(psi,2)); 
    data.prof.psi    =  psi;
    data.prof.dpsidt = zdxdt(psi,data.gene.temps);
    data.cons.flux   = fluxbord;
    groupe.plus =1;
    data.prof.jmoy         = zsample(bile.jmoy,bile.times,xb,data.gene.temps,xz,groupe);
    disp('Utilisation des donnees de l''equilibre de Tprof')
elseif ~isempty(coup.fluxp) & (fjli ~= 1)
    groupe.plus =0;
    data.prof.psi          = zsample(coup.fluxp,coup.times,xb,data.gene.temps,xz,groupe)./2./pi;
    data.prof.dpsidt       = zdxdt(data.prof.psi,data.gene.temps);
    data.cons.flux         = data.prof.psi(1,end)-(data.prof.psi(:,end)-data.prof.psi(1,end));
    groupe.plus =1;
    data.prof.jmoy         = zsample(coup.jmoy.*1e6,coup.times,xb,data.gene.temps,xz,groupe);
    disp('Utilisation des donnees de Coupole')
else
    groupe.plus =0;
    % on recupere le flux au bord du plasma
    [fluxbord,tfluxbord,void,cert] = tsbase(fix(numchoc),'gfluxnoy%4');
    infotsdata = zinfocat(infotsdata,'gfluxnoy','dpolo',cert);
    desctsdata{end+1} = 'dpolo';
    signal.plus      = 0;
    fluxbord         = - zsample(fluxbord,tfluxbord,data.gene.temps,signal) ./ 2 ./ pi;
    fluxbord         = -(fluxbord - fluxbord(1));
    psi              = - zsample(psi_eq,bile.times,xx,data.gene.temps,xz,groupe) ./ 2 ./ pi;
    psi              = psi + (fluxbord - psi(:,end)) * ones(1,size(psi,2)); 
    data.prof.psi    = psi;
    data.prof.dpsidt = zdxdt(data.prof.psi,data.gene.temps);
    data.cons.flux   = fluxbord;
    groupe.plus =1;
    data.prof.jmoy         = zsample(jmoy,bile.times,xx,data.gene.temps,xz,groupe);
    disp('Utilisation des donnees de l''equilibre de analytique')
end 

[nefit,tnefit,xnefit]  = ztsprotect(bile.nefit,bile.times,xb,'nefit');

data.prof.ne           = zsample(nefit,tnefit,xnefit,data.gene.temps,xz,groupe);
ind = find(data.prof.ne<3e17);
if ~isempty(ind)
        data.prof.ne(ind)   = 3e17 .* ones(1,length(ind));
end

[nifit,tnifit,xnifit]  = ztsprotect(abs(bile.nifit),bile.times,xb,'nifit');
ni           = zsample(nifit,tnifit,xnifit,data.gene.temps,xz,groupe);


ind = find(ni<3e17);
if ~isempty(ind)
        ni(ind)   = 3e17 .* ones(1,length(ind));
end
ind = find(ni < (data.prof.ne./10));
if ~isempty(ind)
        ni(ind)   = data.prof.ne(ind) ./10;
end
ind = find(ni > data.prof.ne);
if ~isempty(ind)
        ni(ind)   = data.prof.ne(ind);
end
% ni initiale peut changer dans zineb
data.prof.ni           =  ni;
data.prof.ae           =  data.prof.ni ./ data.prof.ne;

[tefit,ttefit,xtefit]  = ztsprotect(bile.tefit,bile.times,xb,'tefit');
data.prof.te           = facte .* zsample(tefit.*1e3,ttefit,xtefit,data.gene.temps,xz,groupe);

ind = find(data.prof.te<13.6);
if ~isempty(ind)
        data.prof.te(ind)   = 13.6 .* ones(1,length(ind));
end
data.prof.pe           = data.prof.te .* data.prof.ne .* param.phys.e;


% calcul de ti
% if tiop == 0
%    ti_                    = interp1(bile.rhofit,bile.tifit_,param.gene.x,'cubic','extrap')*1e3;
%    ti_(ti_<13.6)          = 13.6;
%    ti                     = ones(size(data.gene.temps)) * ti_;
% elseif tiop == 1
if tiop == 5
	% calcul utilisant la mesure de profile et le flux de neutron
	try
		[newti,fluxn,Win,tp,ecart,val] = calcTi2(fix(numchoc));
		if isempty(newti)
			tiop = 0;
		else
	 		[ti,tti,xti]           = ztsprotect(newti,bile.times,xb,'ti');
    			data.prof.ti           = zsample(ti.*1e3,tti,xti,data.gene.temps,xz,groupe);	
			disp('using one profile Ti measurement and neutron flux for Ti estimation');		
		end
	catch
		tiop = 0;
	end
end
if tiop == 5
	%rien c'est fait
elseif tiop == 1
    ti                     = bile.nefit .* (bile.ti_ne * ones(1,size(bile.rhofit,2)));
	 [ti,tti,xti]           = ztsprotect(ti,bile.times,xb,'ti');
    data.prof.ti           = zsample(ti.*1e3,tti,xti,data.gene.temps,xz,groupe);

elseif tiop == 2
    ti                     = bile.tefit .* (bile.ti_te * ones(1,size(bile.rhofit,2)));
	 [ti,tti,xti]           = ztsprotect(ti,bile.times,xb,'ti');
    data.prof.ti           = zsample(ti.*1e3,tti,xti,data.gene.temps,xz,groupe);

elseif tiop == 3
    ti                     = bile.tefit;
	 [ti,tti,xti]           = ztsprotect(ti,bile.times,xb,'ti');
    data.prof.ti           = zsample(ti.*1e3,tti,xti,data.gene.temps,xz,groupe);

elseif (isempty(tiprof) | (tiop == 4) | (tiop == 0)) & (tiop ~= 5)
         param.from.option.tiop = 4;
	 rm                = zsample(bile.rm,bile.times,xb,data.gene.temps,xz,groupe);
   
	 pfci_ion  = zeros(length(bile.times),1);
    % sceurite
    if isfield(bile.info,'ffci')
      if isfield(bile.info.ffci,'antenne_1')
	      txt = upper(deblank(bile.info.ffci.antenne_1));
		   txt(txt==' ') =[];
	 	   if ~strcmp(txt,'EH') & ~all(txt =='-');
   		   pfci_ion          = pfci_ion + bile.pfci(:,1); 
         end
	   end
      if isfield(bile.info.ffci,'antenne_2')
	      txt = deblank(bile.info.ffci.antenne_2);
		   txt(txt==' ') =[];
	 	   if ~strcmp(txt,'EH') & ~all(txt =='-');
   		    pfci_ion          = pfci_ion + bile.pfci(:,1); 
         end
	   end
      if isfield(bile.info.ffci,'antenne_3')
	      txt = deblank(bile.info.ffci.antenne_3);
		   txt(txt==' ') =[];
	 	   if ~strcmp(txt,'EH') & ~all(txt =='-');
   		   pfci_ion          = pfci_ion + bile.pfci(:,1); 
         end
	   end
    end
    puiss  =  (1/3) .* zsample(pfci_ion,bile.times,data.gene.temps,signal);
	prad   = zsample(min(bile.prad ,0.5 .* bile.ploss),bile.times,data.gene.temps,signal);
    ae0    = zsample(bile.nifit(:,1) ./ max(1e13,bile.nefit(:,1)),bile.times, data.gene.temps,signal);
    if (length(ae0) ~= length(data.gene.temps) ) || (length(prad) ~= length(data.gene.temps))
       prad_  = min(bile.prad ,0.5 .* bile.ploss);
       ae0_   = bile.nifit(:,1) ./ max(1e13,bile.nefit(:,1));
       puiss_ =  (1/3) .* pfci_ion;
       tloc_  = bile.times;
       indbad = find(diff(tloc_) <= 0);
       while ~isempty(indbad)
            prad_(indbad) = [];
            puiss_(indbad) = [];
            ae0_(indbad) = [];
            tloc_(indbad) = [];
            indbad = find(diff(tloc_) <= 0);
       end
       prad  =  interp1(tloc_,prad_,data.gene.temps,'nearest','extrap'); 
       puiss  =  interp1(tloc_,puiss_,data.gene.temps,'nearest','extrap'); 
       ae0   =  interp1(tloc_,ae0_,data.gene.temps,'nearest','extrap'); 
    end
     
	 [data.prof.ti,sbrag] = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,0,puiss - prad,param.phys,ae0);
	 infotsdata = zinfocat(infotsdata,sbrag,'tbragg');
    desctsdata{end+1} = 'tbragg';

	 
	 % calcul de ndne
	 [ndne,main,ndnev,badti,corti,certn,certi] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                       rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);

    infotsdata = zinfocat(infotsdata,'gfluntn','dneutron',certn);
    desctsdata{end+1} = 'dneutron';
    infotsdata = zinfocat(infotsdata,'siht','didn',certi);
    desctsdata{end+1} = 'didn';

	 % 2ieme tentative si badti == 1
	 if badti == 1
		  disp('Incompatibilite entre ti et le flux de neutron -> correction de ti (1)');
	 	  data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,1,puiss - prad,param.phys,ae0);
	 
	 
		  [ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                        rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
		  if badti == 1 
		  		disp('Incompatibilite entre ti et le flux de neutron -> correction de ti (2)');
	 	  		data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,2, ...
				                        puiss - prad,param.phys,ae0);
	 
	 
		  		[ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                        rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
		      
		  
		  		if badti == 1 
		     			disp('Incompatibilite entre ti et le flux de neutron -> correction impossible');
	 					data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,0, ...
						                        puiss - prad,param.phys,ae0);
						 % calcul de ndne
	 					[ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                       rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
						
				end
		  end				
	 end	
										
	 if main ~= param.compo.z(1)
	 		disp('Mauvaise identification du gaz majoritaire : choc en Deuterium')
			% changement de la composition
	      ind=find(param.compo.z==1 & param.compo.a==2);
			oldz=param.compo.z;
			olda=param.compo.a;
			if isempty(ind)
	 		  disp('Deuterium impose en espece majoritaire')
			  param.compo.z(1)   = 1;
			  param.compo.a(1)   = 2;
			else
			  param.compo.z(1)   = oldz(ind);
			  param.compo.a(1)   = olda(ind);
			  param.compo.z(ind) = oldz(1);
			  param.compo.a(ind) = olda(1);
	      end
	 elseif (main == 1) & (ndne <0.3) 
	     disp('Ti incompatible avec la composition du plasma -> correction de ti (-1)');
	 	  data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,-1, ...
		                          puiss - prad,param.phys,ae0);
		  [ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                        rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
		  if ndne < 0.3 
	     		disp('Ti incompatible avec la composition du plasma -> correction de ti (-2)');
	 	  		data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,-2, ....
				                        puiss - prad,param.phys,ae0);
		  		[ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                        rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
				if (ndne <0.3)										 
		     				disp('Incompatibilite entre ti et le flux de neutron -> correction impossible');
	 						data.prof.ti = ztibragg(numchoc,data.gene.temps,param.gene.x,data.prof.ne,data.prof.te,0, ...
						                        	puiss - prad,param.phys,ae0);
						 	% calcul de ndne
	 						[ndne,main,ndnev,badti,corti] = zndne(numchoc,param.compo.z(1),data.gene.temps,param.gene.x, ...
	                                       				rsa,data.prof.ne,data.prof.ti,rm,data.geo.a,puiss);
						
				end
		  end				
			
	 end
	 
	 if ~isempty(corti)
	      vee          = ones(1,size(data.prof.ti,2));
	      tibord       = data.prof.ti(:,end);
			data.prof.ti = data.prof.ti - tibord * vee;
	 		data.prof.ti = data.prof.ti .* (corti * vee);
			data.prof.ti = data.prof.ti + tibord * vee;
	 end		
	 
	 if ~isempty(ndnev)
	        if length(ndnev) == 1
	 	  data.cons.ext = ndnev*ones(size(data.prof.ti,1),1);
		else
		  data.cons.ext = ndnev;
		end
	 end				
	 
	 % calcul de nhnd
	 try
	 	[tnhnd,ti_h,ti_d,nhnd,va,cert]=risotoauto(fix(numchoc),6);
	 	infotsdata = zinfocat(infotsdata,'risotoauto','neutres',cert);
    		desctsdata{end+1} = 'neutres';
	catch
		tnhnd =[];
                ti_h =[];
                ti_d =[];
                nhnd =[];
	end
	 [zmain_spm,nhnd_spm,nonc_spm,void,cert]=zsmasse(fix(numchoc),param.compo.z(1));
	 infotsdata = zinfocat(infotsdata,'zsmasse','dstd2',cert);
         desctsdata{end+1} = 'dstd2';
    
	 if isempty(tnhnd) & ~isempty(nhnd_spm) 
	   if (nhnd_spm > 0.05) & (nhnd_spm <= 0.5)
	   	disp('Utilisation de nhnd provenant des spectro de masse');
	 		tnhnd = data.gene.temps;
			nhnd  = nhnd_spm .* ones(size(tnhnd));
		end
	 end
	 if ~isempty(nonc_spm) & strcmp(param.fonction.impur,'zinebcompo') 
	      disp('Utilisation de rimp provenant des spectro de masse');
	 		param.cons.impur.rimp = nonc_spm;
	 end
	 
	 if real(zmain_spm) ~=  param.compo.z(1)
	 	if real(zmain_spm) == 1
			legspm = 'Deuterium';
		else
			legspm = 'Helium';
		end
	 	if real(param.compo.z(1)) == 1
			legcron = 'Deuterium';
		else
			legcron = 'Helium';
		end
	 	fprintf('Attention : l''espece majoritaire dans cronos est %s alors que sur le spectro de masse c''est du %s\n', ...
		 			legcron,legspm);
	 end
	 
	 
	 [tvoid,consigne,remplissage,ispi,pompage,contenu,gaz,sgaz] = zchimie(fix(numchoc),data.gene.temps);
    noms =fieldnames(sgaz);
    for k = 1:length(noms)
         infotsdata = zinfocat(infotsdata,noms{k},'dgaz',getfield(sgaz,noms{k}));
    end 
    desctsdata{end+1} = 'dgaz';
    
    % recherhe He3
	 if ~isempty(remplissage)
	 	if max(contenu(:,4)) > 1e18
	      if param.compo.z(1) == 1
	      	cminhe3 = mean(contenu(:,4))./ mean(contenu(:,2));
			else
	      	cminhe3 = mean(contenu(:,4))./ mean(contenu(:,5));
			end
			if ~isfinite(cminhe3)
				cminhe3 = 0;
			end
	 		ind = find((param.compo.z == 2 ) & (param.compo.a == 3));
			if isempty(ind)
				disp('Correction composition plasma -> ajout de He3')
	 			ind = find((param.compo.z == 2 ) & (param.compo.a == 4));
				if isempty(ind)
					if param.compo.z(2) > param.compo.z(3)
						param.compo.a(2) =3;
						param.compo.z(2) =2;
						param.cons.impur.cmin1 = cminhe3;
					else
						param.compo.a(3) =3;
						param.compo.z(3) =2;
						param.cons.impur.cmin2 = cminhe3;
					end
				else
					param.compo.a(ind) =3;
					if ind  == 2
						param.cons.impur.cmin1 = cminhe3;
					else
						param.cons.impur.cmin2 = cminhe3;
					end
				end
			else
				if ind  == 2
					param.cons.impur.cmin1 = cminhe3;
				else
					param.cons.impur.cmin2 = cminhe3;
				end
			end
		end
	 end
	 
	 if ~isempty(tnhnd)
	 		 % choix des voies
	 		 lock 	 = sum((real(nhnd) < 0.5) & (real(nhnd) > 0) & (abs(imag(nhnd)) == 0),1) ./ length(tnhnd);
	 		 indok	 = find(lock > (2/3));
	 		 if ~isempty(indok)
			   	indok 	= max(find(lock == max(lock)));
	 		   	nhndv 	= nhnd(:,indok);
	 		 else
	 		   nhndv 	= sum(nhnd(:,indok),2);
	 		 end
	 		 % selection des temps correct
	 		 nhndv  = medfilt1(nhndv,5);
	 		 indok  = find((real(nhndv) < 0.5) & (real(nhndv) > 0) & (imag(nhndv) == 0));
	 		 nhndv  = real(nhndv);
	 		 nhndm  = mean(nhndv(indok));
	 		 indnok = find((real(nhndv) < (0.5 .* nhndm)) | (real(nhndv) > (1.5 .* nhndm))| (imag(nhndv) ~= 0));
	 		 nhndv(indnok) = nhndm;
	 		 nhnd   = interp1(tnhnd,nhndv,data.gene.temps,'nearest');
	 		 ind    = find(~isfinite(nhnd));
	 		 if ~isempty(ind)
	 		   nhnd(ind) = nhndm;
	 		 end
 	       ind1=find(param.compo.z==1 & param.compo.a==1);
			 oldz=param.compo.z;
			 olda=param.compo.a;
			 if isempty(ind1)
	 		   disp('Hydrogene impose en espece minoritaire 2')
			   param.compo.z(3)   = 1;
			   param.compo.a(3)   = 1;
				ind2  = 2;
				ind1  = 3;
	       elseif ind1 == 2
			   ind2  = 3;
			 else
			   ind2  = 2;
			 end
          
          % blindage 
          nhnd   = max(0,min((0.99 - ndnev) ./ ndnev ,nhnd));
          fzi    = param.compo.z(4);
          ndnev  = min(max(0,(fzi - data.cons.zeffm)) ./ 5 ./ (nhnd + 1) ,ndnev);
          
          if strcmp(param.fonction.impur,'zinebcompo') & (param.compo.z(1) == 1)
			   nd     = (ndnev*ones(1,length(param.gene.x))) .* data.prof.ne;
			   nh     = (nhnd*ones(1,length(param.gene.x))) .* nd;
				z1s    = param.compo.z(ind2);
				z1     = param.compo.z(4);
				z2     = param.compo.z(5);
				r      = param.cons.impur.rimp;				
			   deter  = z1s*(z1^2+z2^2*r)-z1s^2*(z1+z2*r);
				nimp   = data.prof.ne-nh-nd;
				zefnim = (data.cons.zeffm*ones(1,length(param.gene.x))).*data.prof.ne-nh-nd;
				nmin2  = (nimp*(z1^2+z2^2*r)-zefnim*(z1+z2*r))/deter;
				nimp1  = (z1s*zefnim-z1s^2*nimp)/deter;
				nimp2  = nimp1*r;
				indok  = find(all((zefnim > 0) & (nimp1 >= 0) & (nimp2 >= 0) & (nmin2 >= 0),2));
				if isempty(indok)
			 		nhnd  = NaN .* data.gene.temps;
					nhndm = NaN;
			   else
	 		      % calcul des concentrations de minoritaire  solution pour le deuterium
					indx  = find(param.gene.x <= 0.8);
	            if ind1 == 2
					   cmin1 = mean(mean(nh(indok,indx))) ./ mean(mean(nd(indok,indx)));
						cmin2 = mean(mean(nmin2(indok,indx))) ./ mean(mean(nd(indok,indx)));
						data.impur.impur = cat(3,nd,nh,nmin2,nimp1,nimp2);
					else
					   cmin2 = mean(mean(nh(indok,indx))) ./ mean(mean(nd(indok,indx)));
						cmin1 = mean(mean(nmin2(indok,indx))) ./ mean(mean(nd(indok,indx)));
						data.impur.impur = cat(3,nd,nmin2,nh,nimp1,nimp2);
						data.impur.impur = data.impur.impur .* (data.impur.impur > 0);
					end
					param.cons.impur.cmin1 = cmin1;
					param.cons.impur.cmin2 = cmin2;
				end
			 elseif strcmp(param.fonction.impur,'zinebcompo') & (param.compo.z(1) == 2)
			   param.compo.z(2) = 1;
			   param.compo.a(2) = 1;
			   param.compo.z(3) = 1;
			   param.compo.a(3) = 2;
				
			   nd     = (ndnev*ones(1,length(param.gene.x))) .* data.prof.ne;
			   nh     = (nhnd*ones(1,length(param.gene.x))) .* nd;
				z1s    = 2;
				z1     = param.compo.z(4);
				z2     = param.compo.z(5);
				r      = param.cons.impur.rimp;				
			   deter  = z1s*(z1^2+z2^2*r)-z1s^2*(z1+z2*r);
				nimp   = data.prof.ne-nh-nd;
				zefnim = (data.cons.zeffm*ones(1,length(param.gene.x))).*data.prof.ne-nh-nd;
				nhe    = (nimp*(z1^2+z2^2*r)-zefnim*(z1+z2*r))/deter;
				nimp1  = (z1s*zefnim-z1s^2*nimp)/deter;
				nimp2  = nimp1*r;
				indok  = find(all((zefnim > 0) & (nimp1 >= 0) & (nimp2 >= 0) & (nhe >= 0),2));
				
				if isempty(indok)
			 		nhnd  = NaN .* data.gene.temps;
					nhndm = NaN;
	         else
	 		      % calcul des concentrations de minoritaire solution pour l'helium
					indx  = find(param.gene.x <= 0.8);
					cmin1 = mean(mean(nh(indok,indx))) ./ mean(mean(nhe(indok,indx)));
					cmin2 = mean(mean(nd(indok,indx))) ./ mean(mean(nhe(indok,indx)));
				   data.impur.impur = cat(3,nhe,nd,nh,nimp1,nimp2);
					param.cons.impur.cmin1 = cmin1;
					param.cons.impur.cmin2 = cmin2;
	         end
			 end 
			 
	 else
	 		nhnd  = NaN .* data.gene.temps;
			nhndm = NaN;
			
			if strcmp(param.fonction.impur,'zinebcompo') & ~isempty(ndnev)
				% seul nd est connu 
				% recherche de D
				indd  = find((param.compo.z == 1) & (param.compo.a == 2));
				if isempty(indd)
					disp('Petit probleme : pas de Deuterium dans le plasma ...');
				elseif indd == 1
					% cas choc D : concentration du minoritaire de charge min
                    if length(ndnev) == 1
				      nd     = (ndnev*ones(length(data.gene.temps),length(param.gene.x))) .* data.prof.ne;  % ndnev est un scalaire, j'ai du ajouter la dimension du temps devant FI 16/01/08
                    else
 				      nd     = (ndnev*ones(1,length(param.gene.x))) .* data.prof.ne;  % ndnev n'est pas un scalaire, j'ai du enleverr la dimension du temps devant VB 16/10/08                 
                    end
                    ne     = data.prof.ne;
					zeff  = data.cons.zeffm * ones(1,length(param.gene.x));
					% recherche du minoritaire a fixe
					nummin = 2;
					if (param.compo.z(2) == 2) & (param.compo.a(2) == 3)
					   % cas He3
						nummin = 2;
					elseif (param.compo.z(3) == 2) & (param.compo.a(3) == 3)
					   % cas He3
						nummin = 3;
					elseif  param.compo.z(2) < param.compo.z(3)
					   nummin = 2;
					else
					   nummin = 3;
					end
					if nummin == 2
						zm1    = param.compo.z(2);
						zm2    = param.compo.z(3);
						c2     = param.cons.impur.cmin2;
					else
						zm1    = param.compo.z(3);
						zm2    = param.compo.z(2);
						c2     = param.cons.impur.cmin1;
					end
					if (c2 == 0) & (zm2 == 1)
						c2 = 0.05;
					elseif (c2 == 0)
						c2 = 0.003;
					end
					z1     = param.compo.z(4);
					z2     = param.compo.z(5);
					r      = param.cons.impur.rimp;

					d      = zm1 .* (z1 .^ 2 +r .* z2 .^ 2) - zm1 .^2 .*  (z1 +r .* z2);
               nmin1  = ((z1 .^ 2 +r .* z2 .^ 2) .* (ne - nd .* (1 + c2 .* zm2)) -  ...
					          (z1 +r .* z2) .* (zeff .* ne - nd .* (1 + c2 .* zm2 .^ 2))) ./ d;

               nimp1  = (ne - nd .* (1+c2 .* zm2) - nmin1 .* zm1) ./ (z1 + r .* z2);
	 				indok = find(all((nd > 0) & (nimp1 >= 0) & (nmin1 >= 0) & isfinite(nd),2));
					if ~isempty(indok)
						indx  = find(param.gene.x <= 0.8);
						cd    = mean(mean(nmin1(indok,indx))) ./ mean(mean(nd(indok,indx))) ;
						if nummin == 2
							param.cons.impur.cmin1  = cd;
							param.cons.impur.cmin2  = c2;
						else
							param.cons.impur.cmin2 = cd;
							param.cons.impur.cmin1 = c2;
						end
					end
				else
			    	if indd == 2
				 		ind2  = 3;
						c     = param.cons.impur.cmin2;
				 	else
				 		ind2  = 2;
						c     = param.cons.impur.cmin1;
				 	end
					zs    = param.compo.z(ind2);
					zm    = param.compo.z(1); 
				   a     = zm + c .* zs;
					am    = zm .^ 2 + c .* zs .^ 2;
					b     = param.compo.z(4) + param.cons.impur.rimp .* param.compo.z(5);
					bm    = param.compo.z(4) .^ 2 + param.cons.impur.rimp .* param.compo.z(5) .^ 2;
					zeff  = data.cons.zeffm * ones(1,length(param.gene.x));
				   nd    = (ndnev*ones(1,length(param.gene.x))) .* data.prof.ne;
    			   d     = am .* b - bm .* a;
	 				nimp1 =   (am .* (data.prof.ne - nd) - a .* (data.prof.ne  .* zeff - nd)) ./ d;
	 				nmain = - (bm .* (data.prof.ne - nd) - b .* (data.prof.ne  .* zeff - nd)) ./ d;
	 				indok = find(all((nmain > 0) & (nimp1 >= 0) & (nd >= 0) & isfinite(nd),2));
					if ~isempty(indok)
						indx  = find(param.gene.x <= 0.8);
						cd    = mean(mean(nd(indok,indx))) ./ mean(mean(nmain(indok,indx))) ;
						if indd == 2
							param.cons.impur.cmin1 = cd;
				 		else
							param.cons.impur.cmin2 =cd;
				 		end
					end
				end
			end
	 end
	 data.cons.nhnd  = nhnd;
	 
	 % cas particulier de zinebcompo
	if strcmp(param.fonction.impur,'zinebcompo')
	   zeff  = data.cons.zeffm * ones(1,size(data.prof.ne,2));
		[ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(data.prof.ne,zeff,param.cons.impur.cmin1, ...
		        param.cons.impur.cmin2, param.cons.impur.rimp, ...
				  param.compo.z(1),param.compo.z(2),param.compo.z(3),param.compo.z(4),param.compo.z(5));
	    ind     = find(~isfinite(ae));
	    indok   = find(isfinite(ae) & (ae > 0) & (ae <= 1));
		 ae(ind) = mean(ae(indok)); 
		 ae   = min(1,max(ae,0.1));
		 data.prof.ae = ae;
		 data.prof.ni = ae .* data.prof.ne;
		 ni           = data.prof.ni;
	end
else
    if mfl >=3
       for k = 1:size(tiprof.Ti,2)
          tiprof.Ti(:,k) = medfilt1(tiprof.Ti(:,k),mfl);
       end
    end
    xTi                    = ones(size(tiprof.t))*linspace(0,1,size(tiprof.Ti,2));
    data.prof.ti           = zsample(tiprof.Ti.*1e3,tiprof.t,xTi,data.gene.temps,xz,groupe);
   desctsdata{end+1} = 'tiprof';
end

% s'il y a un profile dans tprof
if (tiop == 0) && ~isempty(bile.tifit_)
    % extrapolation au plus proche en temps
    if length(sbile.tifit_.temps) == 1
    	ti_ = ones(size(data.gene.temps)) * sbile.tifit_.data * 1e3;
    else
    	ti_ = interp1(sbile.tifit_.temps,sbile.tifit_.data,data.gene.temps,'nearest','extrap') * 1e3;
    end
    xtifit       = zsample(xb,bile.times,(1:size(xb,2)),data.gene.temps,(1:size(xb,2)),groupe);
    % calcul de la loi d'echelle
    tisl_   = 0.6234e3 .* (max(data.prof.ne,[],2) / 1e19) .^ 0.2846 .* (max(data.prof.te,[],2) / 1e3) .^ 0.6538;
    tisl_(~isfinite(tisl_)) = mean(tisl_(isfinite(tisl_)));
    % regression avec les mesures
    tisl_meas = interp1(data.gene.temps,tisl_,sbile.tifit_.temps,'nearest','extrap');
    [ref_,index_]  =  sort(cat(1,0,max(sbile.tifit_.data,[],2) .* 1e3,1e4));
    meas_ = cat(1,0,tisl_meas,1e4);
    meas_ = meas_(index_);
    tisl_ = pchip(meas_,ref_,tisl_);
    % valeur de la loi d'echelle aux points de mesure (modulation en fonction)
    mod_      =  tisl_ ./ max(13.6,ti_(:,1));
    mod_(mod_ < 0.1) = 0.1;
    mod_(mod_ > 10)  = 10;   
    ti_ = ti_ .* (mod_ * ones(1,size(ti_,2)));
    % reechantillonage en x
    data.prof.ti = zsample(ti_,data.gene.temps,xtifit,data.gene.temps,xz,groupe);
    param.from.option.tiop = 0;

end

ind = find(data.prof.ti<13.6);
if ~isempty(ind)
        data.prof.ti(ind)   = 13.6 .* ones(1,length(ind));
end
ind = find(data.prof.ti < (data.prof.te./5));
if ~isempty(ind)
        data.prof.ti(ind)   = data.prof.te(ind) ./5;
end

% correction puissance rayonnees
pradsl  = (min(max(param.compo.z),max(param.compo.z(1),bile.zeff)) - param.compo.z(1)) ./ 7 .*  ...
          (bile.nbar/10) .^ 2 .* (4*pi^2) .* bile.rm(:,end) .* bile.amin;
rapr    = min( 0.5 .* bile.ploss,bile.prad) ./ max(1e-6,pradsl);
indok   = find(isfinite(rapr) & (rapr > 0.1) & (rapr < 10));
if length(indok) > 10
	frad    = mean(rapr(indok));
else
	frad    = 1;
end
% securite sur la composition
if strcmp(param.fonction.impur,'zinebcompo')
   % secutite sur les minoritaire
	if param.cons.impur.cmin1  >= 0.98 
	     param.cons.impur.cmin1 = 0.98;
	end
	if param.cons.impur.cmin2  >= 0.95 
	     param.cons.impur.cmin2 = 0.95;
	end
	% rayonnement de recombinaison
	[crad,lrad,cert] = zbolo(fix(numchoc));
   infotsdata = zinfocat(infotsdata,'gpbolo','tsbolo',cert);
   desctsdata{end+1} = 'tsbolo';
   desctsdata{end+1} = 'dbolo';

	param.cons.impur.crad = crad;
	param.cons.impur.lrad = lrad;
	param.cons.impur.frad = frad;
	
end




% attention ni et ti sont calculer dans zineb et peuvent etre different en sortie

data.prof.pion         =  data.prof.ti .* data.prof.ni .*param.phys.e;
data.prof.ptot         =  data.prof.pion + data.prof.pe;

% donnees experimentales pour le plot
data.exp.ne0 = data.prof.ne(:,1);
data.exp.nea = data.prof.ne(:,end);
data.exp.ni0 = data.prof.ni(:,1);
data.exp.nia = data.prof.ni(:,end);
data.exp.te0 = data.prof.te(:,1);
data.exp.tea = data.prof.te(:,end);
data.exp.ti0 = data.prof.ti(:,1);
data.exp.tia = data.prof.ti(:,end);
data.exp.ip       = data.cons.ip;
data.exp.vloop    = data.cons.vloop;
data.exp.betadia  = zsample(bile.wdia.*1e6,bile.times,data.gene.temps,signal)./ ...
                           ((3/8).*param.phys.mu0.*data.geo.r0.*data.exp.ip .^2);
data.exp.li       = zsample(bile.li,bile.times,data.gene.temps,signal);
data.exp.qa       = zsample(bile.qpsia,bile.times,data.gene.temps,signal);
if isfield(bile,'q0') & ~isfield(coup,'qpsi')
  data.exp.q0  = zsample(bile.q0,bile.times,data.gene.temps,signal);
elseif ~isempty(coup.qpsi)
   data.exp.q0       = zsample(coup.qpsi(:,1),coup.times,data.gene.temps,signal);
end

% les consignes suites
data.cons.ne1 = data.prof.ne(:,end);
data.cons.te1 = data.prof.te(:,end);
data.cons.ti1 = data.prof.ti(:,end);

data.cons.fci=zeros(size(data.cons.fci)).*1e6;
data.cons.hyb=zeros(size(data.cons.hyb)).*1e6;
data.cons.fce=zeros(size(data.cons.fce)).*1e6;
data.cons.idn=zeros(size(data.cons.idn)).*1e6;

signal.plus =1;

% fci
ind = find(any(bile.pfci'>0.2)');
if isempty(ind)
   phase = pi .*ones(1,size(data.cons.fci,2));
else
   if exist('tiprof','var')
     if ~isempty(tiprof)
       phase = tiprof.phase*pi/180;
     else
       [phase,cert] = lecphase(fix(numchoc),bile.times(min(ind)),bile.times(max(ind)), ...
           bile.pfci,bile.times,sbile.pfci.espace);
       phase  = phase .* pi./180;
       infotsdata = zinfocat(infotsdata,'phfci','dfci',cert);
       desctsdata{end+1} = 'dfci';

     end
   else
      [phase,cert] = lecphase(fix(numchoc),bile.times(min(ind)),bile.times(max(ind)), ...
           bile.pfci,bile.times,sbile.pfci.espace);
       phase  = phase .* pi./180
       infotsdata = zinfocat(infotsdata,'phfci','tfci',cert);
       desctsdata{end+1} = 'tfci';
   end
end
ind =find(~isfinite(phase));
if ~isempty(ind)
   phase(ind) = zeros(1,length(ind));
end

for k=1:min(size(data.cons.fci,2),size(bile.pfci,2))
    data.cons.fci(:,k) = zsample(bile.pfci(:,k),bile.times,data.gene.temps,signal) .* ...
                         exp(i.*phase(k)).*1e6;
end
param.cons.fci.frequence = sbile.pfci.espace;
%m1    = bile.info.ffci.antenne_1;
%ind   = find(m1 <= sprintf(' '));
%if ~isempty(ind)
%	m1(ind) =[];
%end
%m2    = bile.info.ffci.antenne_2;
%ind   = find(m2 <= sprintf(' '));
%if ~isempty(ind)
%	m2(ind) =[];
%end
%m3    = bile.info.ffci.antenne_3;
%ind   = find(m3 <= sprintf(' '));
%if ~isempty(ind)
%	m3(ind) =[];
%end
%
% changement des parametres pour zfci
%
%
% determination du scenario de chaque antenne
%
composition = param.compo;
if composition.a(2) ~= 1 & composition.z(2) ~= 1
  disp ('H is not the first minority species inside CRONOS')
  composition.a(2) = param.compo.a(3);
  composition.a(3) = param.compo.a(2);
  composition.z(2) = param.compo.z(3);
  composition.z(3) = param.compo.z(2);
else
  composition = composition;
end
for k=1:nbfci
  freq = sbile.pfci.espace(k);
  pos  = mean(data.geo.r0) + mean(data.geo.a) + 0.02;
  if ~isnan(k)
    [scenar(k),scenstr(k,:)]=scenarFCI(data.geo,param.phys,freq,pos,composition);
  end	  
end
m1    = scenstr(1,:);	
m2    = scenstr(2,:);	
m3    = scenstr(3,:);
param.cons.fci.code      = [scenar(1),scenar(2),scenar(3)];	
param.cons.fci.mode      = {m1,m2,m3};


% puissance hybride
if ~isempty(phyb)
   % spectres avec n// pic simples
   ind=find(diff(tphyb(:,1))==0);
   if ~isempty(ind)
     tphyb(ind,:) = [];
     phyb(ind,:) = [];
   end
   [npar01,npar02,fi1,fi2] = litphashyb(numchoc,tphyb,phyb);
   infotsdata = zinfocat(infotsdata,'gphashyb','thyb',cert);
   desctsdata{end+1} = 'thyb';
   indpar = zeros(size(phyb,1),2);
	indpar(:,1)    = npar01;
   indpar(:,2)    = npar02;
   signal.plus =1;
   for k=1:2
        data.cons.hyb(:,k) = zsample(phyb(:,k),tphyb(:,k),data.gene.temps,signal).*1e6 .*  ...
                             exp(sqrt(-1).*zsample(indpar(:,k),tphyb(:,k),data.gene.temps,signal));
   end

   % spectres continus (traitements GHYBSPEC1 et GHYBSPEC2)
   % ATTENTION, pour l'instant, le temps n'est pas gere (si variation du spectre pendant le choc ?)
   % a definir en fonction de la strategie du traitement GHYBSPEC
   [x,r,t] = tsbase(fix(numchoc),'GHYBSPEC1');
   if ~isempty(x)
      % redimensionnement des tableaux, puis remplissage
      data.exp.hybspec.npar = zeros(param.gene.nbt,length(r),2);
      data.exp.hybspec.pow  = zeros(param.gene.nbt,length(r),2);
      data.exp.hybspec.npar(:,:,1) = ones(length(data.gene.temps),1) * r';
      data.exp.hybspec.pow(:,:,1) = ones(length(data.gene.temps),1) * x(:,1)';
      [x,r,t] = tsbase(fix(numchoc),'GHYBSPEC2');
      data.exp.hybspec.npar(:,:,2) = ones(length(data.gene.temps),1) * r';
      data.exp.hybspec.pow(:,:,2) = ones(length(data.gene.temps),1) * x(:,1)';
   else
      % redimensionnement des tableaux, puis remplissage avec des NaN
      data.exp.hybspec.npar = ones(param.gene.nbt,2,2);
      data.exp.hybspec.pow  = ones(param.gene.nbt,2,2);
   end
end

if ~isempty(xdur)
   data.prof.xdur = xdur;
   desctsdata{end+1} = 'xdur';
end


% fce
data.cons.fce   = zeros(size(data.cons.fce));

% on essaye d'abord le traitement
occ_ecrh   = tsoccur('tfce',fix(numchoc));

if ~isempty(occ_ecrh )
   times = data.gene.temps;
   if any(occ_ecrh == 0)
         [ecrh,secrh]= cgcgettrait(fix(numchoc),'tfce@');
   elseif any((fix(numchoc) + occ_ecrh/10) == numchoc)
        [ecrh,secrh]= cgcgettrait(numchoc,'tfce@');
   else
        [ecrh,secrh]= cgcgettrait(fix(numchoc) + min(occ_ecrh) / 10,'tfce@');
   end  
   infotsdata = zinfocat(infotsdata,secrh,'tfce');
   indok = find(secrh.p.information.cert > -2);
   if length(indok)>3
      indok =indok(1:3);
      disp('Rema peu traiter que trois gyrotron pour le moment');
   end
   
   param.cons.fce.freq_ghz(1:length(indok))     = mean(ecrh.freq(:,indok));
   param.cons.fce.modpolar(1:length(indok))     = 1 + (secrh.modox.data(indok)  == 1);              % mode polarisation : 1 -> X, 2 -> O [2]
   param.cons.fce.harm_min                       = 1 .* ones(1,param.nombre.fce);                                       %  harmonique min prise en compte [1]
   param.cons.fce.harm_max (1:length(indok))    = min(2,1 + secrh.harm.data(indok));                 % harmonique max prise en compte [2]
   param.cons.fce.rant(1:length(indok))         = secrh.r.data(indok);                                % position en R du dernier miroire (m) [TS = 3.53]
   param.cons.fce.zant(1:length(indok))         = secrh.z.data(indok);                               % position en Z du dernier miroire (m) [0]             
   param.cons.fce.synergie                       = 1 .* ones(1,nbfce);                                 % synergie sur le courant
   param.cons.fce.angle_pol(1:length(indok))    = mean(ecrh.polo(:,indok));                          % angle d'injection dans le plan poloidal par rapport a l'horizontal, en degres (Ts = -20 a 20, >0 vers le haut) [0]
   param.cons.fce.angle_tor(1:length(indok))    = mean(ecrh.toro(:,indok));                          % angle d'injection dans le plan toroidal par rapport a la normal a Btor en degres (Ts = -30 a 30, >0 dans le sens de ip) [10]
   param.cons.fce.save                           = 'No';	                                      % sauvegarde des fichiers contextes
   param.cons.fce.angle_var(1:length(indok))     =  1;	                                     
   
   
   % l'angle de la consigne code pour la position du mirroir 
   % attention le codage est multiplexe :
   % angle_toro =fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
   % angle_polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;
   % il sont en degres
   angle_mul    = (360 + rem(ecrh.polo(:,indok),360)) .* 1e-10  + fix((360 + rem(ecrh.toro(:,indok),360)) .* 1e4) .* 1e-7;
   % angle_toro =fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
   % angle_polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;
   
   % la puissance depend du mode
   ox     =  ones(size(ecrh.p,1),1) * secrh.modox.data';
   pfce   = ecrh.p(:,indok) .* (ecrh.fracox(:,indok) .* (ox(:,indok) == 1) +  ...
            (1 - ecrh.fracox(:,indok)) .* (ox(:,indok) == 2));
            
   
   % on regroupe le tout
   pfce_c = pfce .* 1e6 .* exp(sqrt(-1) .* angle_mul);
   data.cons.fce(:,1:length(indok)) = pfce_c;
   
      
   
else
   	% lecture des donnees ecrh directement dans le diagnostic
	shot = fix(numchoc);
	[xika1,tbad]=tsbase(shot,'sika1');  		% gyrotron A1 cathode current
	[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
	if ~isempty(xika1)  & ~isempty(xHAUTTOR)
		[prxA1,tA1]=tsbase(shot,'spia1');   		% gyrotron A1 power
		[prxA2,tA2]=tsbase(shot,'spia2');   		% gyrotron A2 power
		[sonde,tsonde]=tsbase(shot,'sonderf');   	% RF probe on reflectometer
		%[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
		[xHAUTPOL,tHAUTPOL]=tsbase(shot,'SHAUTPOL');		% Poloidal injection anlge top mirror (A1) hysteresis corrected
		[xMILTOR,tMILTOR]=tsbase(shot,'SMILTOR');		% Toroidal injection anlge central mirror (A2) hysteresis corrected
		[xMILPOL,tMILPOL]=tsbase(shot,'SMILPOL');		% Poloidal injection anlge central mirror (A2) hysteresis corrected
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% powers in kW
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		factA1=1; 		% Calibration factor to be applied to  A1
		factA2=1; 		% Calibration factor to be applied to  A2
		PA1=max(0,1e6*prxA1/factA1);  % A1 power in W
		PA2=max(0,1e6*prxA2/factA2);	% A2 power in W
		Ptot=PA1+PA2;		% total power
		
		% position
		rA2  = 3.5300;
		rA1  = 3.5300;
		zA2  = 0;
		zA1  = 0.2000;

		% puissance
		pecrh1      = max(0,interp1(tbad,PA1,data.gene.temps,'nearest',0));
		pecrh1      = pecrh1 .* (pecrh1 >=1e4);
		pecrh2      = max(0,interp1(tbad,PA2,data.gene.temps,'nearest',0));
		pecrh2      = pecrh2 .* (pecrh2 >=1e4);
		% position 
		phi_tor_a1 = interp1(tbad,xHAUTTOR,data.gene.temps,'nearest','extrap');	
		phi_tor_a2 = interp1(tbad,xMILTOR,data.gene.temps,'nearest','extrap');	
		phi_pol_a1 = interp1(tbad,xHAUTPOL,data.gene.temps,'nearest','extrap');	
		phi_pol_a2 = interp1(tbad,xMILPOL,data.gene.temps,'nearest','extrap');	
		
		% remplissage des parametres
		param.cons.fce.freq_ghz(:)    	              = 118;
		param.cons.fce.modpolar(:)    	  	      = 2;              % mode polarisation : 1 -> X, 2 -> O [2], O par defaut
		param.cons.fce.harm_min(:)                    = 1;                                       %  harmonique min prise en compte [1]
		param.cons.fce.harm_max(:)                    = 2;                 % harmonique max prise en compte [2]
		param.cons.fce.rant(1)                        = rA1;                                % position en R du dernier miroire (m) [TS = 3.53]
		param.cons.fce.zant(1)                        = zA1;                               % position en Z du dernier miroire (m) [0]             
		param.cons.fce.rant(2)                        = rA2;                                % position en R du dernier miroire (m) [TS = 3.53]
		param.cons.fce.zant(2)                        = zA2;                               % position en Z du dernier miroire (m) [0]             
		param.cons.fce.synergie(:)                    = 1;                                 % synergie sur le courant
		param.cons.fce.angle_pol(1)                   = mean(phi_pol_a1);                          % angle d'injection dans le plan poloidal par rapport a l'horizontal, en degres (Ts = -20 a 20, >0 vers le haut) [0]
		param.cons.fce.angle_pol(2)                   = mean(phi_pol_a2);                          % angle d'injection dans le plan poloidal par rapport a l'horizontal, en degres (Ts = -20 a 20, >0 vers le haut) [0]
		param.cons.fce.angle_tor(1)   		      = mean(phi_tor_a1);                          % angle d'injection dans le plan toroidal par rapport a la normal a Btor en degres (Ts = -30 a 30, >0 dans le sens de ip) [10]
		param.cons.fce.angle_tor(2)   		      = mean(phi_tor_a2);                          % angle d'injection dans le plan toroidal par rapport a la normal a Btor en degres (Ts = -30 a 30, >0 dans le sens de ip) [10]
		param.cons.fce.save                           = 'No';	                                      % sauvegarde des fichiers contextes
		param.cons.fce.angle_var(:)     	      =  1;	                                     
		
		
		% l'angle de la consigne code pour la position du mirroir 
		% attention le codage est multiplexe :
		% angle_toro =fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
		% angle_polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;
		% il sont en degres
		angle_mul1    = (360 + rem(phi_pol_a1,360)) .* 1e-10  + fix((360 + rem(phi_tor_a1,360)) .* 1e4) .* 1e-7;
		angle_mul2    = (360 + rem(phi_pol_a2,360)) .* 1e-10  + fix((360 + rem(phi_tor_a2,360)) .* 1e4) .* 1e-7;
		% angle_toro =fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
		% angle_polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;
		
		% la puissance depend du mode
		% on regroupe le tout
		data.cons.fce(:,1) = pecrh1 .* exp(sqrt(-1) .* angle_mul1);
		data.cons.fce(:,2) = pecrh2 .* exp(sqrt(-1) .* angle_mul2);
	else


		[pfce,polar,phi_pol,phi_tor,temps,cert] = zecrh(fix(numchoc),data.gene.temps);
		sert.pfce.information = cert;
		infotsdata = zinfocat(infotsdata,sert,'zecrh');
		desctsdata{end+1} = 'zecrh';
		data.cons.fce(isfinite(pfce))      = pfce(isfinite(pfce)) .* 1e6;
		param.cons.fce.modpolar(isfinite(polar)) = polar(isfinite(polar));
		param.cons.fce.angle_pol(isfinite(phi_pol)) = phi_pol(isfinite(phi_pol));
		param.cons.fce.angle_tor(isfinite(phi_tor)) = phi_tor(isfinite(phi_tor));

	end
end
%if ~isempty(pecrh1)
%    data.cons.fce(:,2) = zsample(pecrh1,tecrh1,data.gene.temps,signal) .* 1e6;  % c'est bon 1 -> 2
%end
%if ~isempty(pecrh2)
%    data.cons.fce(:,1) = zsample(pecrh2,tecrh2,data.gene.temps,signal) .* 1e6;
%end

% donnees pour le bord
%[tvoid,consigne,remplissage,ispi,pompage,contenu,gaz] = zchimie(fix(numchoc),data.gene.temps);
data.cons.c      = zeros(size(data.cons.c));
ind    = find((param.compo.z == 1) & (param.compo.a == 1));
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) = consigne(:,1);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,1)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas d''Hydrogene et de l''hydrogene est injecte');
end
ind    = find((param.compo.z == 1) & (param.compo.a == 2));
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) = consigne(:,2);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,2)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas de Deuterium et du Deuterium est injecte');
end
ind    = find((param.compo.z == 2) & (param.compo.a == 3));
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) = consigne(:,4);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,4)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas d''Helium 3 et de l''Helium 3 est injecte');
end
ind    = find((param.compo.z == 2) & (param.compo.a == 4));
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) =  consigne(:,5);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,5)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas d''Helium  et de l''Helium est injecte');
end
ind    = find(param.compo.z == 6);
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) =  consigne(:,6);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,6)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas de Carbone  et du Carbone est injecte');
end
ind    = find(param.compo.z == 7);
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) =  consigne(:,7);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,7)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas d''Azote  et de l''Azote est injecte');
end
ind    = find(param.compo.z == 10);
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) =  consigne(:,8);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,8)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas de Neon  et du Neon est injecte');
end
ind    = find(param.compo.z == 18);
if ~isempty(ind) & ~isempty(consigne)
	data.cons.c(:,ind) =  consigne(:,9);
elseif isempty(contenu)
	% alors rien
elseif max(contenu(:,9)) > 1e18
	disp('Attention : Le plasma (cronos) ne contient pas d''Argon  et de l''Argon est injecte');
end

signal.plus             = 0;
d                = bile.rm - bile.rmaj * ones(size(bile.rhofit));
dp               = pdederive(bile.amin * bile.rhofit,d,0,2,2,1);
vprb             = (4*pi^2) .* (bile.amin * bile.rhofit) .* ((bile.rmaj * ones(size(bile.rhofit))) + d + (bile.amin * bile.rhofit) .* dp ./ 2);
netot            = bile.amin .* trapz(bile.rhofit,zfiltsvd2m(bile.nefit .*vprb,3),2);
dnetotdt         = pdederive(bile.times,netot,2,2,1,1);
dnetotdt         = zsample(dnetotdt,bile.times,data.gene.temps,signal);
if ~isempty(consigne)
	data.cons.pomp   = max(0,sum(consigne,2) - dnetotdt ./ param.compo.z(1));
else
	data.cons.pomp   = max(0,- dnetotdt ./ param.compo.z(1));
end
data.bord.temp   = 480.* ones(size(data.bord.temp));
%data.bord.nb     = 0.01 .* 4 .* pi .^ 2 .* param.phys.avo .*  ...
%                   data.geo.r0 .* data.geo.a .* sqrt(data.geo.e1); 
if ~isempty(remplissage) & ~isempty(contenu)
	data.bord.nb     = sum(gaz.zgaz .* remplissage,2) + sum((ones(size(data.gene.temps)) * gaz.zgaz ).* contenu,2);
else
   data.bord.nb     = 0.01 .* 4 .* pi .^ 2 .* param.phys.avo .*  ...
                     data.geo.r0 .* data.geo.a .* sqrt(data.geo.e1); 
end


% letcure de la puissance idn
[tvoid,pnbi_filled,beam_energy,frac,amass] = read_nbi_power(fix(numchoc),data.gene.temps);
if any(pnbi_filled > 10e3)
   % connection du module
   param.fonction.idn = 'znemo';
   info = znemo(param.nombre.idn);
   param.cons.idn = info.valeur;
   param.cons.idn.energie   = beam_energy .* 1e3;
   param.cons.idn.fraction1 = frac(1);
   param.cons.idn.fraction2 = frac(2);
   param.cons.idn.fraction3 = frac(3);
   param.cons.idn.machine_nom = 'TS';
   param.cons.idn.masse = amass;
   param.cons.idn.fokker = 'FAST';
   param.cons.idn.directivity = 1;               
   param.cons.idn.injector_config = 'On-axis';
   param.cons.idn.align = 0;
   param.cons.idn.type = 3;
   param.cons.idn.charge = 1;
   param.cons.idn.n_out_profiles = 20;
   param.cons.idn.n_output_2d    = 103;
   param.cons.idn.n_pitch_resol  = 101;
   param.cons.idn.spot_rlong     = 21;
   param.cons.idn.spot_nproc     = 1;
   param.cons.idn.spot_ncreated  = 100;
   param.cons.idn.spot_verbose   = 'No';
   param.cons.idn.spot_init     = 'Zero';
   param.cons.idn.spot_save     = 'No';
   param.cons.idn.spot_anomalous     = 'No';
   param.cons.idn.spot_ano_dcoef1     = 0;
   param.cons.idn.spot_ano_dcoef2     = 0;
   param.cons.idn.spot_ano_rhocut     = 0.5000;
   param.cons.idn.spot_ano_vconv     = 0;
   param.cons.idn.save = 'No';


   % la puissance
   data.cons.idn = pnbi_filled(:);
   % passage en mode calcule
   data.mode.idn(:) = 2;
else
   data.mode.idn(:) = 0;
end


% lecture des donnees pour les equilibres a frontiere libre
[nbpol,cablage] = tsnbpol;
%G0,G1h,G2h,G3h,G4h,G4b,G3b,G2b,G1b.
cablage_base = {'A','Bh','Dh','Eh','Fh','Fb','Eb','Db','Bb'};
data.cons.asser.pfcur = NaN .* ones(length(data.gene.temps),nbpol);
% lecture des donnees
[voltage,t_voltage] = tsbase(fix(numchoc),'gtengmes');
if ~isempty(voltage)
  if size(voltage,2) == length(cablage_base)
	voltage = interp1(t_voltage(:,1),voltage,data.gene.temps,'nearest');
	% prise en compte du cablage des generateurs
	voltage(:,2:end) =  voltage(:,2:end) + voltage(:,1) * ones(1,size(voltage,2)-1);
	for k=1:(nbpol-1)
	ind = strmatch(cablage.name{k},cablage_base,'exact');
	data.cons.asser.pfcur(:,k) = real(data.cons.asser.pfcur(:,k)) + sqrt(-1) .* 1e3 .* voltage(:,ind);
	end
  end
end
%A,Bh,Dh,Eh,Fh,Fb,Eb,Db,Bb
cablage_base = {'A','Bh','Dh','Eh','Fh','Fb','Eb','Db','Bb'};
[current,t_current] = tsbase(fix(numchoc),'gcoubob');
if ~isempty(current)
  if size(current,2) == length(cablage_base)
  	currentr = interp1(t_current(:,1),current,data.gene.temps,'nearest');
  	for k=1:(nbpol-1)
      		ind = strmatch(cablage.name{k},cablage_base,'exact');
      		data.cons.asser.pfcur(:,k) = 1e3 .* currentr(:,ind) + sqrt(-1) .* imag(data.cons.asser.pfcur(:,k));
  	end
   end
  %  premagnetisation
   inda    = strmatch('A',cablage_base,'exact');
   indprem = find(current(:,inda) == max(current(:,inda)),1);
   param.cons.equi.prem_time = t_current(indprem);
   premcur = zeros(1,size(data.cons.asser.pfcur,2));
   for k=1:(nbpol-1)
	ind = strmatch(cablage.name{k},cablage_base,'exact');
	premcur(k) = 1e3 .* current(indprem,ind);
   end
   param.cons.equi.prem_iron = mat2str(premcur);

end


% ajout de la poutre comme pseudo bobine 
[current,t_current] = tsbase(fix(numchoc),'GMAGLPT%5');
if ~isempty(current)
      current = interp1(t_current(:,1),current,data.gene.temps,'nearest');
      data.cons.asser.pfcur(:,end) = 1e3 .* current;
end


% lecture de la rotation
% donnees en kilo rad /s
wphi_lec_ok = 0;
if isfield(sbile,'wphifit') && ~isempty(sbile.wphifit.data)
        wphi_lec_ok = 1;
        sbile.wfit_.data   = sbile.wphifit.data;
        sbile.wfit_.temps  = sbile.wphifit.temps;
        sbile.wfit_.espace = sbile.wphifit.espace;
else
   % essai de charger le fichier
   file_ = sprintf('TS_WPHI_CXS_%d.mat',fix(numchoc));
   if exist(file_) == 2
        data_wphi = load(file_);
        sbile.wfit_.data   = data_wphi.wphifit{1};
        sbile.wfit_.temps  = data_wphi.t;
        sbile.wfit_.espace = data_wphi.rhofit;
        wphi_lec_ok = 1;
   end
end
if wphi_lec_ok == 1
	% extrapolation au plus proche en temps
	if length(sbile.wfit_.temps) == 1
		wphi_ = ones(size(data.gene.temps)) * sbile.wfit_.data;
	else
		wphi_ = interp1(sbile.wfit_.temps,sbile.wfit_.data,data.gene.temps,'nearest','extrap');
	end
        groupe_w = groupe; 
        groupe_w.plus = 0; 
	xwfit       = zsample(xb,bile.times,(1:size(xb,2)),data.gene.temps,(1:size(xb,2)),groupe);
	% reechantillonage en x
	wphi_ = zsample(wphi_,data.gene.temps,xwfit,data.gene.temps,xz,groupe_w);
        % grand rayon correspondant (on suppose la mesure dans le plan median)
        vein    = ones(size(bile.rhofit));
        vtin    = ones(size(bile.times));
	shiftin = (bile.d0 * vein) .* ( 1 - (vtin * bile.rhofit) .^ (bile.piqd * vein));
        rwphi_ = bile.rmaj * vein + bile.amin * bile.rhofit + shiftin;
        %
        % data.prof.vtor_exp = wphi_ .* zsample(rwphi_,bile.times,bile.rhofit,data.gene.temps,xz,groupe_w);
        data.prof.vtor_exp = wphi_;
end


% tranport de la matiere
try
	[de,ve,taup] = ztsmatdata(numchoc,data.gene.temps,rsa);
	if ~isempty(taup)
   		data.gene.taune = taup;
	end
	if ~isempty(de)
   		data.coef.nn = ones(size(data.gene.temps)) * de;
   		data.coef.vn = ones(size(data.gene.temps)) * ve;
	end
catch
	warning('TSACCES : no data for density transport analysis');
end



% lecture du bord si 
try
	[sol,brute] = ztssol(fix(numchoc),1,0);
	if isstruct(brute)
    		noms = fieldnames(brute);
		for k = 1:length(noms);
			infotsdata = setfield(infotsdata,noms{k},getfield(brute,noms{k}));
		end
    	end
	desctsdata{end+1} = 'ztssol';
	
	
	xsol = sol.xsol;
	xsol(~isfinite(xsol)) = - inf;
	if ~any(diff(sol.temps))
	  xsol = 1 +  zsample((max(xsol,[],2) - 1),sol.temps,data.gene.temps,signal) * ...
	           linspace(0,1,size(data.bord.x,2));
        else
	  xsol = 1 +  interp1(sol.temps,(max(xsol,[],2) - 1),data.gene.temps) * ...
	           linspace(0,1,size(data.bord.x,2));
	end
	indnok  = find(any(~isfinite(xsol),2));
	xsol(indnok,:)  = ones(length(indnok),1) * linspace(1,1.5,size(xsol,2));
	matsol = ones(size(sol.xsol,1),1) * linspace(-1,0,size(sol.xsol,2));	   
	sol.xsol(~isfinite(sol.xsol)) = matsol(~isfinite(sol.xsol));
	sol.nd(~isfinite(sol.nd)) = 0;
	sol.nz(~isfinite(sol.nz)) = 0;
    if ~any(diff(sol.temps))
	  nd 		= zsample(sol.nd,sol.temps,sol.xsol,data.gene.temps,xsol,groupe);
	  nc 		= zsample(sol.nz,sol.temps,sol.xsol,data.gene.temps,xsol,groupe);
    else
      nxsol     = interp1(sol.temps,sol.xsol,data.gene.temps);
 	  nnd 		= interp1(sol.temps,sol.nd,data.gene.temps);
 	  nnc 		= interp1(sol.temps,sol.nz,data.gene.temps);
      nd        = zeros(size(xsol));
      nc        = zeros(size(xsol));
      for ii=1:size(xsol,1)
        if sum(isnan(nxsol(ii,:))) == 0
 	      nd(ii,:) = interp1(nxsol(ii,:),nnd(ii,:),xsol(ii,:));
        end
        if sum(isnan(nxsol(ii,:))) == 0
 	      nc(ii,:) = interp1(nxsol(ii,:),nnc(ii,:),xsol(ii,:));
        end
      end
    end
    data.bord.impur      = NaN .* ones(param.gene.nbt,size(xsol,2),param.gene.nbg);
	indd            = find((param.compo.z == 1) &  (param.compo.a == 2));
	if ~isempty(indd)
		data.bord.impur(:,:,indd)  = nd;
	end
	indc            = find(param.compo.z == 6);
	if ~isempty(indc)
		data.bord.impur(:,:,indc)  = nc;
	end
	data.bord.nebord     = zsample(sol.nebord,sol.temps,data.gene.temps,signal);
	data.cons.ne1        = data.bord.nebord; 
	data.bord.nibord     = squeeze(data.bord.impur(:,1,:));
	data.bord.tebord     = zsample(sol.tebord,sol.temps,data.gene.temps,signal);
	data.cons.te1        = data.bord.tebord; 
	data.bord.tibord     = zsample(sol.tibord,sol.temps,data.gene.temps,signal);
	data.cons.ti1        = data.bord.tibord; 
	data.bord.fluxgebord = zsample(sol.ge,sol.temps,data.gene.temps,signal);
	data.cons.ge1        = data.bord.fluxgebord; 
	data.bord.fluxplasma = data.bord.fluxgebord;
	data.bord.fluxmur_c  = zsample(sol.gi,sol.temps,data.gene.temps,signal);
	data.bord.fluxqebord = zsample(sol.qe,sol.temps,data.gene.temps,signal);
	data.cons.qe1        = data.bord.fluxqebord; 
	data.bord.fluxqibord = zsample(sol.qi,sol.temps,data.gene.temps,signal);
	data.cons.qi1        = data.bord.fluxqibord; 
	data.bord.fluxplasma = data.bord.fluxgebord;
	sol.te(~isfinite(sol.te)) = 0;
	sol.ti(~isfinite(sol.ti)) = 0;
	sol.ne(~isfinite(sol.ne)) = 0;
	if ~any(diff(sol.temps))
        data.bord.te         = zsample(sol.te,sol.temps,sol.xsol,data.gene.temps,xsol,groupe);
	    data.bord.ti         = zsample(sol.ti,sol.temps,sol.xsol,data.gene.temps,xsol,groupe);
		data.bord.ne         = zsample(sol.ne,sol.temps,sol.xsol,data.gene.temps,xsol,groupe);
        data.gene.taune      = zsample(sol.taup,sol.temps,data.gene.temps,signal);
    else
        nxsol     = interp1(sol.temps,sol.xsol,data.gene.temps);
 	    nnte 	  = interp1(sol.temps,sol.te,data.gene.temps);
 	    nnti      = interp1(sol.temps,sol.ti,data.gene.temps);
 	    nnne 	  = interp1(sol.temps,sol.ne,data.gene.temps);
 	    nntaune   = interp1(sol.temps,sol.taup,data.gene.temps);
        nte       = zeros(size(xsol));
        nti       = zeros(size(xsol));
        nne       = zeros(size(xsol));
        ntaune    = zeros(size(xsol));
        for ii=1:size(xsol,1)
          if sum(isnan(nxsol(ii,:))) == 0
 	        nte(ii,:) = interp1(nxsol(ii,:),nnte(ii,:),xsol(ii,:));
          end
          if sum(isnan(nxsol(ii,:))) == 0
 	        nti(ii,:) = interp1(nxsol(ii,:),nnti(ii,:),xsol(ii,:));
          end
          if sum(isnan(nxsol(ii,:))) == 0
 	        nne(ii,:) = interp1(nxsol(ii,:),nnne(ii,:),xsol(ii,:));
          end
          %if sum(isnan(nxsol(ii,:))) == 0
 	      %  ntaune(ii,:) = interp1(nxsol(ii,:),nntaune(ii,:),xsol(ii,:));
          %end
        end
       data.bord.te = nte;
       data.bord.ti = nti;
       data.bord.ne = nne;
       data.gene.taune = ntaune;
        
    end
    data.bord.ni         = nc + nd;
	data.bord.x          = xsol;
			

end

% securite taune (effet du calcul de la sol
data.gene.taune = data.gene.taune(:,1);

% securite sur prof.rot pour le module neoclass
if isfield(param.cons.neo,'source_vtor')
	switch param.cons.neo.source_vtor 
	case 'internal'
		data.prof.rot(~isfinite(data.prof.rot)) = 0;
	end
end
data.prof.vtor_exp(~isfinite(data.prof.vtor_exp)) = 0;
data.prof.vpol_exp(~isfinite(data.prof.vpol_exp)) = 0;


% les modes
v1   = ones(param.gene.nbt,1);
v0   = zeros(param.gene.nbt,1);
ind  = param.gene.kmin:1:param.gene.kmax;
von  = v0;
voff = v0;            % mis a zeros
von(ind) = v1(ind); 
vlit   = von;         % donnees en entree
vcalc  = 2 .* von;    % calculees
vcopie = 3 .* von;    % recopie du temps precedent
%
data.mode.impur      = vcalc;  

if option.psimode == 1
	data.mode.psi        = vlit    
else
	data.mode.psi        = vcalc;    
end

if option.nelmode == 1
	data.mode.nel        = vlit;                 
else
	data.mode.nel        = vcalc;                 
end

if option.pemode == 1
	data.mode.pe         = vlit;                
else
	data.mode.pe         = vcalc;                
end

if option.pionmode == 1
	data.mode.pion       = vlit;                 
else
	data.mode.pion       = vcalc;                 
end

data.mode.equi       = vcalc;                 
data.mode.neo        = vcalc;                 
data.mode.fluce      = voff;                 
data.mode.flucion    = voff;                 
data.mode.rot        = voff; 

if option.psilim == 0
	data.mode.cons.psi        = v0;           
elseif option.psilim == 1
	data.mode.cons.psi        = v1;           
else 
	data.mode.cons.psi        = 2 .* v1;           
end

data.mode.cons.ne         = v0;            
data.mode.cons.pe         = v0;            
data.mode.cons.pion       = v0;           
data.mode.cons.fluce      = v0;            
data.mode.cons.flucion    = v0;            
data.mode.cons.rot        = v0;            
data.mode.consbord.ne     = v0;            
data.mode.consbord.te     = v0;            
data.mode.consbord.ti     = v0;            
data.mode.consbord.fluxge = v0;            
data.mode.consbord.fluxqe = v0;            
data.mode.consbord.fluxqi = v0;            
data.mode.mhd.dds    = voff;                 
data.mode.mhd.elm    = voff;                 
data.mode.mhd.limite = voff;                 

% (sources)
data.mode.fci        = vcalc;                 
data.mode.fce        = vcalc;                 
data.mode.hyb        = vcalc;                 
%data.mode.idn        = voff;                 
data.mode.n0         = vcalc;                 
data.mode.bord       = vcalc;                 
data.mode.glacon     = voff;                 
data.mode.fus        = vcalc;                 
data.mode.ohm        = vcalc;                 
data.mode.qneo       = vcalc;                 
data.mode.qei        = vcalc;                 
data.mode.prad       = vcalc;                 
data.mode.brem       = vcalc;                 
data.mode.cyclo      = vcalc;                 

if option.modecoef == 1
	data.mode.ee        =  vcalc;                 
	data.mode.ei        =  voff;                 
	data.mode.en        =  voff;                 
	data.mode.ej        =  voff;                 
	data.mode.ve        =  voff;                 
	data.mode.ep        =  voff;                 
	
	data.mode.ie        =  voff;                 
	data.mode.ii        =  vcalc;                 
	data.mode.in        =  voff;                 
	data.mode.ij        =  voff;                 
	data.mode.vi        =  voff;                 
	data.mode.ip        =  voff;                 
	
	data.mode.ne        =  voff;                 
	data.mode.ni        =  voff;                 
	data.mode.nn        =  vcalc;                 
	data.mode.nj        =  voff;                 
	data.mode.vn        =  vcalc;                 
	
	data.mode.fefe      =  voff;                 
	data.mode.fev       =  voff;                 
	
	data.mode.fifi      =  voff;                 
	data.mode.fiv       =  voff;                 
	
	data.mode.rotc      =  vcalc;                 
	data.mode.rotv      =  voff;                 
else
	data.mode.ee        =  vcalc;                 
	data.mode.ei        =  vcalc;                 
	data.mode.en        =  vcalc;                 
	data.mode.ej        =  vcalc;                 
	data.mode.ve        =  vcalc;                 
	data.mode.ep        =  vcalc;                 
	
	data.mode.ie        =  vcalc;                 
	data.mode.ii        =  vcalc;                 
	data.mode.in        =  vcalc;                 
	data.mode.ij        =  vcalc;                 
	data.mode.vi        =  vcalc;                 
	data.mode.ip        =  vcalc;                 
	
	data.mode.ne        =  vcalc;                 
	data.mode.ni        =  vcalc;                 
	data.mode.nn        =  vcalc;                 
	data.mode.nj        =  vcalc;                 
	data.mode.vn        =  vcalc;                 
	
	data.mode.fefe      =  vcalc;                 
	data.mode.fev       =  vcalc;                 
	
	data.mode.fifi      =  vcalc;                 
	data.mode.fiv       =  vcalc;                 
	
	data.mode.rotc       =  vcalc;                 
	data.mode.rotv      =  vcalc;                 
end	

%mode de variables neo
data.mode.eta       = vcalc;    
data.mode.jboot     = vcalc;    

% autres
if option.plotonoff == 1
	data.mode.plot          = von;  
else
	data.mode.plot          = voff;  
end

data.mode.zeff          = vcalc;    
%data.mode.zeffm         = vlit;    
data.mode.ae            = vcalc;    
data.mode.asser         = voff;    

% options par defaut changees a la demande des utilisateurs de TS
param.gene.ti_invar = 1;
param.cons.impur.exposant = 0;

% mide a jour decripteur de donnees
param.from.source.desc = desctsdata;
% inactive pour recherche de bug
param.from.source.info = zinfo2str_i(infotsdata); 

% modification du nom du fichier de sortie
param.gene.file = strcat(param.gene.origine,'_resultat');
param.plot.pause = 0;
[pp,fp,ep]=fileparts(param.gene.file);
param.gene.rapsauve =fullfile(pp,'rapsauve',fp);

% sauvegarde du fichier
if nargin <3
  try
    % compactage des donnees
    data=zreduit(param,data,'compact');
    % sauvegarde
    post=[];
    savedb(param.gene.origine,'param','data','post');
    % compression du fichier
    zgzip(param.gene.origine,'compress');
    %save(param.gene.origine,'param','data');
    data=zreduit(param,data,'uncompact');
  catch
    disp('Probleme lors de la sauvegarde du fichier :')
    disp(lasterr)
    cr =-11005;
    return
  end
end

disp('The End : data ready')


% calcul de la compositions du plasma 
% en entree : 
%    ne      = densite electronique
%    zeff    = profil de zeff
%    c1      = rapport de la densite du 1er minoriatire sur la densite de l'espece principale
%    c2      = rapport de la densite du 2ieme minoriatire sur la densite de l'espece principale
%    rimp    = rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete
%    zion    = numero atomic (ou charge moyenne) de l'espece principale
%    zmin1   = numero atomic (ou charge moyenne) du 1er minoritaire
%    zmin2   = numero atomic (ou charge moyenne) du 2ieme minoritaire
%    zimp1   = numero atomic (ou charge moyenne) de la 1ere impurete
%    zimp2   = numero atomic (ou charge moyenne) de la 2ieme impurete
%    
% en sortie : ae,nion,nmin1,nmin2,nimp1,nimp2 
% 
% equations :
% ae*ne =nion+nmin1+nmin2+nimp1+nimp2
% ni =ae*ne
% nmin1 = c1 *nion
% nmin2 = c2*nion
% nimp2 = rimp*nimp1
% ne = nion*zion +nmin1 *zmin1 +nmin2 *zmin2 +nimp1*zimp1 +nimp2*zimp2
% ne *zeff = nion*zion^2 +nmin1 *zmin1^2 +nmin2 *zmin2^2 +nimp1*zimp1^2 +nimp2*zimp2^2
% 
function [ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(ne,zeff,c1,c2,rimp,zion,zmin1,zmin2,zimp1,zimp2)

% variables de calcul
de =  - zimp1 .^ 2 .* zion - zimp1 .^ 2 .* c1 .* zmin1 - zimp1 .^ 2 .* c2 .* zmin2 - ...
        rimp .* zimp2 .^ 2 .* zion - rimp .* zimp2 .^ 2 .* c1 .* zmin1 -  ...
        rimp .* zimp2 .^ 2 .* c2 .* zmin2 + zion .^ 2 .* zimp1 +  ...
        zion .^ 2 .* rimp .* zimp2 + c1 .* zmin1 .^ 2 .* zimp1 + ...
        c1 .* zmin1 .^ 2 .* rimp .* zimp2 + c2 .* zmin2 .^ 2 .* zimp1 + ...
        c2 .* zmin2 .^ 2 .* rimp .* zimp2;
       
% especes principales            
nion = -ne .* (zimp1 .^ 2 + rimp .* zimp2 .^ 2 - zeff .* zimp1 - zeff .* rimp .* zimp2) ./ de;
nmin1 = c1.*nion;
nmin2 = c2.*nion;

% impuretees:
nimp1 = (zion .^ 2 + c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - ...
         c1 .* zmin1 .* zeff - c2 .* zmin2 .* zeff ) .* ne ./ de;
nimp2 = rimp.*nimp1;

% rapport somme(ni)/ne :
ae = (- zimp1 .^ 2 - rimp .* zimp2 .^ 2 + zeff .* zimp1 + zeff .* rimp .* zimp2 - ...
        zimp1 .^ 2 .* c1 - c1 .* rimp .* zimp2 .^ 2 + c1 .* zeff .* zimp1 + ...
        c1 .* zeff .* rimp .* zimp2 - c2 .* zimp1 .^ 2 - c2 .* rimp .* zimp2 .^ 2 + ...
        c2 .* zeff .* zimp1 + c2 .* zeff .* rimp .* zimp2 + zion .^ 2 + ...
        c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - c1 .* zmin1 .* zeff - ...
        c2 .* zmin2 .* zeff + zion .^ 2 .* rimp + c1 .* zmin1 .^ 2 .* rimp + ...
        c2 .* zmin2 .^ 2 .* rimp - zion .* zeff .* rimp - c1 .* zmin1 .* zeff .* rimp - ...
        c2 .* zmin2 .* zeff .* rimp) ./ de;



function o = zfiltsvd2m(e,n)

[u,s,v]  = svd(e,0);
so       = s;
s        = diag(s);
er       = sqrt(sum(s .^ 2));
indmax   = min(max(find((s./er) >= 0.05)),3);
s((indmax+1):end) = 0;

if nargin >1 
	for k =1:indmax
		[v1,v2,v3,m] = meanfilt1(u(:,k),n); 
		u(:,k)            = m';
	end
end
s=diag(s);
try
    o         = u * diag(s) * v';
catch
    for k=1:size(s,1)
        so(k,k) = s(k,k);
    end
    o         = u * so * v';
end


