function [cr,data,param] = z0dacces(chemin,option,post)
% cette fonction cree un jeu de donnees cronos a partir des donnees du  0d
% declaration des parametres
warning('This function is no more updated and can be inacurrate for present used !');

if nargin <=1 
	valeur.nbrho       = 101;     % nombre de points radiaux [101]
	valeur.nbnbi       = 16;     % nombre d'injecteur
	valeur.nbicrh      = 3;       % nombre de coupleur icrh
	valeur.nbecrh      = 3;      % nombre d'antenne ECRH
	valeur.nblh        = 1;      % nombre d'antennes LH
	valeur.nbline      = 1;      % nombre de lignes d'injection de glacon
	valeur.psimode     = 2;       % mode Psi : 1 -> interpretatif, 2 -> predictif
	valeur.mapequi     = 'off';       % remplit ou non la structure equi
	valeur.nbp_sepa    = 0;       % nombre de poitn dans la separatrice finale (si utilise le points)
	valeur.frac_sepa   = 1;       % 
	
	type.nbrho       = 'integer';    	
	type.nbnbi       = 'integer';
	type.nbicrh      = 'integer';       
	type.nbecrh      = 'integer';       
	type.nblh        = 'integer';       
	type.nbline      = 'integer';       
	type.psimode     = 'integer';
	type.mapequi     = 'integer';
	type.nbp_sepa    = 'integer';
	type.frac_sepa   = 'real';

	borne.nbrho       = {21,51,101,201,501,1001};    
	borne.nbnbi       = [0,30];     % nombre d'injecteur 
	borne.nbicrh      = [0,4];       % nombre de coupleur icrh
	borne.nbecrh       = [0,3];      % nombre d'antenne ECRH
	borne.nblh        = [0,2];      % nombre d'antennes LH
	borne.nbline      = [0,16];      % nombre de lignes d'injection de glacon
	borne.psimode     = {1,2};
	borne.mapequi     = {'off','on','helena'};
	borne.nbp_sepa    = [0,201];     % 
	borne.frac_sepa   = [0.8,1.2];     % 
	
	defaut            = valeur;
	
	info.nbrho       = 'number of radial point [101]';
        info.psimode     = 'if =  1 -> current profile is given, if = 2,  solve the current diffusion equation';
	info.nbnbi       = 'number of pinis for NBI';
	info.nbicrh    = 'number of ICRH antenna';       % nombre de coupleur icrh
	info.nbecrh     = 'number of ECRH mirror';      % nombre d'antenne ECRH
	info.nblh      = 'number of LH launcher';      % nombre d'antennes LH
	info.nbline    = 'number of lines for pellet injection';      % nombre de lignes d'injection de glacon
	info.mapequi    = 'if =''on'' filled the equi structure; if = ''helena'' then recompute equilibrium with HELENA';      
	info.nbp_sepa    = 'number of point in LCFS of CRONOS, if METIS use LCFS given by point; if = 0, use directly METIS LCFS';      % 
	info.frac_sepa   = 'fraction of r/a at which the LCFS is computed for CRONOS, if METIS use LCFS given by point';      % 

        % parametre utile pour helena
        valeur.modep        = 'Ptot';           % facteur de modulation de l'efficacite de generation de courant par le ballooning
        type.modep          = 'list';
	borne.modep         = {'Ptot','Pth','Pth norm'};
        defaut.modep        = 'Ptot';
        info.modep          = 'Pression used in HELNA (Pth norm = Pth * Wtot / Wth)';


	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Module to create Cronos datas from 0D datas';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end

% compatibilite
nbrho        = option.nbrho;

% cr par defaut
cr =0;

% gestion des entrees
if nargin <3
	cr = 1;
	disp('wrong number of input arguments !');
        return
end

% test d'existance des donnees
if ~isfield(post,'zerod') |  ~isfield(post,'z0dinput')
    cr = 2;
    disp('no 0D datas')
    return
elseif isempty(post.zerod.temps)
    cr = 2;
    disp('no 0D datas')
    return
end  


% pour une utilisation plus simple
cons   = post.z0dinput.cons;
geo    = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
op0d   = post.z0dinput.option;
zs     = post.zerod;
temps  = zs.temps;

% securite nbi
if any(imag(zs.pnbi_th))
  option.nbnbi = max(2,option.nbnbi);
end

% creation du nom du fichier
tdebut = min(temps);
tfin   = max(temps);
numchoc = post.z0dinput.shot;
if isempty(numchoc)
   numchoc =1;
end
machine = post.z0dinput.machine;
if isempty(machine)
   machine = 'OD_unknown';
end
if isempty(chemin)
	chemin=strcat(getenv('HOME'),'/zineb/data');
end
fichier = strcat(chemin,'/zineb',int2str(fix(numchoc)),machine,int2str(round(rem(numchoc,1)*10)),'_0D');

% initialisation de la structure
[cr,data,param]=zinit('',post.zerod.temps,option.nbrho,fichier,[],[],option.nbicrh,option.nbecrh,option.nbnbi,option.nblh,option.nbline,5);
if cr ~=0
	return
end
param.edit.currentfile = fichier;

vt = ones(param.gene.nbt,1);
ve = ones(1,param.gene.nbrho);

% parametres generaux
param.gene.modecoef     = 1;   % modele diffussif/convectif - pas de terme non diagonaux
param.gene.lambda       = 5/2;
param.gene.fast         = 0;
param.plot.onoff        = 0;
param.gene.verbose      = 0;
param.gene.rebuilt      = 1;
param.gene.post         = 0;

% module par default
param.fonction.hyb = 'zhybsimple';
param.fonction.fci = 'zfcifile';
param.fonction.fce = 'zremafile';
param.fonction.idn = 'znemo';
param.fonction.fus = 'zfusion2';
param.fonction.cyclo = 'zexatec';
param.fonction.coefa = 'zkiautonew';

% connexion des modules externes
[cr,data,param] = zconnexion(data,param);
if cr ~=0
	return
end

% parametres des modules
param.fonction.coefa = 'zkiautonew';
infocoefa = zkiautonew;
param.cons.coefa = infocoefa.valeur;
param.fonction.hyb = 'zstarwars_LH';
try
    infolh= zstarwars_LH;
catch
    param.fonction.hyb ='zhybsimple';
    infolh = zhybsimple;
end
param.cons.hyb = infolh.valeur;

param.from.machine =  post.z0dinput.machine;
if isempty(param.from.machine)
	param.from.machine = 'undefined';
end
param.from.shot.num = post.z0dinput.shot;
param.from.shot.date = NaN;
param.from.shot.info = '';
param.from.creation.date =clock;
[s,whoami] = unix('whoami');
param.from.creation.user = whoami(whoami > ' ');
param.from.createur = 'zerod';
param.from.option  = op0d;

% parametre d'echantillnage (cf. z0dsample.m)
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
param.from.sample.signal   = signal;      % parametre d'echantillonage  pour un signal simple =s(temps,1)
param.from.sample.groupe   = groupe;      % parametre d'echantillonage  pour un groupe =g(temps,espace)

% la composition du plasma
if isfield(op0d.option,'Sn_fraction') && (op0d.option.Sn_fraction > 0)
    error('CRONOS code is not adapted for tin (Sn) in plasma composition (option.Sn_fraction should be 0)');
end

% on commence par mettre des 0
param.compo.z = zeros(size(param.compo.z));
param.compo.a = zeros(size(param.compo.a));
% palsma de fond
switch op0d.gaz
case 1
   zj = 1;
   aj = 1;
case 2
   zj = 1;
   aj = 2;
case 3
   zj = 1;
   aj = 2;
case 4
   zj = 2;
   aj = 4;
end
param.compo.z(1) = zj;
param.compo.a(1) = aj;

% choix du minoritaire
switch op0d.mino
case 'He3'
   ag = 3;
   zg = 2;
case 'T'
   ag = 3;
   zg = 1;
case 'He4'
   ag = 4;   
   zg = 2;
case 'D'
   ag = 2;   
   zg = 1;
otherwise
   ag = 1;   
   zg = 1;
end
param.compo.z(2) = zg;
param.compo.a(2) = ag;
switch op0d.gaz
case 3
   if ~any((param.compo.z == 1)&(param.compo.a == 3))
        param.compo.z(3) = 1;
        param.compo.a(3) = 3;
        c2 = 1e-16;
   else
        param.compo.z(3) = 2;
        param.compo.a(3) = 4;
        c2 = op0d.frhe0;
   end
otherwise
   if ~any((param.compo.z == 2)&(param.compo.a == 4))
        param.compo.z(3) = 2;
        param.compo.a(3) = 4;
        c2 = mean(zs.nhem) ./ mean(zs.nDm);
   elseif ~any((param.compo.z == 1)&(param.compo.a == 1))
        param.compo.z(3) = 1;
        param.compo.a(3) = 1;
        c2 = 0.05;
   else
        param.compo.z(3) = 2;
        param.compo.a(3) = 3;
        c2 = 0.05;
   end
end
% impuretes
param.compo.z(4) = fix(op0d.zimp);
[void1,void2,void3,void4,void5,Mimp]  = z0yield2(30,30,param.compo.z(4),47);
param.compo.a(4) = round(Mimp);
param.compo.z(5) = fix(op0d.zmax);
[void1,void2,void3,void4,void5,Mimp]  = z0yield2(30,30,param.compo.z(5),47);
param.compo.a(5) = round(Mimp);
rimp             = op0d.rimp;

if strcmp(param.fonction.impur,'zinebcompo')
	param.cons.impur.cmin1    =  op0d.cmin;
	param.cons.impur.cmin2    =  mean(zs.nhem) ./ mean(zs.nDm);
	param.cons.impur.rimp     =  rimp;
	param.cons.impur.zeff     =  1;
end

% 1- preparation des donnees de l'equilibre

% 2 - remplissage des parametres

%
% le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
% lorsque le tokamak est regarde depuis le haut
%
param.gene.signe.ip  = -op0d.signe;              % signe du courant plasma
param.gene.signe.b0  = -1;              % signe du champ toroidal

% 3 - remplissage des signaux simple (pas de rho)

% le temps
data.gene.temps = temps;
% la geometrie
data.geo.r0     = geo.R;
if isfield(geo,'z0')
	data.geo.z0 = geo.z0;
else
	data.geo.z0     = 0 .* geo.R;
end
data.geo.a      = geo.a;
data.geo.e1     = geo.K;
data.geo.b0     = geo.b0;
data.geo.trh1   = geo.d;
data.geo.trb1   = geo.d;
data.geo.ind1   = zeros(size(temps));
data.geo.mode   = 2.*ones(size(temps));

% calcul de la separatrice
if isfield(exp0d,'Rsepa') & isfield(exp0d,'Zsepa')
        if option.nbp_sepa >= 5 
        	[Rsepa,Zsepa] = z0reshape_sepa(geo,zs,exp0d.Rsepa,exp0d.Zsepa - data.geo.z0 * ones(1,size(exp0d.Zsepa,2)), ...
                                               post.z0dinput.option.morphing,option.frac_sepa,option.nbp_sepa);
		data.geo.R = double(Rsepa);
		data.geo.Z = double(Zsepa) +  data.geo.z0 * ones(1,size(Zsepa,2));

		h = findobj(0,'type','figure','tag','z0daccessepa');
		if isempty(h)
			h=figure('tag','z0daccessepa');
		else
			figure(h);
		end
		clf
		set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
		zplotprof(gca,cons.temps,Rsepa,Zsepa,'color','r','marker','o','linestyle','none');
		zplotprof(gca,cons.temps,exp0d.Rsepa,exp0d.Zsepa,'color','b','marker','none','linestyle','-');
		axis('square')
		axis('equal')
                drawnow
	else
		data.geo.R = double(exp0d.Rsepa);
		data.geo.Z = double(exp0d.Zsepa) + data.geo.z0 * ones(1,size(exp0d.Zsepa,2));		
	end

else
	tu  = asin(max(0,min(1,geo.d)));
	u  = linspace(0,2.* pi,size(data.geo.R,2));
	vu = ones(size(u));
	vt = ones(size(tu));
	Rext  = geo.R *vu + (geo.a * vu) .* cos(vt * u + tu * sin(u));
	Zext  = (geo.a .* geo.K) * sin(u);
	%
	data.geo.R      = Rext;
	data.geo.Z      = Zext + data.geo.z0 * vu;
end

% consignes
data.cons.ip      = zs.ip;
data.cons.vloop   = zs.vloop;
vloop = zs.vloop;
vloop(find(~isfinite(vloop))) = 0;
flux = - cumtrapz(temps,vloop,1) ./2 ./ pi;
data.cons.flux    = flux;
data.cons.zeffm   = zs.zeff;
% minoritaire
param.cons.impur.cmin1 = op0d.cmin;
param.cons.impur.cmin2 = c2;
param.cons.impur.rimp = rimp;
% consigne nhnd, utilise pour le T
data.cons.nhnd  = op0d.cmin + sqrt(-1) .* cons.iso;

% profil de courant et de puissance
swlh = abs(zs.plh ./ max(1,zs.pohm)) > 0.2;
pfweh = max(0,min(1,op0d.fwcd ~=0)) .* zs.picrh_th;
picrh = zs.picrh_th - pfweh;
zeffin = zs.zeff;

profli = post.profil0d;
xlao    = profli.xli;
xli    = profli.rmx ./ (profli.rmx(:,end) * ones(1,size(profli.rmx,2)));
tli    = profli.temps;
jfus   = profli.jfus;
sfus   = profli.salf;
pfus   = profli.pfus;

data.prof.jmoy    = z0dsample(profli.jli,tli,xli,temps,param.gene.x,groupe);
data.prof.jeff    = z0dsample(profli.jeff,tli,xli,temps,param.gene.x,groupe);
data.prof.xdur    = z0dsample(profli.jlh,tli,xli,temps,param.gene.x,groupe);
data.prof.q       = z0dsample(profli.qjli,tli,xli,temps,param.gene.x,groupe);
data.prof.psi     = z0dsample(profli.psi,tli,xli,temps,param.gene.x,groupe);
% securite sur psi pour CRONOS
for k=1:size(data.prof.psi,1)
	psik = data.prof.psi(k,:);
	if any(diff(psik) >= 0)
                npsi  = max(1e-38,abs(max(psik) - min(psik)));
                psi_edge = psik(end);
		delta = - npsi ./ length(psik) ./ 10;
		dpsi = diff(psik);
		indbad = find(dpsi >= 0);
		dpsi(indbad) = delta;
                psik = cumsum(cat(2,delta,dpsi));
                data.prof.psi(k,:) = psik ./ max(1e-38,abs(max(psik) - min(psik))) .* npsi; 
		data.prof.psi(k,:) = data.prof.psi(k,:) - data.prof.psi(k,end) + psi_edge;
	end
end
%
data.prof.dpsidt  = z0dsample(profli.dpsidt,tli,xli,temps,param.gene.x,groupe);
data.prof.epar    = z0dsample(profli.epar,tli,xli,temps,param.gene.x,groupe);
data.prof.ej      = z0dsample(profli.ej,tli,xli,temps,param.gene.x,groupe);

%data.prof.ptot   = (data.prof.pion + data.prof.pe) .* ((zs.w ./zs.wth) * ve);
data.prof.ptot   = z0dsample(profli.ptot,tli,xli,temps,param.gene.x,groupe);
data.equi.ptot   = data.prof.ptot;
data.prof.ne     = z0dsample(profli.nep,tli,xli,temps,param.gene.x,groupe);
data.prof.ni     = z0dsample(profli.nip,tli,xli,temps,param.gene.x,groupe);
data.prof.ae     = data.prof.ni ./ data.prof.ne;
data.prof.zeff   = z0dsample(profli.zeff,tli,xli,temps,param.gene.x,groupe);;

data.prof.te     = z0dsample(profli.tep,tli,xli,temps,param.gene.x,groupe);
data.prof.ti     = z0dsample(profli.tip,tli,xli,temps,param.gene.x,groupe);
data.prof.pe     = param.phys.e .* data.prof.ne .* data.prof.te;
data.prof.pion   = param.phys.e .* data.prof.ni .* data.prof.ti;
data.prof.bpol   = z0dsample(profli.bpol,tli,xli,temps,param.gene.x,groupe);

data.coef.ee     = z0dsample(profli.xie,tli,xli,temps,param.gene.x,groupe) .* data.prof.ne;
data.coef.ii     = z0dsample(profli.xii,tli,xli,temps,param.gene.x,groupe) .* data.prof.ni;


data.equi.vpr     =z0dsample(profli.vpr_tor,tli,xli,temps,param.gene.x,groupe);
data.equi.ri      =z0dsample(profli.ri,tli,xli,temps,param.gene.x,groupe);
data.equi.spr     = data.equi.vpr ./ 2 ./ pi .* data.equi.ri;
data.equi.r2i     = z0dsample(profli.r2i,tli,xli,temps,param.gene.x,groupe);
data.equi.c2c     = z0dsample(profli.C2,tli,xli,temps,param.gene.x,groupe);
data.equi.grho2r2 = z0dsample(profli.grho2r2,tli,xli,temps,param.gene.x,groupe);
data.equi.grho2   = z0dsample(profli.grho2,tli,xli,temps,param.gene.x,groupe);
data.equi.grho    = z0dsample(profli.grho,tli,xli,temps,param.gene.x,groupe);
data.equi.raxe    = z0dsample(profli.Raxe,tli,xli,temps,param.gene.x,groupe);
data.equi.F       = z0dsample(profli.fdia,tli,xli,temps,param.gene.x,groupe);
data.equi.ftrap   = z0dsample(profli.ftrap,tli,xli,temps,param.gene.x,groupe);
rmx_              = z0dsample(profli.rmx,tli,xli,temps,param.gene.x,groupe);
data.equi.rhomax  = rmx_(:,end);
data.equi.jmoy    = data.prof.jmoy ;
data.equi.q       = data.prof.q;
data.equi.psi     = data.prof.psi;
data.equi.phi     = z0dsample(profli.phi,tli,xli,temps,param.gene.x,groupe);

% parametre util pour la convergence de HELENA
param.cons.equi.modep = option.modep;

% securite sur phi pour CRONOS
indbad_k = find(any(diff(data.equi.phi,1,2) <= 0,2));
 
switch option.mapequi
case 'off'
      % rien de plus
case 'helena'
    lprev = -1;
    for k =1:length(temps)
        equi_ok = 0;
        l = find(profli.temps >= temps(k),1);
        if l ~= lprev
            fprintf('%d / %d (@ %d) => new time slice\n',k, length(temps),l);
            try
		%[equi,straightline] = z0dhelena(post,temps(k));	    
                datak = zget1t(data,k);
                [equi,mhd_cd,memoire.memoire,straightline] = ...
                     zequi_helena(param.memoire.equi,param.cons.equi,datak.cons.asser.pfcur,datak.geo,...
                                  datak.equi,datak.prof,param.phys,post.zerod.ip(k), ...
                                  param.gene.x,1);
                if equi.fail
		    error('Helena : Non convergence');
		end
	        lprev = l;
                equi_ok = 1;
	    catch
		fprintf('using METIS internal equilibrium as backup solution\n');

		% donnees de metis 
		prof.xli    = profli.xli;
		prof.temps    = profli.temps(l);
		prof.kx   = profli.kx(l,:);     
		prof.dx   = profli.dx(l,:);      
		prof.rmx  = profli.rmx(l,:);     
		prof.Raxe = profli.Raxe(l,:);
		prof.psi  = profli.psi(l,:);
		prof.fdia = profli.fdia(l,:);
		prof.jmoy = profli.jli(l,:);
		prof.ptot = profli.ptot(l,:);
		% interpolation
		prof = z0dinterpx(prof,length(param.gene.x),'standard');
		%
		geo_input.a     = geo.a(k);
		geo_input.R     = geo.R(k);
		geo_input.K     = geo.K(k);
		geo_input.d     = geo.d(k);
		geo_input.b0    = geo.b0(k); 
		geo_input.z0    = geo.z0(k);
		geo_input.sp    = zs.sp(k);
		geo_input.vp    = zs.vp(k);
		geo_input.sext  = zs.sext(k);
		if isfield(profli,'Rsepa') && isfield(profli,'Zsepa')
		    geo_input.Rsepa = profli.Rsepa(l,:);
		    geo_input.Zsepa = profli.Zsepa(l,:);
		    if all(isfinite(geo_input.Rsepa(:))) && all(isfinite(geo_input.Zsepa(:)))
			% that good
		    else
			geo_input.Rsepa = [];
			geo_input.Zsepa = []; 
			op0d.morphing = 0;
		    end 
		else
		    geo_input.Rsepa = [];
		    geo_input.Zsepa = []; 
		    op0d.morphing = 0;
		end
		%
		phys.mu0 = 4 .* pi .* 1e-7;
		[profequi,deuxd,moment,scalaire,facteur] = metis2equi1t(prof,geo_input,phys,size(data.equi.R,3),1,op0d.morphing,0);

		% remplissage de la structure equi
		% passage en coodornnees de flux
		xout = profequi.rho ./ max(profequi.rho);

		% donn??es scalaire
		data.equi.errit(k)          = 0;               % pas de convergence
		data.equi.fail(k)           = 0;               % fonctionne toujours
		data.equi.amix(k)           = NaN;                 % final value of convergence mixing parameter
		data.equi.bnorme(k)         = NaN;
		data.equi.psi0(k)           = profequi.psi(1);
		data.equi.ip(k)             = scalaire.ipout;              % plasma current from equilibrium (A)
		data.equi.li(k)             = scalaire.liout;              % internal inductance
		data.equi.betap(k)          = scalaire.betap;              % poloidal normalised pressure
		data.equi.betat(k)          = scalaire.betat;              % toroidal normalised pressure
		data.equi.betan(k)          = scalaire.betan;              % current normalised pressure
		
		% calcul de rhomax
		%data.equi.rhomax(k)          = profequi.rho(end);
		%data.equi.q(k,:)             = pchip(xout,profequi.q,param.gene.x);
		data.equi.q0(k)               = data.equi.q(k,1);          
		% calcul du shear
		sigma = [length(param.gene.x) ones(size(data.prof.psi(2:end)))];
		%[qout,qd1,qd2] = interpos(param.gene.x,data.prof.q(k,:),-5,[1 1],[0 1e32],sigma);
		qd1   = z0polyderive(param.gene.x,data.prof.q(k,:),5,param.gene.x);
		data.equi.shear(k,:)           = qd1 ./ data.prof.q(k,:) .* param.gene.x;
		data.prof.shear(k,:)           = data.equi.shear(k,:);
		%data.equi.psi(k,:)             = pchip(xout,profequi.psi,param.gene.x); 
		%data.equi.phi(k,:)             = pchip(xout,profequi.phi,param.gene.x); 
		%data.equi.ftrap(k,:)           = pchip(xout,profequi.ftrap,param.gene.x);
		%data.equi.F(k,:)               = pchip(xout,profequi.fdia,param.gene.x);        
		%data.equi.raxe(k,:)            = pchip(xout,moment.raxe',param.gene.x);   
		data.equi.zaxe(k,:)            = pchip(xout,moment.zaxe',param.gene.x);  
		data.equi.d(k,:)               = data.equi.raxe(k,:) - data.equi.raxe(k,end);
		data.equi.a(k,:)               = pchip(xout,moment.a',param.gene.x);
		data.equi.rhog(k,:)            = data.equi.a(k,:) ./ data.equi.a(k,end);
		data.equi.e(k,:)               = pchip(xout,moment.e',param.gene.x);
		data.equi.trh(k,:)             = pchip(xout,moment.trh',param.gene.x);
		data.equi.trl(k,:)             = pchip(xout,moment.trl',param.gene.x);
		data.equi.indh(k,:)            = NaN .* param.gene.x;
		data.equi.indl(k,:)            = NaN .* param.gene.x;
		data.equi.b2i(k,:)             = pchip(xout,profequi.b2i,param.gene.x);  
		data.equi.b2(k,:)              = pchip(xout,profequi.b2,param.gene.x);
		%data.equi.vpr(k,:)             = pchip(xout,profequi.vpr,param.gene.x);  
		data.equi.volume(k,:)          = pchip(xout,profequi.volume,param.gene.x); 
		%data.equi.spr(k,:)             = pchip(xout,profequi.spr,param.gene.x); 
		data.equi.surface(k,:)         = pchip(xout,profequi.surface,param.gene.x);
		%data.equi.grho2r2(k,:)         = pchip(xout,profequi.grho2r2,param.gene.x);
		data.equi.grhor(k,:)           = pchip(xout,profequi.grhor,param.gene.x);        
		%data.equi.ri(k,:)              = pchip(xout,profequi.ri,param.gene.x);
		%data.equi.grho2(k,:)           = pchip(xout,profequi.grho2,param.gene.x);
		%data.equi.c2c(k,:)             = pchip(xout,profequi.C2,param.gene.x);
		%data.equi.grho(k,:)            = pchip(xout,profequi.grho,param.gene.x);
		data.equi.rmoy(k,:)            = pchip(xout,profequi.rmoy,param.gene.x);
		%data.equi.jmoy(k,:)            = pchip(xout,profequi.jmoy,param.gene.x);
		%data.equi.ptot(k,:)            = pchip(xout,profequi.ptot,param.gene.x);
		data.equi.r2(k,:)              = pchip(xout,profequi.r2,param.gene.x);
		%data.equi.r2i(k,:)             = pchip(xout,profequi.r2i,param.gene.x);
		data.equi.grho2b2(k,:)         = pchip(xout,profequi.grho2b2,param.gene.x);    
		data.equi.errcur(k)    	       = 0;                 % final relative error on current
		
		% utilisation des donnnees du mapping pour les indice de temps problematique
		if any(k == indbad_k)
			data.equi.ftrap(k,:)           = pchip(xout,profequi.ftrap,param.gene.x);
			data.equi.phi(k,:)             = pchip(xout,profequi.phi,param.gene.x); 
			data.equi.vpr(k,:)             = pchip(xout,profequi.vpr,param.gene.x); 
			data.equi.spr(k,:)             = pchip(xout,profequi.spr,param.gene.x);
			data.equi.c2c(k,:)             = pchip(xout,profequi.C2,param.gene.x);  
			data.equi.raxe(k,:)            = pchip(xout,moment.raxe',param.gene.x);   
		end
		% plus utiliser dans cronos
		% equi.r3tau3         = NaN .* equi.r2tau2;
		% equi.r3tau          = NaN .* equi.r3tau; 
		% equi.r2tau2         = NaN .* equi.r2tau2;
			
		% le R et le Z
		data.equi.R(k,:,:)            = single(shiftdim(deuxd.R,-1));
		data.equi.Z(k,:,:)            = single(shiftdim(deuxd.Z,-1)) + geo.z0(k);
		data.equi.rhoRZ(k,:)        = single(profequi.rho)';
		data.equi.psiRZ(k,:)        = single(profequi.psi)';
		data.equi.df2RZ(k,:)        = single(2 .* profequi.fdia .* pdederive(profequi.psi,profequi.fdia,2,2,2,1))';
		data.equi.dprRZ(k,:)        = single(pdederive(profequi.psi,profequi.ptot,2,2,2,1))';
		nb_frmode         = NaN;
		data.equi.frmode(k,:)       = NaN .* data.equi.frmode(k,:)  ;
		% La carte de champ
		data.equi.BR(k,:,:)     = -single(shiftdim(deuxd.BR,-1));
		data.equi.BZ(k,:,:)     = -single(shiftdim(deuxd.BZ,-1));
		data.equi.BPHI(k,:,:)   = single(shiftdim(deuxd.BPHI,-1));
		
		% donnees MHD =====> reprendre ici
		rmx        = profequi.rho ;
		dqdrmx     = pdederive(rmx,profequi.q,0,2,2,1);
		dpdrmx     = pdederive(rmx,profequi.ptot,0,2,2,1);
		dimerc     = 2 .* phys.mu0 .*  profequi.q .^ 2 ./ geo.b0(k) .^ 2 ./ max(eps,rmx) .* dpdrmx  .*  ...
				(1 - 1 ./ profequi.q .^ 2) .* (profequi.q ./ max(eps,dqdrmx)) .^ 2;
		indbad          = find((dqdrmx == 0) | (rmx == 0));
		indok           = find((dqdrmx ~= 0) & (rmx ~= 0));
		dimerc(indbad)  = pchip(rmx(indok),dimerc(indok),rmx(indbad));
		
		inte       = rmx .^ 3 ./ profequi.q .^ 2 - 2 .*  phys.mu0 .* moment.raxe' .^ 2 .* max(eps,rmx) ./ geo.b0(k) .^ 2 .* dpdrmx;
		drmerc     = 2 .* phys.mu0 .*  profequi.q .^ 2 ./ geo.b0(k) .^ 2 ./ max(eps,rmx) .* dpdrmx  .* (profequi.q ./ max(eps,dqdrmx)) .^ 2 .* ...
			(1 - 1 ./ profequi.q .^ 2 + profequi.q ./ max(eps,rmx) .^ 3 .* dqdrmx .* cumtrapz(rmx,inte,2));  
		drmerc(indbad)  = pchip(rmx(indok),drmerc(indok),rmx(indbad));
		
		balcrit    =  0.3 .* (profequi.fdia .* profequi.ri) .^ 2 ./ phys.mu0 .* profequi.ri .* rmx ./ profequi.q .^ 3 .* dqdrmx;
		gptot      =  abs(pdederive(rmx./max(rmx),profequi.ptot,2,2,2,1));
		balcrit    = balcrit ./ max(eps,gptot);
		balcrit(1) = balcrit(2);

		data.equi.mhd.ballooning(k,:)     = pchip(xout,(1- balcrit) .* (balcrit ~= 0),param.gene.x);
		data.equi.mhd.mercier(k,:)	= pchip(xout,drmerc,param.gene.x);
		data.equi.mhd.ideal(k,:)  	= pchip(xout,dimerc,param.gene.x);
                % le dernier equilibre n'est pas valide => calcul pour le prochain temps
                lprev = -1;
	    end
	else
            fprintf('%d / %d (@ %d) => previous time slice\n',k, length(temps),l);
            equi_ok = 1;
	end
        if equi_ok == 1
	    equi.drhomaxdt      = NaN;
	    equi.dvprdt         = NaN .*  equi.vpr;
	    equi.dsprdt         = NaN .*  equi.spr;
	    equi.dphidt         = NaN .*  equi.phi;
	    equi.phid1          = NaN .*  equi.phi;
	    equi.phid2          = NaN .*  equi.phi;
            if length(equi.frmode) > size(data.equi.frmode,2)
		equi.frmode         = equi.frmode(1:size(data.equi.frmode,2));
            elseif length(equi.frmode) < size(data.equi.frmode,2)
		equi.frmode(end+1:size(data.equi.frmode,2)) = 0;
            end
	    datak = zget1t(data,k);
 	    datak.equi = equi;
            datak.prof.mhd_cd = mhd_cd;
	    data  = zput1t(data,k,datak);
        end
    end
otherwise
    % boucle de calcul des sorties manquantes
    for k =1: length(temps)
        % profil le plus proche
        d = abs(profli.temps - temps(k));
        l = find(d == min(d),1);
        fprintf('%d / %d (@ %d)\n',k, length(temps),l);
	% donnees de metis 
	prof.xli    = profli.xli;
	prof.temps    = profli.temps(l);
	prof.kx   = profli.kx(l,:);     
	prof.dx   = profli.dx(l,:);      
	prof.rmx  = profli.rmx(l,:);     
	prof.Raxe = profli.Raxe(l,:);
	prof.psi  = profli.psi(l,:);
	prof.fdia = profli.fdia(l,:);
	prof.jmoy = profli.jli(l,:);
	prof.ptot = profli.ptot(l,:);
        % interpolation
        prof = z0dinterpx(prof,length(param.gene.x),'standard');
	%
	geo_input.a     = geo.a(k);
	geo_input.R     = geo.R(k);
	geo_input.K     = geo.K(k);
	geo_input.d     = geo.d(k);
	geo_input.b0    = geo.b0(k); 
	geo_input.z0    = geo.z0(k);
	geo_input.sp    = zs.sp(k);
	geo_input.vp    = zs.vp(k);
	geo_input.sext  = zs.sext(k);
	if isfield(profli,'Rsepa') &&isfield(profli,'Zsepa')
	    geo_input.Rsepa = profli.Rsepa(l,:);
	    geo_input.Zsepa = profli.Zsepa(l,:);
	    if all(isfinite(geo_input.Rsepa(:))) && all(isfinite(geo_input.Zsepa(:)))
		% that good
	    else
		geo_input.Rsepa = [];
		geo_input.Zsepa = []; 
		op0d.morphing = 0;
	    end 
	else
	    geo_input.Rsepa = [];
	    geo_input.Zsepa = []; 
	    op0d.morphing = 0;
	end
	%
	phys.mu0 = 4 .* pi .* 1e-7;
	[profequi,deuxd,moment,scalaire,facteur] = metis2equi1t(prof,geo_input,phys,size(data.equi.R,3),1,op0d.morphing,0);

        % remplissage de la structure equi
	% passage en coodornnees de flux
	xout = profequi.rho ./ max(profequi.rho);

	% donn??es scalaire
	data.equi.errit(k)          = 0;               % pas de convergence
	data.equi.fail(k)           = 0;               % fonctionne toujours
	data.equi.amix(k)           = NaN;                 % final value of convergence mixing parameter
	data.equi.bnorme(k)         = NaN;
	data.equi.psi0(k)           = profequi.psi(1);
	data.equi.ip(k)             = scalaire.ipout;              % plasma current from equilibrium (A)
	data.equi.li(k)             = scalaire.liout;              % internal inductance
	data.equi.betap(k)          = scalaire.betap;              % poloidal normalised pressure
	data.equi.betat(k)          = scalaire.betat;              % toroidal normalised pressure
	data.equi.betan(k)          = scalaire.betan;              % current normalised pressure
	
	% calcul de rhomax
	%data.equi.rhomax(k)          = profequi.rho(end);
	%data.equi.q(k,:)             = pchip(xout,profequi.q,param.gene.x);
	data.equi.q0(k)               = data.equi.q(k,1);          
	% calcul du shear
	sigma = [length(param.gene.x) ones(size(data.prof.psi(2:end)))];
	%[qout,qd1,qd2] = interpos(param.gene.x,data.prof.q(k,:),-5,[1 1],[0 1e32],sigma);
	[qd1,qd2] = z0polyderive(param.gene.x,data.prof.q(k,:),5,param.gene.x);
	data.equi.shear(k,:)           = qd1 ./ data.prof.q(k,:) .* param.gene.x;
	data.prof.shear(k,:)           = data.equi.shear(k,:);
	%data.equi.psi(k,:)             = pchip(xout,profequi.psi,param.gene.x); 
	%data.equi.phi(k,:)             = pchip(xout,profequi.phi,param.gene.x); 
	%data.equi.ftrap(k,:)           = pchip(xout,profequi.ftrap,param.gene.x);
	%data.equi.F(k,:)               = pchip(xout,profequi.fdia,param.gene.x);        
	%data.equi.raxe(k,:)            = pchip(xout,moment.raxe',param.gene.x);   
	data.equi.zaxe(k,:)            = pchip(xout,moment.zaxe',param.gene.x);  
	data.equi.d(k,:)               = data.equi.raxe(k,:) - data.equi.raxe(k,end);
	data.equi.a(k,:)               = pchip(xout,moment.a',param.gene.x);
	data.equi.rhog(k,:)            = data.equi.a(k,:) ./ data.equi.a(k,end);
	data.equi.e(k,:)               = pchip(xout,moment.e',param.gene.x);
	data.equi.trh(k,:)             = pchip(xout,moment.trh',param.gene.x);
	data.equi.trl(k,:)             = pchip(xout,moment.trl',param.gene.x);
	data.equi.indh(k,:)            = NaN .* param.gene.x;
	data.equi.indl(k,:)            = NaN .* param.gene.x;
	data.equi.b2i(k,:)             = pchip(xout,profequi.b2i,param.gene.x);  
	data.equi.b2(k,:)              = pchip(xout,profequi.b2,param.gene.x);
	%data.equi.vpr(k,:)             = pchip(xout,profequi.vpr,param.gene.x);  
	data.equi.volume(k,:)          = pchip(xout,profequi.volume,param.gene.x); 
	%data.equi.spr(k,:)             = pchip(xout,profequi.spr,param.gene.x); 
	data.equi.surface(k,:)         = pchip(xout,profequi.surface,param.gene.x);
	%data.equi.grho2r2(k,:)         = pchip(xout,profequi.grho2r2,param.gene.x);
	data.equi.grhor(k,:)           = pchip(xout,profequi.grhor,param.gene.x);        
	%data.equi.ri(k,:)              = pchip(xout,profequi.ri,param.gene.x);
	%data.equi.grho2(k,:)           = pchip(xout,profequi.grho2,param.gene.x);
	%data.equi.c2c(k,:)             = pchip(xout,profequi.C2,param.gene.x);
	%data.equi.grho(k,:)            = pchip(xout,profequi.grho,param.gene.x);
	data.equi.rmoy(k,:)            = pchip(xout,profequi.rmoy,param.gene.x);
	%data.equi.jmoy(k,:)            = pchip(xout,profequi.jmoy,param.gene.x);
	%data.equi.ptot(k,:)            = pchip(xout,profequi.ptot,param.gene.x);
	data.equi.r2(k,:)              = pchip(xout,profequi.r2,param.gene.x);
	%data.equi.r2i(k,:)             = pchip(xout,profequi.r2i,param.gene.x);
	data.equi.grho2b2(k,:)         = pchip(xout,profequi.grho2b2,param.gene.x);    
	data.equi.errcur(k)    	       = 0;                 % final relative error on current
	
        % utilisation des donnnees du mapping pour les indice de temps problematique
	if any(k == indbad_k)
		data.equi.ftrap(k,:)           = pchip(xout,profequi.ftrap,param.gene.x);
                data.equi.phi(k,:)             = pchip(xout,profequi.phi,param.gene.x); 
                data.equi.vpr(k,:)             = pchip(xout,profequi.vpr,param.gene.x); 
                data.equi.spr(k,:)             = pchip(xout,profequi.spr,param.gene.x);
                data.equi.c2c(k,:)             = pchip(xout,profequi.C2,param.gene.x);  
                data.equi.raxe(k,:)            = pchip(xout,moment.raxe',param.gene.x);   
        end
	% plus utiliser dans cronos
	% equi.r3tau3         = NaN .* equi.r2tau2;
	% equi.r3tau          = NaN .* equi.r3tau; 
	% equi.r2tau2         = NaN .* equi.r2tau2;
		
	% le R et le Z
	data.equi.R(k,:,:)            = single(shiftdim(deuxd.R,-1));
	data.equi.Z(k,:,:)            = single(shiftdim(deuxd.Z,-1)) + geo.z0(k);
	data.equi.rhoRZ(k,:)        = single(profequi.rho)';
	data.equi.psiRZ(k,:)        = single(profequi.psi)';
	data.equi.df2RZ(k,:)        = single(2 .* profequi.fdia .* pdederive(profequi.psi,profequi.fdia,2,2,2,1))';
	data.equi.dprRZ(k,:)        = single(pdederive(profequi.psi,profequi.ptot,2,2,2,1))';
	nb_frmode         = NaN;
	data.equi.frmode(k,:)       = NaN .* data.equi.frmode(k,:)  ;
	% La carte de champ
	data.equi.BR(k,:,:)     = -single(shiftdim(deuxd.BR,-1));
	data.equi.BZ(k,:,:)     = -single(shiftdim(deuxd.BZ,-1));
	data.equi.BPHI(k,:,:)   = single(shiftdim(deuxd.BPHI,-1));
	
	% donnees MHD =====> reprendre ici
	rmx        = profequi.rho ;
	dqdrmx     = pdederive(rmx,profequi.q,0,2,2,1);
	dpdrmx     = pdederive(rmx,profequi.ptot,0,2,2,1);
	dimerc     = 2 .* phys.mu0 .*  profequi.q .^ 2 ./ geo.b0(k) .^ 2 ./ max(eps,rmx) .* dpdrmx  .*  ...
			(1 - 1 ./ profequi.q .^ 2) .* (profequi.q ./ max(eps,dqdrmx)) .^ 2;
	indbad          = find((dqdrmx == 0) | (rmx == 0));
	indok           = find((dqdrmx ~= 0) & (rmx ~= 0));
	dimerc(indbad)  = pchip(rmx(indok),dimerc(indok),rmx(indbad));
	
	inte       = rmx .^ 3 ./ profequi.q .^ 2 - 2 .*  phys.mu0 .* moment.raxe' .^ 2 .* max(eps,rmx) ./ geo.b0(k) .^ 2 .* dpdrmx;
	drmerc     = 2 .* phys.mu0 .*  profequi.q .^ 2 ./ geo.b0(k) .^ 2 ./ max(eps,rmx) .* dpdrmx  .* (profequi.q ./ max(eps,dqdrmx)) .^ 2 .* ...
		(1 - 1 ./ profequi.q .^ 2 + profequi.q ./ max(eps,rmx) .^ 3 .* dqdrmx .* cumtrapz(rmx,inte,2));  
	drmerc(indbad)  = pchip(rmx(indok),drmerc(indok),rmx(indbad));
	
	balcrit    =  0.3 .* (profequi.fdia .* profequi.ri) .^ 2 ./ phys.mu0 .* profequi.ri .* rmx ./ profequi.q .^ 3 .* dqdrmx;
	gptot      =  abs(pdederive(rmx./max(rmx),profequi.ptot,2,2,2,1));
	balcrit    = balcrit ./ max(eps,gptot);
	balcrit(1) = balcrit(2);

	data.equi.mhd.ballooning(k,:)     = pchip(xout,(1- balcrit) .* (balcrit ~= 0),param.gene.x);
	data.equi.mhd.mercier(k,:)	= pchip(xout,drmerc,param.gene.x);
	data.equi.mhd.ideal(k,:)  	= pchip(xout,dimerc,param.gene.x);
    end
end

% calcul des d??riv??es temporelles
switch option.mapequi
case 'off'
      % rien de plus
otherwise
	    data.equi.drhomaxdt      = z0dxdt(data.equi.rhomax,data.gene.temps);
	    data.equi.dvprdt         = zdxdt(data.equi.vpr,data.gene.temps);
	    data.equi.dsprdt         = zdxdt(data.equi.spr,data.gene.temps);
	    data.equi.dphidt         = zdxdt(data.equi.phi,data.gene.temps);
	    % calcul des derivees
	    data.equi.phid1          = pdederive(param.gene.x,data.equi.phi,0,2,2,1);
	    data.equi.phid2          = pdederive(param.gene.x,data.equi.phi,1,2,2,2);
end

% ne,ni,te,ti, pe,pion
% le profil doivent davbord etre calcule sur 21 points (effet peidestal)
vee              = ones(size(xli));
if size(xli,1) > 1
	u21              =  (1 - xli .^ 2);
else
	u21              = vt * (1 - xli .^ 2);
end

% neoclassique
data.coef.eta     =  z0dsample(profli.eta,tli,xli,temps,param.gene.x,groupe);
data.source.jboot =  max(0,z0dsample(profli.jboot,tli,xli,temps,param.gene.x,groupe));

% donnees "experimetale"
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
data.exp.betadia  = zs.w  ./ (3 ./ 8 .* param.phys.mu0 .* geo.R.* zs.ip .^ 2);
data.exp.li       = zs.li;
data.exp.qa       = zs.qa;
data.exp.q0       = zs.q0;

% les consignes suites
data.cons.ne1 = data.prof.ne(:,end);
data.cons.te1 = data.prof.te(:,end);
data.cons.ti1 = data.prof.ti(:,end);

%[wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,fact] = z0rot(zs,op0d,cons,geo);

noms = fieldnames(profli);
for k = 1:length(noms);
	var = profli.(noms{k});
	if size(var,1) > 1
		warning off
		prof2.(noms{k}) = interp1(profli.temps,var,zs.temps,'nearest');
		warning on
	else
		prof2.(noms{k}) = var;
	end
end
op0d.evolution =0;
[wrad,snbi,void_sicrh,void_sfus,sripth,sriplh,sripicrh,sturb,fact,wrot,slh,prof2] = z0rot3(zs,op0d,cons,geo,prof2);


% injection de neutres
if any(imag(cons.pnbi))
    nb1 = fix(option.nbnbi / 2);
    nb2 = option.nbnbi - nb1;
    data.cons.idn(:,1:nb1)        = (real(cons.pnbi) ./ nb1) * ones(1,nb1);
    data.cons.idn(:,(nb1+1):end)  = (imag(cons.pnbi) ./ nb2) * ones(1,nb2);
    param.fonction.idn         = 'znemo';
    param.cons.idn.energie(1:nb1)       = op0d.einj;
    param.cons.idn.energie((nb1+1):end) = op0d.einj2;
    param.cons.idn.align       = zeros(1,option.nbnbi);
    param.cons.idn.fraction1   = ones(1,option.nbnbi);
    param.cons.idn.fraction2   = zeros(1,option.nbnbi);
    param.cons.idn.fraction3   = zeros(1,option.nbnbi);
    param.cons.idn.masse       = 2 .* ones(1,option.nbnbi);
    param.cons.idn.charge      = ones(1,option.nbnbi);
    param.cons.idn.type        = 3;
    param.cons.idn.debut       = 'Avant';
    param.cons.idn.machine     = machine;
else
    data.cons.idn              = (cons.pnbi ./ option.nbnbi) * ones(1,option.nbnbi);
    param.fonction.idn         = 'znemo';
    param.cons.idn.energie     = op0d.einj .* ones(option.nbnbi,1);
    param.cons.idn.align       = zeros(1,option.nbnbi);
    param.cons.idn.fraction1   = ones(1,option.nbnbi);
    param.cons.idn.fraction2   = zeros(1,option.nbnbi);
    param.cons.idn.fraction3   = zeros(1,option.nbnbi);
    param.cons.idn.masse       = 2 .* ones(1,option.nbnbi);
    param.cons.idn.charge      = ones(1,option.nbnbi);
    param.cons.idn.type        = 3;
    param.cons.idn.debut       = 'Avant';
    param.cons.idn.machine     = machine;
end

pv                         = max(0,z0dsample(real(profli.pnbi) + imag(profli.pnbi),tli,xli,temps,param.gene.x,groupe));
pvion                      = max(0,z0dsample(real(profli.pnbi_ion) + imag(profli.pnbi_ion),tli,xli,temps,param.gene.x,groupe));
data.source.idn.el         = pv - pvion;
data.source.idn.ion        = pvion;
data.source.idn.ne         = max(0,z0dsample(real(profli.nbinesource) + imag(profli.nbinesource),tli,xli,temps,param.gene.x,groupe));
data.source.idn.j          = max(0,z0dsample(real(profli.jnbicd) + imag(profli.jnbicd),tli,xli,temps,param.gene.x,groupe));
data.source.idn.psupra     =( 2 ./ 3 ) .* pv .* (((real(zs.esup_nbi) + imag(zs.esup_nbi)) ./  ...
                            max(1,(real(zs.pnbi_th) + imag(zs.pnbi_th)))) * ve);
data.source.idn.paniso     = 0 .* (vt * ve);
data.source.idn.w          = z0dsample(real(profli.rot_nbi) + imag(profli.rot_nbi),tli,xli,temps,param.gene.x,groupe); 
                             %pv .* (((real(snbi) + imag(snbi)) ./ max(1,real(zs.pnbi_th) + imag(zs.pnbi_th))) * ve);
data.gene.paddidn          = real(zs.pnbi_th) + imag(zs.pnbi_th);
data.gene.iidn             = real(zs.inbicd) + imag(zs.inbicd);

% regalge de fci
data.cons.fci            = (cons.picrh ./ option.nbicrh  .* exp(i.*pi)) * ones(1,option.nbicrh);
param.cons.fci.frequence = op0d.freq .* ones(1,option.nbicrh);
switch op0d.rip 
case 0
   param.cons.fci.rip          = 'No';
case 1
   param.cons.fci.rip         = 'Yes';
end
% mode FCI
switch op0d.fwcd
case 0
   switch op0d.mino
   case 'He3'
      [param.cons.fci.mode{1:option.nbicrh}] = deal('HMIN_He3');
   case 'He4'
      [param.cons.fci.mode{1:option.nbicrh}] = deal('HMIN_He ');
   case 'T'
      [param.cons.fci.mode{1:option.nbicrh}] = deal('HMIN_2T ');
   otherwise
      [param.cons.fci.mode{1:option.nbicrh}] = deal('HMIN_H  ');
   end
case 1
      [param.cons.fci.mode{1:option.nbicrh}] = deal('FWCD    ');
case 2
      [param.cons.fci.mode{1:option.nbicrh}] = deal('FWEH    ');
case -1
      [param.cons.fci.mode{1:option.nbicrh}] = deal('FWCCD   ');
end

pvi                        = max(0,z0dsample(profli.picrh,tli,xli,temps,param.gene.x,groupe));
pvion                      = max(0,z0dsample(profli.picrh_ion,tli,xli,temps,param.gene.x,groupe));
pve                        = max(0,z0dsample(profli.pfweh,tli,xli,temps,param.gene.x,groupe)) .* (op0d.fwcd ~= 0);
data.source.fci.el         = (pvi  - pvion) + pve;
data.source.fci.ion        = pvion;
data.source.fci.j          = z0dsample(profli.jfwcd,tli,xli,temps,param.gene.x,groupe);
data.source.fci.psupra     =( 2 ./ 3 ) .* pvi .* ((zs.esup_icrh ./ max(1,zs.pion_icrh)) * ve);
data.source.fci.paniso     = 0 .* (vt * ve);
data.source.fci.w          = pv .* ((sripicrh ./ max(1,zs.picrh_th)) *ve);
data.gene.paddfci          = zs.picrh_th;
data.gene.ifci             = zs.ifwcd;


% reglage ecrh
data.cons.fce     = (cons.pecrh ./ option.nbecrh) * ones(1,option.nbecrh);
param.cons.largeur  = 0.05 .* ones(1,option.nbecrh);
param.cons.centre   = cons.xece;
param.cons.anglepoloidal   = op0d.angle_ece .* ones(1,option.nbecrh);
switch op0d.sens
case -1
      [param.cons.fce.mode{1:option.nbecrh}] = deal('CECCD');
case 1
      [param.cons.fce.mode{1:option.nbecrh}] = deal('ECCD ');
otherwise
      [param.cons.fce.mode{1:option.nbecrh}] = deal('ECRH ');
end
data.source.fce.el         = max(0,z0dsample(profli.pecrh,tli,xli,temps,param.gene.x,groupe));
data.source.fce.ion        = 0 .* (vt * ve);
data.source.fce.j          = z0dsample(profli.jeccd,tli,xli,temps,param.gene.x,groupe);
data.source.fce.psupra     = 0 .* (vt * ve);
data.source.fce.paniso     = 0 .* (vt * ve);
data.source.fce.w          = 0 .* (vt * ve) ;
data.gene.paddfce          = zs.pecrh;
data.gene.ifce             = zs.ieccd;



% reglage hyb
if abs(op0d.etalh)  <= 1
   npar  = (2.01 - abs(op0d.etalh)) ./ 0.63;
   phase = max(0,(npar - 1.8) .* 230);
else
   npar  = 2;
   phase = max(0,(npar - 1.8) .* 230);
end 
data.cons.hyb     = (exp(sqrt(-1) .* phase./ 180  .* pi)  .* cons.plh ./ option.nblh) *ones(1,option.nblh);
param.cons.hyb.centre  = op0d.xlh .* ones(1,option.nblh);
param.cons.hyb.largeur  = op0d.dlh .* ones(1,option.nblh);
param.cons.hyb.etavar   = 'No';
param.cons.hyb.etachaud = 'Yes';
param.cons.hyb.eta = zpmean(temps,cons.plh,zs.etalh0);
param.cons.hyb.ripple = op0d.rip;
param.cons.hyb.xdur   = 1;
switch op0d.lhmode
case {0,1,2} 
   param.cons.hyb.scaling = 'Input';
case 3
   param.cons.hyb.scaling = 'Goniche';
case 4
   param.cons.hyb.scaling = 'SimulTS';
end
data.source.hyb.el         = max(0,z0dsample(profli.plh,tli,xli,temps,param.gene.x,groupe));
data.source.hyb.ion        = 0 .* (vt * ve);
data.source.hyb.j          = z0dsample(profli.jlh,tli,xli,temps,param.gene.x,groupe);
data.source.hyb.psupra     =( 2 ./ 3 ) .* data.source.hyb.el .* ((zs.esup_lh ./ max(1,zs.plh))*ve);
data.source.hyb.paniso     = 0 .* (vt * ve);
data.source.hyb.w          = z0dsample(profli.rot_lh,tli,xli,temps,param.gene.x,groupe) + data.source.hyb.el.* ((sriplh ./ max(1,zs.plh))*ve);
data.gene.paddhyb          = zs.plh_th;
data.gene.ihyb             = zs.ilh;


% sortie pfus
pv                         = max(0,z0dsample(profli.pfus,tli,xli,temps,param.gene.x,groupe));
pvion                      = max(0,z0dsample(profli.pfus_ion,tli,xli,temps,param.gene.x,groupe));
data.source.fus.el         = pv - pvion;
data.source.fus.ion        = pvion;
data.source.fus.j          = z0dsample(jfus,tli,xli,temps,param.gene.x,groupe);
data.source.fus.psupra     =( 2 ./ 3 ) .* pv .* ((zs.esup_fus ./ max(1,zs.pfus_th))*ve);
data.source.fus.paniso     = 0 .* (vt * ve);
data.source.fus.w          = 0 .* (vt * ve);
data.gene.paddfus          = zs.pfus;
data.gene.ifus             = zs.ifus;

% sortie neutral
data.source.n0.el       = - max(0,z0dsample(profli.pioniz,tli,xli,temps,param.gene.x,groupe));
data.source.n0.ion      = 0 .* (vt * ve);
data.source.n0.ne       = max(0,z0dsample(profli.s0+profli.s0m,tli,xli,temps,param.gene.x,groupe));
data.source.n0.j        = 0 .* (vt * ve);
data.source.n0.w        = z0dsample(profli.rot_n0,tli,xli,temps,param.gene.x,groupe);


% sortie ripple
shape = profli.bpol ./ sqrt(profli.tep);
shape = shape ./ max(eps,trapz(profli.xli,shape .* profli.vpr,2) * ones(1,size(shape,2)));
shape = z0dsample(shape,tli,xli,temps,param.gene.x,groupe);
data.source.rip.w       = (sripth * ones(1,size(shape,2))) .* shape;

% source totale
% calcul de l'equipartiton
% constantes
gg   = 0.64e15;
ee   = 0.1602176462e-18;
me   = 9.10938188e-31;
mp   = 1.6726485e-27;
c    =   2.99792458e8;
mu0  =   4*pi*1e-7;
epsi0 = 1./c.^2./mu0;
A = aj;

% excurtion autour de 1
rr  = 10;
telim = 1.2e5;

% calcul du facteur z^2/a
nDm                =  interp1(zs.temps,zs.nDm,tli,'nearest','extrap');
n1m                =  interp1(zs.temps,zs.n1m,tli,'nearest','extrap');
nTm                =  interp1(zs.temps,zs.nTm,tli,'nearest','extrap');
nhem               =  interp1(zs.temps,zs.nhem,tli,'nearest','extrap');
nimpm              =  interp1(zs.temps,zs.nimpm,tli,'nearest','extrap');
z2sa = max(0,n1m - nDm - nTm) + (1/2) .* nDm + (1/3) .* nTm + nhem  + (op0d.zimp./ (7/3)) .* nimpm + op0d.rimp .* (op0d.zmax./ (7/3)) .* nimpm;
z2sa = z2sa ./ (n1m + (1 + op0d.rimp) .* nimpm + nhem);
indnok = find(~isfinite(z2sa));
repli = ( 1 + (A == 4)) .^ 2 ./ A;
z2sa(indnok) = repli(indnok);

%
nDp                =  profli.n1p .* ((nDm ./ n1m) * ones(1,size(xli,2)));
nTp                =  profli.n1p .* ((nTm ./ n1m) * ones(1,size(xli,2)));
nHp                =  max(0,profli.n1p - nDp - nTp);
%
lnei          =  14.9 - 0.5.*log(profli.nep ./ 1e20) + log(profli.tep ./ 1e3);
warning on
ind = find(~isfinite(lnei) | (lnei <10));
if ~isempty(ind)
	lnei(ind) = 10 .* ones(1,length(ind));
end
taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (epsi0 .^ 2 ./ ee .^ 4 .* sqrt(me))  .* ...
		((ee .* profli.tep) .^ (3/2) ./ profli.nep ./ lnei);
qeib0    = 3 .* me ./ mp ./ taues .* (nHp + nDp ./ 2 + nTp ./ 3 + profli.nhep  +  ...
		(op0d.zimp./ (7/3)) .* profli.nzp + rimp .*  (op0d.zmax./ (7/3)) .* profli.nzp + ...
                183.84 .* profli.nwp);
qeib     = qeib0 .* ee .* (profli.tep - profli.tip);


% source totale
data.source.qei     = z0dsample(qeib,tli,xli,temps,param.gene.x,groupe);
data.source.totale.el      = z0dsample(profli.source_el,tli,xli,temps,param.gene.x,groupe) - ...
                             data.source.qei;
data.source.totale.ion      = z0dsample(profli.source_ion,tli,xli,temps,param.gene.x,groupe) + ...
                             data.source.qei;
data.source.totale.j        = z0dsample(profli.jni,tli,xli,temps,param.gene.x,groupe);
data.source.totale.psupra   = data.source.fus.psupra + data.source.hyb.psupra + data.source.fci.psupra + data.source.idn.psupra; 
data.source.totale.paniso   = 0 .* (vt * ve);
data.source.totale.w        = data.source.fus.w + data.source.hyb.w + data.source.fci.w + data.source.idn.w + data.source.n0.w + data.source.rip.w; 
% rayonnement
data.source.prad       =  max(0,z0dsample(profli.prad,tli,xli,temps,param.gene.x,groupe));
data.source.brem       =  max(0,z0dsample(profli.pbrem,tli,xli,temps,param.gene.x,groupe));
data.source.cyclo      =  max(0,z0dsample(profli.pcyclo,tli,xli,temps,param.gene.x,groupe));
data.source.ohm        =  z0dsample(profli.pohm,tli,xli,temps,param.gene.x,groupe);

% ripple
% pas calculer

% n0 ?
% pas calculer

% ext 
% pas de sources externes autres

% impuretes
n1                =  z0dsample(profli.n1p,tli,xli,temps,param.gene.x,groupe);
n1m               =  trapz(param.gene.x,data.equi.vpr .* n1,2) ./ trapz(param.gene.x,data.equi.vpr,2);
nD                =  n1 .* ((zs.nDm ./ zs.n1m) * ve);
nT                =  n1 .* ((zs.nTm ./ zs.n1m) * ve);
nH                =  max(0,n1 - nD - nT);
nHe               =  z0dsample(profli.nhep,tli,xli,temps,param.gene.x,groupe);
nimp              =  z0dsample(profli.nzp,tli,xli,temps,param.gene.x,groupe);

data.impur.impur  =  zeros(size(data.impur.impur));
indD              =  find ((param.compo.a == 2) &(param.compo.z == 1));
if ~isempty(indD)
   data.impur.impur(:,:,indD) = nD;
end
indT              =  find ((param.compo.a == 3) &(param.compo.z == 1));
if ~isempty(indT)
   data.impur.impur(:,:,indT) = nT;
end
indH              =  find ((param.compo.a == 1) &(param.compo.z == 1));
if length(indH) > 1
  error('miss parametrise METIS run: in H , ICRH minority species can''t be H');
end
if ~isempty(indH)
   data.impur.impur(:,:,indH) = nH;
end
indHe              =  find ((param.compo.a == 4) & (param.compo.z == 2));
if ~isempty(indHe)
   data.impur.impur(:,:,indHe) = nHe;
else
   indHe              =  find ((param.compo.a == 3) & (param.compo.z == 2));
   if ~isempty(indHe)
      data.impur.impur(:,:,indHe) = nHe;
   end
end
data.impur.impur(:,:,end-1 ) = nimp;
data.impur.impur(:,:,end)    = rimp .* nimp;
data.impur.zeff   =  data.prof.zeff;
data.impur.ae     =  data.prof.ae;

% bord
data.bord.tebord =  zs.tebord;
data.bord.tibord =  zs.tebord;
data.bord.nebord =  zs.nebord;
data.bord.nibord =  (zs.nebord .* data.impur.ae(:,end)) * ones(1,param.gene.nbg);


% la rotation et  le champ electrique radial
data.prof.rot              = z0dsample(profli.vtor,profli.temps,profli.xli,temps,param.gene.x,groupe);
data.neo.eta               = data.coef.eta;
data.neo.jboot             = data.source.jboot;
data.neo.vtheta(:,:,1)     = z0dsample(profli.vtheta,profli.temps,profli.xli,temps,param.gene.x,groupe);
data.neo.vtor(:,:,1)       = z0dsample(profli.vtor,profli.temps,profli.xli,temps,param.gene.x,groupe);
data.neo.er                = z0dsample(profli.er,profli.temps,profli.xli,temps,param.gene.x,groupe);



% autre donnees

% regalge du split dtmax
taue               = mean(zs.taue(isfinite(zs.taue)));
param.split.dtmax = max(taue ./ 10,param.split.dtmax);


% mise en coherence ptot metis et cronos (difference de forme)
%ptot_cronos   = data.prof.pe + data.prof.pion + data.source.totale.psupra;
%psupra_metis  = (((zs.w - zs.wth) ./ max(eps,zs.w)) * ve) .* z0dsample(profli.ptot,profli.temps,profli.xli,temps,param.gene.x,groupe);
%factor = psupra_metis ./ max(1,data.source.totale.psupra);
ptot_cronos   = data.prof.pe + data.prof.pion + data.source.totale.psupra;
ptot_metis    = data.equi.ptot;
% correction 
delta_ptot    = (ptot_metis - ptot_cronos) ./ max(eps,data.source.totale.psupra);
noms = fieldnames(data.source);
for k = 1: length(noms)
        if isfield(data.source.(noms{k}),'psupra')
		data.source.(noms{k}).psupra = data.source.(noms{k}).psupra + delta_ptot .* data.source.(noms{k}).psupra;
        end
end
% sources particulieres
% qei, prad brem et cyclo ne sont pas calcule ici, disponible uniquement apres le run

% data.gene
% pas dupliquer (disponible uniquement apres le run)
 
% les modes
% les modes
v1   = ones(param.gene.nbt,1);
v0   = zeros(param.gene.nbt,1);
ind  = param.gene.kmin:1:param.gene.kmax;
voff = v0;            % mis a zeros
von = v1; 
vlit   = von;         % donnees en entree
vcalc  = 2 .* von;    % calculees
vcopie = 3 .* von;    % recopie du temps precedent
%
data.mode.impur      = vcalc;  
data.mode.psi        = option.psimode .* von;    
data.mode.nel        = vlit;                 
data.mode.pe         = vlit;                
data.mode.pion       = vlit;                 
data.mode.equi       = vcalc;                 
data.mode.neo        = vcalc;                 
data.mode.fluce      = voff;                 
data.mode.flucion    = voff;                 
data.mode.rot        = voff; 
data.mode.cons.psi        = v0;           
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
data.mode.fce        = vlit;                 
data.mode.hyb        = vcalc;                 
data.mode.idn        = vcalc;                 
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

% mode de variables neo
data.mode.eta       = vcalc;    
data.mode.jboot     = vcalc;    

% autres
data.mode.plot          = voff;  
data.mode.zeff          = vcalc;    
data.mode.ae            = vcalc;    
data.mode.asser         = voff;    

%mode de variables neo
data.mode.eta       = vcalc;    
data.mode.jboot     = vcalc;    

% autres
data.mode.plot          = voff;  

data.mode.zeff          = vcalc;    
data.mode.ae            = vcalc;    
data.mode.asser         = voff;    



% modification specifique du 0d
param.gene.k = param.gene.kmax;
data.gene.ip  = data.cons.ip;


switch option.mapequi
case 'off'
      % rien de plus
otherwise
	% calcul de la structure gene
        fprintf('post process :')
	for k=1:length(data.gene.temps)
                % extraction des donnees a un temps
		datak = zget1t(data,k);
		if k == 1
			datak_old = datak;
			datak.gene.dwdiadt = 0;
			datak.gene.dwbpdt  = 0;
                        param.gene.dt = mean(diff(data.gene.temps));
		else
                	param.gene.dt = datak.gene.temps - datak_old.gene.temps;
		end
                % filtrage deriv?e dans equi
                noms = fieldnames(datak.equi);
		for l = 1:length(noms)
  			nomc = noms{l};
                        if ~isempty(strmatch(strcat(nomc,'d1'),noms)) && ~isempty(strmatch(strcat(nomc,'d2'),noms))
				%[void,datak.equi.(strcat(nomc,'d1')),datak.equi.(strcat(nomc,'d2'))]= ...
                                %interpos(param.gene.x(1:5:end),datak.equi.(nomc)(1:5:end),param.gene.x,-1,[1 1],[0 1e32]);
				[datak.equi.(strcat(nomc,'d1')),datak.equi.(strcat(nomc,'d2'))]= ...
                                z0polyderive(param.gene.x(1:5:end),datak.equi.(nomc)(1:5:end),5,param.gene.x);
			end
		end
 
 		% calcul des divers profiles utiles (Te,Ti,  ...)
		datak = zprofdivers(param.phys,param.gene,datak,param.compo);
                % filtrage deriv?e dans prof
                noms = fieldnames(datak.prof);
		for l = 1:length(noms)
  			nomc = noms{l};
			if ~isempty(strmatch(strcat(nomc,'d1'),noms)) && ~isempty(strmatch(strcat(nomc,'d2'),noms))
				%[void,datak.prof.(strcat(nomc,'d1')),datak.prof.(strcat(nomc,'d2'))]= ...
				%		interpos(param.gene.x(1:5:end),datak.prof.(nomc)(1:5:end),param.gene.x,-1,[1 1],[0 1e32]);
				[datak.prof.(strcat(nomc,'d1')),datak.prof.(strcat(nomc,'d2'))]= ...
						z0polyderive(param.gene.x(1:5:end),datak.prof.(nomc)(1:5:end),5,param.gene.x);
			end
			if  ~isempty(strmatch(strcat('d',nomc,'dt'),noms)) 
				if (k > 1)
					datak.prof.(strcat('d',nomc,'dt')) = (datak.prof.(nomc) - datak_old.prof.(nomc)) ./ param.gene.dt;
				else
					datak.prof.(strcat('d',nomc,'dt'))(:) =  0;                   
				end
				if  ~isempty(strmatch(strcat('d',nomc,'dt3p'),noms))  
					if (k > 2) 
						datak.prof.(strcat('d',nomc,'dt3p')) = (datak.prof.(strcat('d',nomc,'dt')) + datak_old.prof.(strcat('d',nomc,'dt'))) ./ 2;         
					else
						datak.prof.(strcat('d',nomc,'dt3p'))(:) = 0;
					end
				end
			end
		end
		% calcul des flux
		datak = zflux(param.phys,param.gene,param.cons.neomulti,datak);
		% calcul des donnees generiques du plasma
		datak = zcalcdivers(param.phys,param.gene,param.compo,datak,datak_old.gene,param.gene.dt);
		% recolle
		data  = zput1t(data,k,datak);
		datak_old = datak;
		fprintf('.');
	end
	fprintf('\n');
end
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
else
   zuisavenonok;
end



function  s  = zpmean(t,p,e)

indok = find(isfinite(p) & isfinite(e));
if isempty(indok)
   s =1;
else
   s = trapz(t(indok),p(indok) .* e(indok)) ./ trapz(t(indok),eps + p(indok));
end


function out  = z0dsample(val,tin,xin,tout,xout,optvoid)

if size(xin,1) > 1 & size(xin,2) == 1
	inter   = interp1(xin',val',xout','linear')';
	out     = interp1(tin,inter,tout,'linear');
elseif size(xin,1) > 1 & size(xin,2) > 1
	%if size(xout,1) == 1
	%	xout = ones(size(val,1),1) * xout;
	%end
	%inter   = tsplinet(xin,val,xout);
	%out     = interp1(tin,inter,tout,'linear');
	if size(xin,1) == 1
		xin = ones(size(tin,1),1) * xin;
	end
	if size(tin,2) == 1
		tin = tin * ones(1,size(xin,2));
	end
	if size(xout,1) == 1
		xout = ones(size(tout,1),1) * xout;
	end
	if size(tout,2) == 1
		tout = tout * ones(1,size(xout,2));
	end
        % normalisation tin et tout (decorrelation espace et temps)
        delta_t = max(sqrt(eps),max(max(tin(:)),max(tout(:))) - min(min(tin(:)),min(tout(:))));
        delta_x = max(max(xin(:)),max(xout(:))) - min(min(xin(:)),min(xout(:)));
        if (10 .* delta_x) >  delta_t
			fact_t =  10 .* delta_x ./ delta_t; 
                        tin  = fact_t .* tin;
                        tout = fact_t .* tout;
        end
        try
            out = griddata(tin,xin,val,tout,xout,'linear');
        catch
            out = griddata(tin,xin,val,tout,xout,'linear',{'Qt','Qbb','Qc','Qs'});
        end
        if any(~isfinite(out(:)))
            try
                out_alt = griddata(tin,xin,val,tout,xout,'nearest');
            catch
                out_alt = griddata(tin,xin,val,tout,xout,'nearest',{'Qt','Qbb','Qc','Qs'});
            end
            indbad = find(~isfinite(out));
            out(indbad) = out_alt(indbad);
        end
		
else 
	warning off
	inter   = pchip(xin,val,xout);
	warning on
	out     = interp1(tin,inter,tout,'linear','extrap');
end
%out     = interp1(tin,inter,tout,'linear');
