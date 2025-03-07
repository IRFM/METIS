% fonction d'apppel de helena depuis METIS
% post -> metis data
% tc -> temps (s)
% plotonoff -> flag pour le graphes
function [equi,straightline,upar_out,utheta_out] =z0dhelena(post,tc,plotonoff)

if nargin < 3
	plotonoff = 0;
end

% la composition du plasma
if isfield(post.z0dinput.option,'Sn_fraction') && (post.z0dinput.option.Sn_fraction > 0)
    error('NCLASS code interface has not been adapted for tin (Sn) in plasma composition (option.Sn_fraction should be 0)');
end


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


zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
t    = zs.temps;
profli = post.profil0d;


% petit calcul complementaire
d0   = geo.a .^ 2  ./ 2 ./  geo.R  .* (2.*( geo.K .^ 2 + 1)./ (3 .*  geo.K .^ 2 + 1) .* ...
	       ( (zs.w ./ max(1,zs.wth) ) .* zs.betap + zs.li ./ 2) +  ...
               0.5 .* ( geo.K .^ 2 - 1)./ (3 .*  geo.K .^ 2 + 1)) .* ( 1 -  geo.a ./  geo.R);


indcp = min(find(profli.temps >= tc));
indc = min(find(zs.temps >= profli.temps(indcp)));

init    = 0 ;
b       = 0;
df2     = 0 .* ones(1,101);
kkbig   = [];
xrsig   = [99,99];
x101    = linspace(0,1,101);
ptot    = pchip(profli.xli,profli.ptot(indcp,:),x101);
ee      =   1.602176462e-19;
ptot    = ptot ./ 1e3 ./ ee;  % Pa -> keV.*m^-3
ptot_raw = ptot;
%ptot(1) = 1.05 .* max(ptot);
jmoy    = profli.jli(indcp,:);
jmoy(1) = max(jmoy(1:2));
jmoy(jmoy<0) = max(jmoy/3);
jmoy    = pchip(profli.xli,jmoy,x101);
jmoy_raw    = pchip(profli.xli,profli.jli(indcp,:),x101);

psipin  = pchip(profli.xli,profli.psi(indcp,:),x101);
absdpsidx  = abs(pdederive(x101,psipin,0,2,2,1));
ind        = find(absdpsidx == 0);
if ~isempty(ind)
        mindpsidx      = min(absdpsidx(absdpsidx > 0));
	absdpsidx(ind) = mindpsidx;
end
psipin = cumtrapz(x101,absdpsidx);
% changement d'echelle (normalisation)
psipin   =  abs(2 .* pi .* (psipin -psipin(1))) ;

dpdpsi   = pdederive(psipin,ptot,0,2,2,1);
if any(dpdpsi >= 0)
	ind            = find(dpdpsi >= 0);
	indm           = find(dpdpsi <0);
	indin          = min(indm);
	dpdpsiok       = dpdpsi(indm);
	dpdpsimin      = dpdpsiok(min(find(dpdpsi(indm) == max(dpdpsi(indm)))));
	dpdpsi(ind)    = dpdpsimin .* ones(1,length(ind));
	if indin >1
		dpdpsi(1:indin) = dpdpsimin .* ones(1,indin);
	end
end
% petit ajout pour ameliorer la convergence
ptotnew  = cumtrapz(psipin,dpdpsi);
ptot     = ptotnew - ptotnew(end) + ptot(end) .* (ptot(end)>0);
% 2 - borne
dpdpsi    = pdederive(psipin,ptot,2,2,2,1);
mu        = psipin ./ psipin(end);
indmu       = min(find(mu >0.01));
if isempty(ind)
   error('Donnees non valide en entree de l''equilibre : Ptot ');
end
dpdpsimin = mean(-abs(dpdpsi(ind:end)));
%  probleme de signe ?
ind       = find((dpdpsi < dpdpsimin) & (mu < mu(indmu)));
if ~isempty(ind)
     % correction effet d'arrondi d'indice pour DDS
     if length(ind)  >3
            ind(end) = [];
     end
     dpdpsi(ind) = dpdpsimin .* ones(1,length(ind));
     ptotnew  = cumtrapz(psipin,dpdpsi);
     ptot     = ptotnew - ptotnew(end) + ptot(end) .* (ptot(end)>0);
end


% separatrice
if isfield(post.z0dinput.exp0d,'Rsepa')
	Rl    = post.z0dinput.exp0d.Rsepa(indc,:);
	Zl    = post.z0dinput.exp0d.Zsepa(indc,:);
   	ind1 = min(find(Rl == min(Rl)));
   	ind2 = min(find(Rl == max(Rl)));

   	r0_  = 0.5 .* (Rl(ind1) +Rl(ind2));
   	z0_  = Zl(ind2);
   	aa    = atan2(Zl -z0_,Rl - r0_);
   	aa    = aa .* (aa >=0) + (aa + 2*pi) .* (aa< 0);
   	[aa,ind_aa] = sort(aa);
	indnok      = find(diff(aa) <= 0);
	if ~isempty(indnok)
		ind_aa(indnok) =[];
	end
	Rl = Rl(ind_aa);
   	% modification du 06/12/2004
   	Zl  = Zl(ind_aa) - z0_;
	if (Rl(end) ~= Rl(1)) |(Zl(end) ~= Zl(1))
   		Rl(end+1) = Rl(1);
   		Zl(end+1) = Zl(1);
	end
	% calcul de b0 pour helena
	b0 = geo.R(indc) .* geo.b0(indc) ./ r0_; 
        r0 = r0_;

else
	td  =  asin(geo.d(indc));
	th   = linspace(0,2.*pi,201);
	Rl    = geo.R(indc) + geo.a(indc) .* cos(th + (td  .* sin(th)));
	Zl    = geo.K(indc) .* geo.a(indc) .* sin(th);
        b0 = geo.b0(indc);
        r0 = geo.R(indc);
        z0_ = geo.z0(indc);
end
geo.mode = 2;


% petite figure pour comparer les entrees de helena et les donnees de metis
if plotonoff > 0
	h = findobj(0,'type','figure','tag','z0dhelena_1');
	if isempty(h)
	h=figure('tag','z0dhelena_1');
	else
	figure(h);
	end   
	clf
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])
	
	subplot(2,2,1)
	plot(profli.xli,profli.jli(indcp,:),'b',x101,jmoy_raw,'-.r',x101,jmoy,'k--');
	xlabel('r/a')
	ylabel('J (A/m^2)');
	legend('METIS','unmodified','modified');
	
	subplot(2,2,2)
	plot(profli.xli,-(profli.psi(indcp,:) - profli.psi(indcp,1)) .* 2 .* pi,'b',x101,psipin,'-.r');
	xlabel('r/a')
	ylabel('Psi (Wb)');
	legend('METIS','HELENA input');
	
	subplot(2,2,3)
	plot(profli.xli,profli.ptot(indcp,:),x101,ptot_raw .* 1e3 .* ee,'-.r',x101,ptot .* 1e3 .* ee,'k--');
	xlabel('r/a')
	ylabel('Ptot (Pa)');
	legend('METIS','unmodified','modified');

	edition2
end

% 1er appel a helena avec Jmoy non retoucher
init(2) = 1e-4;
init(3) = 1e-7;
[C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor, ...
 	drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
 	r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc, ...
 	drmerc,balcrit,frmode,ifail,nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
 	helmex77(psipin,ptot_raw,jmoy_raw,r0,b0,-zs.ip(indc),2,geo.a(indc),geo.K(indc), ...
 	geo.d(indc),geo.d(indc),Rl,Zl,init,b,df2,kkbig,0,xrsig);

if ifail  ~= 0

	fprintf('Z0DHELENA : HELENA do not converge with METIS current and pressure profiles; restart with modified profiles\n');
	init(2) = 0.1;
	init(3) = 1e-4;
	[C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor, ...
		drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
		r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc, ...
		drmerc,balcrit,frmode,ifail,nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
		helmex77(psipin,ptot,jmoy,r0,b0,zs.ip(indc),2,geo.a(indc),geo.K(indc), ...
		geo.d(indc),geo.d(indc),Rl,Zl,init,b,df2,kkbig,2,xrsig);
        conv = 1;
	amix = 1;
else
	% pour la suite 
	ptot = ptot_raw;
	jmoy = jmoy_raw;
        conv = 2;
	amix = 0.9;
end


if ifail  == 0
	equi.rhomax         = rho(end);
        equi.rho            = rho';
	Vpb         =  Vp;
	Vpb(1)      =  eps;
	grho2r2     =  C2 ./ Vpb;
	grho2r2(1)  = 0;
	grho2r2     =  zcentre(grho2r2);   % probleme
	r2i         =  C3 ./ Vpb;
	r2i         =  zcentre(r2i);
	oor(1)=0;
	oor         =  zcentre(oor);
	drhoor(1)=0;
	drhoor      =  zcentre(drhoor);
	drhoav(1)=0;
	drhoav      =  zcentre(drhoav);
	rav         =  zcentre(rav);
	r2          =  zcentre(r2);
	r2tau2      =  zcentre(r2tau2);
	drho2ob2    =  zcentre(drho2ob2);
	r3otau3     =  zcentre(r3otau3);
	r3otau      =  zcentre(r3otau);
	
	equi.q      = qpsi';
        equi.shear  = rho' ./ equi.q .* pdederive(rho',equi.q,0,2,2,1);
        equi.psi    = - psipout' ./ 2 ./ pi + profli.psi(indcp,1);
        equi.phi    = psiT';
        equi.ftrap  = ftra';
        equi.ftrap(1)       = 0;
        equi.F      = Fdia';
        equi.b2i    = invB2';
        equi.b2     = B2';
        equi.vpr    = Vp';
        equi.spr    = Sp';
        equi.d      = xshift';
	% compatibilite ascendante
	if isfield(equi,'ind')
		equi.ind            = zeros(size(rho'));	
	else
		equi.indh           = NaN .* rho';
		equi.indl           = NaN .* rho';
	end
	equi.trh  = xtriapos';
        equi.trl  = xtrianeg';
  	equi.e    = xell';
  	equi.a    = xrad';
  	equi.grho2r2  = grho2r2';
  	equi.grhor    = drhoor';
  	equi.ri       = oor';
  	equi.grho2    = rho2m';
	% correction au centre de grho2 (23/07/2001)
  	equi.grho2    = zregular(equi.rho./max(equi.rho),equi.grho2);
  	equi.c2c      = C2';
  	equi.grho     = drhoav';
  	equi.rmoy     = rav';
  	equi.jmoy     = jout';
  	equi.ptot     = pout' .* 1e3 .* phys.e;
  	equi.r2       = r2';
  	equi.r2i      = r2i';
  	equi.r2tau2   = r2tau2';
  	equi.grho2b2  = drho2ob2';
  	equi.r3tau3   = r3otau3';
  	equi.r3tau    = r3otau';
	% fin du code modifie

	% ip et li
  	equi.ip = ipout .* max(xrad) .* b0 ./ phys.mu0;
  	equi.li = li;
	equi.c2c       = equi.vpr .* equi.grho2r2;
  	equi.c2c(1)        = eps;
  	equi.R     = single(shiftdim(R,-1));
	% z est compatible avec la separatrice
	% correction du 06/12/2004
	if geo.mode == 2
	    equi.Z     = single(shiftdim(Z,-1) + z0_);
	else
	  equi.Z     = single(shiftdim(Z,-1) + z0);
	end
	equi.rhoRZ = single(shiftdim(rho,-1));
	equi.psiRZ = single(shiftdim(-psipout./ 2 ./ pi + equi.psi(1),-1));
	equi.df2RZ = single(shiftdim(df2,-1));
	equi.dprRZ = single(shiftdim(dpr,-1));
	equi.frmode= single(shiftdim(frmode,-1));

	% La carte de champ
	BPHI       = (Fdia * ones(1,size(R,2))) ./ R;
	equi.BR    = single(shiftdim(BPR,-1));
	equi.BZ    = single(shiftdim(BPZ,-1));
	equi.BPHI  = single(shiftdim(BPHI,-1));


	% calcul de rhog
	% ajout du 10/07/2001
	equi.rhog = equi.a ./ equi.a(end);
	% ajout du 12/07/2001
	raxe      = (min(R') + max(R')) ./ 2;
	equi.raxe = raxe;
	equi.zaxe = NaN .* equi.raxe;

	% donnees MHD
	equi.mhd.ballooning     = (1- balcrit') .* (balcrit' ~= 0);
	equi.mhd.mercier	= - drmerc';
	equi.mhd.ideal  	= 0.25 - dimerc';

	% geodesique (nchi,cpsurf,radius,gem11,gem12,gem33)
	straightline.nchi   = nchi;
	straightline.cpsurf = cpsurf;
	straightline.radius = radius;
	straightline.psi    = psipout;
	straightline.rho    = rho;
	straightline.q      = qpsi;
	straightline.Fdia   = Fdia;
	
	straightline.gem11  = reshape(gem11,nchi,length(psipout))';
	straightline.gem11(:,end+1) = straightline.gem11(:,1);
	
	straightline.gem12  = reshape(gem12,nchi,length(psipout))';
	straightline.gem12(:,end+1) = straightline.gem12(:,1);
	
	straightline.gem33  = reshape(gem33,nchi,length(psipout))';
	straightline.gem33(:,end+1) = straightline.gem33(:,1);
	
	straightline.rspot = reshape(rspot,nchi,length(psipout))';
	straightline.rspot(:,end+1) = straightline.rspot(:,1);


	if geo.mode == 2
		straightline.zspot = reshape(zspot,nchi,length(psipout))' + z0_;
	else
		straightline.zspot = reshape(zspot,nchi,length(psipout))' + z0;
	end
	straightline.zspot(:,end+1) = straightline.zspot(:,1);

	straightline.brspot = reshape(brspot,nchi,length(psipout))';
	straightline.brspot(:,end+1) = straightline.brspot(:,1);

	straightline.bzspot = reshape(bzspot,nchi,length(psipout))';
	straightline.bzspot(:,end+1) = straightline.bzspot(:,1);


	% si le calcul de deltaprime est demande
	m_deltap = NaN .* ones(1,length(equi.rho));
	deltap = NaN .* ones(1,length(equi.rho));
	% ecriture dans la structure equi.mhd
	equi.mhd.deltap   = deltap;
	equi.mhd.m_deltap = m_deltap;


        % donnees complementaires si disponible
	equi.drhomaxdt      = NaN;
	equi.dvprdt         = NaN .*  equi.vpr;
	equi.dsprdt         = NaN .*  equi.spr;
	equi.dphidt         = NaN .*  equi.phi;
	equi.phid1          = pdederive(equi.rho,equi.phi,0,2,2,1) ./ max(equi.rho);
	equi.phid2          = pdederive(equi.rho,equi.phi,0,2,2,2) ./ max(equi.rho);
	equi.psi0           = equi.psi(1);
        equi.betap          = NaN;
        equi.betat          = NaN;
        equi.betan           = NaN;
        equi.q0              = equi.q(1);
        equi.volume          = cumtrapz(equi.rho,equi.vpr);
        equi.surface         = cumtrapz(equi.rho,equi.spr);
        equi.bnorme          = b;
        equi.conv            = conv;
        equi.oscil           = 0;
        equi.fail	     = ifail;
        equi.errcur          = init(2);
        equi.errit           = init(3);
        equi.amix            = amix;
        equi.mhd.gamma1      = NaN;
        equi.mhd.gamma2      = NaN;
        equi.mhd.gamma3      = NaN;
        equi.mhd.vper1       = NaN .*  equi.vpr;
        equi.mhd.vper2       = NaN .*  equi.vpr;
        equi.mhd.vper3       = NaN .*  equi.vpr;
        equi.mhd.fail        = NaN;
        equi.mhd.conv1       = NaN;
        equi.mhd.conv2       = NaN;
        equi.mhd.conv3       = NaN;
        equi.mhd.erreur1     = NaN;
        equi.mhd.erreur2     = NaN;
        equi.mhd.erreur3     = NaN;
        equi.free.pfcur      = NaN .*  equi.vpr;
        equi.free.vext       = NaN .*  equi.vpr;
        equi.free.psi_plasma = NaN;
        equi.free.psi_ext    = NaN;
        equi.free.L_plasma   = NaN;
        equi.free.L_int      = NaN;
        equi.free.L_ext      = NaN;
        equi.free.probe.br   = NaN .* ones(1,1024);
        equi.free.probe.bz   = NaN .* ones(1,1024);
        equi.free.probe.psi  = NaN .* ones(1,1024);
        equi.free.probe.gaps = NaN .*  equi.vpr;
        equi.free.strength.bmax = NaN .*  equi.vpr;
        equi.free.strength.fr = NaN .*  equi.vpr;
        equi.free.strength.fz = NaN .*  equi.vpr;
        equi.free.strength.flux = NaN .*  equi.vpr;



        % si les graphes sont demandes
        if plotonoff > 0

		psiloc   =profli.psi(indcp,:);
		psiloc   =  abs(2 .* pi .* (psiloc -psiloc(1))) ;
		ipout    = ipout .* max(xrad) .* geo.b0(indc) ./ (4.*pi.*1e-7);
	
		volume   = trapz(rho',Vp',2);
		xout     = xrad ./ max(eps,max(xrad));
	
		% le plot
		hz =findobj(0,'type','figure','tag','equi0d1');
		if isempty(hz)
			hz=figure('tag','equi0d1','name','Metis : Helena run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
		subplot(2,2,1)
		plot(profli.xli,psiloc,'om',xout,psipout,'r',profli.xli,profli.phi(indcp,:),'oc',xout,psiT,'b');
		title('psi & phi')	
                xlabel('r/a')
	        ylabel('Psi (Wb)');
	        legend('METIS Psi','HELENA output Psi','METIS Phi','HELENA output Phi');

		subplot(2,2,2)
		plot(profli.xli,geo.d(indc) .* profli.xli .^ 3 ,'oc',xout,xtriapos,'r',xout,xtrianeg,'m');
		title('triapos & trianeg')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
		subplot(2,2,3)
		plot(profli.xli,profli.Raxe(indcp,:) - profli.Raxe(indcp,end),'oc',xout,xshift,'r',0,d0(indc),'*g');
		title('shift')
                xlabel('r/a')
	        ylabel('m');
	        legend('METIS','HELENA output');
		subplot(2,2,4)
		plot(profli.xli,profli.kx(indcp,:),'oc',xout,xell,'r');
		title('K(x)')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
	
		hz =findobj(0,'type','figure','tag','equi0d2');
		if isempty(hz)
			hz=figure('tag','equi0d2','name','Metis : Helena run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
		subplot(2,2,1)
		plot(x101,ptot .* 1e3 .* ee,'b',profli.xli,profli.ptot(indcp,:),'oc',xout,pout .* 1e3 .* ee,'r');
		title('Ptot')
                xlabel('r/a')
	        ylabel('Pa');
	        legend('Helena input','METIS','HELENA output');
		subplot(2,2,2)
		plot(x101,jmoy,'b',profli.xli,profli.jli(indcp,:),'oc',xout,jout,'r');
		title('Jmoy')
                xlabel('r/a')
	        ylabel('A/m^2');
	        legend('Helena input','METIS','HELENA output');
		subplot(2,2,3) 
	%  	pp = polyfit(psiloc,profli.phi(indcp,:),7);
	%  	ppd = pp(1:end-1) .* (5:-1:1);
	%  	qpp = polyval(ppd,psiloc)
	%  	plot(profli.xli,profli.qjli(indcp,:),'oc',xout,qpsi,'r',1,zs.qeff(indc),'*g', ...
	%       xout,pdederive(psipout,psiT,0,2,1,1),'+',profli.xli,qpp,'x');
		plot(profli.xli,profli.qjli(indcp,:),'oc',xout,qpsi,'r',1,zs.qeff(indc),'*g');
		title('q')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
		subplot(2,2,4)
		plot(profli.xli,profli.fdia(indcp,:),'oc',xout,Fdia,'r',1,r0.*b0,'g*');
		title('Fdia')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
	
		hz =findobj(0,'type','figure','tag','equi0d3');
		if isempty(hz)
			hz=figure('tag','equi0d3','name','Metis : Helena run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
		subplot(3,2,1)
		plot(profli.xli,profli.vpr_tor(indcp,:),'oc',xout,Vp,'r');
		title('Vpr')
                xlabel('r/a')
	        ylabel('m^2');
	        legend('METIS','HELENA output');
		subplot(3,2,2)
		plot(profli.xli,profli.C2(indcp,:),'om',xout,C2 ,'r');
		hold on
		plot(profli.xli,profli.C3(indcp,:) ,'oc',xout,C3,'b');
		hold off
		title('C2 & C3')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
		subplot(3,2,3)
		plot(profli.xli,profli.grho2r2(indcp,:),'om',xout,C2./max(eps,Vp),'r');
		hold on
		plot(profli.xli,profli.r2i(indcp,:),'oc',xout,C3./max(eps,Vp),'b');
		plot(profli.xli,profli.ri(indcp,:),'og',xout,oor,'g');
		hold off
		title('grho2r2, r2i & ri')
                xlabel('r/a')
	        ylabel('m^{-2}');
	        legend('METIS','HELENA output');
		subplot(3,2,4)
		plot(profli.xli,profli.ftrap(indcp,:),'oc',xout,ftra,'r');
		title('ftrap')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
	% equi.grho2          = interp1(rho,rho2m,rhoin,'spline')';
	% equi.grho           = interp1(rho,drhoav,rhoin,'spline')';
		subplot(3,2,5)
		plot(profli.xli,profli.grho(indcp,:),'om',xout,drhoav,'r');
		title('grho')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
		subplot(3,2,6)
		plot(profli.xli,profli.grho2(indcp,:),'oc',xout,rho2m,'b');
		title('grho2')
                xlabel('r/a')
	        ylabel('');
	        legend('METIS','HELENA output');
	
		hz =findobj(0,'type','figure','tag','equi0d4');
		if isempty(hz)
			hz=figure('tag','equi0d4','name','Metis : Helena run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
		plot(R',Z','b',Rl,Zl,'r',geo.R(indc),0,'+r')
		title('Flux surface')
		xlabel('R (m)')
		ylabel('Z (m)')
		set(gca,'xlim',[0,max(Rl)]);
		axis('equal')
	
	
		% appel de nclass
		clear p1 p2 p3 p4 p5 p6 p7 p8 p9
		info = zneo;
		parametre = info.valeur;

                if post.z0dinput.option.W_effect

			% charge
			zW = trapz(profli.xli,z0wavez(profli.tep(indcp,:)) .* profli.vpr(indcp,:),2)'  ./ max(eps,trapz(profli.xli,profli.vpr(indcp,:),2))';
			% e H D T He Imp
			p1    = [-1;1;1;1;2;op0d.zimp;op0d.zmax;zW];
			% nombre de masse
			p2    = [1;1;2;3;4;round(7/3 .* op0d.zimp);round(7/3 .* op0d.zmax);183.84];
			% champ electrique radial et sa derivees ou indices du gaz injecte
			p3=[0;0];
			% parametres de la fonction
			p9  = [parametre.bootmodel;parametre.nbeq;parametre.banane;1;1;length(profli.xli)];
		
			% correction epar <> 0
			epar    = profli.epar(indcp,:);
			ind     = find(~isfinite(epar));
			if ~isempty(ind)
			epar(ind) = 0;
			end
			eparnz  = 0.01 ./ 2 ./ pi ./ geo.R(indc);
			epar    = eparnz .* (epar == 0) + epar .* (epar ~= 0);
		
			% temperatures
			p4  = [profli.tep(indcp,:)' profli.tip(indcp,:)' * ones(1,7)]  .* 1e-3;
		
			% densites
			nDp = max(1e13,zs.nDm(indc) ./ zs.n1m(indc) .* profli.n1p(indcp,:));
			nTp = max(1e13,zs.nTm(indc) ./ zs.n1m(indc) .* profli.n1p(indcp,:));
			nHp = max(1e13,profli.n1p(indcp,:) - nDp - nTp);
			p5  = [max(1e13,profli.nep(indcp,:))'  nHp' nDp' nTp' max(1e13,profli.nhep(indcp,:))' max(1e11,profli.nzp(indcp,:))' op0d.rimp .* max(1e11,profli.nzp(indcp,:))',max(1e9,profli.nwp(indcp,:))'];


		else

			% charge
			% e H D T He Imp
			p1    = [-1;1;1;1;2;op0d.zimp;op0d.zmax];
			% nombre de masse
			p2    = [1;1;2;3;4;round(7/3 .* op0d.zimp);round(7/3 .* op0d.zmax)];
			% champ electrique radial et sa derivees ou indices du gaz injecte
			p3=[0;0];
			% parametres de la fonction
			p9  = [parametre.bootmodel;parametre.nbeq;parametre.banane;1;1;length(profli.xli)];
		
			% correction epar <> 0
			epar    = profli.epar(indcp,:);
			ind     = find(~isfinite(epar));
			if ~isempty(ind)
			epar(ind) = 0;
			end
			eparnz  = 0.01 ./ 2 ./ pi ./ geo.R(indc);
			epar    = eparnz .* (epar == 0) + epar .* (epar ~= 0);
		
			% temperatures
			p4  = [profli.tep(indcp,:)' profli.tip(indcp,:)' * ones(1,6)]  .* 1e-3;
		
			% densites
			nDp = max(1e13,zs.nDm(indc) ./ zs.n1m(indc) .* profli.n1p(indcp,:));
			nTp = max(1e13,zs.nTm(indc) ./ zs.n1m(indc) .* profli.n1p(indcp,:));
			nHp = max(1e13,profli.n1p(indcp,:) - nDp - nTp);
			p5  = [max(1e13,profli.nep(indcp,:))'  nHp' nDp' nTp' max(1e13,profli.nhep(indcp,:))' max(1e11,profli.nzp(indcp,:))' op0d.rimp .* max(1e11,profli.nzp(indcp,:))'];

	        end
		% derivees des temperatures
		p6  = pdederive(profli.xli',p4,0,2,1,1);
		p7  = pdederive(profli.xli',p5,0,2,1,1);
		%dx  = mean(diff(profli.xli));
		%p6(end-1,:)  = (p4(end,:) - p4(end-1,:)) ./ dx;
		%p7(end-1,:)  = (p5(end,:) - p5(end-1,:)) ./ dx;
	
		% calcul de b2
		b2      = pchip(xout,B2,profli.xli);
		b2i     = pchip(xout,invB2,profli.xli);
		drho2ob2(1) = drho2ob2(2);
		grho2b2 = pchip(xout,drho2ob2,profli.xli);
		ftrap = pchip(xout,ftra,profli.xli);
		fdia  = pchip(xout,Fdia,profli.xli);
		oor(1) = oor(2);
		ri  = pchip(xout,oor,profli.xli);
		oor2(1) = oor2(2);
		ri2  = pchip(xout,oor2,profli.xli);
		qjli  = pchip(xout,qpsi,profli.xli);
		Rext  = pchip(xout,max(R,[],2),profli.xli);
		Raxe  = pchip(xout,max(R,[],2) + min(R,[],2),profli.xli) ./ 2;
		vpr  = pchip(xout,Vp,profli.xli);
		spr  = pchip(xout,Sp,profli.xli);
		% structure complexe
		ven =ones(size(profli.xli));
		p8  = [Rext;     ...
		profli.xli;                               ...
		qjli;                          ...
		ftrap;                      ...
		fdia  .* ri;                      ...
		vpr  .* max(rho);        ...
		max(rho)*ven;  ...
		profli.Raxe(indcp,end) *ven;       ...
		b2 ;                        ...
		b2i ;                       ...
		qjli(1)*ven;    ...
		ri2 .* vpr.* max(rho);        ...
		grho2b2;                    ...
		xell(1) .* ven;    ...
		epar .* geo.b0(indc)];
		%
		% dimension interne profil -> 200
		%
		ndim          = 200;
		p4(ndim,6)    = 0;
		p5(ndim,6)    = 0;
		p6(ndim,6)    = 0;
		p7(ndim,6)    = 0;
		p8(15,ndim)   = 0;
		p24  =p8(15,:)';
		rpts                  = length(profli.xli);
		p25 = zeros(size(p24));
		p26 = zeros(size(p24));
		p27 = zeros(size(p24));
		p28 = zeros(size(p24));
		p29 = zeros(size(p24));
		pold        = [p25 p26 p27 p28 p29]';	
		tab         = bootHoulbergMat6(p1,p2,p3,p4,p5,p6,p7,p8,p9,pold);
%  		% nouvel appel pour nclass_1
%  		%
%  		% decopue de p8 en 15 composants
%  		%       
%  		p10  =p8(1,:)';
%  		p11  =p8(2,:)';
%  		p12  =p8(3,:)';
%  		p13  =p8(4,:)';
%  		p14  =p8(5,:)';
%  		p15  =p8(6,:)';
%  		p16  =p8(7,:)';
%  		p17  =p8(8,:)';
%  		p18  =p8(9,:)';
%  		p19  =p8(10,:)';
%  		p20  =p8(11,:)';
%  		p21  =p8(12,:)';
%  		p22  =p8(13,:)';
%  		p23  =p8(14,:)';
%  		p24  =p8(15,:)';
%  
%  		memoire.data  =[];
%                  [tab,memoire]    = neocall(p1,p2,p3,p4,p5,p6,p7,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,memoire);
%  
		neojboot             = tab(1,1:rpts) ./ geo.b0(indc);    % valide
		neoeta               = tab(9,1:rpts);              % valide
	
		% separation
		nspec      = prod(size(p1));
		indmem     = 26;
		upar_out   = tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons
                indmem     = indmem+nspec;
		utheta_out = - post.z0dinput.option.signe .* tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons

	
		% valeur de sauter locale
		dpsidx = -pchip(xout,pdederive(xout',psipout'./2./pi ,0,2,2,1),profli.xli);
		Raxe   = (min(R') + max(R')) ./ 2;
		Raxe   = pchip(xout,Raxe,profli.xli);
		epsi   = pchip(xout,xrad,profli.xli)./Raxe;
		rmx    = pchip(xout,rho,profli.xli);
		r2i    = pchip(xout,r2i,profli.xli);

		jboot2  = zsauter0d(profli.xli,profli.tep(indcp,:),profli.tip(indcp,:),profli.nep(indcp,:),profli.nip(indcp,:), ...
				qjli,profli.zeff(indcp,:),(op0d.gaz < 4) +1,Raxe,ftrap,epsi,dpsidx,fdia,1)./ b0;
	
		%jboot3  = zsauter0d_Koh(profli.xli,profli.tep(indcp,:),profli.tip(indcp,:),profli.nep(indcp,:),profli.nip(indcp,:), ...
		%		        qjli,profli.zeff(indcp,:),(op0d.gaz < 4) +1,Raxe,ftrap,epsi,dpsidx,fdia,1)./ b0;
				        
		jboot4  = zsauter0d_HC(profli.xli,profli.tep(indcp,:),profli.tip(indcp,:),profli.nep(indcp,:),profli.nip(indcp,:), ...
				        qjli,profli.zeff(indcp,:),(op0d.gaz < 4) +1,Raxe,ftrap,epsi,dpsidx,fdia,1,1,1,rmx,r2i,grho2r2,0)./ b0;
	
		iboot_nclass_helena = equi.rhomax .* trapz(profli.xli,neojboot .* spr,2);
		iboot_sauter_helena = equi.rhomax .* trapz(profli.xli,jboot2 .* spr,2);
		%iboot_sauter_koh    = equi.rhomax .* trapz(profli.xli,jboot3 .* spr,2);
		iboot_sauter_HC     = equi.rhomax .* trapz(profli.xli,jboot4 .* spr,2);
		iboot_sauter_metis  = zs.iboot(indc);

		hz =findobj(0,'type','figure','tag','neo0d1');
		if isempty(hz)
			hz=figure('tag','neo0d1','name','Metis : NClass run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
	
		subplot(2,1,1)
		semilogy(profli.xli,profli.eta(indcp,:),'or',profli.xli,neoeta,'b');
		title('eta')
                xlabel('r/a')
	        ylabel('ohm * m');
	        legend('METIS','NCLASS (with HELENA equi.)');
	
		subplot(2,1,2)
		plot(profli.xli,profli.jboot(indcp,:),'or',profli.xli,neojboot,'b',profli.xli,jboot2,'k',profli.xli,jboot4,'g');
		title('jboot');
                xlabel('r/a')
	        ylabel('A/m^2');
	        legend(sprintf('METIS (Sauter with METIS equi @ %g A)',iboot_sauter_metis), ...
                       sprintf('NCLASS (with HELENA equi. @ %g A)',iboot_nclass_helena), ...
                      sprintf('Sauter (with HELENA equi. @ %g A)',iboot_sauter_helena), ...
                      sprintf('Sauter modified H-C (with HELENA equi. @ %g A)',iboot_sauter_HC));


		% vitesse poloidale
		hz =findobj(0,'type','figure','tag','neo0d2');
		if isempty(hz)
			hz=figure('tag','neo0d2','name','Metis : NClass run');
		else
			figure(hz);
		end
		clf
		set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
	
		plot(profli.xli,profli.utheta(indcp,:),'ok',profli.xli,utheta_out(2:end,:));
		title('Utheta')
                xlabel('r/a')
	        ylabel('m/s/T');
		if post.z0dinput.option.W_effect
	        	legend('METIS','H','D','T','He','Zimp','Zmax','W');
		else
	        	legend('METIS','H','D','T','He','Zimp','Zmax');
		end
	end

	

else
	error('Helena : Non convergence');
	figure
	subplot(2,2,1)
	plot(psipin,ptot)
	subplot(2,2,2)
	plot(psipin,jmoy)
	subplot(2,2,3)
	plot(Rl,Zl)
end

% integralle volumique 
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax   
function s=zintvol(e,x,vpr,rhomax)   

  s = rhomax.*trapz(x,vpr .* e,2);
  


% interpolation dans helena
function yout = helena_interp(xin,yin,xout)

yout = pchip(xin,yin,xout);

