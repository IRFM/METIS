% calcul le courant d'electron decouple
% modele these  F.Sourd
function [irun,vl1] = z0irun(cons,geo,zs,profil,option)


phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)



% correction des pieges
if isfield(profil,'ftrap')
	ntrap = trapz(profil.xli,profil.nep .* profil.ftrap .* profil.vpr,2) ./  ...
	        max(1,trapz(profil.xli,profil.nep .* profil.vpr,2));
	% expression Rosenbluth : correction du champ critique
	corp = 1 + ntrap ./ zs.nem;
else
	corp = 2 .* ones(size(zs.nem));
	ntrap = 0.5 .* zs.nem;
end
nem    = zs.nem; 
tem    = zs.tem; 
zeff   = zs.zeff; 
vloop  = zs.vloop;

% utilisation des infomation de profile
if isfield(profil,'ftrap') && isfield(profil,'tep') && isfield(profil,'nep') && isfield(profil,'zeff') && isfield(profil,'epar') && (option.runaway >= 4)
      %jrun = zjrun(profil.epar,profil.nep,profil.tep,profil.zeff,profil.ftrap,zs.tauj,1);
      %jrun = max(jrun,eps);
      jrun = (1 - profil.ftrap) ./ profil.eta;
      jrun = jrun ./ (max(jrun,[],2) * ones(1,size(jrun,2)));
      %figure(21);plot(profil.xli,jrun);drawnow
      nem = trapz(profil.xli,profil.nep .* jrun .* profil.vpr,2) ./  ...
	       max(eps,trapz(profil.xli,jrun .* profil.vpr,2));
      if option.runaway == 5
	  n0m = trapz(profil.xli,profil.n0m .* jrun .* profil.vpr,2) ./  ...
		      max(eps,trapz(profil.xli,jrun .* profil.vpr,2));
      else
	  n0m   = zeros(size(nem));            
      end
      tem = trapz(profil.xli,profil.tep .* jrun .* profil.spr,2) ./  ...
	       max(eps,trapz(profil.xli,jrun .* profil.spr,2));
      zeff = trapz(profil.xli,profil.zeff .* jrun .* profil.spr,2) ./  ...
	       max(eps,trapz(profil.xli,jrun .* profil.spr,2));
      epar = trapz(profil.xli,profil.epar .* jrun .* profil.spr,2) ./  ...
	       max(eps,trapz(profil.xli,jrun .* profil.spr,2));
      ntrap = trapz(profil.xli,profil.nep .* profil.ftrap .* jrun .* profil.vpr,2) ./  ...
	       max(eps,trapz(profil.xli,jrun .* profil.vpr,2));
      corp = 1 + ntrap ./ nem;
      
      vloop = epar .* 2 .* pi .* geo.R;
      epar   = abs(epar);
else
      % pour une initialisation propre
      epar  = max(eps,abs(vloop) ./ 2 ./ pi ./ geo.R);
      n0m   = zeros(size(nem));       
end


% valeur de confinement avec limite + mecanisme pour disruption
fh      = max(eps,min(1,cons.hmore)); 
fh(cons.hmore > 0.1) = 1;

% ref R. Jaspers NF vol 33 1993 p 1775
% ref R. Jayakumar et al, Physics letters A, vol 172 ,1993 pp 447-451
% ref M.N. Rosenbluth NF vol 37 pp 1355- , 1997
% Tokamaks, Wesson p 72-74
% ref : H. Smith et al, Physics of Plasma 13 (2006) 102502
lng         = 14.9 - 0.5 .* log(nem ./ 1e20) + log(tem ./ 1e3);
vth         = min(phys.c,sqrt(tem .* phys.e ./ phys.me));
irun   = zs.irun;
%
switch option.runaway
case {4,5} 
        % deja calcule
case {2,3}
	if isfield(profil,'epar')
		epar 	    = max(profil.epar,[],2);
	else
		epar        = max(eps,abs(vloop) ./ 2 ./ pi ./ geo.R);
	end
	% anticipation de la variation de P_LH car il n'y a pas de temps de thermalisation sur les electrons rapides.
	epar  = cat(1,epar(2:end),epar(end));
otherwise
	epar        = max(eps,abs(vloop) ./ 2 ./ pi ./ geo.R);
end

%
% probleme de stabilite numerique
vcrit0      = min(phys.c,sqrt(nem .* phys.e .^ 3 .* lng ./ 4 ./ pi ./ phys.epsi0 .^ 2 ./ phys.me ./ max(1e-6,epar) .* (2 + zeff)));
% correction relativiste
vcrit  = min(phys.c,(vcrit0 ./ phys.c) .* sqrt(sqrt(vcrit0 .^ 4 + 4 .* phys.c .^ 4 ) - vcrit0 .^ 2)  ./ sqrt(2));
%
tau         = 4 .* pi .* phys.epsi0 .^ 2 .* phys.me .^ 2 .* vth .^ 3 ./ (nem+n0m) ./ phys.e .^ 4 ./ lng; 
% critical electric field
er          = corp .* nem .* phys.e .^ 3 .* lng ./ 4 ./ pi ./ phys.epsi0 .^ 2 ./ phys.me ./ phys.c .^ 2;
switch option.runaway
case {2,3,4}
    if option.lhmode < 5
	vnpar = phys.c ./ option.npar0;
	erlh  = corp .* nem .* phys.e .^ 3 .* lng ./ 4 ./ pi ./ phys.epsi0 .^ 2 ./ phys.me ./ vnpar .^ 2;
	flh   = tanh(zs.plh  ./ max(1,zs.pohm));
	%figure(21);clf;subplot(2,1,1);plot(cons.temps,flh);subplot(2,1,2);plot(cons.temps,er,'b',cons.temps,erlh,'r',cons.temps,er .* (1 - flh) + erlh .*flh,'g');drawnow
	er    = er .* (1 - flh) + erlh .*flh;
    end
end
% Dreicer electric field
ec          = corp .* nem .* phys.e .^ 3 .* lng ./ 4 ./ pi ./ phys.epsi0 .^ 2 ./ phys.me ./ vth .^ 2;

if (option.evolution == 1) && isfield(profil,'epar')
     % rien dans ce cas on n'est pas au debut de la simulation.
elseif option.breakdown < 0
	epar(1)     = max(eps,abs(option.breakdown) ./ 2 ./ pi ./ geo.R(1));
else
	epar(1)     = option.breakdown .* ec(1);
end
vl1         = epar(1) .*  2  .* pi .* geo.R;

% 1ere generation
kzeff       = pchip([1 2 10],[0.43 0.32 0.14],zeff);
se          = kzeff .* (epar ./ ec) .^( -3 .* (zeff +1) / 16 ) .* exp( - ec ./ epar ./ 4 - sqrt( (zeff + 1) .* ec ./epar));
source      = zs.sp .* sign(vloop) .* max(0,nem - ntrap) .* max(0,min(1,se ./ tau)) .* phys.e .* vcrit .* (epar > er);
% secondaire
azeff        = (2 + zeff) ./ 3;
tauc         = 4 .* pi .* phys.epsi0 .^ 2 .* phys.me .^ 2 .* phys.c .^ 3 ./ (nem + n0m)  ./ phys.e .^ 4 ./ lng;
%tauc2        = phys.me .* phys.c ./ (phys.e .* er);
ga           = 1 ./ max(eps,1 + 1.46 .* sqrt(zs.rm ./ geo.R) + 1.72 .* zs.rm ./ geo.R);
djdtsj       = sqrt(pi .* ga ./ 3 ./ (5 + zeff)) ./ tauc ./ lng .* max(0,epar ./ er  - 1) ./  ...
		sqrt(1 - er ./ max(eps,epar) + 4 .* pi .* (zs.zeff + 1) .^ 2./ 3 ./ ga ./ (zeff + 5 ) ./ ...
		(epar .^ 2 ./ er .^ 2 + 4 ./ ga .^ 2 - 1)); 	
	
% solution
% temps de collision (valeur minimale de tau)
taus    = (4/3) .* pi .* phys.epsi0 .^ 2 .* phys.me .^ 2 .*  vcrit .^ 3 ./ (nem + n0m)  ./ phys.e .^ 4 ./ lng;
% temps d'etablissement d'un runaway
taur    =  3 .* lng .* phys.me .* phys.c ./ (phys.e .* epar);
% limite de la source
source  = sign(source) .* min(abs(source),zs.ip ./ taur);

% valeur de confinement avec limite + mecanisme pour disruption
switch option.runaway
case 4
    tauconf = max(zs.taue,max(taus,taur));  
case 3
    tauconf = max(taus,taur);  
otherwise
    tauconf = max(taus,zs.taue ./ fh .^ 2);
end
taueff  = 1 ./ (1 ./ max(eps,tauconf) - max(eps,djdtsj));

% graine permanente
switch option.runaway
case {2,3,4}
	if option.wlh > 0
		vnlh      = phys.c ./ (option.npar0 - phys.c ./ (option.freqlh .* 1e9) ./ option.wlh);
	else
		vnlh      = phys.c ./ option.npar0;
	end
	flh       = max(0,vnlh - vcrit) ./ max(eps,vnlh - vth) .* (max(vnlh,vcrit) > vth);
	source_lh = flh .* zs.ilh ./ taueff;
	%figure(21);clf;
        %semilogy(cons.temps,eps+source,'b',cons.temps,eps+source_lh,'r',cons.temps,eps+source + source_lh,'g');
	%drawnow
	source = source  + source_lh;

otherwise

	iseed = zs.ip .* (1 - fh);      
	source = source + max(0,- iseed ./ taueff);
end
		
% LH supprine les runaway
switch option.runaway
case 1
	taueff(zs.plh > (1e3 .* zs.vp)) = 1e-6;
end
if option.evolution == 1
	iini = zs.irun(1);
else
	iini = max(0,source(1) .* zs.taue(1));
end
[trep,jed] = z0ode(cons.temps,source,taueff,iini,-1);

% courant + stabilisation numerique
if all(isfinite(zs.ipar))  && all(zs.ipar ~= 0)
    imax   =  0.99 .* abs(zs.ipar - zs.iboot); 
else
    imax   =  0.99 .* abs(zs.ip); 
end
switch option.runaway
case 4
    imax   =  min(imax, abs(zs.ip)); 
end    
irun = 0.7 .* irun + 0.3 .* max(-imax,min(imax,jed));

%     figure(22);clf
%  %  % plot(cons.temps,tauc,cons.temps,tauc2)
%      subplot(5,1,1);plot(cons.temps,abs(source .* tauconf),'b',cons.temps,abs(jed),'r');set(gca,'ylim',[0,2.*max(zs.ip)]);
%      subplot(5,1,2);plot(cons.temps,djdtsj,'r',cons.temps,1./tauconf,'b');
%      subplot(5,1,3);semilogy(cons.temps,zs.taue,'r',cons.temps,taus,'b',cons.temps,taueff,'g');
%      subplot(5,1,4);semilogy(cons.temps,vcrit,cons.temps,vth);
%      subplot(5,1,5);plot(cons.temps,zs.vloop);
%      drawnow
%  %  
 
if isappdata(0,'RUNAWAY_EXP')
	runexp = getappdata(0,'RUNAWAY_EXP');
	ipref  = interp1_ex(runexp.temps,runexp.ip,cons.temps,'nearest','extrap');
	irun   = interp1_ex(runexp.temps,runexp.irun,cons.temps,'nearest','extrap');
	irun   = irun ./ ipref .* zs.ip;
	if all(isfinite(zs.ipar))  && all(zs.ipar ~= 0)
	    imax   = abs(pi .* zs.ipar);
	else 
	    imax   =  0.99 .* abs(zs.ip);  % for fisrt loop
	end
	irun   = max(-imax,min(imax,irun));
end





