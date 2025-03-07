% calcul le flux de neutron tt
function [neutron_total,neutron_th,neutron_nbi_th,neutron_nbi_nbi,pttfus,proton_tt,picrh_nbi,einj] = ...
	             z0neutron_tt(option,cons,zs,profli)

% constante physique (phys)
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

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);


% % case of Hydrogen NBI in DT plasma
% if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
%    gas_nbi = -1; 
% else
%    gas_nbi = option.gaz; 
% end
% 
% 
% data preparation
x = profli.xli;
temps = cons.temps;
ux  = 1 - x .^ 2;
ve  = ones(size(x));
vt  = ones(size(temps));
if isfield(profli,'temps') && (length(profli.temps) ~= length(zs.temps))
  noms = fieldnames(profli);
  tp = profli.temps;
  for k = 1:length(noms)
      if ~isempty(profli.(noms{k})) && all(all(isfinite(profli.(noms{k})))) && (size(profli.(noms{k}),1) == length(tp))
%            var = NaN * ones(length(temps),size(profli.(noms{k}),2));
%            for l=1:size(profli.(noms{k}),2)
%  	      var(:,l) = interp1(tp,profli.(noms{k})(:,l),temps,'pchip','extrap');
%            end
	  profli.(noms{k}) = interp1(tp,profli.(noms{k}),temps,'pchip','extrap'); 
      end
  end
end
spr = profli.spr;
vpr = profli.vpr;
nep = profli.nep;
tep = profli.tep;
tip = profli.tip;
nip = profli.nip;
n1p = profli.n1p;
nD  = profli.n1p .* ((zs.nDm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
nT  = profli.n1p .* ((zs.nTm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
pnbi = profli.pnbi;
zeffp = profli.zeff;


% output initialisation
neutron_total   = 0 .* vt;
neutron_th      = 0 .* vt;
neutron_nbi_th  = 0 .* vt;
neutron_nbi_nbi = 0 .* vt;
pttfus          = 0 .* vt;
proton_tt       = 0 .* vt;
picrh_nbi       = 0 .* vt;
einj            = 0 .* vt;

% palsma de fond
switch option.gaz
case 3
   zj = 1;
   aj = (2  + 3 .* cons.iso) ./  (1 + cons.iso);
   %     Wdt  = (phys.e .* (3/2) .* zs.vp .* zs.fracmino .* zs.nDm .* zs.tem .* zs.tite + ...
   %  	  (1 - cons.ftnbi) .* zs.esup_nbi) ./ zs.vp;
otherwise
   % No tritium
   % case gaz == 5 -> can be added later if neeeded  
   return
end


% gestion de la puissance NBI
if isfield(option,'nb_nbi')
    nb_nbi = option.nb_nbi;
else
    nb_nbi = 1;
end
% only T part of NBI is to be used
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    ftnbi = zeros(size(cons.ftnbi));
else
    ftnbi = real(cons.ftnbi);
end
pnbi_th = real(zs.pnbi_th);
pnbi  = profli.pnbi ./ (max(1,trapz(x,vpr .* real(profli.pnbi),2)) * ve) .* (real(zs.pnbi_th) * ve);
if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        ftnbi2 = zeros(size(cons.ftnbi));
    else
        ftnbi2 = imag(cons.ftnbi);
    end
    pnbi_th2 = imag(zs.pnbi_th);
    pnbi2  = profli.pnbi ./ (max(1,trapz(x,vpr .* imag(profli.pnbi),2)) * ve) .* (imag(zs.pnbi_th) * ve);
    pnbi   = pnbi + sqrt(-1) .* pnbi2;
    pnbi_th = pnbi_th + sqrt(-1) .* pnbi_th2;
    ftnbi   = ftnbi + sqrt(-1) .* ftnbi2;
end

% choix du minoritaire
switch option.mino
    case 'T'
        ag = 3;
        zg = 1;
        lg = 7.92e-3;
        %correction de einj pour tenir compte de l'interaction avec icrh
        %fdicrh    = 1 .* zs.nem; %?
        fdicrh    = 0.25e-18 .* zs.nem.^2; %?
        if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
            Wtt  = 0 .* real(zs.esup_nbi) ./ zs.vp;
        else
            Wtt  = real(cons.ftnbi) .* real(zs.esup_nbi) ./ zs.vp;
        end
        wctt      = 15.2e6 .* zj ./ aj .* profli.fdia(:,1) ./ profli.Raxe(:,1);
        if option.cmin > 0
            ficrh_nbi = fdicrh .* Wtt ./ zs.nmino ./ wctt;
        else
            ficrh_nbi = 0 .* vt;
        end
        ficrh_nbi = max(0,min(10,ficrh_nbi));
        picrh_nbi = ficrh_nbi ./ (1 + ficrh_nbi) .* zs.picrh;
        einj      = option.einj .* max(1,min(10,(picrh_nbi + real(pnbi_th)) ./ max(1,real(pnbi_th))));
        pnbi_th   = picrh_nbi + pnbi_th;
        
        if nb_nbi == 2
            if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
                Wtt2  = 0 .* real(zs.esup_nbi) ./ zs.vp;
            else
                Wtt2  = imag(cons.ftnbi) .* imag(zs.esup_nbi) ./ zs.vp;
            end
            if option.cmin > 0
                ficrh_nbi2 = fdicrh .* Wtt2 ./ zs.nmino ./ wctt;
            else
                ficrh_nbi2 = 0 .* vt;
            end
            ficrh_nbi2 = max(0,min(10,ficrh_nbi2));
            picrh_nbi2 = ficrh_nbi2 ./ (1 + ficrh_nbi2) .* zs.picrh;
            einj2      = option.einj2 .* max(1,min(10,(picrh_nbi2 + imag(pnbi_th)) ./ max(1,imag(pnbi_th))));
            pnbi_th    = pnbi_th + sqrt(-1) .* picrh_nbi2;
            ficrh_nbi  = ficrh_nbi + sqrt(-1) .* ficrh_nbi2;
            picrh_nbi  = picrh_nbi + sqrt(-1) .* picrh_nbi2;
            einj       = einj + sqrt(-1) .* einj2;
        end
    otherwise
        ficrh_nbi  = 0 .* vt;
        einj       =  option.einj + sqrt(-1) .* option.einj2;
end
%  figure(21);clf
%  subplot(3,1,1)
%  plot(cons.temps,ficrh_nbi);
%  subplot(3,1,2)
%  plot(cons.temps,einj);
%  subplot(3,1,3)
%  plot(cons.temps,picrh_nbi,'b',cons.temps,pnbi_th,'r');
%  drawnow
%keyboard

%  figure(21)
%  clf;
%  plot(cons.temps,zs.nmino)
%  drawnow

% fusion TT du plasma thermique (correction ,doit etre petit, l'enrichissement du milieu en T, H et He3 n'est pas prise en compte)
% 2 neutrons par reaction
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(max(tip(:),13.6));
sn       = 0.5 .* nT .^ 2 .* reshape(tt.sv,size(tip)) .* (tip >= dd_n.timin) .*  (tip <= dd_n.timax);
neutron_th    = trapz(x,vpr .* sn,2);

% calcul de l'interraction faisceau-plasma
pmod        = max(1,real(pnbi) + imag(pnbi));
nTi         = max(1e13,trapz(x,nT .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
tii         = max(30,trapz(x,tip .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
tei         = max(30,trapz(x,tep .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
nei         = max(1e13,trapz(x,nep .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
%taus_nbi    = 6.27e8 .* 2 .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
fact        = (zs.nDm ./2 + zs.nTm ./ 3 + (zs.n1m - zs.nTm - zs.nDm) + zs.nhem  + zs.nimpm ./2) ./ zs.nem; 
%ecrit_nbi   = max(30,14.8 .* tei .* (2 .^ (3/2)  .* fact) .^ (2/3));   % c'est l'energie liee a vc
% c'est l'energie liee a vc
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    ecrit_nbi = max(30,14.8 .* tei .* (1.^ (3/2)  .* fact) .^ (2/3));
    taus_nbi       = 6.27e8 .* 1 .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
else
    ecrit_nbi = max(30,14.8 .* tei .* ((2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
    taus_nbi       = 6.27e8 .* (2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
end
if nb_nbi == 2
    % c'est l'energie liee a vc
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        ecrit_nbi2 = max(30,14.8 .* (tei .* 1.^ (3/2)  .* fact) .^ (2/3));
        taus_nbi2       = 6.27e8 .* 1 .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
    else
        ecrit_nbi2 = max(30,14.8 .* tei .* ((2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
        taus_nbi2       = 6.27e8 .* (2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
    end
    taus_nbi   = taus_nbi  + sqrt(-1) .* taus_nbi2;
end

if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    [neutron_nbi_th,emean] = znbi_tt(nTi,tii,real(pnbi_th) .* 0,real(taus_nbi),real(einj),real(ecrit_nbi));
else
    [neutron_nbi_th,emean] = znbi_tt(nTi,tii,real(pnbi_th) .* real(ftnbi),real(taus_nbi),real(einj),real(ecrit_nbi));
end
emean                  = max(zs.tem,emean);
if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        [neutron_nbi_th2,emean2] = znbi_tt(nTi,tii,imag(pnbi_th) .* 0,imag(taus_nbi),imag(einj),imag(ecrit_nbi));
    else
        [neutron_nbi_th2,emean2] = znbi_tt(nTi,tii,imag(pnbi_th) .* imag(ftnbi),imag(taus_nbi),imag(einj),imag(ecrit_nbi));
    end
    emean2                   = max(zs.tem,emean2);
    neutron_nbi_th           = neutron_nbi_th + neutron_nbi_th2;
    emean                    = emean + sqrt(-1) .* emean2;
end

% approximation faisceau-faisceau avec la meme formule (c'est faux)
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    [esupranbi,pTth,tauseff] = zsupra0(temps,real(pnbi_th) .* 0 ,real(taus_nbi),1e-6*ones(size(taus_nbi)),real(ecrit_nbi),real(einj),2);
else
    [esupranbi,pTth,tauseff] = zsupra0(temps,real(pnbi_th) .* real(ftnbi),real(taus_nbi),1e-6*ones(size(taus_nbi)),real(ecrit_nbi),real(einj),2);
end
fon  = ((real(pnbi_th) - real(picrh_nbi)) ./ zs.vp) > 1e3;
if any(fon)
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        ctr  = 0 .* real(tauseff) ./ 2 ./ zs.vp ./ (1.602176462e-19 .* real(emean));
    else
        ctr  = real(ftnbi) .* real(tauseff) ./ 2 ./ zs.vp ./ (1.602176462e-19 .* real(emean));
    end
    ctr  = fon .* ctr + (~fon) .* mean(ctr(find(fon)));
	ttr  = fon .* real(taus_nbi) + (~fon) .* mean(real(taus_nbi(find(fon))));
	neutron_nbi_nbi   = 0.5 .* znbi_tt(real(pnbi_th) .* ctr, real(emean),pTth,ttr,real(einj),real(ecrit_nbi));
else
	neutron_nbi_nbi = zeros(size(cons.temps));
end

neutron_nbi_nbi_1 = neutron_nbi_nbi;
if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        [esupranbi2,pTth2,tauseff2] = zsupra0(temps,imag(pnbi_th) .* 0,imag(taus_nbi),1e-6*ones(size(taus_nbi)),imag(ecrit_nbi),imag(einj),2);
    else
        [esupranbi2,pTth2,tauseff2] = zsupra0(temps,imag(pnbi_th) .* imag(ftnbi),imag(taus_nbi),1e-6*ones(size(taus_nbi)),imag(ecrit_nbi),imag(einj),2);
    end
    fon2  = ((imag(pnbi_th) - imag(picrh_nbi)) ./ zs.vp) > 1e3;
    if any(fon2)
        if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
            ctr2  = 0 .* tauseff2 ./ 2 ./ zs.vp ./ (1.602176462e-19 .* imag(emean));
        else
            ctr2  = imag(ftnbi).* tauseff2 ./ 2 ./ zs.vp ./ (1.602176462e-19 .* imag(emean));
        end
        ctr2  = fon2 .* ctr2 + (~fon2) .* mean(ctr2(find(fon2)));
        ttr2  = fon2 .* imag(taus_nbi) + (~fon2) .* mean(imag(taus_nbi(find(fon2))));
    end
    % 2 sur 2
    if any(fon2)
        neutron_nbi_nbi_2   = 0.5 .* znbi_tt(imag(pnbi_th) .* ctr2,imag(emean),pTth2,ttr2,imag(einj),imag(ecrit_nbi));
    else
         neutron_nbi_nbi_2   = 0 * neutron_nbi_nbi_1;
      end 
      
      if any(fon2) && any(fon)
	  % 1 sur 2
	  neutron_nbi_nbi_3   = 0.5 .* znbi_tt(real(pnbi_th) .* ctr,real(emean),pTth2,ttr2,imag(einj),imag(ecrit_nbi));
	  % 2 sur 1
	  neutron_nbi_nbi_4   = 0.5 .* znbi_tt(imag(pnbi_th) .* ctr2,imag(emean),pTth,ttr,real(einj),real(ecrit_nbi));
      else
         neutron_nbi_nbi_3   = 0 * neutron_nbi_nbi_1;
         neutron_nbi_nbi_4   = 0 * neutron_nbi_nbi_1;
      end
      esupranbi    = esupranbi + sqrt(-1) .* esupranbi2;
      pTth         = pTth + sqrt(-1) .* pTth2;
      tauseff      = tauseff  + sqrt(-1) .* tauseff2;
      neutron_nbi_nbi = neutron_nbi_nbi_1 + neutron_nbi_nbi_2 + neutron_nbi_nbi_3 + neutron_nbi_nbi_4;
end

%  figure(21);clf
%  subplot(3,1,1)
%  plot(cons.temps,ctr);
%  subplot(3,1,2)
%  plot(cons.temps,ttr);
%  subplot(3,1,3)
%  plot(cons.temps,neutron_nbi_nbi);
%  drawnow




% autre calcul a partir de J.D Strachan NF 33,7 (1993) p 991-	
% normalisation beam beam sur beam-plasma		    
%sbt  = 0.167e-9 .* nDi .* pnbi_th.* (1 -ftnbi) .* taus_nbi;
%sbb  = 1.3e13 .* (pnbi_th.* (1 -ftnbi)) .^ 2 .* taus_nbi .^ 2 ./ option.einj .^2;
%disp([mean(neutron_nbi_th),mean(neutron_nbi_nbi)]);
%disp([mean(sbt),mean(sbb)]);
%disp([mean(sbb ./ sbt .* neutron_nbi_th),mean(sbb),mean(neutron_nbi_nbi)]);
%neutron_nbi_nbi = sbb ./ sbt .* neutron_nbi_th; 
			    

% nombre de neutron total
% 2 neutrons par reaction
neutron_th        = 2 .* max(0,neutron_th);
neutron_nbi_th    = 2 .* max(0,neutron_nbi_th);
neutron_nbi_nbi   = 2 .* max(0,neutron_nbi_nbi);
neutron_total     = neutron_th + neutron_nbi_th +  neutron_nbi_nbi;

% puissance de fusion complementaire
% only one He4 for 2 neutrons
% E_He4 ~ 1.6 MeV
pttfus    = phys.e .* 1.6e6 .* neutron_total ./ 2;

%    figure(21);clf
%    plot(cons.temps,neutron_th,'g',cons.temps,neutron_nbi_th,'b',...
%    cons.temps,neutron_nbi_nbi,'m', cons.temps,neutron_total,'c', ...
%    evalin('base','z0dinput.cons.temps'),evalin('base','z0dinput.exp0d.ndd'),'r');
%    drawnow
%keyboard



function  s  = zpmv(t,p,e)

%  figure(21);clf
%  subplot(2,1,1)
%  plot(t,p)
%  subplot(2,1,2);
%  plot(t,e);
%  drawnow
ve = ones(1,size(e,2));

indok = find(isfinite(p) & all(isfinite(e),2));
if length(indok) < 3
	indok = find(all(isfinite(e),2));
	if isempty(indok)
		s = NaN*ve;
		disp('zpmean : invalid data ...');
	else
		s = mean(e(indok,:),1);
	end
else
	t = t(indok);
	p = p(indok);
	e = e(indok,:);
	indp = find(p > 0);
	%disp(length(indp))
	if length(indp) > 2
		s = trapz(t(indp),(p(indp) * ve) .* e(indp,:),1) ./ (trapz(t(indp),eps + p(indp),1) * ve);
	else
		s = mean(e,1);
	end
end


function  s  = zpmean(t,p,e)

%  figure(21);clf
%  subplot(2,1,1)
%  plot(t,p)
%  subplot(2,1,2);
%  plot(t,e);
%  drawnow

indok = find(isfinite(p) & isfinite(e));
if length(indok) < 3
	indok = find(isfinite(e));
	if isempty(indok)
		s = NaN;
		disp('zpmean : invalid data ...');
	else
		s = mean(e(indok));
	end
else
	t = t(indok);
	p = p(indok);
	e = e(indok);
	indp = find(p > 0);
	%disp(length(indp))
	if length(indp) > 2
		s = trapz(t(indp),p(indp) .* e(indp)) ./ trapz(t(indp),eps + p(indp));
	else
		s = mean(e);
	end
end



