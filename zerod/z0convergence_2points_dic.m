% appel du model a 2 points pour METIS + convergence
% test : [tebord,nelim,telim,qpl,err,nb,indbad,fmom,delta_qpl_rad,delta_qpl_neutral,qpl_in_loc,pl,zeff_div,gamma,mach_target] =  ...
% z0convergence_2points_dic(post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,post.zerod,post.profil0d);
function [tebord,nelim,telim,qpl,err,nb,indbad,fmom,delta_qpl_rad,delta_qpl_neutral,qpl_in_loc,pl, ...
          zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond,profli] = ...
          z0convergence_2points_dic(option,cons,geo,zs,profli,mode_conv)

% test de la convergence
% mode_conv = 0 ou non fourni -> dicotomie par defaut
% mode_conv = 1 -> depart par la valeur max
% mode_conv = 2 -> simple boucle a pas fixe
% mode_conv = 3 -> depart proche de 0
% mode_conv = 4 -> tirage aleatoire avec selection de la meilleur solution.
% mode_conv = 5 -> genre recuit
% test des entrees
if nargin < 6
	mode_conv = 0;
end


% recycling control
if option.Recycling_target == 0
      % use global recycling coeficient
      Recycling = option.Recycling;
else
      % use taget recycling coeficient
      Recycling = option.Recycling_target;

end

% Sn
if ~isfield(option,'Sn_fraction')
    Sn_fraction = 0;
else
    Sn_fraction = option.Sn_fraction;
end

if ~isfield(option,'te_max')
    te_max = 105;
else
    te_max = max(option.te_max /1e3 +5);
end

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
    

% compatibilite mode fast pour plot
if isfield(profli,'temps') && (length(cons.temps) ~= length(profli.temps))
        times = cons.temps;
        temps = profli.temps;
	% modification des donnees
	noms = fieldnames(profli);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(profli,nomc);
		if all(size(val) > 1)
			warning off
			valn  = interp1(temps,val,times,'linear');
			warning on
			indbad      = find(any(~isfinite(valn),2));
			if ~isempty(indbad)
				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
			end
			profli = setfield(profli,nomc,valn);
		end
	end
	profli.temps = times;
	
end

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


% variables locales pour la convergence (first guess)
tebord     = zs.tebord;
if isfield(zs,'nelim')
	nelim      = zs.nelim;
else
	nelim      = 2 .* zs.nebord;
end
if isfield(zs,'telim')
	telim      = zs.telim;
else
	telim      = zs.tebord;
end

if isfield(profli,'vpr_tor')
    switch option.sol_model
    case '2_points'
        %  on fait les calculs
    otherwise
	return
    end
else
    return
end

% puissance conduite a la separatrice
if option.W_effect == 1 && ~isfield(option,'plot2points')
	pl        = max(zs.pin .* sqrt(eps),(zs.pin - zs.pioniz - zs.pcyclo));
else
	pl        = max(zs.pin .* sqrt(eps),(zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz - zs.pradsol));
end
% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* option.lcx .* zs.qa) .^ 2);  
switch option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lcpol;
case 3
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lclim;
otherwise
	lc  = lcx;
end

if isfield(option,'Sq') && (option.Sq > 0)
  f_Sq   = 1.64 .* option.Sq ./ (zs.dsol + 1.64 .* option.Sq);
  fpower = option.fpower.* f_Sq;
elseif isfield(option,'Sq') && (option.Sq < 0)
  %A. Scarabosio paper (Journal of Nuclear Materials 463 (2015) 49-54)
  Sq_mm = abs(option.Sq) .* (double(zs.modeh) .* (0.12 .* (zs.plim ./ 1e6) .^ 0.21 .*  max(1e13/1e19,zs.nebord ./ 1e19) .^ -0.02) .* profli.bpol(:,end) .^ -0.82 .* geo.R .^ 0.71 +  ...
	  (1- double(zs.modeh)) .* (0.13 .* max(1e13/1e19,zs.nebord ./ 1e19) .^ 1.1) .* profli.bpol(:,end) .^ -0.84);
  %figure(21);plot(cons.temps,Sq_mm);drawnow
  f_Sq   = 1.64 .* Sq_mm ./ (zs.dsol + 1.64 .* Sq_mm);
  fpower = option.fpower .* f_Sq;
else
  fpower = option.fpower; 
end


% flux //  (formula 5.64 et 5.75 Stangeby)
rext         = profli.Raxe + geo.a * profli.xli;
btor         = profli.fdia ./ rext;
grho         = abs((profli.rmx(:,end) * ones(size(profli.xli)))./ max(eps,pdederive(profli.xli,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(profli.xli,profli.psi,0,2,2,1) ./ rext .* grho ./ (profli.rmx(:,end) * ones(size(profli.xli)));
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));
% ut = atan(profli.bpol(:,end) ./  profli.fdia(:,end) .* profli.Raxe(:,end));
%qpl = pl .* (1 - fesol) ./ (4 .* pi  .* rext(:,end) .* zs.dsol .* sin(ut));
%qpl_in = fpower .* pl  ./ (4 .* pi  .* profli.Raxe(:,end) .* zs.dsol .* sin(ut));
% le flux n'est pas forcement divise en 2 :
% ref / Stangeby, NF 51 (2011) 063001
qpl_in = fpower .* pl  ./ (2 .* pi  .* profli.Raxe(:,end) .* zs.dsol .* sin(ut));
qpl_in = max(1,qpl_in);

% initialisation qpl
switch mode_conv
case 1
	qpl      = qpl_in;
	qpl_step = qpl;
	loop = 0;
case 2 
	qpl      = qpl_in;
	qpl_step = qpl;
	loop = 1;
case 3
	qpl      = qpl_in .* sqrt(eps);
	qpl_step = qpl_in;
	loop = 0;
case 4 
	qpl      = qpl_in;
	qpl_step = qpl;
	loop = 2;
case 5 
	qpl      = 0.5 .* qpl_in;
	qpl_step = qpl_in;
	loop = 3;
otherwise
	qpl      = 0.5 .* qpl_in;
	qpl_step = 0.5 .* qpl_in;
	loop = 0;
end

% calcul de la contribution du tungsten si besoin
switch option.imp_div 
case 'auto'
    zeff_div = max(1,profli.zeff(:,end));
    %figure(21);clf;plot(cons.temps,zeff_div);hold on
    for k=1:7
      fzmax_div = z0fleak_zmax(geo,zs,profli,option.zmax,option.rimp,zeff_div,lc,option.residence_time);
      fact_ne = ones(size(zs.nebord));
      if option.W_effect == 1
	      % dans le modele initial, seul la temperature de bord intervient; la majeur partie de la ligne de champs a une temperature proche de tebord
	      ya0_w =  profli.nwp(:,end) ./ profli.nep(:,end);
	      zwave  = (1 - Sn_fraction) .* z0wavez(zs.tebord) + Sn_fraction .* z0snavez(zs.tebord);
	      zeff_div = max(1,profli.zeff(:,end)) + 0.01 .* abs(fzmax_div) .* option.zmax .* (option.zmax -1) + ya0_w .* zwave .* max(1,zwave -1);
	      fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax + ya0_w .* zwave;
      else
	      zeff_div = max(1,profli.zeff(:,end)) + 0.01 .* abs(fzmax_div) .* option.zmax .* (option.zmax -1);
	      fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax;
      end
     %figure(21);plot(cons.temps,zeff_div);drawnow
  end
   % convergence must be checked
otherwise
    % fixed
    fzmax_div = option.fzmax_div;
    fact_ne = ones(size(zs.nebord));
    if option.W_effect == 1
	    % dans le modele initial, seul la temperature de bord intervient; la majeur partie de la ligne de champs a une temperature proche de tebord
	    ya0_w =  profli.nwp(:,end) ./ profli.nep(:,end);
	    zwave  = (1 - Sn_fraction) .* z0wavez(zs.tebord) + Sn_fraction .* z0snavez(zs.tebord);
	    zeff_div = max(1,profli.zeff(:,end)) + 0.01 .* abs(fzmax_div) .* option.zmax .* (option.zmax -1) + ya0_w .* zwave .* max(1,zwave -1);
	    fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax + ya0_w .* zwave;
    else
	    zeff_div = max(1,profli.zeff(:,end)) + 0.01 .* abs(fzmax_div) .* option.zmax .* (option.zmax -1);
	    fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax;
    end

end

% carbon feed back on plate 
if ~isfield(option,'carbonblow')
    option.carbonblow = 0;
elseif  (option.carbonblow ~= 0) && (option.zimp == 6)
   [ya0_C,zavez_C]    = z0carbonblow(option,zs,profli);
   zeff_div = zeff_div + min(1,abs(option.carbonblow) .* ya0_C) .* zavez_C .* max(1,zavez_C - 1);
   fact_ne = fact_ne + min(1,abs(option.carbonblow) .* ya0_C) .* zavez_C;
end


% les donn?es, ici viennent du transport de coeur
if isfield(profli,'tip')
	fie = 1 + zs.tibord ./ zs.tebord;
	tite_loc = zs.tibord ./ zs.tebord;
	if isfield(profli,'qe') && isfield(profli,'qi')
		fpe = min(1,max(0.1,profli.qe(:,end) ./ max(eps,profli.qe(:,end) + profli.qi(:,end))));
	else
		fpe = min(1,max(0.1,zs.pel ./ (zs.pel + zs.pion)));
	end
	%figure(21);clf;plot(cons.temps,fie,cons.temps,fpe);drawnow
else
	fie = 4;
	fpe = 0.5;
	tite_loc = zs.tite;
end

% radial position of the target effect
% ref: T;W. Petrie et al, Nuclear Fusion 53 (2013) 113024.
if isfield(option,'fR_target') && option.fR_target < 0
      f_R = abs(option.fR_target) .* geo.R ./ rext(:,end);
      warning off 
      corr_f_R = log(f_R) ./ (f_R - 1);
      warning on
      corr_f_R(~isfinite(corr_f_R)) = 1;
      fpe = corr_f_R .* fpe;
      %figure(21);clf;plot(cons.temps,corr_f_R);drawnow
end

% sheath factor (include ion and electron)
% STANGEBY p 649  25.46
gamma      = 2.5 .* tite_loc + 2 ./ (1 - option.de)  -  ...
             0.5 .* log((2 .* pi.* phys.me ./ phys.mp) .* (1 + tite_loc) .* (1 - option.de) .^ -2); % 25.46
% we consider only  electron for some correction term
gamma_e      = 2 ./ (1 - option.de)  -  0.5 .* log((2 .* pi.* phys.me ./ phys.mp) .* (1 + tite_loc) .* (1 - option.de) .^ -2); % 25.46

% alpha_e pour la correction cinetique
%alpha_e = 0.82; % et non pas 1.5 a cause de tanh dans la formule

%  initilisation fmom
if option.fmom == 0
	fmom = ones(size(qpl_in));
else
	fmom = option.fmom .* ones(size(qpl_in));
end
fmom_step = 0.5 .* ones(size(qpl_in));

% 1ere valeur 
if ~isfield(option,'plot2points')
	[tebord,telim,nelim,indbad,noconv] = z0twopoints(qpl_in,fact_ne .* zs.nebord,lc,zs.zeff,zs.meff,abs(option.fcond) .* fpe,fmom .* fie,tite_loc,gamma);
end

% preparation de la convergence
nbmax = 301;
nb = nbmax;
err_best  = Inf .* ones(size(qpl_in));
qpl_best  = qpl_in;
qpl_in_loc_best = qpl_in;
fmom_best = ones(size(qpl_in));
zeff_div_best= ones(size(qpl_in));
delta_qpl_rad_best = zeros(size(qpl_in));
delta_qpl_neutral_best = zeros(size(qpl_in));
err_stop = sqrt(eps);
zs_loc   = zs;
sgn = 1;
fact_ne_ext = fact_ne;
if ~isfield(option,'plot2points')
  fact_ne     = ones(size(zs.nebord));
end

% debut de la boucle de convergence
while (nb > 0) && (max(qpl_step ./ max(1 , qpl_in)) > sqrt(eps))
	
        % memorisation
	telim_mem   = telim;
        tebord_mem  = tebord;
        nelim_mem   = nelim;
        fact_ne_mem = fact_ne;
        % dans la dichotomie, le model a 2 point est appele en premier
	
	% boucle sur fmom
	%  initilisation fmom
	if option.fmom == 0
		fmom = ones(size(qpl_in));
	else
		fmom = option.fmom .* ones(size(qpl_in));
	end
        %figure(27);clf
	for kmom=1:11
		fmom_mem = fmom;
		% appel model a 2 points 
		if ~isfield(option,'plot2points')
			[tebord,telim,nelim,indbad,noconv] = z0twopoints(qpl,fact_ne .* zs.nebord,lc,zs.zeff,zs.meff,abs(option.fcond) .* fpe .*  qpl_in ./ max(1,qpl),fmom .* fie,tite_loc,gamma);
			
			% calcul de la correction cinetique
			% ref : Stangeby ch 26 + 
                        % W. Fundamenski, topical review, PPCF 47 (2005) p R163-R208, formule 10
                        % alpha_e ~ 1.5 recommande mais avec tanh, on choisi alpha_e = 0.82
			if option.fcond < 0
				% point fixe
%				figure(23);clf; 
				for lcond=1:5
					qelim = option.alpha_e .* sqrt(phys.e .* telim ./ phys.me) .*  phys.e .* telim .* nelim;
					fcond = abs(option.fcond) .* fpe ./ (1 - min(0.9,tanh(qpl ./ qelim))).*  qpl_in ./ max(1,qpl);	
					[tebord,telim,nelim,indbad,noconv] = z0twopoints(qpl,fact_ne .* zs.nebord,lc,zs.zeff,zs.meff,fcond,fmom .* fie,tite_loc,gamma);
					
%  					figure(23);
%  					subplot(3,1,1);
%  					semilogy(cons.temps,qelim/1e6,cons.temps,qpl/1e6);
%  					hold on
%  					subplot(3,1,2);
%  					plot(cons.temps,fcond);
%  					hold on
%  					subplot(3,1,3);
%  					semilogy(cons.temps,tebord,cons.temps,telim);
%  					hold on
%  					%drawnow
				end
%					drawnow
                        else
				fcond = abs(option.fcond) .* fpe .*  qpl_in ./ max(1,qpl);
		        end
		else
			qelim = option.alpha_e .* sqrt(phys.e .* telim ./ phys.me) .*  phys.e .* telim .* nelim;
			fcond = abs(option.fcond) .* fpe ./ (1 - min(0.9,tanh(qpl ./ qelim))).*  qpl_in ./ max(1,qpl);	
	
		end
		%figure(21); plot( qpl_in ./ max(1,qpl));drawnow		
		if option.mach_corr == 1
			% nombre de mach (model simple pour detachement)
			% ref JET-P(64) 05 , P. H. Harbour and A. Loarte
                        if  option.eioniz > 0
				Kt          = min(0.9,max(-0.9,0.5 .* (1 - Recycling .* option.eioniz ./ gamma_e ./ max(eps,telim))));
			else
				Kt          = min(0.9,max(-0.9,0.5 .* (1 - Recycling .* z0eioniz_div(telim,nelim) ./ gamma_e ./ max(eps,telim))));
			end
			mach_target = sqrt((1 - Kt) ./ (1 + Kt));
			%figure(21);plot(mach_target);drawnow
			
			if ~isfield(option,'plot2points')
				% formule Stangeby 14.4-14.11
				ftm = (gamma_e .* (1 + mach_target .^ 2) ./ (2 .* abs(mach_target) .* (gamma_e - 1 + mach_target .^ 2))) .^ 2;
				telim = ftm .* telim;
				fnm   = 2 ./ (ftm .* (1 + mach_target .^ 2));
				nelim = fnm .* nelim;
				% tebord est inchange
				% le flux sur la plaque de divertor est change !
			end
		else	
			mach_target = ones(size(qpl));
		end
 		% calcul de fmom d'apres Pitcher, Self et Ewald dans Plasma edge physics for tokamaks de Ralf Schneider.
                % attention notation de Stangeby, ici fmom = fm du papier original fm = (1-fmom) !!!
		if option.fmom == 0
			% longueur de ionisation pres des plaques
			[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(telim,telim .* (zs.tibord ./ zs.tebord));  % attention ici le rapport ti/te est calculer a l'exterieur, tebord ne doit pas etre mis a jour
			%% equilibre entre 1s et 2s pour le neutres de centre 
			alphas = nelim .* sss ./ (nelim .* sss + Ass);
			% etat d'equilibre 1s/2s
			sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
			alpha  = (sviss  + sii)  ./ (sviss  + sii+ svcx);
			fmom_corr    = max(0.1,min(1,2 .* (alpha ./ (alpha + 1)) .^ ((alpha + 1) ./ 2)));
		else
			break;		
		end
		if kmom < 3
			fmom       = fmom_corr; 
		elseif kmom < 7
			fmom       = 0.5 .* fmom + 0.5 .* fmom_corr; 
		else
			fmom       = 0.7 .* fmom + 0.3 .* fmom_corr; 
		end
		fmom_step  = fmom_mem - fmom_corr;
%  		figure(27);
%  		plot(cons.temps,fmom,'b',cons.temps,fmom_corr,'r',cons.temps,fmom_step,'g');
%  		hold on
%  		drawnow;
		if all(abs(fmom_step) < err_stop)
			break;
		end
	end
	%figure(27);clf;hist(fmom_step);title(sprintf('kmom = %d',kmom));drawnow
    % Si W, Prad depend tres fortement de la concentration, il faut recalculer qpl_in
    if option.W_effect == 1 && ~isfield(option,'plot2points')
        
        % dans le modele initial, seul la temperature de bord intervient; la majeur partie de la ligne de champs a une temperature proche de tebord
        zs_loc.tebord = tebord;
        zs_loc.telim  = telim;
        zs_loc.nelim  = nelim;
        zs_loc.tibord = tebord ./ max(eps,zs.tebord) .* zs.tibord; % conservation du rapport calculer en dehors de cette fonction
        [nwp,nwm,zu1w,zu2w] = z0acctungsten(option,geo,cons,zs_loc,profli);
        profli.nwp = nwp;
        zs_loc.nwm = nwm;
        
        switch option.gaz
            case 1
                zgaz = 1;
            case 2
                zgaz = 1;
            case 3
                zgaz = 1;
            case 4
                zgaz = 2;
            case 5
                zgaz = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
            case 11
                zgaz = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        end
        % selon la configuration
        switch option.configuration
            case {2,3,4}
                telim       = zs.xpoint .* telim;
                nesol_ave   = (1 - option.fnesol) .* zs.nebord + option.fnesol .* nelim;
            otherwise
                telim       = 0 .* telim;
                nesol_ave   = zs.nebord;
        end
        switch option.gaz
            case 11
                nboronm = zs_loc.nTm;
            otherwise
                nboronm = zeros(size(zs_loc.nem));
        end

        [zs_loc.prad,zs_loc.pbrem,void_loc.pcyclo,zs_loc.pradsol,profli_void] = zrad0(zs_loc.nem,zs_loc.tem, ...
            zs_loc.zeff ./ zs_loc.zmszl,option.zmax,option.zimp,option.rimp,zgaz, ...
            geo.R,geo.a,geo.K,geo.b0,zs_loc.ane,zs_loc.ate,sqrt(-1),zs.vp,zs.sext, ...
            (option.modeh > 0),option.frad,zs.te0,zs.nbar,option.matthews, ...
            telim,nesol_ave,zs_loc.taue,zs_loc.dsol,0,option.z_prad,option.sol_rad, ...
            option.gaunt,option.noncoronal == 1,cons.temps,profli, ...
            option.Sn_fraction,option.te_max,option.brem_rel_cor,nboronm);
        % puissance conduite a la separatrice
        pl_loc     = max(zs.pin  .* sqrt(eps),(zs_loc.pin - zs_loc.prad - zs_loc.pbrem - zs.pcyclo - zs.pioniz - zs_loc.pradsol));
        qpl_in_loc = fpower .* pl_loc  ./ (2 .* pi  .* profli.Raxe(:,end) .* zs.dsol .* sin(ut));
        qpl_in_loc = max(1,qpl_in_loc);
        prad_loc  = zs_loc.prad;
        pradsol_loc  = zs_loc.pradsol;
        
    else
        nwp = profli.nwp;
        qpl_in_loc = qpl_in;
        prad_loc  = zs.prad;
        pradsol_loc  = zs.pradsol;
    end

        % nouvelle evaluation de la puissance compte tenu du rayonnment et des neutrees
        % avec les dernieres valeurs de temperature et densite
	% correction du rayonnement dans le divertor
	% ref : R. Clark et al, Journal of Nuclear Materials, vol 220-222 (1995) p 1028-1032
	% fzmax_div est en %
	if any(option.fzmax_div ~= 0) || ((option.carbonblow ~= 0) && (option.zimp == 6)) || (strcmp(option.imp_div,'auto')) 

		switch option.imp_div 
		case 'auto'
		    zs_loc.tebord = tebord;
		    zs_loc.telim  = telim;
		    zs_loc.nelim  = nelim;
		    zs_loc.tibord = tebord ./ max(eps,zs.tebord) .* zs.tibord; % conservation du rapport calculer en dehors de cette fonction
		    fzmax_div = z0fleak_zmax(geo,zs_loc,profli,option.zmax,option.rimp,zeff_div,lc,option.residence_time);
		otherwise
		    % fixed
		    fzmax_div = option.fzmax_div;
		end

		% calcul de la contribution du tungsten si besoin
        fact_ne = ones(size(zs.nebord));
        if option.W_effect == 1
            ya0_w =  nwp(:,end) ./ profli.nep(:,end);
            zwave  = (1 - Sn_fraction) .* z0wavez(zs.tebord) + Sn_fraction .* z0snavez(zs.tebord);
            zeff_div = max(1,profli.zeff(:,end)) + 0.01 .* abs(fzmax_div) .* option.zmax .* (option.zmax -1) + ya0_w .* zwave .* (zwave -1);
            fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax + ya0_w .* zwave;
        else
            fact_ne = fact_ne + 0.01 .* abs(fzmax_div) .* option.zmax;
        end
        
        % carbon feed back on plate
        if  (option.carbonblow ~= 0) && (option.zimp == 6)
            zs_loc.tebord = tebord;
            zs_loc.telim  = telim;
            zs_loc.nelim  = nelim;
            zs_loc.tibord = tebord ./ max(eps,zs.tebord) .* zs.tibord; % conservation du rapport calculer en dehors de cette fonction
            [ya0_C,zavez_C]    = z0carbonblow(option,zs_loc,profli);
            %figure(21);plot(cons.temps,ya0_C);drawnow
            zeff_div = zeff_div + min(1,abs(option.carbonblow) .* ya0_C) .* zavez_C .* (zavez_C -1);
            fact_ne = fact_ne + min(1,abs(option.carbonblow) .* ya0_C) .* zavez_C;
            if option.carbonblow < 0
                fact_ne = 0.3 .* fact_ne + 0.7 .* fact_ne_mem;
            end
        else
            ya0_C    = zeros(size(telim));
            zavez_C    = zeros(size(telim));
        end
        %figure(21);plot(cons.temps,fact_ne,cons.temps,fact_ne_ext);drawnow
        % vecteur temperature
        tepar = linspace(0,1,101);
        if option.fzmax_div < 0
            teloc  = min(telim , tebord) * ones(size(tepar)) + abs(telim - tebord) * tepar;
            tescale = abs(telim - tebord);
        else
            teloc = tebord * tepar;
            tescale =  tebord;
        end
        teloc = max(sqrt(eps),teloc);
        
        % residence time for impurities:
        if option.noncoronal == -1
            if option.residence_time == 0
                cs_max = trapz(tepar,sqrt(phys.e .* teloc .* ((zeff_div + tite_loc) * ones(size(tepar))) ./ phys.ua ./ getafromz_imas(option.zmax)),2);
                tau_resid_max = lc ./ cs_max;
                cs_imp = trapz(tepar,sqrt(phys.e .* teloc .* ((zeff_div + tite_loc) * ones(size(tepar))) ./ phys.ua ./ getafromz_imas(option.zimp)),2);
                tau_resid_imp = lc ./ cs_imp;
                cs_he = trapz(tepar,sqrt(phys.e .* teloc .* ((zeff_div + tite_loc) * ones(size(tepar))) ./ phys.ua ./ 4),2);
                tau_resid_he = lc ./ cs_he;
                nenc = (fact_ne .* zs.nebord)  * ones(size(tepar));
                %figure(21);clf;plot(cons.temps,tau_resid_max,cons.temps,tau_resid_imp,cons.temps,tau_resid_he);drawnow
            else
                tau_resid_max = option.residence_time * ones(size(cons.temps));
                tau_resid_imp = option.residence_time * ones(size(cons.temps));
                tau_resid_he  = option.residence_time * ones(size(cons.temps));
                nenc          = (fact_ne .* zs.nebord)  * ones(size(tepar));
            end
        end
        
        
        % zmax
        [ua,uz,post]=z0coefpost;
        dz     = abs(uz - option.zmax);
        indimp =  min(find(dz == min(dz)));
        lzimp    = post(indimp).lz;
        tzimp    = post(indimp).te;
        lzimp    = noncoronal(option,tzimp ,lzimp);
        lzimp_coef    = cat(2,1e-38,1e-37,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-27);
        tzimp    = cat(2,1e-37,1e-4,tzimp,te_max);
        fz_bord  = option.rimp .* max(0,profli.nzp(:,end)) ./ profli.nep(:,end) .* 100;
        fz  = (abs(fzmax_div) + fz_bord) ./ zeff_div;
        if option.noncoronal == -1
            lzimp  = zlightznoncoronal(teloc,nenc,tau_resid_max * ones(size(tepar)),option.zmax) .* 1e13;
            if any(~isfinite(lzimp(:)))
                lzimp_p   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
                %figure(21);subplot(2,2,1);loglog(teloc(:),lzimp(:),'or',teloc(:),lzimp_p(:),'b.');
                lzimp(~isfinite(lzimp)) = lzimp_p(~isfinite(lzimp));
            end
        else
            lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        end
        lzimp(~isfinite(lzimp(:))) = 0;
        inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        delta_qpl_rad = 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        
        
        % zimp
        dz     = abs(uz - option.zimp);
        indimp =  min(find(dz == min(dz)));
        lzimp    = post(indimp).lz;
        tzimp    = post(indimp).te;
        lzimp    = noncoronal(option,tzimp ,lzimp);
        lzimp_coef    = cat(2,1e-38,1e-37,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-27);
        tzimp    = cat(2,1e-37,1e-4,tzimp,te_max);
        fz_bord  = max(0,profli.nzp(:,end)) ./ profli.nep(:,end) .* 100 + min(1,abs(option.carbonblow) .* ya0_C) .* 100;
        fz  = fz_bord ./ zeff_div;
        if option.noncoronal == -1
            lzimp  = zlightznoncoronal(teloc,nenc,tau_resid_imp * ones(size(tepar)),option.zimp) .* 1e13;
            if any(~isfinite(lzimp(:)))
                lzimp_p   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
                %figure(21);subplot(2,2,2);loglog(teloc(:),lzimp(:),'or',teloc(:),lzimp_p(:),'b.');
                lzimp(~isfinite(lzimp)) = lzimp_p(~isfinite(lzimp));
            end
        else
            lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        end
        lzimp(~isfinite(lzimp(:))) = 0;
        inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        
        % he
        dz     = abs(uz - 2);
        indimp =  min(find(dz == min(dz)));
        lzimp    = post(indimp).lz;
        tzimp    = post(indimp).te;
        lzimp    = noncoronal(option,tzimp ,lzimp);
        lzimp_coef    = cat(2,1e-38,1e-37,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-27);
        tzimp    = cat(2,1e-37,1e-4,tzimp,te_max);
        switch option.gaz
            case 5
                     fz_bord  = (option.frhe0 + max(0,profli.nhep(:,end)) ./ profli.nep(:,end)) .* 100;
            otherwise
                     fz_bord  = max(0,profli.nhep(:,end)) ./ profli.nep(:,end) .* 100;
        end
        fz  = fz_bord ./ zeff_div;
        if option.noncoronal == -1
            lzimp  = zlightznoncoronal(teloc,nenc,tau_resid_he * ones(size(tepar)),2) .* 1e13;
            if any(~isfinite(lzimp(:)))
                lzimp_p   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
                %figure(21);subplot(2,2,3);loglog(teloc(:),lzimp(:),'or',teloc(:),lzimp_p(:),'b.');
                lzimp(~isfinite(lzimp)) = lzimp_p(~isfinite(lzimp));
            end
        else
            lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        end
        lzimp(~isfinite(lzimp(:))) = 0;
        inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        %void = tau_resid_he * ones(size(tepar));
        %figure(21);subplot(2,2,4);loglog(teloc(:),nenc(:).* void(:),'or',teloc(:),1e16,'g',teloc(:),1e19,'g');drawnow
        % HDT
        %  		tzimp         = [1e-37,1e-4,	0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
        %  		lzimp_coef    = [1e-38,1e-37,	7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 13.6 .* 1.6022e-19 .* 10;
        %  		fz_bord  =profli.n1p(:,end) ./ profli.nep(:,end) .* 100;
        %  		fz  = fz_bord ./ zeff_div;
        %  		lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        %  		lzimp(~isfinite(lzimp(:))) = 0;
        %  		inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        %  		delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        %  tzimp         = [1e-37,1e-4,	0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
        %  lzimp_coef    = [1e-38,1e-37,	7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 13.6 .* 1.6022e-19 .* 10;
        %  fz_bord  =profli.n1p(:,end) ./ profli.nep(:,end) .* 100;
        %  fz  = fz_bord ./ zeff_div;
        %  lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        %  lzimp(~isfinite(lzimp(:))) = 0;
        %  inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        %  old_value    = 1e9 .* 2.5e5 .* (zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        
        % calcul avec la nouvelle version de Lz provenant de la base de donnes ADAS
        dz     = abs(uz - 1);
        indimp =  min(find(dz == min(dz)));
        lzimp    = post(indimp).lz;
        tzimp    = post(indimp).te;
        %pour HDT -> valide ?
        lzimp    = noncoronal(option,tzimp ,lzimp);
        lzimp_coef    = cat(2,1e-38,1e-37,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-27);
        tzimp    = cat(2,1e-37,1e-4,tzimp,te_max);
        fz_bord  = max(0,profli.n1p(:,end)) ./ profli.nep(:,end) .* 100;
        fz  = fz_bord ./ zeff_div;
        lzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp_coef),log10(teloc(:) ./ 1e3)),size(teloc));
        lzimp(~isfinite(lzimp(:))) = 0;
        inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzimp,2);
        delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);
        
        switch option.gaz
            case 11
                dz     = abs(uz - 5);
                indB =  min(find(dz == min(dz)));
                lzB    = post(indB).lz;
                tzB    = post(indB).te;
                %pour HDT -> valide ?
                lzB    = noncoronal(option,tzB ,lzB);
                lzB_coef    = cat(2,1e-38,1e-37,lzB,lzB(end) + 5.355e3 .* uz(indB)  .*(sqrt(te_max)-sqrt(tzB(end))) .* 1e-27);
                tzB    = cat(2,1e-37,1e-4,tzB,te_max);
                fz_bord  = (max(0,zs.nTm)./ max(1e13,zs.n1m)) .* max(0,profli.n1p(:,end)) ./ profli.nep(:,end) .* 100;
                fz  = fz_bord ./ zeff_div;
                lzB   = reshape(10 .^ pchip(log10(tzB),log10(lzB_coef),log10(teloc(:) ./ 1e3)),size(teloc));
                lzB(~isfinite(lzB(:))) = 0;
                inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzB,2);
                delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div);               
        end
        
        %figure(21);clf;plot(cons.temps,1e9 .* 2.5e5 .* (zs.nebord ./ 1e20) .* sqrt(fz) .* sqrt(tebord .^ 2 .* inte_rad_div),'r',cons.temps,old_value,'b');drawnow
        
        % calcul de la contribution du tungsten si besoin
        if option.W_effect == 1
            
            indw  = find(uz == 74,1);
            lzw    = post(indw).lz;
            tzw    = post(indw).te;
            lzw    = noncoronal(option,tzw,lzw);
            lzw_coef    = cat(2,1e-37,lzw,lzw(end) + 5.355e3 .* uz(indw)  .*(sqrt(te_max)-sqrt(tzw(end))) .* 1e-27);
            tzw    = cat(2,1e-4,tzw,te_max);
            fz_w   = (1 -Sn_fraction) .* ya0_w .* 100 ./ zeff_div;
            lzw    = reshape(10 .^ pchip(log10(tzw),log10(lzw_coef),log10(teloc(:) ./ 1e3)),size(teloc));
            lzw(~isfinite(lzw(:))) = 0;
            inte_rad_div = tescale .* trapz(tepar,sqrt(teloc) .* lzw,2);
            delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz_w) .* sqrt(tebord .^ 2 .* inte_rad_div);
            if Sn_fraction > 0
                indsn  = find(uz == 50,1);
                lzsn    = post(indsn).lz;
                tzsn   = post(indsn).te;
                lzsn    = noncoronal(option,tzsn,lzsn);
                lzsn_coef    = cat(2,1e-37,lzsn,lzsn(end) + 5.355e3 .* uz(indsn)  .*(sqrt(te_max)-sqrt(tzsn(end))) .* 1e-27);
                tzsn    = cat(2,1e-4,tzsn,te_max);
                fz_sn   = Sn_fraction .* ya0_w .* 100 ./ zeff_div;
                lzsn   = reshape(10 .^ pchip(log10(tzsn),log10(lzsn_coef),log10(teloc(:) ./ 1e3)),size(teloc));
                lzsn(~isfinite(lzsn(:))) = 0;
                inte_rad_div = (1 -Sn_fraction) .* tescale .* trapz(tepar,sqrt(teloc) .* lzsn,2);
                delta_qpl_rad = delta_qpl_rad + 1e9 .* 2.5e5 .* (fact_ne .* zs.nebord ./ 1e20) .* sqrt(fz_sn) .* sqrt(tebord .^ 2 .* inte_rad_div);                
            end
        end
    else
        delta_qpl_rad = 0 .* zs.nebord;
    end
    
    % selon le recyclage au niveau des plaques de divertor
    if  Recycling == 0
        qpl_corr   = qpl_in_loc - delta_qpl_rad;
        delta_qpl_neutral = 0 .* zs.nebord;
    else
        if option.detach ~= 0
            % adding Apiwat/Pegourie aditionnal term
            c_st   = sqrt(2 .* phys.e .* telim ./ zs.meff ./ phys.ua);
            qpl_cx = phys.e .* gamma_e .* c_st .* nelim .* telim .*  ...
                ((zs.nebord .* tebord ./ 2  - nelim .* telim) ./max(1,nelim .* telim)) .^ option.detach;
            qpl_cx  = max(0,real(qpl_cx));
            qpl_cx = min(max(0,qpl_in_loc - delta_qpl_rad),qpl_cx);
            %figure(37);plot(cons.temps,qpl_cx,cons.temps,qpl_in_loc - delta_qpl_rad);hold on; drawnow
        else
            qpl_cx = 0 .* zs.nebord;
        end
        
        % recyclage au niveau des plaques de diverteur
        % reference : P.C. Stangeby NF 51 (2011) 063001
        %sheath factor (include ion and electron)
        if option.eioniz > 0
            qpl_corr = min((qpl_in_loc - delta_qpl_rad) ,(qpl_in_loc - delta_qpl_rad - qpl_cx)./  ...
                (1 + Recycling .* option.eioniz ./ gamma_e ./ max(0.1,telim)));
        else
            qpl_corr = min((qpl_in_loc - delta_qpl_rad) ,(qpl_in_loc - delta_qpl_rad - qpl_cx)./  ...
                (1 + Recycling .* z0eioniz_div(telim,nelim) ./ gamma_e ./ max(0.1,telim)));
        end
        delta_qpl_neutral = qpl_in_loc - min(qpl_in_loc,qpl_corr + delta_qpl_rad);
        
    end
    
    
    % graphe de controle
    if 0
        figure(21);clf
        subplot(5,1,1)
        semilogy(cons.temps,telim,'r',cons.temps,telim_mem,'b.');
        subplot(5,1,2)
        semilogy(cons.temps,nelim,'r',cons.temps,nelim_mem,'b.');
        subplot(5,1,3)
        semilogy(cons.temps,tebord,'r',cons.temps,tebord_mem,'b.');
        subplot(5,1,4)
        plot(cons.temps,qpl,'r',cons.temps,max(0,qpl_corr),'b.',cons.temps,qpl_step,'g',cons.temps,qpl_best,'c');
        subplot(5,1,5)
        %plot(cons.temps,fmom,'r',cons.temps,fmom_corr,'b.',cons.temps,fmom_step,'g',cons.temps,fmom_best,'c');
        drawnow
    end
    
    err_loc = (qpl - qpl_corr) .^ 2 ./ max(1,qpl_in) .^ 2;
    
    indbest = find(err_loc < err_best);
    if ~isempty(indbest)
        err_best(indbest)  = err_loc(indbest);
        qpl_best(indbest)  = qpl(indbest);
        fmom_best(indbest) = fmom(indbest);
        delta_qpl_rad_best(indbest) = delta_qpl_rad(indbest);
        delta_qpl_neutral_best(indbest) = delta_qpl_neutral(indbest);
        qpl_in_loc_best(indbest) =  qpl_in_loc(indbest);
        zeff_div_best(indbest) = zeff_div(indbest);
    end
    
    if  (max(qpl_step ./ max(1 , qpl_in)) <= sqrt(eps))
        break
    elseif isfield(option,'plot2points')
        break
    end
    
    % evolution
    if loop == 0
        qpl = min(qpl_in,max(qpl_in .* sqrt(eps),qpl + qpl_step .* sign(qpl_corr - qpl)));
        qpl_step = qpl_step ./ sqrt(2);
    elseif loop == 1
        qpl = min(qpl_in,max(qpl_in .* sqrt(eps),qpl - qpl_step / nbmax));
    elseif loop == 2
        qpl = min(qpl_in,max(qpl_in .* sqrt(eps),qpl_in .* rand(size(qpl_in))));
    else
        indbest = find(isfinite(err_best));
        qpl(indbest) = qpl_best(indbest);
        qpl = min(qpl_in,max(qpl_in .* sqrt(eps),qpl + qpl_step .* (rand(size(qpl_in)) - (0.5 - eps))));
        qpl_step = qpl_step .* 0.95;
    end
    % comptage
        nb    = nb - 1;
%         rep = whos;
%         if any(cat(1,rep(:).complex))
%             keyboard
%         end
	
end		   

err = max(sqrt(err_loc));
indnon_best = find((err_loc > err_stop) & (err_loc > err_best));
if ~isempty(indnon_best) && ~isfield(option,'plot2points')

	qpl(indnon_best)  = qpl_best(indnon_best);
	fmom(indnon_best) = fmom_best(indnon_best);
	delta_qpl_rad(indnon_best) = delta_qpl_rad_best(indnon_best);
	delta_qpl_neutral(indnon_best) = delta_qpl_neutral_best(indnon_best);
        qpl_in_loc(indnon_best) = qpl_in_loc_best(indnon_best);
        zeff_div(indnon_best) = zeff_div_best(indnon_best);

	% recalcul final
 	if ~isfield(option,'plot2points')
       		[tebord,telim,nelim,indbad,noconv] = z0twopoints(qpl,fact_ne .* zs.nebord,lc,zs.zeff,zs.meff,abs(option.fcond) .* fpe .*  qpl_in ./ max(1,qpl),fmom .* fie,tite_loc,gamma);
		% calcul de la correction cinetique
		% ref : Stangeby ch 26 + 
		% W. Fundamenski, topical review, PPCF 47 (2005) p R163-R208, formule 10
		% alpha_e ~ 1.5 recommande mais avec tanh, on choisi alpha_e = 0.82
		if option.fcond < 0
			% point fixe
			for lcond=1:5
				qelim = option.alpha_e .* sqrt(phys.e .* telim ./ phys.me) .*  phys.e .* telim .* nelim;
				fcond = abs(option.fcond) .* fpe ./ (1 - min(0.9,tanh(qpl ./ qelim))) .*  qpl_in ./ max(1,qpl);	
				[tebord,telim,nelim,indbad,noconv] = z0twopoints(qpl,fact_ne .* zs.nebord,lc,zs.zeff,zs.meff,fcond,fmom .* fie,tite_loc,gamma);					
			end
                else
			fcond = abs(option.fcond) .* fpe .*  qpl_in ./ max(1,qpl);
		end
	else
		qelim = option.alpha_e .* sqrt(phys.e .* telim ./ phys.me) .*  phys.e .* telim .* nelim;
		fcond = abs(option.fcond) .* fpe ./ (1 - min(0.9,tanh(qpl ./ qelim))) .*  qpl_in ./ max(1,qpl);	
	
	end
	if option.mach_corr == 1
		% nombre de mach (model simple pour detachement)
		% ref JET-P(64) 05 , P. H. Harbour and A. Loarte
                if  option.eioniz > 0
			Kt          = min(0.5,max(-0.5,0.5 .* (1 - 2 .* tanh(Recycling .* option.eioniz ./ gamma_e ./ max(eps,telim)))));
                else
			Kt          = min(0.5,max(-0.5,0.5 .* (1 - 2 .* tanh(Recycling .* z0eioniz_div(telim,nelim) ./ gamma_e ./ max(eps,telim)))));
		end
		mach_target = sqrt(3 .* (1 - Kt) ./ (1 + Kt));
		
 		if ~isfield(option,'plot2points')
			% formule Stangeby 14.4-14.11
			ftm = (gamma_e .* (1 + mach_target .^ 2) ./ (2 .* abs(mach_target) .* (gamma_e - 1 + mach_target .^ 2))) .^ 2;
			telim = ftm .* telim;
			fnm   = 2 ./ (ftm .* (1 + mach_target .^ 2));
			nelim = fnm .* nelim;
			% tebord est inchange
			% le flux sur le target est change !
		end
	end
else
	indbad = [];
        noconv = 0;
end
if ~isempty(indbad)
	tebord(indbad)= 30;
	nelim(indbad) = 2 .* zs.nebord(indbad);
	telim(indbad) = 13.6;
end


% graphe de controle
if 0
	figure(22);clf
	subplot(5,1,1)
	plot(cons.temps,telim,'r',cons.temps,telim_mem,'b.');
	set(gca,'ylim',[0.1,Inf]);
	title(sprintf('nb = %d',nbmax -nb));
	subplot(5,1,2)
	plot(cons.temps,nelim,'r',cons.temps,nelim_mem,'b.');
	set(gca,'ylim',[0,1e22]);
	subplot(5,1,3)
	plot(cons.temps,tebord,'r',cons.temps,tebord_mem,'b.');
	set(gca,'ylim',[0.1,Inf]);
	subplot(5,1,4)
	plot(cons.temps,qpl,'r',cons.temps,max(0,qpl_corr),'b.',cons.temps,qpl_step,'g',cons.temps,qpl_best,'c',cons.temps,qpl_in,'k:',cons.temps,qpl_in_loc_best,'k-.');
	subplot(5,1,5)
	%plot(cons.temps,fmom,'r',cons.temps,fmom_corr,'b.',cons.temps,fmom_step,'g',cons.temps,fmom_best,'c');
	drawnow
end


%fprintf('nb = %d, err = %g, nb_bad  = %d\n',nbmax - nb, err,length(indbad));

if (noconv == 1) || ~isempty(indnon_best)
	fprintf(':');
end 
%figure(21);clf;plot(cons.temps,telim,cons.temps,tebord);hold on;drawnow
%figure(21);plot(cons.temps,telim_x,'r');hold on;drawnow


%  % diagnostic sur neu
%  neu_verif         =  (4 .* telim .* nelim) ./ (option.fmom .* fie) ./ tebord;
%  
%  if option.mach_corr == 1
%  	neu_verif = neu_verif ./ ftm ./ fnm;
%  end
%  
%  figure(21);
%  plot(cons.temps,zs.nebord,'or',cons.temps,neu_verif,'.b')
%  drawnow

% correction non coronal 
function lz = noncoronal(option,tz,lz)

if option.noncoronal > 0
    switch option.noncoronal
    case 3
      frac = 0.5;
    case 4
      frac = 0.2;
    case 5
      frac = 0.1;
    otherwise
	frac = 1;
    end
    % temperature du maximum de lz
    tz_bar = mean(tz(lz == max(lz)));
    lz_bar = mean(lz(lz == max(lz)));
    lz_mem = lz;
    lz   = lz .* (tz <= tz_bar) + (lz_bar .* frac + (1 - frac) .* lz_mem) .* (tz > tz_bar);
%      figure(21);clf
%      loglog(tz,lz_mem,'b.',tz,lz,'r',tz_bar,lz_bar,'*k');
%      pause(1)
end


% function computing concentration of impurities in divertor, knowing concentration in core plasma 
% with simple screening function
function fzmax_div = z0fleak_zmax(geo,zerod,profil0d,zmax,rimp,zeff_div,lc,residence_time)

persistent tabmat
% nouvelle donnees ADAS
if isempty(tabmat)
  try 
    load('Lz_zave.mat')   
  end
end


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

% offline
if ~isfield(zerod,'temps') || ~isfield(profil0d,'temps')
  % nothing
elseif length(zerod.temps) ~= length(profil0d.temps)
    noms = fieldnames(profil0d);
    tp   = profil0d.temps;
    for k=1:length(noms)
      var = profil0d.(noms{k});
      if size(var,1) == length(tp)
	profil0d.(noms{k}) = interp1(tp,var,zerod.temps,'nearest','extrap');
      end
    end
  
end


x      = profil0d.xli;
vt     = ones(size(profil0d.tep,1),1);
ve     = ones(size(x));
rext         = profil0d.Raxe + geo.a * x;
btor         = profil0d.fdia ./ rext;
grho         = abs((profil0d.rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profil0d.psi,0,2,2,1) ./ rext .* grho ./ (profil0d.rmx(:,end) * ve);
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));


% local ion charge
if residence_time == 0
    cs_max = sqrt(phys.e .* zerod.telim .* (zeff_div + zerod.tibord ./ zerod.tebord) ./ phys.ua ./ getafromz_imas(zmax));
    tau_resid = lc ./ cs_max;
else
    tau_resid_max = residence_time;
end
[Lz,Z_ave] = zlightznoncoronal(zerod.telim,zerod.nelim,tau_resid,zmax);
if any(~isfinite(Z_ave))
  [A,label] = getafromz_imas(zmax);
  if isfield(tabmat,label)
      tetab = tabmat.(label).data(:,1) * 1e3;
      Zavetab = tabmat.(label).data(:,3);
      Zave_alt = interp1(tetab,Zavetab,zerod.telim,'pchip',NaN);
      Zave_alt(zerod.telim < min(tetab)) = eps;
      Z_ave(~isfinite(Z_ave)) = Zave_alt(~isfinite(Z_ave));      
      %figure(21);plot(zerod.telim,Zave_alt,'.r',tetab,Zavetab,'b');drawnow
  end 
end
Z_ave(~isfinite(Z_ave)) = zmax;

% longueur de ionisation pres des plaques
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(zerod.telim,zerod.telim .* zerod.tibord ./ zerod.tebord);
%% equilibre entre 1s et 2s pour le neutres de centre 
alphas = zerod.nelim .* sss ./ ( zerod.nelim .* sss + Ass);
% etat d'equilibre 1s/2s
sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
s3  = sviss  + sii + svcx;
delta_s = sqrt(2 .* phys.e .* zerod.telim  ./ zerod.meff ./ phys.mp) ./ zerod.nelim ./  s3 ./ ut;
% zmax may be a problem in SOL
tleak      = 3.8e-9 .* Z_ave .* sqrt(zerod.nelim .* delta_s);
% inverse of leaking
fleak         = min(1,exp(-(tleak ./ zerod.telim) .^ 2));
fz_bord       = rimp .* profil0d.nzp(:,end) ./ profil0d.nep(:,end) .* 100;
fzmax_div     = min(100 ./ zmax ./ 2,fz_bord ./ fleak);

%tv = 1:length(fzmax_div);
%figure(21);plot(tv,fzmax_div,'r',tv,fz_bord,'b');drawnow
