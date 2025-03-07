%
% same species order, but #9 is electron
%
function [psupra,psupra_para,psupra_perp,nfast] = nfast4imas(zerod,profil0d,option,frhe3,geo,cons)

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


% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = option.gaz; 
end


if (length(zerod.temps) ~= length(profil0d.temps)) || ~all(zerod.temps == profil0d.temps)

	% resample of zerod on profil0d time slices (if needed)
	noms = fieldnames(zerod);
        temps = zerod.temps;
	for k=1:length(noms)
		var = zerod.(noms{k});
		if length(var) == length(temps)
            if all(isfinite(var))
                zerod.(noms{k})= interp1_imas(temps,var,profil0d.temps,'pchip','extrap');
            else
                zerod.(noms{k})= interp1_imas(temps,var,profil0d.temps,'nearest','extrap');                
            end
		end	
	end
end

if (length(cons.temps) ~= length(profil0d.temps)) || ~all(cons.temps == profil0d.temps)


	% resample of cons on profil0d time slices (if needed)
    temps = cons.temps;
	noms = fieldnames(cons);
	for k=1:length(noms)
		var = cons.(noms{k});
		if length(var) == length(temps)
			cons.(noms{k})= interp1_imas(temps,var,profil0d.temps,'pchip','extrap');
		end	
	end

	% resample of geo on profil0d time slices (if needed)
	noms = fieldnames(geo);
	for k=1:length(noms)
		var = geo.(noms{k});
		if length(var) == length(temps)
			geo.(noms{k})= interp1_imas(temps,var,profil0d.temps,'pchip','extrap');
		end	
	end

end

% vector unity
ve = ones(size(profil0d.xli));

% ATTENTION : nromalisation des mass sur le gaz principal 
% density of fast ions
% [n_fast,source_fast,emoy_fast,tau_nfast,ecrit] =zerod_fast_ions_density(nep,tep,zeff,meff,Afast,Zfast,Einj,p_supra)
% fast alpha
psup_alpha  = profil0d.pfus;
psup_alpha  = psup_alpha ./ (max(1, trapz(profil0d.xli,psup_alpha .* profil0d.vpr,2)) * ve) .*  (zerod.esup_fus * ve) ./ (3/2);
nfast_alpha = zerod_fast_ions_density(profil0d.nep,profil0d.tep,profil0d.zeff,zerod.meff * ve,4,2,3.56e6,psup_alpha);
%talpha      = psup_alpha ./ max(1,nfast_alpha) ./  phys.e;
psup_perp_alpha = 2/3 .* psup_alpha;
psup_para_alpha = 1/3 .* psup_alpha;


% NBI
% le (2/3) vient des observation sur JET , la pente dans le tanh doit etre ajustee plus finement
frpar_nbi   = (1/3) + (2/3) .* abs(real(zerod.mu0_nbi)) .*  max(0,tanh(option.einj ./ max(eps,real(zerod.ecrit_nbi)))) + ...
              sqrt(-1) .* ((1/3) + (2/3) .* abs(imag(zerod.mu0_nbi)) .*  max(0,tanh(option.einj2 ./ max(eps,imag(zerod.ecrit_nbi)))));
frpar_nbi   = frpar_nbi * ve;

% first injecteur
switch gas_nbi
    case -1
        minj =  1 .* real(cons.ftnbi) * ve  + sqrt(-1) .* (1 .* imag(cons.ftnbi) * ve);
    case 11
        minj = 1 .* (1-cons.ftnbi) + 11 .* cons.ftnbi;
    case {3,5}
        minj = 2 .* (1-real(cons.ftnbi) * ve ) + 3 .* real(cons.ftnbi) * ve  + ...
            sqrt(-1) .* (2 .* (1-imag(cons.ftnbi) * ve ) + 3 .* imag(cons.ftnbi) * ve);
    otherwise
        minj = 2 .* (1-real(cons.ftnbi) * ve ) + 1 .* real(cons.ftnbi) * ve  + ...
            sqrt(-1) .* (2 .* (1-imag(cons.ftnbi) * ve ) + 1 .* imag(cons.ftnbi) * ve);
end
%
psup_nbi1  = real(profil0d.nbinesource);
psup_nbi1  = psup_nbi1 ./ (max(1, trapz(profil0d.xli,psup_nbi1 .* profil0d.vpr,2)) * ve) .*  (real(zerod.esup_nbi) * ve) ./ (3/2);
nfast_nbi1 = zerod_fast_ions_density(profil0d.nep,profil0d.tep,profil0d.zeff,zerod.meff * ve,real(minj),1,option.einj,psup_nbi1);
%tnbi1      = psup_nbi1  ./ max(1,nfast_nbi1) ./  phys.e;
psup_perp_nbi1 = (1 - real(frpar_nbi)) .* psup_nbi1;
psup_para_nbi1 = real(frpar_nbi) .* psup_nbi1;
%
psup_nbi2  = imag(profil0d.nbinesource);
psup_nbi2  = psup_nbi2 ./ (max(1, trapz(profil0d.xli,psup_nbi2 .* profil0d.vpr,2)) * ve) .*  (imag(zerod.esup_nbi) * ve) ./ (3/2);
nfast_nbi2 = zerod_fast_ions_density(profil0d.nep,profil0d.tep,profil0d.zeff,zerod.meff * ve,imag(minj),1,option.einj2,psup_nbi2);
%tnbi2      = psup_nbi2 ./ max(1,nfast_nbi2) ./  phys.e;
psup_perp_nbi2 = (1 - imag(frpar_nbi)) .* psup_nbi2;
psup_para_nbi2 = imag(frpar_nbi) .* psup_nbi2;

% choix du minoritaire
switch option.mino
    case 'He3'
        ag = 3;
        zg = 2;
        lg = 7.92e-3;
    case 'T'
        ag = 3;
        zg = 1;
        lg = 7.92e-3;
    case 'He4'
        ag = 4;
        zg = 2;
        lg = 4.55e-3;
    case 'D'
        ag = 2;
        zg = 1;
        lg = 6.46e-3;
    case 'B'
        ag = 11;
        zg = 5;
        lg = 4.576e-3 .* sqrt(11) / 5;
    otherwise
        ag = 1;
        zg = 1;
        lg = 4.576e-3;
end
psup_icrh   = profil0d.picrh;
psup_icrh   = psup_icrh ./ (max(1, trapz(profil0d.xli,psup_icrh .* profil0d.vpr,2)) * ve) .*  (zerod.esup_icrh * ve ) ./ (3/2);
nfast_icrh  = zerod_fast_ions_density(profil0d.nep,profil0d.tep,profil0d.zeff,zerod.meff * ve ,ag,zg,zerod.einj_icrh * ve ,psup_icrh);
%ticrh       = psup_icrh ./ max(1,nfast_icrh) ./  phys.e;
% this law is made to mimic PION results for JET [ref? Lars ?]
frpar_icrh  = (1/3) + (1/3) .*  max(0,tanh(zerod.einj_icrh ./ zerod.ecrit_icrh)) * ve;
psup_perp_icrh = (1 - frpar_icrh) .* psup_icrh;
psup_para_icrh = frpar_icrh .* psup_icrh;

warning off
% pas de meilleurs formulation pour le moment
if option.lhmode == 5
	nfast_lh  = zeros(size(profil0d.plh));
else
	nfast_lh  = (max(0,zerod.esup_lh) ./ (zerod.einj_lh .* phys.e) .* 2 ./ zerod.vp ) * ve;
end
nfast_lh(~isfinite(nfast_lh)) = 0;
warning on
psup_lh   = profil0d.plh;
nfast_lh  = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr,2)) * ve) .*  nfast_lh;
psup_lh   = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr,2)) * ve) .*  (zerod.esup_lh * ve ) ./ (3/2);
%tlh       = psup_lh ./ max(1,nfast_lh) ./  phys.e;
% le (1/3) pour LH provient d'un calcul DKE fait par J.Decker pour TS.	  
psup_perp_lh = (1/3) .* psup_lh;
psup_para_lh = (2/3) .* psup_lh;

% fast electron due to runaway electrons addedt to LH terms
if option.runaway ~= 0
    % assuming v// ~ c
    nfast_lh     = nfast_lh  + profil0d.jrun ./ phys.c ./ phys.e;
    psup_para_lh = psup_para_lh + profil0d.jrun .* phys.me .* phys.c ./ phys.e;
    psup_lh      = psup_lh + profil0d.jrun .* phys.me .* phys.c ./ phys.e;
end

% hydrogenoid density
nDp   = max(0,profil0d.n1p .* (zerod.nDm * ve ) ./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2) * ve) .* (trapz(profil0d.xli,profil0d.vpr,2) * ve));
nTp   = max(0,profil0d.n1p .* (zerod.nTm * ve ) ./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2) * ve) .* (trapz(profil0d.xli,profil0d.vpr,2) * ve));
switch option.gaz
    case 11
        nbp   = nTp;
        nTp   = zeros(size(nTp));
        nHp   = max(0,profil0d.n1p - nDp);
    otherwise
        nbp   = zeros(size(nTp));
        nHp   = max(0,profil0d.n1p - nTp - nDp);
end
% helium
switch option.gaz
    case 5
        nHep3  = profil0d.nhep;
        nHep   = option.frhe0 .* profil0d.nep;
    otherwise
        nHep  = max(1e13,profil0d.nhep .* (1 - frhe3));
        nHep3 = max(1e13,profil0d.nhep .* frhe3);
end
%nHep  = profil0d.nhep;

% compute corrected density depending on minority scheme
switch option.mino
    case 'H'
        nHp_mem = nHp;
        nHp = max(0, nHp - nfast_icrh);
        % security
        nfast_icrh = nHp_mem - nHp;
        
    case 'D'
        nDp_mem = nDp;
        nDp = max(0, nDp - nfast_icrh);
        % security
        nfast_icrh = nDp_mem - nDp;
        
    case 'T'
        nTp_mem = nTp;
        nTp = max(0, nTp - nfast_icrh);
        % security
        nfast_icrh = nTp_mem - nTp;
        
    case 'He3'
        nHep3_mem = nHep3;
        nHe3 = max(0, nHep3 - nfast_icrh);
        % security
        nfast_icrh = nHep3_mem - nHep3;
        
    case 'He4'
        nHep_mem = nHep;
        nHep = max(0, nHep - nfast_icrh);
        % security
        nfast_icrh = nHep_mem - nHep;
        
    case 'B'
    otherwise
        error(sprintf('METIS4IMAS call: minority species %s not yet implemanted',option.mino));
end

% correction of fast NBI ions
switch gas_nbi
    case -1
        nHp_mem = nHp;
        nDp_mem = nDp;
        nHp = max(0, nHp - nfast_nbi1);
        nDp = max(0, nDp);
        % security
        nfast_nbi1 = nHp_mem - nHp + nDp_mem - nDp;
        nHp_mem = nHp;
        nDp_mem = nDp;
        nHp = max(0, nHp - nfast_nbi2);
        nDp = max(0, nDp);
        % security
        nfast_nbi2 = nHp_mem - nHp + nDp_mem - nDp;
        
    case 3
        nTp_mem = nTp;
        nDp_mem = nDp;
        nTp = max(0, nTp - (real(cons.ftnbi) * ve ) .* nfast_nbi1);
        nDp = max(0, nDp - (1 - real(cons.ftnbi) * ve ) .* nfast_nbi1);
        % security
        nfast_nbi1 = nTp_mem - nTp + nDp_mem - nDp;
        nTp_mem = nTp;
        nDp_mem = nDp;
        nTp = max(0, nTp - (imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        nDp = max(0, nDp - (1 - imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        % security
        nfast_nbi2 = nTp_mem - nTp + nDp_mem - nDp;
        
    case 5
        nHep3_mem = nHep3;
        nDp_mem   = nDp;
        nHep3 = max(0, nHep3 - (real(cons.ftnbi) * ve ) .* nfast_nbi1);
        nDp = max(0, nDp - (1 - real(cons.ftnbi) * ve ) .* nfast_nbi1);
        % security
        nfast_nbi1 = nHep3_mem - nHep3 + nDp_mem - nDp;
        nHep3_mem = nHep3;
        nDp_mem = nDp;
        nHep3 = max(0, nHep3 - (imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        nDp = max(0, nDp - (1 - imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        % security
        nfast_nbi2 = nHep3_mem - nHep3 + nDp_mem - nDp;       
    case 11
        nHp_mem = nHp;
        nBp_mem = nBp;
        nBp = max(0, nBp - (real(cons.ftnbi) * ve ) .* nfast_nbi1);
        nHp = max(0, nHp - (1 - real(cons.ftnbi) * ve ) .* nfast_nbi1);
        % security
        nfast_nbi1 = nHp_mem - nHp + nBp_mem - nBp;
        nHp_mem = nHp;
        nBp_mem = nBp;
        nBp = max(0, nBp - (imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        nHp = max(0, nHp - (1 - imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        % security
        nfast_nbi2 = nHp_mem - nHp + nBp_mem - nBp;
        
    otherwise
        
        nHp_mem = nHp;
        nDp_mem = nDp;
        nHp = max(0, nHp - (real(cons.ftnbi) * ve ) .* nfast_nbi1);
        nDp = max(0, nDp - (1 - real(cons.ftnbi) * ve ) .* nfast_nbi1);
        % security
        nfast_nbi1 = nHp_mem - nHp + nDp_mem - nDp;
        nHp_mem = nHp;
        nDp_mem = nDp;
        nHp = max(0, nHp - (imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        nDp = max(0, nDp - (1 - imag(cons.ftnbi) * ve ) .* nfast_nbi2);
        % security
        nfast_nbi2 = nHp_mem - nHp + nDp_mem - nDp;
        
end

% correction of fast alpha
nHep_mem = nHep;
nHep = max(0, nHep - nfast_alpha);
% security
nfast_alpha = nHep_mem - nHep;

% correction of fast electron
nfast_lh = profil0d.nep - max(0,profil0d.nep - nfast_lh);

% ============================
% now all densities are correct.
% ============================

% fast ions NBI
switch gas_nbi
    case -1
        nfast_T = zeros(size(profil0d.nep));
        nfast_H = nfast_nbi1  + nfast_nbi2;
        nfast_D = zeros(size(profil0d.nep));
        nfast_B = zeros(size(profil0d.nep));
        nfast_He3 = zeros(size(profil0d.nep));
        
        psup_T = zeros(size(profil0d.nep));
        psup_H = psup_nbi1 + psup_nbi2 ;
        psup_D = zeros(size(profil0d.nep));
        psup_B = zeros(size(profil0d.nep));
        psup_He3 = zeros(size(profil0d.nep));
       
        psup_para_T = zeros(size(profil0d.nep));
        psup_para_H = psup_para_nbi1 + psup_para_nbi2;
        psup_para_D = zeros(size(profil0d.nep));
        psup_para_B = zeros(size(profil0d.nep));
        psup_para_He3 = zeros(size(profil0d.nep));
        
        psup_perp_T = zeros(size(profil0d.nep));
        psup_perp_H = psup_perp_nbi1 + psup_perp_nbi2;
        psup_perp_D = zeros(size(profil0d.nep));
        psup_perp_B = zeros(size(profil0d.nep));
        psup_perp_He3 = zeros(size(profil0d.nep));
        
    case 3
        nfast_H = zeros(size(profil0d.nep));
        nfast_T = nfast_nbi1  .* real(cons.ftnbi * ve) + nfast_nbi2 .* imag(cons.ftnbi * ve);
        nfast_D = nfast_nbi1  .* (1 - real(cons.ftnbi * ve))  + nfast_nbi2.* (1 - imag(cons.ftnbi * ve ));
        nfast_B = zeros(size(profil0d.nep));
        nfast_He3 = zeros(size(profil0d.nep));
        
        psup_H = zeros(size(profil0d.nep));
        psup_T = psup_nbi1  .* real(cons.ftnbi * ve) + psup_nbi2 .* imag(cons.ftnbi * ve);
        psup_D = psup_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_B = zeros(size(profil0d.nep));
        psup_He3 = zeros(size(profil0d.nep));
        
        psup_para_H = zeros(size(profil0d.nep));
        psup_para_T = psup_para_nbi1  .* real(cons.ftnbi * ve) + psup_para_nbi2 .* imag(cons.ftnbi * ve);
        psup_para_D = psup_para_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_para_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_para_B = zeros(size(profil0d.nep));
        psup_para_He3 = zeros(size(profil0d.nep));
        
        psup_perp_H = zeros(size(profil0d.nep));
        psup_perp_T = psup_perp_nbi1  .* real(cons.ftnbi * ve) + psup_perp_nbi2 .* imag(cons.ftnbi * ve);
        psup_perp_D = psup_perp_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_perp_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_perp_B = zeros(size(profil0d.nep));
        psup_perp_He3 = zeros(size(profil0d.nep));
        
    case 5
        nfast_H = zeros(size(profil0d.nep));
        nfast_He3 = nfast_nbi1  .* real(cons.ftnbi * ve) + nfast_nbi2 .* imag(cons.ftnbi * ve);
        nfast_D = nfast_nbi1  .* (1 - real(cons.ftnbi * ve))  + nfast_nbi2.* (1 - imag(cons.ftnbi * ve ));
        nfast_B = zeros(size(profil0d.nep));
        nfast_T = zeros(size(profil0d.nep));
        
        psup_H = zeros(size(profil0d.nep));
        psup_He3 = psup_nbi1  .* real(cons.ftnbi * ve) + psup_nbi2 .* imag(cons.ftnbi * ve);
        psup_D = psup_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_B = zeros(size(profil0d.nep));
        psup_T = zeros(size(profil0d.nep));
        
        psup_para_H = zeros(size(profil0d.nep));
        psup_para_He3 = psup_para_nbi1  .* real(cons.ftnbi * ve) + psup_para_nbi2 .* imag(cons.ftnbi * ve);
        psup_para_D = psup_para_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_para_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_para_B = zeros(size(profil0d.nep));
        psup_para_T = zeros(size(profil0d.nep));
        
        psup_perp_H = zeros(size(profil0d.nep));
        psup_perp_He3 = psup_perp_nbi1  .* real(cons.ftnbi * ve) + psup_perp_nbi2 .* imag(cons.ftnbi * ve);
        psup_perp_D = psup_perp_nbi1  .* (1 - real(cons.ftnbi * ve))  + psup_perp_nbi2.* (1 - imag(cons.ftnbi * ve ));
        psup_perp_B = zeros(size(profil0d.nep));
        psup_perp_T = zeros(size(profil0d.nep));
        
    case 11
        nfast_T = zeros(size(profil0d.nep));
        nfast_B = nfast_nbi1  .* real(cons.ftnbi * ve) + nfast_nbi2 .* imag(cons.ftnbi * ve) ;
        nfast_H = nfast_nbi1  .* (1 - real(cons.ftnbi * ve) ) + nfast_nbi2 .* (1 - imag(cons.ftnbi * ve));
        nfast_D = zeros(size(profil0d.nep));
        nfast_He3 = zeros(size(profil0d.nep));
        
        psup_T = zeros(size(profil0d.nep));
        psup_B = psup_nbi1  .* real(cons.ftnbi * ve) + psup_nbi2 .* imag(cons.ftnbi * ve);
        psup_H = psup_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_D = zeros(size(profil0d.nep));
        psup_He3 = zeros(size(profil0d.nep));
        
        psup_para_T = zeros(size(profil0d.nep));
        psup_para_B = psup_para_nbi1  .* real(cons.ftnbi * ve) + psup_para_nbi2 .* imag(cons.ftnbi * ve);
        psup_para_H = psup_para_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_para_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_para_D = zeros(size(profil0d.nep));
        psup_para_He3 = zeros(size(profil0d.nep));
        
        psup_perp_T = zeros(size(profil0d.nep));
        psup_perp_B = psup_perp_nbi1  .* real(cons.ftnbi * ve) + psup_perp_nbi2 .* imag(cons.ftnbi * ve);
        psup_perp_H = psup_perp_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_perp_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_perp_D = zeros(size(profil0d.nep));
        psup_perp_He3 = zeros(size(profil0d.nep));

        
    otherwise
        nfast_T = zeros(size(profil0d.nep));
        nfast_H = nfast_nbi1  .* real(cons.ftnbi * ve) + nfast_nbi2 .* imag(cons.ftnbi * ve) ;
        nfast_D = nfast_nbi1  .* (1 - real(cons.ftnbi * ve) ) + nfast_nbi2 .* (1 - imag(cons.ftnbi * ve));
        nfast_B = zeros(size(profil0d.nep));
        nfast_He3 = zeros(size(profil0d.nep));
        
        psup_T = zeros(size(profil0d.nep));
        psup_H = psup_nbi1  .* real(cons.ftnbi * ve) + psup_nbi2 .* imag(cons.ftnbi * ve);
        psup_D = psup_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_B = zeros(size(profil0d.nep));
        psup_He3 = zeros(size(profil0d.nep));
        
        psup_para_T = zeros(size(profil0d.nep));
        psup_para_H = psup_para_nbi1  .* real(cons.ftnbi * ve) + psup_para_nbi2 .* imag(cons.ftnbi * ve);
        psup_para_D = psup_para_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_para_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_para_B = zeros(size(profil0d.nep));
        psup_para_He3 = zeros(size(profil0d.nep));
        
        psup_perp_T = zeros(size(profil0d.nep));
        psup_perp_H = psup_perp_nbi1  .* real(cons.ftnbi * ve) + psup_perp_nbi2 .* imag(cons.ftnbi * ve);
        psup_perp_D = psup_perp_nbi1  .* (1 - real(cons.ftnbi * ve) ) + psup_perp_nbi2 .* (1 - imag(cons.ftnbi * ve));
        psup_perp_B = zeros(size(profil0d.nep));
        psup_perp_He3 = zeros(size(profil0d.nep));
        
end

% ICRH contribution
%nfast_He3     =  zeros(size(nfast_H));
%psup_He3      =  zeros(size(nfast_H));
%psup_para_He3 =  zeros(size(nfast_H));
%psup_perp_He3 =  zeros(size(nfast_H));
switch option.mino
case 'H'
  nfast_H = nfast_H + nfast_icrh;  
  psup_H = psup_H + psup_icrh;  
  psup_para_H = psup_para_H + psup_para_icrh;  
  psup_perp_H = psup_perp_H + psup_perp_icrh;  
case 'D'
  nfast_D = nfast_D + nfast_icrh;  
  psup_D = psup_D + psup_icrh;  
  psup_para_D = psup_para_D + psup_para_icrh;  
  psup_perp_D = psup_perp_D + psup_perp_icrh;  
case 'T'
  nfast_T = nfast_T + nfast_icrh;  
  psup_T = psup_T + psup_icrh;  
  psup_para_T = psup_para_T + psup_para_icrh;  
  psup_perp_T = psup_perp_T + psup_perp_icrh;  
case 'He3'
  nfast_He3 = nfast_He3 + nfast_icrh; 
  psup_He3  = psup_He3 + psup_icrh;  
  psup_para_He3 = psup_para_He3 + psup_para_icrh;  
  psup_perp_He3 = psup_perp_He3 + psup_perp_icrh;  
case 'He4'
  nfast_alpha = nfast_alpha + nfast_icrh;
  psup_alpha = psup_alpha + psup_alpha;  
  psup_para_alpha = psup_para_alpha + psup_para_alpha;  
  psup_perp_alpha = psup_perp_alpha + psup_perp_alpha;  
case 'B'
  nfast_B = nfast_B + nfast_icrh;  
  psup_B = psup_B + psup_icrh;  
  psup_para_B = psup_para_B + psup_para_icrh;  
  psup_perp_B = psup_perp_B + psup_perp_icrh;  
otherwise
  error(sprintf('METIS4IMAS call: minority species %s not yet implemanted',option.mino));
end

% remplissage des tableaux
psupra        = zeros(length(profil0d.temps),length(profil0d.xli),10);
psupra_para   = zeros(length(profil0d.temps),length(profil0d.xli),10);
psupra_perp   = zeros(length(profil0d.temps),length(profil0d.xli),10);
nfast         = zeros(length(profil0d.temps),length(profil0d.xli),10);

% order = {'H','D','T','He3','He4',Imp_1,Imp_2,'W','e-'}
psupra(:,:,1) = psup_H;
psupra(:,:,2) = psup_D;
psupra(:,:,3) = psup_T;
psupra(:,:,4) = psup_He3;
psupra(:,:,5) = psup_alpha;
psupra(:,:,9) = psup_lh;
psupra(:,:,10) = psup_B;
%
psupra_para(:,:,1) = psup_para_H;
psupra_para(:,:,2) = psup_para_D;
psupra_para(:,:,3) = psup_para_T;
psupra_para(:,:,4) = psup_para_He3;
psupra_para(:,:,5) = psup_para_alpha;
psupra_para(:,:,9) = psup_para_lh;
psupra_para(:,:,10) = psup_para_B;
%
psupra_perp(:,:,1) = psup_perp_H;
psupra_perp(:,:,2) = psup_perp_D;
psupra_perp(:,:,3) = psup_perp_T;
psupra_perp(:,:,4) = psup_perp_He3;
psupra_perp(:,:,5) = psup_perp_alpha;
psupra_perp(:,:,9) = psup_perp_lh;
psupra_perp(:,:,10) = psup_perp_B;
%
nfast(:,:,1) = nfast_H;
nfast(:,:,2) = nfast_D;
nfast(:,:,3) = nfast_T;
nfast(:,:,4) = nfast_He3;
nfast(:,:,5) = nfast_alpha;
nfast(:,:,9) = nfast_lh;
nfast(:,:,10) = nfast_B;

