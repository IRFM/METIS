function mat = metis_gaz_imas(z0dstruct,data_zerod,profil0d)


% compatibilite
if ~isfield(z0dstruct,'profil')
    z0dstruct.profil =  z0dstruct.profil0d;
end
% script d'estimation des flux de matiere
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
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3

% compatibilite ascendante
if ~isfield(z0dstruct.zerod,'nbar_nat') || all(~isfinite(z0dstruct.zerod.nbar_nat))
    ulh = 0.25;
    fh  = z0dstruct.zerod.modeh;
    z0dstruct.zerod.nbar_nat   = min(z0dstruct.zerod.negr,max(1e13, 1e20 .* (z0dstruct.z0dinput.geo.b0 ./ z0dstruct.zerod.q95 ./ z0dstruct.z0dinput.geo.R) .^ 0.6  .*  (ulh + (1 - ulh) .* fh)));
    data_zerod.nbar_nat = inpterp1(z0dstruct.zerod.temps,z0dstruct.zerod.nbar_nat,data_zerod.temps);
    z0dstruct.zerod.frac_pellet = z0dstruct.zerod.frac_pellet .* (z0dstruct.zerod.frac_pellet > 1e-2);
    data_zerod.frac_pellet = inpterp1(z0dstruct.zerod.temps,z0dstruct.zerod.frac_pellet,data_zerod.temps);
end
pellet_fraction = data_zerod.frac_pellet;
flag_compute = 'new';

switch flag_compute
    
    case 'new'
        
        % choix du coefficiet de recyclage
        eta_p = 0.5;
        filter_width = 11;
        Recycling = z0dstruct.z0dinput.option.Recycling;
        eta_g_default = 0.1;
        
        % temps
        mat.temps        = data_zerod.temps;
        % cette formule n'est pas precise numeriquement
        % flux de matiere sortant
        mat.output      = profil0d.ge(:,end) .* profil0d.grho2(:,end) .* profil0d.vpr_tor(:,end) ./ (4.41e-4 .* phys.avo);
        if length(data_zerod.temps) > 1
            mat.output      = interp1_imas(profil0d.temps,mat.output,data_zerod.temps,'linear','extrap');
        end
        % bilan contenu du plasma
        ntot              = trapz(z0dstruct.profil.xli,z0dstruct.profil.vpr .* z0dstruct.profil.nep,2) ./ (4.41e-4 .* phys.avo);
        mat.contents      = interp1_imas(z0dstruct.profil.temps,ntot,mat.temps,'linear','extrap');
        mat.plasmanetflux = interp1_imas(z0dstruct.profil.temps,z0dxdt(ntot,z0dstruct.profil.temps),data_zerod.temps,'linear','extrap'); % flux net qui sort du plasam
        
        % source de recylcage dans le plasma (a partir de la source de neutre qui entre dans le plasam)
        mat.s0_in       = data_zerod.n0a ./ (4.41e-4 .* phys.avo);  % recyclage  + gas puff
        
        % effect a recycling fraction (at LCFS not in the divertor)
        % sorte de minimum pour le flux de gaz qui circule dans la SOL
        switch z0dstruct.z0dinput.option.configuration
            case {0,1}
                mat.s0       = mat.s0_in  ./ max(eps,z0dstruct.z0dinput.option.fn0a);
            case {2,3}
                mat.s0       = mat.s0_in  ./ max(eps,data_zerod.xpoint .* z0dstruct.z0dinput.option.fn0a_div + (~data_zerod.xpoint) .* z0dstruct.z0dinput.option.fn0a);
            otherwise
                mat.s0       = mat.s0_in  ./ max(eps,z0dstruct.z0dinput.option.fn0a_div);
        end
        % s0_star est la source qui sort du plasma
        % computation of suggested fn0a
        fn0a_div = mat.output(data_zerod.xpoint ~= 0) ./ mat.s0(data_zerod.xpoint ~= 0);
        fn0a_div = mean(fn0a_div(isfinite(fn0a_div) & (fn0a_div>0)));
        fn0a     = mat.output(data_zerod.xpoint == 0) ./ mat.s0(data_zerod.xpoint == 0);
        fn0a     = mean(fn0a(isfinite(fn0a) & (fn0a>0)));
        try
            if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
                mat.s0          = medfilt1(mat.s0,filter_width);
                mat.s0_in       = medfilt1(mat.s0_in,filter_width);
                mat.output      = medfilt1(mat.output ,filter_width);
            end
        catch
            disp('No signal toolbox licence available : filtering removed in metis_gaz_imas');
        end
        if length(data_zerod.temps) > 1
            mat.bilan_s0    = cumtrapz(data_zerod.temps,mat.s0);
            mat.bilan_s0_in = cumtrapz(data_zerod.temps,mat.s0_in);
        else
            mat.bilan_s0    = 0;
            mat.bilan_s0_in = 0;
        end
        % security
        output_mem       = mat.output;
        mat.output       = max(mat.output,1.1 .* mat.s0);
        if length(data_zerod.temps) > 1
            mat.bilan_output = cumtrapz(data_zerod.temps,mat.output);
        else
            mat.bilan_output = 0;
        end
        % source due au glacon (dans le plasma)
        mat.pellet_star   = trapz(profil0d.xli,profil0d.spellet .* profil0d.vpr,2)  ./ (4.41e-4 .* phys.avo);
        mat.pellet_star(~isfinite(mat.pellet_star)) = 0;
        if (filter_width > 1) &&  (z0dstruct.z0dinput.option.pif < 1) && (length(data_zerod.temps) > filter_width)
            try
                mat.pellet_star          = medfilt1(mat.pellet_star,filter_width);
            end
        end
        if length(data_zerod.temps) > 1
            mat.pellet_star          = interp1_imas(profil0d.temps,mat.pellet_star,data_zerod.temps,'linear','extrap');
        end
        if length(data_zerod.temps) > 1
            mat.bilan_pellet_star  = max(eps,cumtrapz(data_zerod.temps,mat.pellet_star));
        else
            mat.bilan_pellet_star  = 0;
        end
        % flux de  matiere du a l'injection de gla�on
        mat.pellet = mat.pellet_star ./ eta_p;
        mat.bilan_pellet = mat.bilan_pellet_star ./ eta_p;
        
        % source due a l'injection de neutres
        mat.snbi_star        = real(data_zerod.pnbi) ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo) + ...
            imag(data_zerod.pnbi) ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo);
        if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
            try
                mat.snbi_star          = medfilt1(mat.snbi_star,filter_width);
            end
        end
        if length(data_zerod.temps) > 1
            mat.bilan_nbi_star   = max(eps,cumtrapz(data_zerod.temps,mat.snbi_star));
        else
            mat.bilan_nbi_star   = eps;
        end
        % flux de matiere du a l'injection de neutres
        mat.snbi             = real(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(data_zerod.frnbi)) + ...
            imag(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(data_zerod.frnbi));
        if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
            try
                mat.snbi          = medfilt1(mat.snbi,filter_width);
            end
        end
        mat.perte_snbi       = real(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(data_zerod.frnbi)) .* (1 - real(data_zerod.frnbi)) + ...
            imag(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(data_zerod.frnbi)) .* (1 - imag(data_zerod.frnbi));
        if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
            try
                mat.perte_snbi          = medfilt1(mat.perte_snbi,filter_width);
            end
        end
        if length(data_zerod.temps) > 1
            mat.bilan_nbi   = max(eps,cumtrapz(data_zerod.temps,mat.snbi));
        else
            mat.bilan_nbi   = eps;
        end
        
        % calcul gaz puff (a partir du scaling pour la densite naturelle)
        fact_nbar = data_zerod.nem ./ max(1,data_zerod.nbar);
        if (z0dstruct.z0dinput.option.tauhemul < 0) && (z0dstruct.z0dinput.option.Recycling < 1)
            % modele qui prend en compte le confinement reel et le recyclage dans le divertor
            % ref Stangeby section 6.7
            % ref originale : D. Reiter et al, PPCF vol 33 (1991) p 1579-1600
            tau_ref    = data_zerod.tauhe - z0dstruct.z0dinput.option.Recycling ./ (1 - z0dstruct.z0dinput.option.Recycling) .* data_zerod.taup;
        else
            tau_ref    = data_zerod.tauhe;
        end
        % gaz fuelling effciency
        mat.eta_g_low  =  data_zerod.taup ./ tau_ref;
        % source correspondant a la difference entre la densit� et la densit� naturelle du plasma
        mat.dfuelling_bilan_dt  = fact_nbar .* (data_zerod.nbar - data_zerod.nbar_nat) ./ tau_ref .* data_zerod.vp ./ (4.41e-4 .* phys.avo) - ...
            (mat.snbi_star + mat.pellet_star);
        if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
            try
                mat.dfuelling_bilan_dt     = medfilt1(mat.dfuelling_bilan_dt,filter_width);
            end
        end
        if length(data_zerod.temps) > 1
            mat.fuelling_bilan         = cumtrapz(data_zerod.temps,mat.dfuelling_bilan_dt);
        else
            mat.fuelling_bilan        = 0;
        end
        % source de gaz puff dans le plasma si efficacite de 1
        mat.gas_puff          =  mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt > 0);
        mat.wall_pumping      = -mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt < 0);
        if ~isfinite(fn0a)
            fn0a = z0dstruct.z0dinput.option.fn0a;
        end
        if ~isfinite(fn0a_div)
            fn0a_div = z0dstruct.z0dinput.option.fn0a_div;
        end
        switch z0dstruct.z0dinput.option.configuration
            case {0,1}
                mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a);
            case {2,3}
                mat.gas_puff      = mat.gas_puff ./ max(eps,data_zerod.xpoint .* fn0a_div + (~data_zerod.xpoint) .* fn0a);
            otherwise
                mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a_div);
        end
        
        % source totale dans le plasma
        mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
        % bilan = input - losses = 0 + changement de densite du plasma
        mat.bilan = mat.input - mat.plasmanetflux;
        % correction de gaz puff
        mat.gas_puff = mat.gas_puff - 1.1 .* mat.bilan .* (mat.bilan < 0);
        %mise a jour  bilan = input - losses = 0 + changement de densite du plasma
        mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
        mat.bilan = mat.input - mat.plasmanetflux;
        % pumping
        mat.pumping = max(mat.bilan,mat.wall_pumping);
        if length(data_zerod.temps) > 1
            mat.bilan_pumping = cumtrapz(data_zerod.temps,mat.pumping);
        else
            mat.bilan_pumping =  0;
        end
        % correction to gas puff
        puff_correc  = (mat.pumping - mat.bilan);
        mat.gas_puff = mat.gas_puff + puff_correc .* (puff_correc> 0);
        if length(data_zerod.temps) > 1
            mat.bilan_gas_puff = cumtrapz(data_zerod.temps,mat.gas_puff);
        else
            mat.bilan_gas_puff = 0;
        end
        mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
        % bilan final
        mat.bilan = mat.input - mat.plasmanetflux  -mat.pumping;
        if length(data_zerod.temps) > 1
            mat.bilan_input   = cumtrapz(data_zerod.temps,mat.input);
        else
            mat.bilan_input   = 0;
        end
        
        
        % flux to divertor/limiter
        mat.flux_divlim_total = mat.output + mat.gas_puff  + (1 - eta_p) .* mat.pellet + mat.perte_snbi - mat.s0_in;
        mat.recycle           = mat.flux_divlim_total - mat.pumping;
        mat.Reff              = mat.recycle ./ mat.flux_divlim_total;
        if length(data_zerod.temps) > 1
            mat.bilan_recycle     = cumtrapz(data_zerod.temps,mat.recycle);
        else
            mat.bilan_recycle     = 0;
        end
        mat.eta_g_high        = mat.s0_in ./ max(eps,mat.output + mat.gas_puff + mat.recycle);
        
        
        
        
    otherwise
        % flux de matiere sortant
        mat.temps        = z0dstruct.zerod.temps;
        mat.output       = interp1_imas(z0dstruct.profil.temps,z0dstruct.profil.ge(:,end) .* z0dstruct.profil.grho2(:,end) .* z0dstruct.profil.vpr_tor(:,end), ...
            z0dstruct.zerod.temps,'pchip','extrap') ./ (4.41e-4 .* phys.avo);
        mat.recycle      = z0dstruct.zerod.n0a ./ (4.41e-4 .* phys.avo);  % recyclage  + gas puff
        mat.pellet       = interp1_imas(z0dstruct.profil.temps,trapz(z0dstruct.profil.xli,z0dstruct.profil.spellet .* z0dstruct.profil.vpr,2), ...
            z0dstruct.zerod.temps,'pchip','extrap')  ./ (4.41e-4 .* phys.avo);
        mat.pellet(~isfinite(mat.pellet)) = 0;
        mat.snbi         = z0dstruct.zerod.pnbi  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo);
        mat.input        = mat.recycle + mat.pellet + mat.snbi;
        
        mat.bilan_input   = cumtrapz(z0dstruct.zerod.temps,mat.input);
        mat.bilan_recycle = cumtrapz(z0dstruct.zerod.temps,mat.recycle);
        mat.bilan_pellet  = max(eps,cumtrapz(z0dstruct.zerod.temps,mat.pellet));
        mat.bilan_nbi     = max(eps,cumtrapz(z0dstruct.zerod.temps,mat.snbi));
        mat.bilan_output  = cumtrapz(z0dstruct.zerod.temps,mat.output);
        
        ntot                 = trapz(z0dstruct.profil.xli,z0dstruct.profil.vpr .* z0dstruct.profil.nep,2) ./ (4.41e-4 .* phys.avo);
        mat.contents         = interp1_imas(z0dstruct.profil.temps,ntot,z0dstruct.zerod.temps,'pchip','extrap');
        mat.bilan_io         = mat.bilan_input - mat.bilan_output;
        % la precision numerique du calcul est limite (la source est mieux connu que le flux sortant)
        mat.bilan_error      = mat.contents - mat.bilan_io -mat.contents(1);
        
        %mat.fuelling_pumping = z0dxdt(mat.bilan_io,z0dstruct.zerod.temps);
        mat.plasmanetflux    = z0dxdt(mat.contents,z0dstruct.zerod.temps); % flux net
        mat.flux_error       = z0dxdt(mat.bilan_error,z0dstruct.zerod.temps); % flux convecte et/ou diffuse
        % separation
        % correction flux sortant
        mat.output        = mat.output - mat.flux_error ;
        mat.bilan_output  = mat.bilan_output - mat.bilan_error;
        mat.bilan_io      = mat.bilan_input - mat.bilan_output;
        
        
        % calcul gaz puff (a partir du scaling pour la densite naturelle)
        % scaling
        fact_nbar = z0dstruct.zerod.nem ./ max(1,z0dstruct.zerod.nbar);
        tauref = min(1e3,max(min(z0dstruct.zerod.tauhe,z0dstruct.zerod.taup) ./ max(eps,1 - z0dstruct.z0dinput.option.Recycling),1e-6));
        snem   = fact_nbar .* z0dstruct.zerod.nbar_nat./ tauref .* z0dstruct.zerod.vp + z0dstruct.zerod.pnbi_th ./ z0dstruct.z0dinput.option.einj;
        [tntot,mat.ntot_nat] = z0ode(z0dstruct.zerod.temps,snem,tauref,fact_nbar(1) .* z0dstruct.zerod.nbar_nat(1) .* z0dstruct.zerod.vp(1));
        %else
        %  mat.ntot_nat = fact_nbar .* z0dstruct.zerod.nbar_nat .* z0dstruct.zerod.vp;
        %end
        % facteur de conversion
        mat.fuelling_bilan         = z0dstruct.zerod.vp .* z0dstruct.zerod.nem  - mat.ntot_nat - (mat.bilan_pellet + mat.bilan_nbi) .* (4.41e-4 .* phys.avo);
        mat.dfuelling_bilan_dt     = z0dxdt(mat.fuelling_bilan,z0dstruct.zerod.temps);
        mat.gas_puff               = mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt> 0) ./ (4.41e-4 .* phys.avo);
        mat.pumping                = -mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt < 0) ./ (4.41e-4 .* phys.avo);
        mat.recycle                = max(0,mat.recycle - mat.gas_puff);
        mat.coef_recycle           = mat.recycle  ./ max(mat.input,mat.output); % securite a cause du bruit sur les donnees + convection
        mat.bilan_recycle          = cumtrapz(z0dstruct.zerod.temps,mat.recycle);
        mat.bilan_gas_puff         = cumtrapz(z0dstruct.zerod.temps,mat.gas_puff);
        
        % extraction du temps d'interet
        if length(data_zerod.temps) == 1
            noms = fieldnames(mat);
            for k=1:length(noms)
                mat.(noms{k}) = interp1_imas(z0dstruct.zerod.temps,mat.(noms{k}),data_zerod.temps,'nearest','extrap');
            end
        end
end
