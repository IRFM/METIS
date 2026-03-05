% This function extracts neoclassical fluxes and transport coefficients for three species
% associated to QLKZ, assuming linearity in the fluxes

function [G,D,V] = facit_std(facitparams, option, z0dinput, profil, tmp, itime, fac)

% Initialize output arrays
irok = 1:length(profil.xli);
TE = interp1(profil.temps, profil.tep(:,irok), profil.temps(itime));
Q = interp1(profil.temps, profil.qjli(:,irok), profil.temps(itime));
insize = size(Q);

% Initialize output arrays with proper dimensions
D = zeros([insize, 6]);
V = zeros([insize, 6]);
G = zeros([insize, 6]);

% dzeff = max((option.zimp.^2 - option.zimp) .* profil.nzp(itime,:) ./ profil.nep(itime,:));
Dzeff = zeros(6, 1);

% Define species to process
species = {
    struct('name', 'Z1', 'A', 2*option.zimp, 'Z', option.zimp, ...
            'Nfunc', @() interp1(profil.temps, profil.nzp(:,irok), profil.temps(itime)) * fac, 'idx', 2),...
    struct('name', 'Z2', 'A', 2*option.zmax, 'Z', option.zmax, ...
            'Nfunc', @() interp1(profil.temps, profil.nzp(:,irok) .* option.rimp, profil.temps(itime))* fac, 'idx', 3),...
    struct('name', 'Z3', 'A', 184, 'Zfunc', @() round((1 - SnFraction(option)) .* z0wavez(TE) + ...
            SnFraction(option) .* z0snavez(TE)), ...
            'Nfunc', @() interp1(profil.temps, profil.nwp(:,irok), profil.temps(itime)) * fac, 'idx', 4),...
    struct('name', 'H', 'A', 1, 'Z', 1, ...
            'Nfunc', @() interp1(profil.temps, profil.nip(:,irok) .* option.cmin, profil.temps(itime))* fac, 'idx', 5)
};

% Computes dzeff and ratio
for i = 1:length(species)
    s = species{i};

    % Handle Z value (can be function or scalar)
    if isfield(s, 'Zfunc')
        Zimp = s.Zfunc();
    else
        Zimp = s.Z;
    end

    Nimp = s.Nfunc();
    % disp(mean(Nimp ./ interp1(profil.temps, profil.nip(:,irok), profil.temps(itime))))
    ratio_imp = Nimp ./ interp1(profil.temps, profil.nip(:,irok), profil.temps(itime));           
    Zfac      = Zimp.^2 - Zimp;                   
    dzeff_imp_max = max(ratio_imp .* Zfac);
    ratio_max(i)=max(ratio_imp);
    Dzeff(i) =  dzeff_imp_max;
end
    
% Process each species
for i = 1:length(species)
    s = species{i};
    fprintf('%s computing\n', s.name);
    
    % Handle Z value (can be function or scalar)
    if isfield(s, 'Zfunc')
        Zimp = s.Zfunc();
    else
        Zimp = s.Z;
    end

    Nimp = s.Nfunc();
    % Compute coefficients
    yfacit_out = facit_coeff(Dzeff, ratio_max,facitparams, option, z0dinput, profil, itime, s.A, Zimp, Nimp, fac);
    
    % Store results
    D(:,:,s.idx) = reshape(yfacit_out.Da_neo, insize);
    V(:,:,s.idx) = reshape(yfacit_out.Va_neo, insize);
    G(:,:,s.idx) = reshape(yfacit_out.fluxImp, insize);
end

function sf = SnFraction(option)
    if isfield(option, 'Sn_fraction')
        sf = option.Sn_fraction;
    else
        sf = 0;
    end
end

function [yfacit_out] = facit_coeff(dzeff, ratio, facitparams, option, z0dinput, profil, itime, Aimp, Zimp, Nimp, fac)
    iplot = 0;
    yfacit_out = struct('fluxImp', 0, 'Da_neo', 0, 'Va_neo', 0);
    
    if sum(dzeff, 'all') > 0.1 || facitparams.MULTIS
        collision_outputs = struct();
        
        % Define collision types
        collision_types = {
            struct('name', 'H', 'condition', Aimp ~= 1, 'multis', 4),...
            struct('name', 'Z2', 'condition', Aimp ~= 2*option.zmax, 'multis', 3),...
            struct('name', 'Z1', 'condition', Aimp ~= 2*option.zimp, 'multis', 2)
        };
%            struct('name', 'Z3', 'condition', Aimp ~= 184, 'multis', 4),...
	
        for j = 1:length(collision_types)
            ct = collision_types{j};
            if ct.condition
                fprintf('Coll with %s\n', ct.name);
                facitparams.MULTIS = ct.multis;
                collision_outputs.(ct.name) = run_facit(facitparams, option, z0dinput, profil, profil.temps(itime), iplot, Aimp, Zimp, Nimp, fac);		
            else
                collision_outputs.(ct.name) = struct('fluxImp', 0, 'Da_neo', 0, 'Va_neo', 0);
            end
        end
        
        % Main ion collision
        facitparams.MULTIS = 0;
        disp('Coll with main ions');
        yfacit_out_i = run_facit(facitparams, option, z0dinput, profil, profil.temps(itime), iplot, Aimp, Zimp, Nimp, fac);

        % Initialize outputs
        yfacit_out.fluxImp = 0;
        yfacit_out.Da_neo = 0;
        yfacit_out.Va_neo = 0;

        dzeff_threshold = 0.1;
        ratio_threshold = 0.05;	

        if dzeff(1) ./ sum(dzeff, 'all')>= dzeff_threshold || ratio(1)>ratio_threshold
            yfacit_out.fluxImp = yfacit_out.fluxImp + collision_outputs.Z1.fluxImp;
            yfacit_out.Da_neo = yfacit_out.Da_neo + collision_outputs.Z1.Da_neo;
            yfacit_out.Va_neo = yfacit_out.Va_neo + collision_outputs.Z1.Va_neo;
        end

        if dzeff(2) ./ sum(dzeff, 'all') >= dzeff_threshold || ratio(2)>ratio_threshold
            yfacit_out.fluxImp = yfacit_out.fluxImp + collision_outputs.Z2.fluxImp;
            yfacit_out.Da_neo = yfacit_out.Da_neo + collision_outputs.Z2.Da_neo;
            yfacit_out.Va_neo = yfacit_out.Va_neo + collision_outputs.Z2.Va_neo;
        end

%        if dzeff(3) ./ sum(dzeff, 'all') >= dzeff_threshold
%            yfacit_out.fluxImp = yfacit_out.fluxImp + collision_outputs.Z3.fluxImp;
%            yfacit_out.Da_neo = yfacit_out.Da_neo + collision_outputs.Z3.Da_neo;
%            yfacit_out.Va_neo = yfacit_out.Va_neo + collision_outputs.Z3.Va_neo;
%        end

        if dzeff(4) ./ sum(dzeff, 'all') >= 0  || ratio(4)>ratio_threshold
            yfacit_out.fluxImp = yfacit_out.fluxImp + collision_outputs.H.fluxImp;
            yfacit_out.Da_neo = yfacit_out.Da_neo + collision_outputs.H.Da_neo;
            yfacit_out.Va_neo = yfacit_out.Va_neo + collision_outputs.H.Va_neo;
        end

        % On ajoute toujours la contribution des ions principaux
        yfacit_out.fluxImp = yfacit_out.fluxImp + yfacit_out_i.fluxImp;
        yfacit_out.Da_neo = yfacit_out.Da_neo + yfacit_out_i.Da_neo;
        yfacit_out.Va_neo = yfacit_out.Va_neo + yfacit_out_i.Va_neo;

    else
        yfacit_out = run_facit(facitparams, option, z0dinput, profil, profil.temps(itime), iplot, Aimp, Zimp, Nimp, fac);
    end
end
end
