function sawtooth_data = metis4muscle3_mimic_sawtooth(internal_state)

% backport from JT-60SA simulator
% get METIS parameter
qdds_metis     = internal_state.metis.output.z0dstruct.z0dinput.option.qdds;
q0_dds_trig    = internal_state.metis.output.z0dstruct.z0dinput.option.q0_dds_trig;
s1crit         = internal_state.metis.output.z0dstruct.z0dinput.option.s1crit;
betap1crit     = internal_state.metis.output.z0dstruct.z0dinput.option.betap1crit;
qddsmode       = internal_state.metis.output.z0dstruct.z0dinput.option.ddsmode;
w1             = internal_state.metis.output.z0dstruct.z0dinput.option.w1;
epsq           = internal_state.metis.output.z0dstruct.z0dinput.option.epsq;
% METIS coordinate
xli            = internal_state.metis.output.z0dstruct.profil0d.xli;

% get data from external equilibrium
% select which port is used for equilibrium
% switch internal_state.workflow.used_ports
%     case 'input_s'
%         equilibrium = internal_state.actor.io.data.equilibrium_input_s.data.equilibrium;
%     otherwise
%         if isfield(internal_state.actor.io.data.equilibrium_f_init.data,'equilibrium') && ...
%                 ~isempty(internal_state.actor.io.data.equilibrium_f_init.data.equilibrium)
%             equilibrium = internal_state.actor.io.data.equilibrium_f_init.data.equilibrium;
%         else
%             equilibrium = internal_state.actor.io.data.equilibrium_input_s.data.equilibrium;
%         end
% end
%
[~,equilibrium] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'equilibrium',internal_state.workflow.used_ports);
%
q                = abs(equilibrium.time_slice{end}.profiles_1d.q(:)');
rho_tor_norm     = equilibrium.time_slice{end}.profiles_1d.rho_tor_norm(:)';
rho_tor          = equilibrium.time_slice{end}.profiles_1d.rho_tor(:)';
psi              = equilibrium.time_slice{end}.profiles_1d.psi(:)';
phi              = equilibrium.time_slice{end}.profiles_1d.phi(:)';
gm1              = equilibrium.time_slice{end}.profiles_1d.gm1(:)';
gm2              = equilibrium.time_slice{end}.profiles_1d.gm2(:)';
dvolume_drho_tor = equilibrium.time_slice{end}.profiles_1d.dvolume_drho_tor(:)';
darea_dpsi       = equilibrium.time_slice{end}.profiles_1d.darea_dpsi(:)';
pressure         = equilibrium.time_slice{end}.profiles_1d.pressure(:)';
f                = equilibrium.time_slice{end}.profiles_1d.f(:)';
j_parallel       = equilibrium.time_slice{end}.profiles_1d.j_parallel(:)';
% need Lao coordinate from NICE to convert index
R_out = equilibrium.time_slice{end}.profiles_1d.r_outboard;
R_in = equilibrium.time_slice{end}.profiles_1d.r_inboard;
x_ext = abs(R_out - R_in) ./ 2;
x_ext = x_ext(:)' ./ max(x_ext);
% COCOS section
% psi should be decresing in METIS
if psi(1) < psi(end)
    sign_psi = -1;
else
    sign_psi = 1;
end
% 2pi factor
q_test = abs((phi(end) - phi(1)) ./ (psi(end) - psi(1)));
if q_test <  (3 * mean(q))
      factor_2pi = 2 .* pi;
else
      factor_2pi = 1;  
end
psi = sign_psi .* psi ./ factor_2pi;
sign_b0          = sign(mean(f));
f = sign_b0 .* f;
sign_ip = sign(mean(j_parallel));
j_parallel = sign_ip .* j_parallel;
darea_dpsi =  sign_psi .* factor_2pi .* darea_dpsi;

% mechanism to mimick sawtooth
% swap data 
sawtooth_data = internal_state.workflow.sawtooth;
% swap data
sawtooth_data.previous = sawtooth_data.next;
sawtooth_data.next.delta_jni = [];
sawtooth_data.next.delta_eta = [];
sawtooth_data.next.x_ext     = x_ext;

% start of the model
mu0 = 4*pi*1e-7;
indice_inv = 0;
if (qdds_metis ~= 0)
    if qdds_metis > 0
        indcor    = max(find(q <= qdds_metis));
        if ~isempty(indcor) &&(indcor > 4) && (indcor < 87)
            indice_inv = indcor;
        end
    elseif qdds_metis < 0
        indtrig   = max(find(q <= abs(qdds_metis)));
        % case with additional trigger on critical magnetic shear
        if ~isempty(indtrig) && (s1crit > 0) && (indtrig > 4) && (indtrig < 87)
            shear_loc = abs(rho_tor_norm ./ q .*  pdederive(rho_tor_norm,q,2,2,2,1));
            %
            % betap1
            psid1k    = 2 .* pi .* pdederive(rho_tor,psi,0,2,2,1);
            bpol1k  = sqrt(gm2 .* abs(psid1k));
            emp1k   = trapz(rho_tor(1:indtrig),bpol1k(1:indtrig) .^ 2 .* dvolume_drho_tor(1:indtrig),2) ./ 2 ./ mu0;
            if emp1k(1) > 0
                betap1k  = trapz(rho_tor(1:indtrig),pressure(1:indtrig) .* dvolume_drho_tor(1:indtrig),2) ./ emp1k;
            else
                betap1k = 0;
            end
            %disp([indtrig,shear_loc(indtrig), s1crit,TokaMesh.equilibrium.profiles_1d.q(2),q0_dds_trig]);
            if (shear_loc(indtrig) > s1crit)
                % ST
                indice_inv = indtrig;
            elseif (betap1k > betap1crit)
                % ST
                indice_inv = indtrig;
            elseif (TokaMesh.equilibrium.profiles_1d.q(2) < q0_dds_trig)
                % ST
                indice_inv = indtrig;
            end
        end
        % time suppression for ST rebound
        if (internal_state.workflow.clock.t_cur - sawtooth_data.LAST_ST_TIME) <= (internal_state.workflow.clock.dt * internal_state.actor.settings.value.metis.sawtooth.nbtmin)
            indice_inv = 0;
        end
        
    end
    
    % ind_inv should be put on METIS grid
    if indice_inv > 0
        x_inv = x_ext(indice_inv);
        indice_inv = min(length(xli),max(1,round(interp1(xli,1:length(xli),x_inv,'linear','extrap'))));
    end
    
    % compute perturbation to be applied on eta and jni to mimick crash of
    % sawtooth in  NICE
    if indice_inv > 0
        % compute equivalent change in J_ni
        if qdds_metis > 0
            qnew = q;
            qnew(1:indice_inv) = max(qnew(1:indice_inv),qdds_metis);
            indmix = indice_inv;
        else
            [qnew,indmix] = z0dds(rho_tor_norm,psi,q,qddsmode,w1,epsq);
        end
        %Cq            = - abs(TokaMesh.equilibrium.profiles_1d.q) .* TokaMesh.equilibrium.profiles_1d.dpsi_drho_tor /2/pi;
        Cq            = f .* dvolume_drho_tor ./ (4*pi^2) .* gm1;
        Cq(1) = Cq(2);
        delta         = Cq ./ qnew - Cq ./ q;
        variation_psi = cumtrapz(rho_tor,abs(delta),2);
        variation_psi(1) = variation_psi(2);
        inte          = dvolume_drho_tor .* gm2 ./ f .* delta;
        inte          = sgolayfilt(medfilt1(inte,3),1,5); % anti aliasing filter
        delta_jnib0   =   f .^ 2 ./ dvolume_drho_tor ./ mu0 .* pdederive(rho_tor,inte,2,2,2,1);
        delta_jnib0(1) =  delta_jnib0(2);
        fact_fjnib     = internal_state.metis.output.profil.eta ./ f ./ gm1 .* abs(j_parallel);
        fjnib_vect     = variation_psi ./ internal_state.workflow.clock.dt .* ...
                         cumtrapz(psi,darea_dpsi,2) ./ cumtrapz(psi,darea_dpsi .* fact_fjnib,2);
        fjnib_vect(1)  = fjnib_vect(2);
        fjnib_vect     = max(1,min(max(abs(j_parallel)) ./ max(abs(delta_jnib0)),fjnib_vect));
        % reconnexion on one time step
        if qdds_metis < 0
            fjnib = max(1,interp1(1:length(fjnib_vect),fjnib_vect,indmix,'linear','extrap'));
        else
            fjnib = 1;
        end        
        % in any case, delta_jnib0 is = 0 if no real
        % reconnexion
        % should be added to jni before IDS creation
        sawtooth_data.next.delta_jni = fjnib .* delta_jnib0;
        if sawtooth_data.next.delta_jni(2) > 0
            sawtooth_data.next.delta_jni = - sawtooth_data.next.delta_jni;
        end
        
        if qdds_metis > 0
            % clamping
            indice_eta_change = 1:ceil(indice_inv);
            sawtooth_data.next.delta_eta(indice_eta_change) = max(internal_state.metis.output.profil.eta(indice_eta_change))- ...
                                                           internal_state.metis.output.profil.eta(indice_eta_change);
            sawtooth_data.LAST_ST_TIME  = internal_state.workflow.clock.t_cur;
        elseif (qdds_metis < 0) && ~isempty(indmix) && (indmix > 1)
            sawtooth_data.LAST_ST_TIME = internal_state.workflow.clock.t_cur;
            indice_inv = indmix;
            switch qddsmode
                case 1
                    % Porcelli model
                    disp('Porcelli model for reconnexion is not yet implemented');
                    disp('Using Kadomtsev incomplete model instead');
                    % Porcelli incomplete
                    indice_eta_change = find(delta_jnib0 < 0);
                    if ~isempty(indice_eta_change)
                        sawtooth_data.next.delta_eta = max(internal_state.metis.output.profil.eta(indice_eta_change))- ...
                                                           internal_state.metis.output.profil.eta(indice_eta_change);
                    end
                    
                case 2
                    % Kadomtsev (full radius)
                    indice_eta_change = 1:ceil(indice_inv);
                    sawtooth_data.next.delta_eta = max(internal_state.metis.output.profil.eta(indice_eta_change))- ...
                                                       internal_state.metis.output.profil.eta(indice_eta_change);
                    
                otherwise
                    % Kadomtsev incomplete
                    indice_eta_change = 1:max(1,ceil((1 - w1) * indice_inv));
                    sawtooth_data.next.delta_eta = max(internal_state.metis.output.profil.eta(indice_eta_change))- ...
                                                       internal_state.metis.output.profil.eta(indice_eta_change);
            end
        else
            indice_inv = 0; % no real ST crash
        end
    end
end

% update data
sawtooth_data.INDICE_INV_ST = indice_inv;
sawtooth_data.MINIMAL_TAU_ST = internal_state.workflow.clock.dt .* internal_state.actor.settings.value.metis.sawtooth.nbtmin;
      




