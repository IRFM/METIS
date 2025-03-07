if isappdata(0,'userpath')
    rmappdata(0,'userpath');
end
if isempty(which('zineb_path'))
    addpath ../../../../../SVN_ITER/metis/
    if isempty(which('zero1t'))
        zineb_path
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% here insert FBE initialisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METIS initailisation           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% equilibrium specific input for metis
create_data4metis.time	 = 604;        % current time  [1] (s)
create_data4metis.Ip     = 15.0E+6;     % plasma current [1] (A)
create_data4metis.Bp_error_vacuum  = 0; % error magnetic field, i.e. vaccum poloidal magnetic field for breakdown computation  [1] (T)

% separatrice
%create_data4metis.R_LCFS  =  some R vector;
%create_data4metis.Z_LCFS  =  some Z vector;

% psi values an 21 iso contour values, radius of magentic axis and elongation 
create_data4metis.Psi	  = Psi; %poloidal flux vector Psi from magnetic axis to LCFS (Wb / 2pi, maximum value on magnetic axis) [N]
create_data4metis.Raxe	  = Raxe;%radius of magnetic axis versus PSI[n] (m)
create_data4metis.K	  = K;   %fux surface elongation K = (max(Z) - min(Z)) ./ (max(R) - min(R)) versus psi [N]

%% metis call INITIALISATION
[metis_data4create,internal_state,equi_geo] = metis_create_link([],create_data4metis,[],parameters);

% adapt B0;
fbe.B0 = metis_data4create.fdia(end)/fbe.R0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial FBE inverse problem for equilibrium initailisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convergence loop on initial solution
for i = 1:3
    
    %% we solve an inverse problem, with the profile and boundary data from metis
    prof_data = [(metis_data4create.psi'-metis_data4create.psi(1))/...
        (metis_data4create.psi(end)-metis_data4create.psi(1)),...
        metis_data4create.dptotdpsi', 1/2*metis_data4create.df2dpsi'];
    prof_data(1,2) = prof_data(2,2);
     
    % call of FBE    
    % inverse problem 
    
    %recall of METIS with new LCFS
    % separatrice
    [P_LCFS,~,~,~] = get_PlasmaBnd_Points(Mesh,psi);
    B = [mean(P_LCFS(:,1)),mean(P_LCFS(:,2))];
    Angle = angle((P_LCFS(:,1)-B(1)+1i*(P_LCFS(:,2)-B(2))));
    [Angle,id] = sort(Angle); P_LCFS = P_LCFS(id,:);
    %figure; plot_geometry(plist.mesh_file, 'f'); hold on;
    %plot(P_LCFS(:,1),P_LCFS(:,2),'.-');
    create_data4metis.R_LCFS  = ;
    create_data4metis.Z_LCFS  = ;
    
    % psi values an 21 iso contour values, radius of magentic axis and elongation 
    create_data4metis.Psi	  = Psi; %poloidal flux vector Psi from magnetic axis to LCFS (Wb / 2pi, maximum value on magnetic axis) [N]
    create_data4metis.Raxe	  = Raxe;%radius of magnetic axis versus PSI[n] (m)
    create_data4metis.K	      = K;   %fux surface elongation K = (max(Z) - min(Z)) ./ (max(R) - min(R)) versus psi [N] 
     
    %% metis call
    metis_data4create2 = metis_data4create;
    [metis_data4create,internal_state] = metis_create_link(internal_state,create_data4metis,[],parameters);
    zplotstruct(metis_data4create,metis_data4create2)
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct call of FBE for checking    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time initialisation  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for evolution (nominal value = 3)
% METIS current diffusion equation boundary condition mode
%parameters.ipfluxmode = 2; 
parameters.ipfluxmode = 3; 
%parameters.ipfluxmode = 1;

% set psi offset
internal_state.psi_offset = -internal_state.profil0d.psi(end)+psi(boundary);


% FBE evolution initialisation


% loop on time
for k = 1:n_tsteps
    
    % stop VS at k = 100
    if k > 100
      ctr_on_off = 0;
    end
    
    prof_data_old = prof_data;
    
    %% equilibrium specific input for metis
    create_data4metis.time	 = ?      % current time  [1] (s)
    time(k) = create_data4metis.time;
    create_data4metis.Ip     = ?;     % plasma current [1] (A)
    create_data4metis.Bp_error_vacuum  = 0;  % error magnetic field, i.e. vaccum poloidal magnetic field for breakdown computation  [1] (T)
   
    fprintf('#### timestep %d and time (in s) %d ####\n \n',k,create_data4metis.time);
   
    % separatrice
    create_data4metis.R_LCFS  = ;
    create_data4metis.Z_LCFS  = ;
    
    % psi values an 21 iso contour values, radius of magentic axis and elongation 
    create_data4metis.Psi	  = Psi; %poloidal flux vector Psi from magnetic axis to LCFS (Wb / 2pi, maximum value on magnetic axis) [N]
    create_data4metis.Raxe	  = Raxe;%radius of magnetic axis versus PSI[n] (m)
    create_data4metis.K	      = K;   %fux surface elongation K = (max(Z) - min(Z)) ./ (max(R) - min(R)) versus psi [N] 
     
    %% metis call
    metis_data4create2 = metis_data4create;
    [metis_data4create,internal_state] = metis_create_link(internal_state,create_data4metis,[],parameters);
    FBE.Ip_evo = metis_data4create.ip;
    fprintf('Ip %d \n', FBE.Ip_evo);
    %% the profile data from metis
    prof_data = [(metis_data4create.psi'    -metis_data4create.psi(1))/...
                 (metis_data4create.psi(end)-metis_data4create.psi(1)),...
                  metis_data4create.dptotdpsi', 1/2*metis_data4create.df2dpsi'];
    
    
    %% update voltage from controler 
    FBE.Supplydata.voltages = voltages_init + ctr_on_off .* ctr.voltages';

    % for current computation
    psi_mem = psi;
    
    % call of FBE in evolution mode 
    call_fbe_evol
    
    %% for controler, compute some useful quantities from fbe solution
    %% compute obs
    obs.ip(end+1)      = ;
    obs.Rcj(end+1)     = ;
    obs.Zcj(end+1)     = ;
    obs.time(end+1)    = create_data4metis.time;
    %% call controler
    ctr = simple_vs_ctr(obs,ctr)

    %display FBE solution
    
end
post.zerod = internal_state.z0dstruct.zs;
post.profil0d = internal_state.z0dstruct.profil;
post.z0dinput = internal_state.z0dstruct.z0dinput;




