waitfor(zuicreefunform('gui_ImpTrans','ImpTrans_parameters',1,0,'','visible'));
if getappdata(0,'ZGUI_CANCEL')
    return;
end

path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')));
if isempty(getappdata(0,'local'))
  path_to_METIS_FACIT_data = [path_to_METIS_coef '/FACIT_TOOLS/data'];
  path_to_METIS_data = [path_to_METIS_coef '/data'];
  path_to_METIS_NEO_data = [path_to_METIS_coef '/NEO_TOOLS/data'];
else
  path_to_METIS_FACIT_data = [getappdata(0,'local') '/data_facit'];
  path_to_METIS_data = [getappdata(0,'local') '/data'];
  path_to_METIS_NEO_data = [getappdata(0,'local') '/data_neo'];
end  

unix(sprintf(['mkdir -p ' path_to_METIS_FACIT_data]));
unix(sprintf(['mkdir -p ' path_to_METIS_data]));
unix(sprintf(['mkdir -p ' path_to_METIS_NEO_data]));

% ---------------------------------------------------------------
% Useful flags
% ---------------------------------------------------------------
useNeo  = ImpTrans_parameters.neoclass_transport ~= "0";
useTurb = ImpTrans_parameters.turb_transport    ~= "0";

% Output arrays (initialised to zero)
irok = find(post.profil0d.xli > 0);
irok = 1:length(post.profil0d.xli);
G_n = nan(1,length(irok),6);  D_n = G_n;  V_n = G_n;   % neoclassical
G_t = G_n;                     D_t = G_n;  V_t = G_n;   % turbulent

% -------------------------------------------------------------------------
%  AUTO_CONS = 0  →  RAW DATA ONLY
% -------------------------------------------------------------------------
if ImpTrans_parameters.auto_cons == 0
% Determine transport case based on user selection
nbloop_max = 1;
qlkzparams.ext_nbloop = 1;
switch true
    % ----------------------------------------------------------------------
    % Case 1: Both turbulent and neoclassical contributions
    % ----------------------------------------------------------------------
    case (useTurb && useNeo)
        % ---------- QLKZ ----------
        switch ImpTrans_parameters.turb_transport
            case "QLKZ_e"
                disp("Option not available")
            case "QLKZ_imp"
            qlkzparams.imp_trans = 1;
            disp('[ImpTrans] QLKZ_imp in RAW mode (imp_trans = 1)');
            clear QLK_nn_data
            z0working_point_qlkz_std;
            [tmp itime] = min(abs(post.profil0d.temps-twp));
            G_t = G_loop(1,irok,:);
            D_t = D_loop(1,irok,:);
            V_t = V_loop(1,irok,:);
            otherwise
                warning('[ImpTrans] Unknown turbulent transport option.');
        end
        
    % ---------- NEO or FACIT ----------
    switch ImpTrans_parameters.neoclass_transport
        case "FACIT"
            waitfor(zuicreefunform('gui_facit','facitparams',1,0,'','visible'));
            if getappdata(0,'ZGUI_CANCEL'), return; end 
            [G_n,D_n,V_n] = facit_std(facitparams,post.z0dinput.option,post.z0dinput, post.profil0d, tmp, itime,1);
            disp('ici')
        case "NEO"
            waitfor(zuicreefunform('gui_neo','neoparams',1,0,'','visible'));
            if getappdata(0,'ZGUI_CANCEL'), return; end
            yneo_out = neo_std(neoparams,post.z0dinput.option,post.z0dinput, post.profil0d,itime,twp);
            G_n = yneo_out.G; D_n = yneo_out.D; V_n = yneo_out.V;
        otherwise
            warning('[ImpTrans] Unknown neoclassical transport option.');
    end
    flnm=[path_to_METIS_data '/QLKZimp_neocl_data.mat'];
    disp(['=========== DATA SAVED TO ' flnm '==============='])
    save(flnm,'G_n','G_t','D_n','D_t','V_n','V_t');

    % ----------------------------------------------------------------------
    % Case 2: Turbulent only (QLKZ)
    % ----------------------------------------------------------------------
    case (useTurb && ~useNeo)
        switch ImpTrans_parameters.turb_transport
            case "QLKZ_e"
                disp("Option not available")
            case "QLKZ_imp"
                qlkzparams.imp_trans = 1;
                disp('[ImpTrans] QLKZ_imp in RAW mode (imp_trans = 1)');
                clear QLK_nn_data
                z0working_point_qlkz_std;
                G_t = G_loop(1,irok,:);
                D_t = D_loop(1,irok,:);
                V_t = V_loop(1,irok,:);
            otherwise
                warning('[ImpTrans] Unknown turbulent transport option.');
        end

    % ----------------------------------------------------------------------
    % Case 3: Neoclassical only (FACIT or NEO)
    % ----------------------------------------------------------------------
    case (useNeo && ~useTurb)
        % Select time
        clear QLK_nn_data
        z0plot_reference;
        title('choose a time')
        drawnow	
        [twp,void] = ginput(1);

        % Ask user to enter the time manually
        prompt = {'Check the time value:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {num2str(twp)};
        answer = inputdlg(prompt, dlgtitle, dims, definput);

        if isempty(answer)
          return; 	% Exit if the user cancels the input dialog
        end

        twp = str2double(answer{1});
        close(gcf);
        drawnow

        [tmp itime] = min(abs(post.profil0d.temps-twp));
        switch ImpTrans_parameters.neoclass_transport
            case "FACIT"
                waitfor(zuicreefunform('gui_facit','facitparams',1,0,'','visible'));
                if getappdata(0,'ZGUI_CANCEL'), return; end
                [G_n,D_n,V_n] = facit_std(facitparams,post.z0dinput.option,post.z0dinput, post.profil0d, tmp, itime,1);
                flnm=[path_to_METIS_FACIT_data '/FACIT_data.mat'];
                disp(['=========== DATA SAVED TO ' flnm '==============='])
                save(flnm);
                
            case "NEO"
                waitfor(zuicreefunform('gui_neo','neoparams',1,0,'','visible'));
                if getappdata(0,'ZGUI_CANCEL'), return; end
                yneo_out = neo_std(neoparams,post.z0dinput.option,post.z0dinput, post.profil0d,itime,twp);
                G_n = yneo_out.G; D_n = yneo_out.D; V_n = yneo_out.V;
                flnm=[path_to_METIS_NEO_data '/NEO_data.mat'];
                disp(['=========== DATA SAVED TO ' flnm '==============='])
                save(flnm);
            otherwise
                warning('[ImpTrans] Unknown neoclassical transport option.');
        end

    % ----------------------------------------------------------------------
    % Case 4: No contribution at all
    % ----------------------------------------------------------------------
    otherwise
        disp('[ImpTrans] No transport model selected (neither turbulent nor neoclassical).');
end

  % -------------------------------------------------------------------------
  %  AUTO_CONS = 1  →  SELF-CONSISTENT SOLUTIONS
  % -------------------------------------------------------------------------
else
    
    nbloop_max = ImpTrans_parameters.inter_nbloop_max;
    qlkzparams.ext_nbloop = ImpTrans_parameters.ext_nbloop_max;
      
    % ------------------------------------------------------------------
    % Determine QLKZ imp_trans flag
    % ------------------------------------------------------------------
    if  useTurb && ~useNeo
        switch ImpTrans_parameters.turb_transport
            case "QLKZ_e"
                mode = 0;
                G_loop = zeros(1,21,1);
            case "QLKZ_imp"
                mode = 2;
            otherwise
                warning('[ImpTrans] Unknown turbulent transport option.');
        end
        qlkzparams.imp_trans = mode;
        clear QLK_nn_data
        disp(['[ImpTrans] Self-consistent QLKZ (imp_trans = ' num2str(mode) ')']);
        z0working_point_qlkz_std;
        [tmp itime] = min(abs(post.profil0d.temps-twp));

        G_t = G_loop(1,irok,:);
        D_t = D_loop(1,irok,:);
        V_t = V_loop(1,irok,:);

        if mode == 2
            flnm=[path_to_METIS_data '/QLKZimp_data.mat'];
            disp(['=========== DATA SAVED TO ' flnm '==============='])
            save(flnm);
        else
        end
            
    elseif useTurb && useNeo && ImpTrans_parameters.neoclass_transport == "FACIT"
        mode = 3;   % QLKZ + FACIT


        qlkzparams.imp_trans = mode;
        clear QLK_nn_data
        disp(['[ImpTrans] Self-consistent QLKZ (imp_trans = ' num2str(mode) ')']);
        z0working_point_qlkz_std;
        [tmp itime] = min(abs(post.profil0d.temps-twp));

        G_t = G_loop(1,irok,:);
        D_t = D_loop(1,irok,:);
        V_t = V_loop(1,irok,:);

        G_n = G_loop_n(1,irok,:);
        D_n = D_loop_n(1,irok,:);
        V_n = V_loop_n(1,irok,:);

        flnm=[path_to_METIS_data '/QLKZimp_FACIT_data.mat'];
        disp(['=========== DATA SAVED TO ' flnm '==============='])
        save(flnm,'G_n','G_t','D_n','D_t','V_n','V_t');

    elseif useTurb && useNeo && ImpTrans_parameters.neoclass_transport == "NEO"
        mode = 4;   % QLKZ + NEO
        nbloop_max = ImpTrans_parameters.inter_nbloop_max;

        qlkzparams.imp_trans = mode;
        clear QLK_nn_data
        disp(['[ImpTrans] Self-consistent QLKZ (imp_trans = ' num2str(mode) ')']);
        z0working_point_qlkz_std;
        [tmp itime] = min(abs(post.profil0d.temps-twp));
        
        G_t = G_loop(1,irok,:);
        D_t = D_loop(1,irok,:);
        V_t = V_loop(1,irok,:);

        G_n = G_loop_n(1,irok,:);
        D_n = D_loop_n(1,irok,:);
        V_n = V_loop_n(1,irok,:);

        flnm=[path_to_METIS_data '/QLKZimp_NEO_data.mat'];
        disp(['=========== DATA SAVED TO ' flnm '==============='])
        save(flnm,'G_n','G_t','D_n','D_t','V_n','V_t');
    else
        error(['[ImpTrans] Unsupported self-consistent combination; ' ...
               'pure neoclassical auto-consistency not implemented yet.']);
    end
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=4;
if  useTurb && ~useNeo
    switch ImpTrans_parameters.turb_transport
        case "QLKZ_e"
            idx = 1;
        case "QLKZ_imp"
        otherwise
            warning('[ImpTrans] Unknown neoclassical transport option.');
    end
else
end

% === FIGURE 1 : Flux ===
figure;
hold on;
box on;
grid on;

species_names = {'$\Gamma_{neo}$', '$\Gamma_{turb}$'};

colors = lines(6);  % Une seule couleur pour z3
markers = {'s'};  % Un seul marqueur pour z3

h_n = gobjects(1, 1);  % 1 courbe pour NEO (z3)
h_t = gobjects(1, 1);  % 1 courbe pour FACIT (z3)


% Plot neoclassical contribution
h_n = plot(post.profil0d.xli(:, irok), G_n(:,irok,idx), ...
    '-', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

% Plot turbulent contribution
h_t = plot(post.profil0d.xli(:, irok), G_t(:,:,idx), ...
    ':', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

% set(gca, 'YScale', 'log');
xlabel('r/a', 'FontSize', 14);
ylabel('\Gamma (m^{-2}/s)', 'Interpreter', 'tex', 'FontSize', 14);
% title('R_N=10^{-4}', 'Interpreter', 'tex', 'FontSize', 14');
legend([h_n, h_t], species_names,  'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
title('W')
	     
% === FIGURE 2 : Coefficients V ===
figure;
hold on;
box on;
grid on;

species_names = {'$V_{neo}$', '$V_{turb}$'};

colors = lines(6);  % Une seule couleur pour z3
markers = {'s'};  % Un seul marqueur pour z3

h_n = gobjects(1, 1);  % 1 courbe pour NEO (z3)
h_t = gobjects(1, 1);  % 1 courbe pour FACIT (z3)
 

% Plot neoclassical contribution
h_n = plot(post.profil0d.xli(:, irok), V_n(:,irok,idx), ...
    '-', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

% Plot turbulent contribution
h_t = plot(post.profil0d.xli(:, irok), V_t(:,:,idx), ...
    ':', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

% set(gca, 'YScale', 'log');
xlabel('r/a', 'FontSize', 14);
ylabel('V (m/s)', 'Interpreter', 'tex', 'FontSize', 14);
% title('R_N=10^{-4}', 'Interpreter', 'tex', 'FontSize', 14');
legend([h_n, h_t], species_names,  'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
title('W')
	     
% === FIGURE 3 : Coefficients de diffusion D ===
figure;
hold on;
box on;
grid on;

species_names = {'$D_{neo}$', '$D_{turb}$'};

colors = lines(6);  % Une seule couleur pour z3
markers = {'s'};  % Un seul marqueur pour z3

h_n = gobjects(1, 1);  % 1 courbe pour NEO (z3)
h_t = gobjects(1, 1);  % 1 courbe pour FACIT (z3) 

% Plot neoclassical contribution
h_n = plot(post.profil0d.xli(:, irok), D_n(:,irok,idx), ...
    '-', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

% Plot turbulent contribution
h_t = plot(post.profil0d.xli(:, irok), D_t(:,:,idx), ...
    ':', 'Color', colors(idx, :), 'LineWidth', 1.5, ...
    'Marker', markers{1}, 'MarkerSize', 5);

xlabel('r/a', 'FontSize', 14);
ylabel('D (m^{2}/s)', 'Interpreter', 'tex', 'FontSize', 14);
% title('R_N=10^{-4}', 'Interpreter', 'tex', 'FontSize', 14');
legend([h_n, h_t], species_names,  'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
if useTurb && useNeo
    set(gca, 'YScale', 'log');
else 
    set(gca, 'FontSize', 12);
end
title('W')








