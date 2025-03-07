function option_feeqs_eqdsk = compute_eqdsk_feeqs_based(option_feeqs_eqdsk)



% GUI generation
if nargin < 1
    
    valeur.eqksk_file_root_name	      = '';
    type.eqksk_file_root_name         = 'string';
    borne.eqksk_file_root_name        = '';
    defaut.eqksk_file_root_name       = '';
    info.eqksk_file_root_name         = 'Root name for EQDSK files; leave it empty for automatic name generation based on METIS information (machine name and shot name)';
    
    valeur.first_time     = 0;
    type.first_time       = 'float';
    borne.first_time      = [0,Inf];
    defaut.first_time     = 0;
    info.first_time      = 'First time of the METIS simulation at which EQDSK file would be computed (if both first_time and last_time =0, graphical time selection will be proposed at the next step)';
    
    valeur.last_time     = 0;
    type.last_time       = 'float';
    borne.last_time      = [0,Inf];
    defaut.last_time     = 0;
    info.last_time      = 'Last time of the METIS simulation for which EQDSK file would be computed (if equal to first_time, only one EQDSK file will be created)';

    valeur.grid_mode    = 'automatic';
    type.grid_mode      = 'string';
    borne.grid_mode     = {'automatic','manual'};
    defaut.grid_mode    = 'automatic';
    info.grid_mode      = 'if = automatic, grid definition will be automatically generated; otherwise manual radial and vertical number of point will be used';
    
    valeur.nbr_points   = 51;
    type.nbr_points     = 'integer';
    borne.nbr_points    = [5,501];
    defaut.nbr_points   = 51;
    info.nbr_points     = 'if grid_mode = manual, number of radial points in the grid';

    valeur.nbz_points   = 71;
    type.nbz_points     = 'integer';
    borne.nbz_points    = [7,701];
    defaut.nbz_points   = 71;
    info.nbz_points     = 'if grid_mode = manual, number of vertical points in the grid';

    valeur.precision    = 'high';
    type.precision      = 'string';
    borne.precision     = {'low','medium','high'};
    defaut.precision    = 'high';
    info.precision      = 'precision of equilibrium computation and interpolation (computational time increases with precision)';
    
    valeur.extrapolation    = 'G-S polynomial';
    type.extrapolation      = 'string';
    borne.extrapolation     = {'Interpolation','G-S polynomial'};
    defaut.extrapolation    = 'G-S polynomial';
    info.extrapolation      = 'Method used to extrapolate outside the LCFS:\n1/Interpolation = use simple interpolation method;\n2/G-S polynomial = fit a polynomial solution of G-S constained with flux and magnetic field on each point of LCFS';

    valeur.plotonoff    = 1;
    type.plotonoff     = 'integer';
    borne.plotonoff     = {0,1};
    defaut.plotonoff    = 1;
    info.plotonoff      = sprintf('Set level of displayed graphs:  0 = no graphs; 1 = results only\n(for multi time slices computation, it is recommended to switch off graohics display)');

 
    
    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
    
    option_feeqs_eqdsk.valeur     = valeur;
    option_feeqs_eqdsk.type       = type;
    option_feeqs_eqdsk.borne      = borne;
    option_feeqs_eqdsk.defaut     = defaut;
    option_feeqs_eqdsk.info       = info;
    option_feeqs_eqdsk.interface  = interface;
    
    option_feeqs_eqdsk.description = sprintf('Computation of EQDSK files using fast fixe boundary equilibrium solver of FEEQS.M;\nextrapolation outside the lCFS will not included real x-point and is valid only far from external current sources\n(FEEQS.M is not part of METIS and should be installed separtely)');   % description (une ligne) de la fonction
    
    option_feeqs_eqdsk.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    option_feeqs_eqdsk.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    option_feeqs_eqdsk.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
    
    % end of GUI form declaration
    return
    
end

% get METIS data
post = evalin('base','post');
time_metis = post.z0dinput.cons.temps;

% time selection
if (option_feeqs_eqdsk.first_time == 0) && (option_feeqs_eqdsk.last_time == 0)
   % selection of the time slice
    evalin('base','z0plot_reference_post;');
    subplot(3,1,1);
    title('choose 1 time slice for EQDSK computation')
    drawnow	
    [time_eqdsk,~] = ginput(1);
    close(gcf);
    drawnow
    option_feeqs_eqdsk.first_time = time_eqdsk;
    option_feeqs_eqdsk.last_time  = time_eqdsk;
elseif option_feeqs_eqdsk.first_time == option_feeqs_eqdsk.last_time 
    indt = find(time_metis >= option_feeqs_eqdsk.first_time,1);
    if isempty(indt)
        time_eqdsk = time_metis(end);
    else
        time_eqdsk =time_metis(indt);
    end
else
    time_eqdsk = time_metis((time_metis >= option_feeqs_eqdsk.first_time) & ...
                            (time_metis <= option_feeqs_eqdsk.last_time));
    if isempty(time_eqdsk)
        indt = find(time_metis >= ((option_feeqs_eqdsk.first_time + option_feeqs_eqdsk.last_time) / 2),1);
        if isempty(indt)
            time_eqdsk = time_metis(end);
        else
            time_eqdsk =time_metis(indt);
        end        
    end
end

% grid 
switch option_feeqs_eqdsk.grid_mode
    case 'manual'
        ngrid = [option_feeqs_eqdsk.nbr_points,option_feeqs_eqdsk.nbz_points];
    otherwise
        ngrid = [];
end

disp(' ');
% call FEESQ.M
hf_void = figure;   % if some degug graph appears
try
    [profiles_2d,FPSI,FBR,FBZ,FBPHI,plist_out] =   ...
        make_eqdsk_with_feeqs(post,time_eqdsk,ngrid, ...
        option_feeqs_eqdsk.eqksk_file_root_name, ...
        option_feeqs_eqdsk.precision,option_feeqs_eqdsk.plotonoff,option_feeqs_eqdsk.extrapolation);
catch
    error(lasterror);
    keyboard
end
if ishandle(hf_void)
    close(hf_void);
end
if length(time_eqdsk) == 1
    option_feeqs_eqdsk.profiles_2d = profiles_2d;
    option_feeqs_eqdsk.FPSI        = FPSI;
    option_feeqs_eqdsk.FBR         = FBR;
    option_feeqs_eqdsk.FBZ         = FBZ;
    option_feeqs_eqdsk.FBPHI       = FBPHI;
    option_feeqs_eqdsk.plist_out   = plist_out;
end





