function z0dinput = zerod_init_compass(mode_exp,shot,gaz,temps,z0dinput)

    langue                 =  'anglais';
    % cas donnees TS
    z0dinput.exp0d         = [];
    % 1 - numero du choc
    if isempty(gaz) || isempty(shot)

        try 
            numshot = fix(evalin('base','param.from.shot.num'));
        catch
            numshot = 6023;
        end 

        switch langue 
        case 'francais'
            prompt={'Numero du choc :','Charge du gaz principal :'};
            def={sprintf('%d',numshot),'1'};
            dlgTitle='Lecture des donnees de COMPASS';
        otherwise
            prompt={'shot number :','charge of main gas:' };
            def={sprintf('%d',numshot),'1'};
            dlgTitle='Access to COMAPSS data';
        end
        lineNo=1;
        answer=zinputdlg(prompt,dlgTitle,lineNo,def);
        if isempty(answer)
	    disp('Operation cancelled');
	    z0dinput = [];
            return
        end
        shot  = str2num(answer{1});
        gaz  = str2num(answer{2});

    end

    % lecture des donnees
    % connect to CDB
    cdb = cdb_client();
    % get Ip
    Ip_efit = cdb.get_signal(compass_str_id('R_bound/EFIT',shot));
    t_efit = Ip_efit.time_axis.data * 1e-3;
    sig = cdb.get_signal(compass_str_id('I_plasma',shot));
    signe_ip = sign(mean(sig.data));
    ip = abs(sig.data);
    tip = sig.time_axis.data * 1e-3;
    if isempty(ip)
        disp('No plasma')
        return
    end
    % Ip limit
    % ip_mask = tip >=0 & ip > 1.1*ip(end);
    % use EFIT times
    ip_mask = tip >= t_efit(1) & tip <= t_efit(end);
    tt = tip(ip_mask);
    if isempty(tt)
        disp('No plasma current');
        return
    end


    if ~isempty(temps)
        % on prend le vecteur donnee en entree
    elseif gaz < 0
        gaz = abs(gaz);
        temps = tt;
        temps(find(diff(temps)<=0)) =[];
    else
        t_max = min(max(tt), t_efit(end));
        t_min = max(min(tt), t_efit(1));
        dt    = max(mean(diff(tt)), 0.001);
        temps = (t_min:dt:t_max)';

        if length(temps) > 1001
            temps = linspace(t_min,t_max,1001)';
        end 
    end

    % ip_filt = medfilt1(tt, ip(ip_mask), floor(length(tt) / length(temps) / 2), 1024);
    ip_filt = sgolayfilt(ip(ip_mask), 1, 2*floor(length(tt) / length(temps) / 2) + 1);    
    z0dinput.cons.temps    = temps;
    z0dinput.cons.ip       = max(1,interp10d(tt,ip_filt,temps,'nearest'));
    
    % default values
    vt    = ones(size(temps));
    z0dinput.cons.flux    = 0      .* vt;
    z0dinput.cons.nbar    = 5e19 .* vt;
    z0dinput.cons.picrh   = 0    .* vt;
    z0dinput.cons.plh     = 0    .* vt;
    z0dinput.cons.pnbi    = 0    .* vt;
    z0dinput.cons.pecrh   = 0  .* vt;
    z0dinput.cons.zeff    = 2      .* vt;
    z0dinput.cons.xece    = 0      .* vt;
    z0dinput.option.li    = 1;
    z0dinput.cons.hmore   = 1      .* vt;

    z0dinput.geo.a        = 0.2   .* vt;
    z0dinput.geo.R        = 0.65    .* vt; 
    z0dinput.geo.K        = 1      .* vt;  
    z0dinput.geo.d        = 0      .* vt;
    z0dinput.geo.z0       = 0      .* vt;
    z0dinput.geo.b0       = 1.2    .* vt;
    z0dinput.geo.vp       = [];
    z0dinput.geo.sp       = [];
    z0dinput.geo.sext     = [];


    % z0dinput.cons.flux = signe_ip .* zechan(tfluxbord,fluxbord,temps,'nearest')./2./pi; 

    % equilibrium
    % sig = cdb.get_signal(compass_str_id('minor_radius/EFIT',shot))
    % tplasma = sig.time_axis.data * 1e-3;
    % z0dinput.geo.a = interp10d(tplasma,sig.data,temps,'nearest');
    % sig = cdb.get_signal(compass_str_id('R_mag_axis/EFIT',shot))
    % z0dinput.geo.R = interp10d(tplasma,sig.data,temps,'nearest');
    % sig = cdb.get_signal(compass_str_id('Z_mag_axis/EFIT',shot))
    % z0dinput.geo.z0 = interp10d(tplasma,sig.data,temps,'nearest');
    % sig = cdb.get_signal(compass_str_id('elongation_lcfs/EFIT',shot))
    % z0dinput.geo.K = interp10d(tplasma,sig.data,temps,'nearest');
    % triangularity
    % sig = cdb.get_signal(compass_str_id('triangularity_lower_lcfs/EFIT',shot))
    % z0dinput.geo.d = interp10d(tplasma,sig.data,temps,'nearest');
    % average upper + lower
    % sig = cdb.get_signal(compass_str_id('triangularity_upper_lcfs/EFIT',shot))
    % z0dinput.geo.d = 0.5 * (z0dinput.geo.d + interp10d(tplasma,sig.data,temps,'nearest'));
    % separatrix
    sig = cdb.get_signal(compass_str_id('R_bound/EFIT',shot));
    Rsepa = sig.data;
    t_efit = sig.time_axis.data * 1e-3;
    sig = cdb.get_signal(compass_str_id('Z_bound/EFIT',shot));
    Zsepa = sig.data;

    z0dinput = z0dsepageo(z0dinput, Rsepa, Zsepa, t_efit);
    
    % B vacuum
    sig = cdb.get_signal(compass_str_id('B_vac_R_geom/EFIT',shot));
    rb0 = sig.data;
    sig = cdb.get_signal(compass_str_id('R_geom_axis/EFIT',shot));
    rb0 = abs(rb0 .* sig.data);
    z0dinput.geo.b0 = interp10d(t_efit, rb0, temps, 'nearest') ./ z0dinput.geo.R ;
    sig = cdb.get_signal(compass_str_id('li/EFIT',shot));
    li_mask = isfinite(sig.data);
    z0dinput.exp0d.li     = interp10d(t_efit(li_mask), sig.data(li_mask), temps, 'nearest');
    z0dinput.option.li    = mean(z0dinput.exp0d.li);
    sig = cdb.get_signal(compass_str_id('psi_lcfs/EFIT',shot));
    z0dinput.cons.flux   = signe_ip .* interp10d(t_efit,sig.data,temps,'nearest');
    % TODO
    z0dinput.geo.vp       = [];
    z0dinput.geo.sp       = [];
    z0dinput.geo.sext     = [];
    % TODO 
    filter_ne = true;
    try
        % elongation + non-linearity corrections included
        sig = cdb.get_signal(compass_str_id('MARTE_NODE.PositionDataCollection2.CorrectedDensity',shot));
        nl = sig.data * 1e19;
        nl_fact = 1;
    catch
        try
            sig = cdb.get_signal(compass_str_id('MWA3_Interferometer_Unambiguous_PD_455kHz_with_6xAD8302',shot));
            % sig = cdb.get_signal(compass_str_id('MWB1_Interferometer_Unambiguous_PD_133MHz_triggered',shot));
            % this signal is multiplied by 0.4 m to get volume averaged ne
            % http://wiki.tok.ipp.cas.cz/mediawiki/index.php/Microwave_interferometer
            nl = sig.data * 0.4 * 1e19;
            nl_fact = 1 ./ z0dinput.geo.a ./ 2 ./ z0dinput.geo.K;
        catch
            sig = cdb.get_signal(compass_str_id('ne_avg/THOMSON',shot));
            nl = sig.data;
            nl_fact = 1;
            filter_ne = false;
            warning('using TS ne_avg')
        end
    end
    tnl = sig.time_axis.data * 1e-3;
    if filter_ne
        nl = sgolayfilt(nl, 1, 2*floor(length(nl) / length(temps) / 2) + 1);
        nl = nl - nl(1);
    end
    % remove the offset
    if isempty(nl)
        warning('No density measurment available')
    elseif all(nl(:) <= 3e17)
        warning('Density too low')
    else
        ind = find(diff(tnl) <=0);
        while(~isempty(ind) && ~isempty(tnl))
          tnl(ind) = [];
          nl(ind,:)=[];
          ind = find(diff(tnl) <=0);
        end
    end
    nbar                   = interp10d(tnl, nl, temps,'nearest') .* nl_fact;
    % securite anti Nan
    nbarm                  = mean(nbar(isfinite(nbar)));
    nbar(~isfinite(nbar))  = nbarm;
    nbar                   = max(1e17,nbar);
    z0dinput.cons.nbar     = nbar;
    % Zeff
    disp('No Zeff data, using Zeff = 2')
    zeffm = 2 .* ones(size(temps));
    if gaz == 2
        z0dinput.cons.zeff    = max(2,min(16,zeffm));
    else
        z0dinput.cons.zeff    = max(1,min(7,zeffm));
    end
        
    % TODO NBI

    if gaz == 1
        z0dinput.option.gaz = 2;
    else
        z0dinput.option.gaz = 4;
    end
    z0dinput.option.zmax  = 8;


    % donnees experimentales
    z0dinput.exp0d.temps = temps;
    disp('Reading W_dia')
    try
        sig = cdb.get_signal(compass_str_id('diamagnet_PP_EnergyBT/MAGNETICS',shot));
    catch
        sig = cdb.get_signal(compass_str_id('W/EFIT',shot));
    end
    sig.data = sgolayfilt(sig.data, 1, 2*floor(length(sig.data) / length(z0dinput.exp0d.temps) / 2) + 1);
    z0dinput.exp0d.w = interp10d(sig.time_axis.data*1e-3, sig.data, z0dinput.exp0d.temps, 'nearest');
    z0dinput.exp0d.dwdt  = pdederive(temps,z0dinput.exp0d.w,2,2,1,1); 
    
    try
       sig = cdb.get_signal(compass_str_id('Prad_insep',shot));
       sig.data = sgolayfilt(sig.data, 1, 2*floor(length(sig.data) / length(temps) / 2) + 1);
       z0dinput.exp0d.prad = interp10d(sig.time_axis.data*1e-3, sig.data, temps, 'nearest');
    catch
        disp('Prad_insep signal not foundm trying imrisek')
        prad_file = ['/compass/Shared/Exchange/imrisek/MATLAB/Tomography/Prad_', num2str(shot), '.mat'];
        if exist(prad_file, 'file')
            disp(['Reading Prad from', prad_file])
            Pr = load(prad_file);
            z0dinput.exp0d.prad = sgolayfilt(Pr.Prad_InSep, 1, 2*floor(length(Pr.Prad_InSep) / length(z0dinput.exp0d.temps) / 2) + 1);
            z0dinput.exp0d.prad = interp10d(Pr.tvec, z0dinput.exp0d.prad, z0dinput.exp0d.temps, 'nearest');
        else
            disp('No Prad data')
        end
    end
    try
       sig = cdb.get_signal(compass_str_id('P_sep',shot));
       sig.data = sgolayfilt(sig.data, 1, 2*floor(length(sig.data) / length(temps) / 2) + 1);
       z0dinput.exp0d.ploss = interp10d(sig.time_axis.data*1e-3, sig.data, temps, 'nearest');
    catch
        disp('No P_sep')
    end
    % TODO Thomson                  
    disp('Reading U_loop')
    try
        sig = cdb.get_signal(compass_str_id('U_loop_01',shot));
        disp('Got U_loop_01')
    catch
        sig = cdb.get_signal(compass_str_id('U_loop_04',shot));
        disp('Got U_loop_04')
    end
    sig.data = sgolayfilt(sig.data, 1, 2*floor(length(sig.data) / length(temps) / 2) + 1);
    z0dinput.exp0d.vloop = interp10d(sig.time_axis.data*1e-3, sig.data, temps,'nearest');

    z0dinput.exp0d.ip    = z0dinput.cons.ip;

    z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
    
    try
        sig = cdb.get_signal(compass_str_id('ne/THOMSON',shot));
        % TODO use a better centreal density
        z0dinput.exp0d.ne0 = interp10d(sig.time_axis.data*1e-3, sig.data(:, 1), temps, 'nearest');
    catch
        disp('no ne from Thomson')
    end
    try
        sig = cdb.get_signal(compass_str_id('Te/THOMSON',shot));
        % TODO use a better centreal density
        z0dinput.exp0d.te0 = interp10d(sig.time_axis.data*1e-3, sig.data(:, 1), temps, 'nearest');
    catch
        disp('no Te from Thomson')
    end

    z0dinput.exp0d.edgeflux    = z0dinput.cons.flux;
    z0dinput.machine     = 'COMPASS';
    z0dinput.shot        = shot;

end


function str_id = compass_str_id(signal_name, shot)
    str_id = [signal_name,':',num2str(shot)];
end


function yi = interp10d(x,y,xi,methode)

    if (size(x,1)>1) && (size(x,2)==1)
        indnok = 0;
        nb     = 100;
        while (~isempty(indnok)) && (nb >0)
            indnok = find(diff(x)<=0);
            if ~isempty(indnok)
                x(indnok) = [];
                y(indnok,:) = [];
        
            end
            nb = nb - 1;
        end
    end
    yi = interp1(x,y,xi,methode);

end

function data = decimate_data(t, data, t_out)

    % fast resample data using a running average filter

%    interp1(data_out.U_loop.time_axis.data(t_resample/2:t_resample:end-t_resample/2)/1e3,compass_resample_data(data_out.U_loop.data,t_resample),shotime,'linear');%

    % enforce column vector
    % TODO adjust automatically for row/column vectors
    data = reshape(data, [], 1);
    ww = floor(min(diff(t_out)) / min(diff(t)));
    % cut the last data points
    % TODO properly treat the last data points
    n_mod = mod(length(data), ww);
    if n_mod > 0
        last_data = mean(data(end-n_mod:end));
    else
        last_data = [];
    end
    data = data(1:end-n_mod);
    % reshape the data and calculate the mean
    data = reshape(data, ww, []);
    data = mean(data, 1);

end
