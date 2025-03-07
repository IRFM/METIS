%============================================================
% acces au donnees de WEST
%============================================================
function z0dinput = zerod_init_west(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
% cas donnees TS
z0dinput.exp0d         = [];

% test access to tools_dc
if exist('run_gui_dc') ~= 2
        warndlg('WEST data access requires IMAS framework and tools_dc functions. You must load the module tools_dc before launching METIS. This feature is not yet available under Windows or Mac OS/X','Access to data not avialable');
        z0dinput = [];
        return
end


% 1 - numero du choc
if isempty(gaz) | isempty(shot)
    
    try
        numshot = fix(evalin('base','param.from.shot.num'));
    catch
        try
            numshot = tsdernier_choc;
            rignitron = tsbase(numshot,'rignitron');
            if isempty(rignitron)
                numshot = numshot - 1;
            end
        catch
            numshot = 50123;
        end
    end
    
    prompt={'shot number :'};
    def={sprintf('%d',numshot)};
    dlgTitle='Access to WEST data';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
        z0dinput = [];
        return
    end
    shot  = str2num(answer{1});
    %gaz  = str2num(answer{2});
    if (shot < 50000) && (shot > 0);
        warndlg('Shot number < 50000 is Tore Supra shot','This is not a WEST shot');
        z0dinput = [];
        return
    end
    
end

% lecture des donnees
t_igni = tsbase(abs(shot),'rignitron');
[ip,tip] = tsbase(abs(shot),'smag_ip');
if isempty(ip)
    disp('No plasma')
    z0dinput = [];
    return
elseif any(ip > 1e7)
    ip = ip ./ 1e3;
end
if shot < 0
    tt = tip(tip >=0);
else
    tt = tip(tip >=0 & ip > 3);
end
if isempty(tt)
    disp('No plasma current');
    z0dinput = [];
    return
end
% gaz (version for  C1
[gt,tgt] = tsbase(abs(shot),'GTENSION');
tgt = tgt(:,1);
indok = (tgt> 10);
if any(gt(indok,7) > 30) || any(gt(indok,7) > 30)
    gaz = 2
else
    gaz = 1
end


if ~isempty(temps)
    % on prend le vecteur donnee en entree
    if shot < 0
        temps = tt((tt >= min(temps)) & (tt <= max(temps)));
    else
        temps = temps(:);
    end
elseif shot < 0
    gaz = abs(gaz);
    temps = tt;
    temps(find(diff(temps)<=0)) =[];
else
    dt    = max(mean(diff(tt)),0.1);
    temps = (min(tt):dt:max(tt))';
    
    if length(temps) > 1001
        temps = linspace(min(tt),max(tt),1001)';
    elseif length(temps) < 4
        temps = linspace(min(tt),max(tt),11)';
    end
end
shot = abs(shot);

%
% probleme tip
%
indtip = find(diff(tip) <= 0);
if ~isempty(indtip)
    tip(indtip) = [];
    ip(indtip) = [];
end
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = max(1,interp10d(tip,ip,temps,'nearest') .* 1e3);
ip   = z0dinput.cons.ip;

% flux au bord
% must be checked !
[fluxbord,tfluxbord] = tsbase(shot,'gmag_flxcs%3');
z0dinput.cons.flux    = interp10d(tfluxbord,fluxbord,temps,'nearest') ./ 2 ./ pi;


[gplasma,tplasma]       = tsbase(shot,'gmag_geom');
tplasma                = tplasma(:,1);
%
% probleme tplasma
%
indtplasma = find(diff(tplasma) <= 0);
if ~isempty(indtplasma)
    tplasma(indtplasma) = [];
    gplasma(indtplasma,:) = [];
end

% conversion to m from mm
gplasma = gplasma ./ 1000; 

z0dinput.geo.a         = interp10d(tplasma,gplasma(:,3),temps,'nearest');
z0dinput.geo.R         = interp10d(tplasma,gplasma(:,1),temps,'nearest');
z0dinput.geo.z0        = interp10d(tplasma,gplasma(:,2),temps,'nearest');
z0dinput.geo.K         = interp10d(tplasma,gplasma(:,5),temps,'nearest');
z0dinput.geo.d         = interp10d(tplasma,mean(gplasma(:,6:7),2),temps,'nearest');

if any(z0dinput.geo.R < 1)
    [gbary,tbary]        = tsbase(shot,'gmag_bary');
    tbary                = tbary(:,1);
    % parois
    [rp,zp] = west_limiter;
    if (rp(end) == rp(1)) &&(zp(end) == zp(1))
        %rien
    else
        % on ferme
        rp(end+1) =rp(1);
        zp(end+1) =zp(1);
    end
    
    lrz = cumtrapz(sqrt(diff(rp).^ 2 + diff(zp).^ 2));
    lrz(end+1) = lrz(end) + sqrt((rp(end) - rp(2)) .^ 2  + (zp(end) - zp(2)) .^ 2);
    lint = union(linspace(min(lrz),max(lrz),1001)',lrz);
    rpint = interp1(lrz,rp,lint,'linear','extrap')';
    zpint = interp1(lrz,zp,lint,'linear','extrap')';
    
    z0dinput.geo.R         = interp10d(tbary,gbary(:,1),temps,'nearest');
    z0dinput.geo.z0        = interp10d(tbary,gbary(:,2),temps,'nearest');
    z0dinput.geo.K         = ones(size(temps));
    z0dinput.geo.d         = zeros(size(temps));
    
    ve = ones(size(rpint));
    vt = ones(size(temps));
    dd     = sqrt((vt * rpint -  z0dinput.geo.R * ve) .^ 2 + (vt * zpint -  z0dinput.geo.z0 * ve) .^ 2);
    mask   = double(dd == (min(dd,[],2) * ve));
    aminor = sum(dd .* mask,2) ./ max(1,sum(mask,2));
    
    z0dinput.geo.a = aminor;
    
end


% separatrice
[grho,tgrho]   = tsbase(shot,'gmag_bnd');
if ~isempty(grho)
    tgrho          = tgrho(:,1);
    indnok   = find(all(grho <= 0,2));
    grho(indnok,:) = [];
    tgrho(indnok) = [];
end
with_equinox = 0;
if ~isempty(grho)
    
    grho     = interp1(tgrho,grho,temps,'nearest','extrap');
    Rext             = grho(:,1:2:end-1);
    Zext             = grho(:,2:2:end);
    Rext(:,end)    = Rext(:,1);
    Zext(:,end)    = Zext(:,1);32.2764
    
    % recalcul des parametres pour verification
    ve    = ones(1,size(Rext,2));
    rmin  = min(Rext,[],2);
    rmax  = max(Rext,[],2);
    ra    = max(1,0.5 .* (rmin + rmax));
    a     = max(0.01,0.5 .* (rmax - rmin));
    zmin  = min(Zext,[],2);
    zmax  = max(Zext,[],2);
    za    = (zmin + zmax) ./ 2;
    b     = 0.5 .* (zmax - zmin);
    k     = max(0.5,b ./ a);
    mask1 = (Zext == (max(Zext,[],2)*ve));
    mask2 = (Zext == (min(Zext,[],2)*ve));
    
    rzmax = max(Rext .* mask1,[],2);
    rzmin = max(Rext .* mask2,[],2);
    cl    = ra - rzmin;
    cu    = ra - rzmax;
    d     = (cl+cu) ./2 ./ a;
    z0dinput.exp0d.Rsepa = Rext;
    z0dinput.exp0d.Zsepa = Zext - za * ones(1,size(Zext,2));
    z0dinput.geo.a       = a;
    z0dinput.geo.R       = ra;
    z0dinput.geo.z0      = za;
    z0dinput.geo.K       = k;
    z0dinput.geo.d       = d;
else
    % try with equinox data
    [Rsepa_eqx,Zsepa_eqx,ip_eqx,li_eqx,betap_eqx,qax_eqx,q95_eqx,wmhd_eqx,R0_eqx,Z0_eqx, ...
        Rax_eqx,Zax_eqx,aminor_eqx,K_eqx,d_eqx,dlow_eqx,dup_eqx,vloop_eqx,flux_edge_eqx]=read_equinox_data(shot,temps);
    z0dinput.exp0d.Rsepa = Rsepa_eqx;
    z0dinput.exp0d.Zsepa = Zsepa_eqx - Z0_eqx * ones(1,size(Zsepa_eqx,2));
    z0dinput.geo.a       = aminor_eqx;
    z0dinput.geo.R       = R0_eqx;
    z0dinput.geo.z0      = Z0_eqx;
    z0dinput.geo.K       = K_eqx;
    z0dinput.geo.d       = d_eqx;
    with_equinox = 1;
    
    % benchmark
    %figure;plot(z0dinput.cons.temps,z0dinput.cons.ip / 1e3,'b',temps,ip_eqx / 1e3,'r');
    %legend('Vacth','Equinox');
    %xlabel('time (s)')
    %ylabel('I_p (kA)');
    
end

% if available use NICE with polarimetry data
[Rsepa_eqx,Zsepa_eqx,ip_eqx,li_eqx,betap_eqx,qax_eqx,q95_eqx,wmhd_eqx,R0_eqx,Z0_eqx, ...
        Rax_eqx,Zax_eqx,aminor_eqx,K_eqx,d_eqx,dlow_eqx,dup_eqx,vloop_eqx,flux_edge_eqx]=read_nice_data(shot,temps);
% replace data for time where NICE data are available
if ~isempty(Rsepa_eqx) && ~isempty(Zsepa_eqx)
    Rsepa_vacth = z0dinput.exp0d.Rsepa;
    Zsepa_vacth = z0dinput.exp0d.Zsepa;
    z0dinput.exp0d.Rsepa = NaN * ones(length(temps),size(Rsepa_eqx,2));
    z0dinput.exp0d.Zsepa = NaN * ones(length(temps),size(Zsepa_eqx,2));    
    fprintf('LCFS from NICE:');
    for k=1:length(temps)
        if all(isfinite(Rsepa_eqx(k,:))) && all(isfinite(Zsepa_eqx(k,:))) && ...
                all(Rsepa_eqx(k,:) > 1.5) && all(abs(Zsepa_eqx(k,:)) < 1)
            z0dinput.exp0d.Rsepa(k,:) = Rsepa_eqx(k,:);
            z0dinput.exp0d.Zsepa(k,:) = Zsepa_eqx(k,:) - Z0_eqx(k) * ones(1,size(Zsepa_eqx,2));
            z0dinput.geo.a(k)       = aminor_eqx(k);
            z0dinput.geo.R(k)       = R0_eqx(k);
            z0dinput.geo.z0(k)      = Z0_eqx(k);
            z0dinput.geo.K(k)       = K_eqx(k);
            z0dinput.geo.d(k)       = d_eqx(k);
            with_equinox = 1;
            fprintf('o');
        else
            [z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:)] =  ...
                reshape_LCFS(Rsepa_vacth(k,:),Zsepa_vacth(k,:),size(Rsepa_eqx,2));
            
            fprintf('x');
        end
        %figure(21);plot(Rsepa_vacth(k,:),Zsepa_vacth(k,:),'.r',z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:),'b');pause(0.1);
    end
    fprintf('\n');
end


[gitor,titor]          = tsbase(shot,'gmag_itor');
itor                   = gitor(:,1);
ditor                  = gitor(:,2);
titor                  = titor(:,1);
ind        	       = find(titor(:,1) >=0);
itor                   = itor(ind);
ditor                  = ditor(ind);
titor                  = titor(ind);
[b,a]                  = butter(11,0.1);
itor                   = filtfilt(b,a,itor);
ditor                   = filtfilt(b,a,ditor);
rb0                    = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi;
rb0                    = interp10d(titor,rb0,temps,'nearest');
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
itor                   = interp10d(titor,itor,temps,'nearest');

% energy mag
Ltor     = 2 .* 600e6 ./ 1400 .^ 2;
phi_tor  = Ltor .* (ditor  - ditor(1));
phi_tor  = interp10d(titor,phi_tor,temps,'nearest');

% density
% density
try
    [tmix,lidmix,lidbar,avdensity,netot,cert,delta,phidcn,phih2o,tfjdcn,fjdcn,tfjh2o,fjh2o,struct] = TINTimasinterfero(shot,0);
%      cert(2) = -2;
%      cert(1) = -2;
%      cert(12) = -2;
%      cert(11) = -2;
    cert(4:9) = -2;
    cert(1:2) = -2;
    cert(14:19) = -2;
    cert(11:12) = -2;
catch
    lidmix = [];
    avdensity = NaN * tmag;
end
if ~isempty(lidmix) && ~isempty(lidbar) && any(cert(1:10) > -2) && any(cert(11:20) > -2)
    Z2500=[-0.78 -0.8 0.03 0.3 0.20 0.25 0.13 -0.2 -0.37 -0.04];
    %t_igni = tsbase(shot,'rignitron');
    trt = tmix;
    indok = find(cert(1:10) > -2);
    ncorde = max(lidmix(:,indok),[],2);
    indok = find(cert(11:20) > -2);
    nl     =  max(lidbar(:,indok),[],2).* 1e19;
    ovni.tne = tmix';
    indok = find(cert(1:10) > -2);
    ovni.ne  = sum(lidmix(:,indok),2);
    ovni.decl = double(sqrt(sgolayfilt((ovni.ne-sgolayfilt(ovni.ne,1,5)).^2,1,5)) > 0.05);
    ovni.decl(ovni.decl == 0) = NaN;
    if isempty(avdensity)
        nemoy = NaN * temps;
        ane  =  NaN * temps;
    else
        avdensitym                  = mean(avdensity(isfinite(avdensity)));
        avdensity(~isfinite(avdensity))  = avdensitym;
        avdensity                   = max(1e-3,avdensity);
        nemoy = interp1(trt,avdensity,temps,'nearest',NaN) .* 1e19;
        ncordem                  = mean(ncorde(isfinite(ncorde)));
        ncorde(~isfinite(ncorde))  = ncordem;
        ncorde                   = max(1e-3,ncorde);
        ane       =  interp1(trt,ncorde,temps,'nearest',NaN) .* 1e19 ./ 2 ./ z0dinput.geo.a ./ nemoy;
    end
    nlm                  = mean(nl(isfinite(nl)));
    nl(~isfinite(nl))  = nlm;
    nl                   = max(1e17,nl);
    nbar    = interp10d(trt,nl,temps,'nearest');
    z0dinput.cons.nbar     = nbar;
    ne0     = nemoy .* ane;
    
else
    avdensity = NaN * temps;
    [LIDRT,trt,rrt]=tsbase(shot,'GINTLIDRT');
    if isempty(LIDRT)
        disp('No density measurment available')
        z0dinput = [];
        return
    end
    tnl = trt;
    CHANNEL_ON=tsmat(shot,'DINTRPOL;LIDRT;CHANNEL');
    if isempty(CHANNEL_ON)==1
        CHANNEL_ON=ones(1,10);
    end
    stdnl = std(LIDRT,0,1);
    mendnl = mean(LIDRT((end-10):end,:),1);
    CHANNEL_ON = CHANNEL_ON .* (abs(mendnl) < stdnl);
    if all(CHANNEL_ON == 0)
        CHANNEL_ON=tsmat(shot,'DINTRPOL;LIDRT;CHANNEL');
    end
    if isempty(CHANNEL_ON)==1
        CHANNEL_ON=ones(1,10);
    end
    LIDRT = sgolayfilt(LIDRT,1,5) .* 1e19;
    indok = find(CHANNEL_ON);
    nl = LIDRT(:,indok);
    
    
    if isempty(nl)
        disp('No density measurment available')
        z0dinput = [];
        return
    elseif all(nl(:) <= 3e17)
        disp('No density measurment available')
        z0dinput = [];
        return
    else
        tnl=tnl(:,1);
        ind = find(diff(tnl) <=0);
        while(~isempty(ind) & ~isempty(tnl))
            tnl(ind) = [];
            nl(ind,:)=[];
            ind = find(diff(tnl) <=0);
        end
    end
    nl = nl - mean(nl((tnl <0) & (tnl > -2)));
    nbar                   = interp10d(tnl,max(nl,[],2),temps,'nearest')./ 2 ./ z0dinput.geo.a;
    % securite anti Nan
    nbarm                  = mean(nbar(isfinite(nbar)));
    nbar(~isfinite(nbar))  = nbarm;
    nbar                   = max(1e17,nbar);
    z0dinput.cons.nbar     = nbar;
    
    [nemoy,tne]            = tsbase(shot,'snmoy');
    if ~isempty(nemoy)
        nemoy                   = interp10d(tne,nemoy,temps,'nearest');
        piqne                  = max(1.1,min(5,interp10d(tne,piqne,temps,'nearest')));
        ane                    = piqne - 1;
        ne0                    = nemoy .* piqne;
        avdensity              = nemoy;
    else
        [nemoy,tne]            = tsbase(shot,'strnmoy');
        if ~isempty(nemoy)
            nemoy                   = interp10d(tne,nemoy,temps,'nearest');
            piqne                  = 1.5 .* ones(size(itor));
            ane                    = piqne - 1;
            ne0                    = nbar .*  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi));
            avdensity              = nemoy;
        else
            nemoy = nbar;
            piqne = 1.5 .* ones(size(itor));
            ane                    = piqne - 1;
            ne0                    = nbar .*  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi));
            avdensity              = nemoy;
        end
    end
end

%GMAG_BILAN  	: VACTH: Energy and power. 1-6: ohmic, LHCD, ICRH, ECRH, Total, radiated power; 7: diamagnetic energy; 8: diamagnetic time constant; 9: magnetic energy; 10: MHD energy; 11: smoothed total power; 12: spare
[gbilan,tbilan]         = tsbase(shot,'gmag_bilan');
if ~isempty(gbilan) && ~all(gbilan(:) == 0);
    tbilan                 = tbilan(:,1);
    %
    % probleme tbilan
    %
    indtbilan = find(diff(tbilan) <= 0);
    if ~isempty(indtbilan)
        tbilan(indtbilan) = [];
        gbilan(indtbilan,:) = [];
    end
    z0dinput.cons.picrh    =  max(0,interp10d(tbilan,max(0,gbilan(:,3)),temps,'nearest')) .* 1e6;
    z0dinput.cons.plh      = max(0,interp10d(tbilan,gbilan(:,2),temps,'nearest')) .* 1e6;
    z0dinput.cons.pnbi     = zeros(size(temps));
    pecrh                  = zecrh(shot,temps);
    z0dinput.cons.pecrh    = max(0,sum(pecrh,2).*1e6);
    z0dinput.cons.pecrh    = z0dinput.cons.pecrh .* (z0dinput.cons.pecrh >=1e4);
    z0dinput.cons.hmore    = ones(size(temps));
    ptot                   =  max(0,interp10d(tbilan,gbilan(:,5),temps,'nearest')) .* 1e6;
    pohm                   =  interp10d(tbilan,gbilan(:,1),temps,'nearest') .* 1e6;
    wdia                   = interp10d(tbilan,gbilan(:,7),temps,'nearest');
    prad                   = interp10d(tbilan,gbilan(:,6),temps,'nearest');
    w                     = interp10d(tbilan,gbilan(:,9),temps,'nearest');
    
    
else
    [pfci,tpfci]=tsbase(shot,'spuiss');
    if isempty(pfci)
        [pfci,tpfci]=tsbase(shot,'gptra');
        if isempty(pfci)
            [pfci,tpfci]=tsbase(shot,'gpfci');
        end
        if ~isempty(pfci)
            tpfci = tpfci(:,1);
            pfci  = sum(pfci .* (pfci >=0),2);
        end
    end
    if ~isempty(pfci)
        z0dinput.cons.picrh    =  max(0,interp10d(tpfci,max(0,pfci),temps,'nearest')) .* 1e6;
    else
        z0dinput.cons.picrh    =  zeros(size(temps));
    end
    %[phyb,tphyb]=tsbase(shot,'GPHYB%3');
    [phyb,tphyb] = read_phyb(shot);
     %if isempty(phyb)
    %    [phyb,tphyb]=tsbase(shot,'GHYB%3');
    %end
    if ~isempty(phyb)
        tphyb = tphyb(:,1);
        phyb  = phyb(:,3);
        z0dinput.cons.plh      = max(0,interp10d(tphyb,phyb,temps,'nearest')) .* 1e6;
    else
        z0dinput.cons.plh      = zeros(size(temps));
    end
    z0dinput.cons.pnbi     = zeros(size(temps));
    z0dinput.cons.pecrh    = zeros(size(temps));
    z0dinput.cons.hmore    = ones(size(temps));
    [pohm,tpohm]=tsbase(shot,'spohm');
    if ~isempty(pohm)
        pohm      = max(0,interp10d(tpohm,pohm,temps,'nearest')) .* 1e6;
    else
        pohm      = zeros(size(temps));
    end
    ptot                   =  pohm + z0dinput.cons.picrh + z0dinput.cons.plh;
    
    if  with_equinox == 1
         wdia = wmhd_eqx;
         w    = wmhd_eqx;
    else
        wdia  = NaN * temps;
        w     = NaN * temps;
    end
    
end
%  try
%      [Prad,Pbulk,Pdiv,Pchan,trad]=radpwr(shot,0);
%      prad                   = interp10d(trad,Prad,temps,'nearest');
%  catch
%      prad  = NaN * temps;
%  end
% Prad
% try to read IMAS data first
bolometer = imas_west_get(shot,'bolometer');
if isempty(bolometer) || isempty(bolometer.time)
    try
	%[Prad,Pbulk,Pdiv,~,trad]=radpwr(shot,0);
	[Prad,Pbulk,Pdivb,Pdivh,~,trad]=radpwr_link(shot,0);
        prad                   = interp10d(trad,Pbulk,temps,'nearest');
        pradsol                = max(0,interp10d(trad,Prad - Pbulk,temps,'nearest'));
    catch
	prad     = NaN * time;
	pradsol  = NaN * time;
    end
else
	disp('using IMAS bolometer IDS from west imas_public run 0');
	prad     = interp1(bolometer.time - t_igni,bolometer.power_radiated_inside_lcfs,temps,'nearest') / 1e6; 
	pradsol  = max(0,interp1(bolometer.time - t_igni,bolometer.power_radiated_total - bolometer.power_radiated_inside_lcfs,temps,'nearest') / 1e6); 
end

% waiting for available Zeff
z0dinput.cons.zeff    = zeffscaling(z0dinput.cons.nbar./ 1e19,ptot ./ 1e6,ip ./ 1e6, ...
                         z0dinput.geo.a,z0dinput.geo.R,itor,gaz);


% try to read Zeff from treatement
summary   = imas_west_get(shot,'summary');
if ~isempty(summary) && ~isempty(summary.time) && ~isempty(summary.line_average.zeff.value)
    stime 	  = summary.time - t_igni;
    zeff_value = summary.line_average.zeff.value;
    zeff_value(zeff_value < 0) = NaN;
    if sum(isfinite(zeff_value)) >3
        indbad     = find(diff(stime) <= 0);
        if ~isempty(indbad)
            while  ~isempty(indbad)
                zeff_value(indbad) = [];
                stime(indbad) = [];
                indbad     = find(diff(stime) <= 0);
            end
        end
        if length(stime) > 3
            z0dinput.cons.zeff  	 = interp1(stime,zeff_value,temps,'nearest',NaN);
            indbad = find(~isfinite(z0dinput.cons.zeff) | (z0dinput.cons.zeff <= 1) |(z0dinput.cons.zeff  > 11 ));
            if ~isempty(indbad)
                if gaz == 2
                    z0dinput.cons.zeff(indbad)    = 3;
                else
                    z0dinput.cons.zeff(indbad)    = 1.5;
                end
            end
        end
    end
end



%GMAG_BELI   	: VACTH: Diamagnetic parameters. 1: beta dia; 2: li; 3: beta MHD; 4: beta N; 5: lik; 6: li dia
d0                    = 0.1 .* ones(size(temps));
[gqbeli,tgqbeli]      = tsbase(shot,'GMAG_BELI');
if ~isempty(gqbeli) && ~all(gqbeli(:) == 0)
    tgqbeli               = tgqbeli(:,1);
    li                    = interp10d(tgqbeli,abs(gqbeli(:,2)),temps,'nearest');
    z0dinput.option.li    = li(1);
    beta                  = interp10d(tgqbeli,abs(gqbeli(:,3)),temps,'nearest');
    
elseif with_equinox
        li                    = li_eqx;
        z0dinput.option.li    = li(find(li>0,1));
else
	li = NaN * temps;
end
if  with_equinox
      d0 = Rax_eqx - R0_eqx;
end
% donnee calculee dans le zerod
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];


% les parametres
z0dinput.option = map_option_west(z0dinput.option);

% fci
[pant,tpant]          = tsbase(shot,'gpuifci');
if ~isempty(pant)
    indant              = find(max(pant(:,1:3),[],1) > 0.3);
    if isempty(indant)
        indant            = 1;
    end
    frequence           = tsbase(shot,'sfreqfci');
    if ~isempty(frequence)
        z0dinput.option.freq = max(1,mean(frequence(indant)));
    end
end
z0dinput.option.mino     = 'H';
geom.a                   = z0dinput.geo.a;
geom.r0                  = z0dinput.geo.R;
geom.b0                  = z0dinput.geo.b0;
geom.d0                  = d0;
pos                      = geom.r0 + geom.a + 0.02;
[scenar,scenstr,Rres]    = scenarzerodfci(geom,z0dinput.option.freq,pos,z0dinput.option.mino);
if isfinite(Rres)
    z0dinput.option.fwcd = 0;
else
    z0dinput.option.fwcd = 2;
end

%  % rapport iso (does not work for WEST)
%  [rgaz,rdcx] = lit_nhnd_ts(shot);
%  if ~isempty(rdcx)
%      z0dinput.option.cmin = min(0.3,max(0.01,rdcx ./ (1 - rdcx)));
%      z0dinput.option.mino = 'H';
%  elseif ~isempty(rgaz)
%      z0dinput.option.cmin = min(0.3,max(0.01,rgaz ./ (1 - rgaz)));
%      z0dinput.option.mino = 'H';
%  end
nhonhpnd = min(1-eps,max(0,isotope_ratio_west(shot,temps,z0dinput.cons.picrh)));
z0dinput.option.cmin = nhonhpnd ./ (1 - nhonhpnd);
z0dinput.option.mino = 'H';


% reading hard x ray inversion from IMAS
% read HXR is available
try
    hxr=imas_west_get(shot,'hard_x_rays');
catch
    hxr = [];
end
hxr_lh = [];
if ~isempty(hxr) && (length(hxr.emissivity_profile_1d) > 0)
    for k=1:length(hxr.emissivity_profile_1d)
        if (hxr.emissivity_profile_1d{1}.lower_bound == 60) && (hxr.emissivity_profile_1d{1}.upper_bound == 80)
            hxr_lh = hxr.emissivity_profile_1d{1};
        end
    end
end



% si presence hybride
% read Hybrid power
%[phyb,thyb]=tsbase(shot,'GPHYB');
[phyb,thyb] = read_phyb(shot);
if isempty(phyb)
    [pfw1,tpfw1] = tsbase(shot,'SHYBPFORW1');
    [pfw2,tpfw2] = tsbase(shot,'SHYBPFORW2');
    [pref1,tpref1] = tsbase(shot,'SHYBPREFL1');
    [pref2,tpref2] = tsbase(shot,'SHYBPREFL2');
    if ~isempty(pfw1)
        phyb = 1e-3 .* cat(2,max(0,interp10d(tpfw1,sgolayfilt(pfw1-pref1,1,5),temps,'nearest')),max(0,interp10d(tpfw2,sgolayfilt(pfw2-pref2,1,5),temps,'nearest')));
        phyb = cat(2, phyb,sum(phyb,2));
        phyb(phyb <0) = 0;
        thyb = temps;
    else
        phyb = [];
        thyb = [];
    end
else
    thyb = thyb(:,1);
end
if ~isempty(phyb) && all(z0dinput.cons.plh == 0)
    z0dinput.cons.plh = interp10d(thyb,max(0,phyb(:,3)) .* 1e6,temps,'nearest');
end

% reading request n//0
[npar,tnpar] = tsbase(shot,'GNPARDEM');
if ~isempty(phyb) && ~isempty(npar)
    phyb_loc = interp10d(thyb(:,1),phyb(:,1:2),temps,'nearest');
    npar = interp10d(tnpar(:,1),npar,temps,'nearest');
    npar0          =  npar(:,1) .* phyb_loc(:,1) +  npar(:,2) .* phyb_loc(:,2);
    phybt          =  phyb_loc(:,1) +  phyb_loc(:,2);
    indok          =  find(isfinite(npar0));
    if isempty(indok)
        npar0 = 1.8;
    else
        npar0 = sum(npar0(indok)) ./ sum(phybt(indok));
    end
else
    npar0 = 1.8;
end
z0dinput.option.npar0 =npar0;
z0dinput.option.wlh = 0.58;
z0dinput.option.lhmode = 0;
z0dinput.option.etalh  = min(1,max(0.1,2.01 - 0.63 .* npar0));

if ~isempty(hxr_lh) && ~isempty(hxr_lh.emissivity)
    g3 = hxr_lh.emissivity';
    g3(abs(g3) > 1e40) = 0;
    g3(g3 < 0) = 0;
    t3 = hxr_lh.time(:) - t_igni;
    x3 = hxr_lh.rho_tor_norm(:)';
    mask_fusion = ones(size(t3)) *max(0,min(1,(0.6 - x3)*5)) ;
    g3 = mask_fusion .* g3;
    mask = all(g3 == 0,2);
    g3(mask,:)      = [];
    t3(mask) = [];
    indbad = find((diff(t3) <=0));
    while (~isempty(indbad) & ~isempty(t3))
        g3(indbad,:) = [];
        t3(indbad)   = [];
        indbad = find((diff(t3) <=0));
    end
    
    if ~isempty(g3)
        % change coordinate from rho_tor_norm to Lao
        % read NICE data
        equi = imas_west_get(shot,'equilibrium',0,1);
        if isempty(equi)
            disp('NICE + Polarimetry not available: try NICE with magnetic alone');
            equi = imas_west_get(shot,'equilibrium');
        end
        if ~isempty(equi)
            % equi Lao coordinate
            roa  = (equi.profiles_1d.r_outboard - equi.profiles_1d.r_inboard) / 2;
            roa  = roa ./ (max(roa,[],2) * ones(1,size(roa,2)));
            troa = (equi.time - t_igni)' * ones(1,size(roa,2)); 
            % interpolant from rho_tor_norm to Lao
            %figure;plot(equi.profiles_1d.rho_tor_norm',roa',[0,1],[0,1],'k:');drawnow
            Froa = scatteredInterpolant(troa(:),equi.profiles_1d.rho_tor_norm(:),roa(:),'natural','linear');
            % resmapling 
            x3   = min(1,max(0,Froa(t3 * ones(1,length(x3)),ones(length(t3),1) * x3)));
            %
            XDUR.v = g3;
            XDUR.t = t3;
            XDUR.x = x3;
            %
            % use mean coordinate for dlh & xlh computation
            x3m  = median(x3,1);
            g3s  = sum(g3,1);
            %xlh  = trapz(d3,(1-d3) .* d3.* g3s,2) ./trapz(d3,(1-d3) .* g3s,2);
            xlh  = x3m(max(find(g3s ==max(g3s))));
            dlh  = sqrt(trapz(x3m,x3m .^ 2 .* g3s,2) ./max(eps,trapz(x3m,g3s,2)) - xlh.^ 2);
            z0dinput.option.xlh = xlh;
            z0dinput.option.dlh = dlh;
            %
            % stored in experimental data
            %
            z0dinput.option.lhmode = 3;
            z0dinput.exp0d.XDURt = t3;
            z0dinput.exp0d.XDURx = x3;
            z0dinput.exp0d.XDURv = g3;
        end
    end
end


% lecture de ICRH WEST
[picrh_tot, tpicrh]  = tsbase(shot,'SICHPTOT');
if ~isempty(picrh_tot)
    z0dinput.cons.picrh = interp10d(tpicrh,max(0,picrh_tot)*1e3,temps,'nearest');
    freq_par              = tsmat(shot, 'DFCI;PILOTAGE;ICHFREQ');
    if ~ isempty(freq_par)
        [pq1, tq1]  = tsbase(shot,'SICHPQ1');
        [pq2, tq2]  = tsbase(shot,'SICHPQ2');
        [pq4, tq4]  = tsbase(shot,'SICHPQ2');
        freq = 0;
        pref = 0;
        if ~isempty(pq1)
            freq = freq + trapz(tq1,max(0,pq1) .* freq_par(1));
            pref = pref + trapz(tq1,max(0,pq1));
        end
        if ~isempty(pq2)
            freq = freq + trapz(tq2,max(0,pq2) .* freq_par(2));
            pref = pref + trapz(tq2,max(0,pq2));
        end
        if ~isempty(pq4) & 0    % for the momemnt there is a noise on the signal.
            freq = freq + trapz(tq4,max(0,pq4) .* freq_par(3));
            pref = pref + trapz(tq4,max(0,pq4));
        end
        freq = freq ./ max(eps,pref);
        if (freq > 10) && (freq < 100)
            z0dinput.option.freq = freq;
        end
    end
end

% lecture des donnees ecrh directement dans le diagnostic
[xika1,tbad]=tsbase(shot,'sika1');  		% gyrotron A1 cathode current
[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
if ~isempty(xika1)  & ~isempty(xHAUTTOR)
    [prxA1,tA1]=tsbase(shot,'spia1');   		% gyrotron A1 power
    [prxA2,tA2]=tsbase(shot,'spia2');   		% gyrotron A2 power
    [sonde,tsonde]=tsbase(shot,'sonderf');   	% RF probe on reflectometer
    %[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
    [xHAUTPOL,tHAUTPOL]=tsbase(shot,'SHAUTPOL');		% Poloidal injection anlge top mirror (A1) hysteresis corrected
    [xMILTOR,tMILTOR]=tsbase(shot,'SMILTOR');		% Toroidal injection anlge central mirror (A2) hysteresis corrected
    [xMILPOL,tMILPOL]=tsbase(shot,'SMILPOL');		% Poloidal injection anlge central mirror (A2) hysteresis corrected
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % powers in kW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    factA1=1; 		% Calibration factor to be applied to  A1
    factA2=1; 		% Calibration factor to be applied to  A2
    PA1=max(0,1e6*prxA1/factA1);  % A1 power in W
    PA2=max(0,1e6*prxA2/factA2);	% A2 power in W
    Ptot=PA1+PA2;		% total power
    
    % position
    rA2  = 3.5300;
    rA1  = 3.5300;
    zA2  = 0;
    zA1  = 0.2000;
    
    % puissance
    z0dinput.cons.pecrh      = max(0,interp10d(tbad,Ptot,temps,'nearest'));
    z0dinput.cons.pecrh      = z0dinput.cons.pecrh .* (z0dinput.cons.pecrh >=1e4);
    % position
    phi_tor_a1 = interp10d(tbad,xHAUTTOR,temps,'nearest')./180*pi;
    phi_tor_a2 = interp10d(tbad,xMILTOR,temps,'nearest')./180*pi;
    phi_pol_a1 = interp10d(tbad,xHAUTPOL,temps,'nearest')./180*pi;
    phi_pol_a2 = interp10d(tbad,xMILPOL,temps,'nearest')./180*pi;
    
    % calcul de la position
    [xece1,nh1,apol1] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,d0,z0dinput.geo.a,rA1,zA1,phi_pol_a1,phi_tor_a1,z0dinput.geo.b0,118);
    [xece2,nh2,apol2] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,d0,z0dinput.geo.a,rA2,zA2,phi_pol_a2,phi_tor_a2,z0dinput.geo.b0,118);
    %[xece,nh,apol] = r2xece(r0,z0,d0,a,rant,zant,phi_pol,phi_tor,b0,freq_ghz)
    
    % moyenne
    pa1 = interp10d(tbad,PA1,temps,'nearest');
    pa2 = interp10d(tbad,PA2,temps,'nearest');
    z0dinput.cons.xece   = max(0,min(1,(xece1 .* pa1 + xece2 .* pa2) ./ max(1,pa1 + pa2)));
    z0dinput.option.sens = sign(trapz(temps,phi_tor_a1 .* pa1 + phi_tor_a2 .* pa2)) ;
    z0dinput.option.angle_ece = trapz(temps,abs(apol1) .* pa1 + abs(apol2) .* pa2) ./ max(1,trapz(temps,pa1 +pa2));
end

% try to see if shot is in helium
[heshot,nHeon1]= get_helium_west(abs(shot));
if heshot
   z0dinput.option.gaz = 4;
   z0dinput.option.frhe0 = 0;
else
   z0dinput.option.gaz = 2;
   z0dinput.option.frhe0 = nHeon1 ./ (1 + 2 .* nHeon1);
end
z0dinput.option.zmax  = 8;

% read gaz puff reference
[sdeb,tdeb]           = tsbase(shot,'sdeb');
try
  [void,flux_electrons] = debit_west_gas(shot,temps);
catch
  flux_electrons = [];
end 
if ~isempty(flux_electrons)
    z0dinput.cons.nbar  = z0dinput.cons.nbar + sqrt(-1) .* flux_electrons;
elseif ~isempty(tdeb)
    if z0dinput.option.gaz  == 4
        sdeb = (4.41e-4 .* 6.02214199e23) .* sdeb;   % conversion to electron/s from  Pa.m^3/s
    else
        sdeb = 2 .* (4.41e-4 .* 6.02214199e23) .* sdeb;   % conversion to electron/s from  Pa.m^3/s
    end
    sdeb = sgolayfilt(sdeb,1,51);
    sdeb(sdeb < 0) = 0;
    z0dinput.cons.nbar  = z0dinput.cons.nbar + sqrt(-1) .* interp10d(tdeb,sdeb,temps,'nearest');
end
% prefill pressure
[gptore,tptore]       = tsbase(shot,'gptore');
if ~isempty(gptore)
    pressure   = gptore(:,1) .* 10 .^ gptore(:,2);
    indpress = find((tptore(:,1) < 0) & (tptore(:,1) >= -0.2));
    z0dinput.option.p_prefill = mean(pressure(indpress));
end
% breakdown data
% z0dinput.option.berror = to be read from freebie_id


% donnees experimentales
z0dinput.exp0d.temps = temps;
z0dinput.exp0d.pin   = z0dinput.cons.picrh + z0dinput.cons.plh + z0dinput.cons.pecrh + pohm;

z0dinput.exp0d.ploss = ptot;
z0dinput.exp0d.zeff  = z0dinput.cons.zeff;
z0dinput.exp0d.ane  = ane;
z0dinput.exp0d.nem  = nemoy;
z0dinput.exp0d.ne0  = ne0;
z0dinput.exp0d.w     = w;
z0dinput.exp0d.wdia   = wdia;
z0dinput.exp0d.dwdt  = pdederive(temps,z0dinput.exp0d.w,2,2,1,1);
if ~isempty(gbilan) && ~all(gbilan(:) == 0)
    z0dinput.exp0d.taue  = interp10d(tbilan,gbilan(:,7),temps,'nearest');
else
    z0dinput.exp0d.taue  = NaN .* ones(size(temps));
end
z0dinput.exp0d.pw    = ptot;

try
    ece_ids = imas_west_get(shot,'ece');
catch
    try
      ece_ids = imas_west_get_old(shot,'ece');
    catch
      ece_ids = [];
    end
end
if ~isempty(ece_ids)
  time_ece = ece_ids.time - t_igni;
  te_ece = NaN * ones(length(time_ece),length(ece_ids.channel));
  R_ece  = NaN * ones(length(time_ece),length(ece_ids.channel));
  for k=1:length(ece_ids.channel)
      if ~isempty(ece_ids.channel{k}.t_e.data)
	  te_ece(:,k) = ece_ids.channel{k}.t_e.data(:);  
      end
      if  ~isempty(ece_ids.channel{k}.position.r.data)
 	  R_ece(:,k) = ece_ids.channel{k}.position.r.data(:);  
     end
      if  ~isempty(ece_ids.channel{k}.position.phi.data)
 	  phi_ece(:,k) = ece_ids.channel{k}.position.r.data(:);  
     end
  end
  tsh             = time_ece;
  gshte           = te_ece ./ 1000;
  gshte(~isfinite(gshte)) = 0;
  gshr           = R_ece;
  
end


if ~isempty(gshte)
    indtok             = find(tsh > max(min(temps),2) & tsh < (min(max(temps),max(tsh)) -1));
    indok              = find(all(gshte(indtok,:)>=0,1));
    gshr               = mean(gshr(:,indok),1);
    gshte              = gshte(:,indok);
    if ~isempty(gshte)
        indok              = find(gshr > 2.35 & gshr <= 2.45);
        if isempty(indok)
            d     = abs(gshr - 2.4);
            indok = max(find(d == min(d)));
        end
        if ~isempty(indok)
            te0                = max(gshte(:,indok),[],2) .* 1e3;
            tshc = tsh;
            te0c = te0;
            indr = find(diff(tsh)<=0);
            if ~isempty(indr)
                tshc(indr) = [];
                te0c(indr) = [];
            end

            z0dinput.exp0d.te0   = interp10d(tshc,te0c,temps,'nearest');
        end
    end
else
    [tethom,tthom,zthom]      = tsbase(shot,'GTETHOM');
    if isempty(tethom)
        [temic,tmic] = tsbase(shot,'gmictenv');
        if isempty(temic)
            [temic,tmic] = tsbase(shot,'gmicte');
            if isempty(temic)
                disp('No temperature measurement available')
            end
        end
        if ~isempty(temic)
            [rmic,tmic,voies,cert] = tsbase(shot,'gmicrnv');
            indtok             = find(tmic >2 & tmic < (max(tsh) -1));
            indok              = find(all(temic(indtok,:)>=0,1));
            rmic               = mean(rmic(:,indok),1);
            temic              = gshte(:,indok);
            indok              = find(rmic > 2.35 & rmic <= 2.45);
            if isempty(indok)
                d     = abs(rmic - 2.4);
                indok = max(find(d == min(d)));
            end
            te0                = mean(temic(:,indok),2) .* 1e3;
            tshc = tmic;
            te0c = te0;
            indr = find(diff(tsh)<=0);
            if ~isempty(indr)
                tshc(indr) = [];
                te0c(indr) = [];
            end
            z0dinput.exp0d.te0   = interp10d(tshc,te0c,temps,'nearest');
        else
            [tefab,tfab] = tsbase(shot,'gtefab');
            if isempty(tefab)
                disp('No temperature measurement available')
            else
                z0dinput.exp0d.te0   = interp10d(tfab,max(tefab .* (tefab >0),[],2),temps,'nearest').* 1e3;
            end
            
            
        end
    else
        
        indok              = find(abs(zthom) == min(abs(zthom)));
        te0                = mean(tethom(:,indok),2) .* 1e3;
        z0dinput.exp0d.te0   = interp10d(tthom,te0,temps,'nearest');
    end
end 


% te0 from eceIDS is available
if ~isempty(ece_ids) && isfield(ece_ids,'t_e_central') && ~isempty(ece_ids.t_e_central.data)
    if sum(ece_ids.t_e_central.data > 0) > 2
        indok = find(ece_ids.t_e_central.data > 0);
        z0dinput.exp0d.te0   = interp10d(ece_ids.time(indok),ece_ids.t_e_central.data(indok),temps + t_igni,'nearest');
    end
end


% loop voltage
[fv,tfv] = tsbase(shot,'GMAG_FV');
if all(fv(:)  == 0) || (shot > 52800)
    % direct computation
    [flux,t_flux] = tsbase(shot,'gmag_flxloop');
    [vloop_mag,t_vloop_mag] = tsbase(shot,'gmag_vloop');
    for k=1:size(flux,2)
        % correction offset
        indo = find((t_flux(:,k) >= 22) & (t_flux(:,k) <= 25));
        flux(:,k) = flux(:,k) - mean(flux(indo,k));
        flux(:,k) = sgolayfilt(flux(:,k),1,5);
        vloop(:,k) = - z0dxdt(flux(:,k),t_flux(:,k));
    end
    %figure(21);clf;subplot(2,1,1);plot(t_vloop_mag,vloop_mag,'r',t_flux,vloop,'b');
    %subplot(2,1,2);plot(t_flux,flux);drawnow
    ind_flux_ok = [1,3,4,6];
    flux_ok    =  mean(flux(:,ind_flux_ok),2);
    vloop_ok   =  mean(vloop(:,ind_flux_ok),2);
%     figure(21);clf;subplot(3,1,1);plot(t_vloop_mag(:,1),vloop_mag(:,ind_flux_ok),'r',t_flux(:,1),vloop(:,ind_flux_ok),'b',t_flux(:,4),vloop_ok,'.k',t_flux(:,4),vloop_ok,'k');
%     flux_from_vloop = -cumtrapz(t_flux(:,4),vloop_ok);
%     flux_from_vloop = flux_from_vloop- mean(flux_from_vloop) + mean(flux_ok);
%     subplot(3,1,2);plot(t_flux(:,1),flux(:,ind_flux_ok),'r',t_flux(:,4),flux_ok,'k',t_flux(:,4),flux_ok,'.k',t_flux(:,4),9.8  - (-7.82),'g');
%     subplot(3,1,3);
%     plot(t_flux(:,4),flux_ok,t_flux(:,4),flux_from_vloop);
%     keyboard
%     drawnow
tfv        = t_flux(:,ind_flux_ok(1));
if any(diff(tfv) <= 0)
    indbad = min(length(tfv),find(diff(tfv) <= 0) + 1);
    while ~isempty(indbad)
        tfv(indbad) = [];
        vloop_ok(indbad) = [];
        flux_ok(indbad) = [];
        indbad = min(length(tfv),find(diff(tfv) <= 0) + 1);
    end
end

    z0dinput.exp0d.vloop = interp10d(tfv,sgolayfilt(vloop_ok,1,5),temps,'nearest');
    z0dinput.exp0d.edgeflux = interp10d(tfv,flux_ok,temps,'nearest');
    z0dinput.cons.flux = z0dinput.exp0d.edgeflux / 2 /pi; % ?
    
    if false
        figure(22);
        subplot(2,1,1)
        plot(temps,z0dinput.exp0d.vloop,'r',temps,vloop_eqx,'b');
        ylabel('Vloop_{LCFS} (V)');
        legend('in METIS (experimental from DMAG)','EQUINOX');
        title(sprintf('shot WEST #%d',shot));
        subplot(2,1,2)
        plot(temps,z0dinput.exp0d.edgeflux,'r',temps,-flux_edge_eqx,'b');
        ylabel('Psi_{LCFS} (Wb)');
        xlabel('time (s)');
        edition2
        drawnow
    end
else
    z0dinput.exp0d.vloop = interp10d(tfv(:,2),sgolayfilt(fv(:,2),1,5),temps,'nearest');
    z0dinput.exp0d.edgeflux = interp10d(tfv(:,1),sgolayfilt(fv(:,1),1,5),temps,'nearest');
    z0dinput.cons.flux = z0dinput.exp0d.edgeflux / 2 /pi; % ?
end

z0dinput.exp0d.pohm  = pohm;

% alternative Pohm computation if available using NICE data
try
  [time_pohm,pohm_equi2d,pohm_equi1d,pohm_bilan,vres_equi2d,vres_equi1d,vres_bilan,vloop_bilan,ip_bilan,ip_equi1d,ip_equi2d,vloop_LCFS,flux_LCFS,flux_loop] = compute_pohm(shot);
  if any(isfinite(pohm_equi2d))
    % replace NaN by bilan data
    pohm_equi2d(~isfinite(pohm_equi2d)) = pohm_bilan(~isfinite(pohm_equi2d));
    %vres_equi2d(~isfinite(vres_equi2d)) = vres_bilan(~isfinite(vres_equi2d));
    vloop_LCFS(~isfinite(vloop_LCFS)) = vloop_bilan(~isfinite(vloop_LCFS));
    flux_LCFS(~isfinite(flux_LCFS)) = flux_loop(~isfinite(flux_LCFS));
    flux_LCFS(flux_LCFS == 0) = flux_loop(flux_LCFS == 0);
    % time interpolation
    z0dinput.exp0d.pohm = interp10d(time_pohm,pohm_equi2d,temps,'nearest');
    z0dinput.exp0d.vloop = interp10d(time_pohm,vloop_LCFS,temps,'nearest');
    z0dinput.exp0d.vmes = interp10d(time_pohm,vloop_bilan,temps,'nearest');
    z0dinput.exp0d.edgeflux = interp10d(time_pohm,flux_LCFS,temps,'nearest');
    z0dinput.cons.flux = z0dinput.exp0d.edgeflux / 2 /pi; % ?
    disp('usind 2D NICE grid for computing Pohm and Vloop');
  end
catch
  disp('Not able to use NICE data to compute Pohm and vloop');
end



if with_equinox
    z0dinput.exp0d.qa    = qax_eqx;
end
%	GMAG_BELI   	: VACTH: Diamagnetic parameters. 1: beta dia; 2: li; 3: beta MHD; 4: beta N; 5: lik; 6: li dia
if ~isempty(gqbeli) & ~all(gqbeli(:) == 0)
    z0dinput.exp0d.betap = interp10d(tgqbeli(:,1),gqbeli(:,3),temps,'nearest');
elseif with_equinox
    z0dinput.exp0d.betap = betap_eqx;
end
z0dinput.exp0d.ip    = z0dinput.cons.ip;
if ~isempty(prad)
    z0dinput.exp0d.prad  = prad .* 1e6;
    z0dinput.exp0d.pradsol  = pradsol .* 1e6;
end
z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
z0dinput.exp0d.li    = li;
z0dinput.exp0d.picrh = z0dinput.cons.picrh;
z0dinput.exp0d.plh   = z0dinput.cons.plh;
z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;

z0dinput.exp0d.edgeflux    = z0dinput.cons.flux;
z0dinput.machine     = 'WEST';
z0dinput.shot        = shot;
z0dinput.option.first_wall = '';
z0dinput.option.available_flux =  9.8  - (-7.82);


% [fn,tfn] = tsbase(shot,'GFLUNTN');
% % protection against problem in data acquisition
% if any(diff(tfn(:,1)) <= 0) ||any(diff(tfn(:,2)) <= 0)
%     indbad = find((diff(tfn(:,1)) <= 0) | (diff(tfn(:,2)) <= 0));
%     while ~isempty(indbad)
%         tfn(indbad,:) = [];
%         fn(indbad,:)  = [];
%         indbad        = find((diff(tfn(:,1)) <= 0) | (diff(tfn(:,2)) <= 0));
%     end
% end
% if ~isempty(fn)
%   neutron = interp1(tfn(:,2),sgolayfilt(10 .^ fn(:,2),1,5),temps,'nearest',NaN);
%   if any(neutron == 0)
% 	neutron_alt = interp1(tfn(:,1),sgolayfilt(10 .^ fn(:,1),1,5),temps,'nearest',NaN);
%         neutron(neutron == 0) = neutron_alt(neutron == 0);
%   end
% else
%   neutron = NaN .* temps;
% end
% z0dinput.exp0d.ndd = neutron;

[tfn,neutron,neutron1,neutron2] = process_neutrons(shot);
if ~isempty(neutron)
    neutron = interp1(tfn,neutron,temps,'nearest',NaN);
else
   neutron = NaN .* temps;   
end
z0dinput.exp0d.ndd = neutron;


function yi = interp10d(x,y,xi,methode)

if (size(x,1)>1) & (size(x,2)==1)
    indnok = 0;
    nb     = 100;
    while (~isempty(indnok)) & (nb >0)
        indnok = find(diff(x)<=0);
        if ~isempty(indnok)
            x(indnok) = [];
            y(indnok,:) = [];
            
        end
        nb = nb - 1;
    end
end
yi = interp1(x,y,xi,methode);


function [Rsepa,Zsepa,a,ra,za,k,d] = reshape_LCFS(Rext,Zext,nbp)


Rsepa = [];
Zsepa = [];
a = [];
ra = [];
za = [];
k = [];
d = [];

indbad = 0;
indok  = find(isfinite(Rext) & isfinite(Zext));
if length(indok) < 5
    indbad = 1;
else
    r = Rext(indok);
    r = r(:);
    z = Zext(indok);
    z = z(:);		    
%      KH = sort(unique(convhull(r,z)));
%      if (length(KH) ~= length(r))
%  	    r = r(KH);
%  	    z = z(KH);
%  	    r = r(:);
%  	    z = z(:);
%      end
    if length(r) < 5
	indbad = 1;
    else
	r0   = (min(r) + max(r)) ./ 2;
	z0   = (min(z) + max(z)) ./ 2;
	cc   = (r - r0) + sqrt(-1) .* (z - z0);
	thc  = unwrap(angle(cc));
	thc(thc <0) = thc(thc<0) + 2 .* pi;
	rhoc = abs(cc);
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(1,rhoc,rhoc,rhoc);
	thc = cat(1,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	while (length(thc) > 3) && ~isempty(indnok)
	  thc(indnok) =[];
	  rhoc(indnok)  = [];
	  indnok = find(diff(thc)<=0);
	end  
	teta  = linspace(0,2*pi,nbp)';
	rho = spline(thc,rhoc,teta);
	R  = r0 + rho .* cos(teta);
	Z  = z0 + rho .* sin(teta);
	R(end) = (R(1)+ R(end)) ./ 2;
	R(1)   = R(end);
	Z(end) = (Z(1)+ Z(end)) ./ 2;
	Z(1)   = Z(end);
	
	if any(~isfinite(R)) || any(~isfinite(Z))		    
	  indbad = 1;
	end
	%figure(21);clf;plot(Rext,Zext,'b',R,Z,'.r');drawnow
	%pause(0.1);
    end
end

if indbad == 1
  return
end

Rext = R';
Zext = Z';		
if any(~isfinite(Rext(:))) || any(~isfinite(Zext(:)))
  keyboard
end 


% recalcul des parametres pour verification
ve    = ones(1,size(Rext,2));
rmin  = min(Rext,[],2);
rmax  = max(Rext,[],2);
ra    = max(0.1,0.5 .* (rmin + rmax));
a     = max(0.01,0.5 .* (rmax - rmin));
zmin  = min(Zext,[],2);
zmax  = max(Zext,[],2);
za    = (zmin + zmax) ./ 2;
b     = 0.5 .* (zmax - zmin);
k     = max(0.5,b ./ a);
mask1 = (Zext == (max(Zext,[],2)*ve));
mask2 = (Zext == (min(Zext,[],2)*ve));

rzmax = max(Rext .* mask1,[],2);
rzmin = max(Rext .* mask2,[],2);
cl    = ra - rzmin;
cu    = ra - rzmax;
d     = (cl+cu) ./2 ./ a;
Rsepa = Rext;
%Zsepa = Zext - za * ones(1,size(Zext,2));
Zsepa = Zext;



