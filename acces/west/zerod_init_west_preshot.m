%============================================================
% acces au donnees de WEST avant choc (preparation)
%============================================================
function z0dinput = zerod_init_west_preshot(mode_exp,shot,gaz,temps,z0dinput)

% constants
phys = cphys;
phys.pam3        =   (4.41e-4 .* phys.avo);

langue                 =  'anglais';
% cas donnees TS preparation
z0dinput.exp0d         = [];

% 1 - numero du choc
if isempty(gaz) & isempty(shot)

    numshot = 0;
    refshot = 0;
    try 
        refshot = evalin('base','jeux1.post.z0dinput.shot');
    catch
        
       try
            refshot = tsdernier_choc;
            rignitron = tsbase(refshot,'rignitron');
            if isempty(rignitron)
                refshot = refshot - 1;
            end
        catch
            refshot = 0;
        end
    end 

    prompt={'shot number: (0= shot in preparation)','reference shot number: (0= no reference shot)'};
    def={sprintf('%d',numshot),sprintf('%d',refshot)};
    dlgTitle='Access to WEST data';

    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
            z0dinput = [];
	    return
    end
    shot     = str2num(answer{1});
    refshot  = str2num(answer{2});
    if (shot < 50000) && (shot > 0);
	warndlg('Shot number < 50000 is Tore Supra shot','This is not a WEST shot');
        z0dinput = [];
	return	
    end
    if (refshot < 50000) && (refshot > 0);
	warndlg('Shot number < 50000 is Tore Supra shot','This is not a WEST shot');
        z0dinput = [];
	return	
    end
else
    refshot = gaz;
end

%function checkup
% Check-up du pilote avant choc

% *** lecture des consignes principales
fprintf('Reading PCS data for shot %d',shot);
wf = pulse_settings(shot);
disp('...done');

if isempty(wf) || isempty(wf.WEST_PCS.Plasma.Ip.waveform.value) || all(wf.WEST_PCS.Plasma.Ip.waveform.value == 0)
	disp('No plasma') 
        z0dinput = [];
	return
else
  fprintf('Preparing METIS simulation with pre-shot PCD data scheduled for shot number %d',shot);
end

% common time 
t = time_waveforms(wf);
t = t(:);
%t = cat(1,t,max(t) + 1);
t = sort(t);
while (any(diff(t)<=0))
	t(find(diff(t)==0)+1)=t(find(diff(t)==0)+1) + 1e-6;
end



% end shot when ip goes back to 0
ind_last= max(find(wf.WEST_PCS.Plasma.Ip.waveform.value > 0)) + 1;

if ind_last <= length(wf.WEST_PCS.Plasma.Ip.waveform.value)
  tmax = wf.WEST_PCS.Plasma.Ip.waveform.time(ind_last);
  t = t(t<=tmax);
end
ind_first = find(wf.WEST_PCS.Plasma.Ip.waveform.value > 0,1) - 1;
if ind_first > 0
   tmin = wf.WEST_PCS.Plasma.Ip.waveform.time(ind_first);
   t = t(t>=tmin);   
end
% resampling of t
wt = data_waveforms(wf,t);
wt.time = t;
wt.time = t - 40;


% first time 
ind_first = max(1,find(wt.Plasma_Ip > 0,1) - 1);
ind_end   = min(length(wt.time),max(find(wt.Plasma_Ip > 0)) + 1);
ind_ip    = ind_first:ind_end;
noms = fieldnames(wt);
for k = 1:length(noms)
    if ~isempty(wt.(noms{k}))
	wt.(noms{k}) = wt.(noms{k})(ind_ip);
    end
end
wt.time = wt.time - wt.time(1);
if nargin < 4
    temps = wt.time;
elseif isempty(temps)
    temps = wt.time;   
else
    noms = fieldnames(wt);
    t = wt.time
    for k = 1:length(noms)
        if ~isempty(wt.(noms{k})) && (length(wt.(noms{k})) == length(t))
            wt.(noms{k}) = interp1(wt.time,wt.(noms{k}),temps,'linear',0);
        end
    end
    
end
    

% if ~isempty(temps)
%     % we keep the set of tile slices given in input
% else
%     % creation de la base temps
%     ind_ipmax = find(wt.ip == max(wt.ip),1);
%     if (ind_ipmax + 1) < length(wt.ip)
%         temps = cat(2,wt.time(1):0.01:wt.time(ind_ipmax),wt.time(ind_ipmax + 1):0.1:wt.time(end))';
%     else
%         temps = (wt.time(1):0.01:wt.time(ind_ipmax))';
%     end
%     temps = union(temps,wt.time);
%     temps(diff(temps)< 1e-4) = [];
%     if isempty(temps)
%         dt    = max(mean(diff(wt.time)),0.1);
%         temps = (min(wt.time):dt:max(wt.time))';
%         if length(temps) > 1001
%             temps = linspace(min(wt.time),max(wt.time),1001)';
%         end
%     end
% end
% noms = fieldnames(wt);
% t    = wt.time;
% for k = 1:length(noms)
%     if ~isempty(wt.(noms{k}))
% 	wt.(noms{k}) = interp1(t,wt.(noms{k}),temps,'linear',0);
%     end
% end

% backup value
wt.Plasma_Rgeo(wt.Plasma_Rgeo <= 0) = min(wt.Plasma_Rgeo(wt.Plasma_Rgeo > 0));

% computation of minor radius without X points
% real time dependant section is not yet available
if shot == 0
    [fw.R,fw.Z] = west_limiter(refshot);
else
    [fw.R,fw.Z] = west_limiter(shot);
end
wt.Plasma_Rgeo(wt.Plasma_Rgeo <= min(fw.R + 1e-2)) = min(fw.R + 1e-2);
vt = ones(size(wt.time));
vfw = ones(1,length(fw.R));
dd  = sqrt((wt.Plasma_Rgeo * vfw - vt * fw.R(:)') .^ 2 + (wt.Plasma_Zgeo * vfw - vt * fw.Z(:)') .^ 2);
wt.amin = min(dd,[],2);
wt.Rext = wt.Plasma_Rgeo + wt.amin;
wt.amin = wt.amin - wt.Plasma_EROG;

% default value
wt.dl = zeros(size(wt.time));  
wt.Kl = ones(size(wt.time));
wt.dh = zeros(size(wt.time));  
wt.Kh = ones(size(wt.time));

% xpoint position
tau_x = 0.2;
if ~isempty(wt.Plasma_dXLow)
  wt.Plasma_dXLow(wt.Plasma_dXLow <= 0) = NaN;
  Rl = 2.218;Zl = -0.705;
  Rh = 2.341;Zh = -0.522;
  VR = Rh - Rl;
  VZ = Zh - Zl;
  NV = sqrt(VR .^ 2 + VZ .^ 2);
  VR = VR / NV;
  VZ = VZ /NV;
  wt.Rxl = Rl + VR .* wt.Plasma_dXLow;
  wt.Zxl = Zl + VZ .* wt.Plasma_dXLow;
  wt.dl(isfinite(wt.Plasma_dXLow)) = (wt.Plasma_Rgeo(isfinite(wt.Plasma_dXLow)) - wt.Rxl(isfinite(wt.Plasma_dXLow))) ./ wt.amin(isfinite(wt.Plasma_dXLow));
  [tntot,wt.dl] = z0ode(wt.time,wt.dl ./ tau_x,tau_x .* ones(size(wt.time)),wt.dl(1));
  wt.Kl(isfinite(wt.Plasma_dXLow)) = (wt.Plasma_Zgeo(isfinite(wt.Plasma_dXLow)) - wt.Zxl(isfinite(wt.Plasma_dXLow))) ./ wt.amin(isfinite(wt.Plasma_dXLow));
  [tntot,wt.Kl] = z0ode(wt.time,wt.Kl ./ tau_x,tau_x .* ones(size(wt.time)),wt.Kl(1));
end
if ~isempty(wt.Plasma_dXUp)
  wt.Plasma_dXUp(wt.Plasma_dXUp <= 0) = NaN;
  Rl = 2.218;Zl = 0.705;
  Rh = 2.341;Zh = 0.522;
  VR = Rh - Rl;
  VZ = Zh - Zl;
  NV = sqrt(VR .^ 2 + VZ .^ 2);
  VR = VR / NV;
  VZ = VZ /NV;
  wt.Rxh = Rl + VR .* wt.Plasma_dXUp;
  wt.Zxh = Zl + VZ .* wt.Plasma_dXUp;
  wt.dh(isfinite(wt.Plasma_dXUp)) = (wt.Plasma_Rgeo(isfinite(wt.Plasma_dXUp)) - wt.Rxh(isfinite(wt.Plasma_dXUp))) ./ wt.amin(isfinite(wt.Plasma_dXUp));
  [tntot,wt.dh] = z0ode(wt.time,wt.dh ./ tau_x,tau_x.* ones(size(wt.time)),wt.dh(1));
  wt.Kh(isfinite(wt.Plasma_dXUp)) = (wt.Plasma_Zgeo(isfinite(wt.Plasma_dXUp)) - wt.Zxh(isfinite(wt.Plasma_dXUp))) ./ wt.amin(isfinite(wt.Plasma_dXUp));
  [tntot,wt.Kh] = z0ode(wt.time,wt.Kh ./ tau_x,tau_x.* ones(size(wt.time)),wt.Kh(1));
end

% gas puff
wt.gas_puff = zeros(size(wt.time));
fact_gas    = ones(size(wt.time));
for k=1:21
    nom = sprintf('Actuators_Gas_valve%2.2d',k);
    if ~isempty(wt.(nom)) && ~all(wt.(nom) == 0)
        wt.gas_puff = wt.gas_puff + fact_gas(k) .* wt.(nom);
    end
end

% nbar
fact_nbar    = [1 1 1];
wt.nbar      = zeros(size(wt.time));
for k=1:3
    nom = sprintf('Actuators_Gas_REF%d',k);
    if ~isempty(wt.(nom)) && ~all(wt.(nom) == 0)
        wt.nbar = max(wt.nbar , fact_nbar(k) .* wt.(nom));
    end
end

% calcul de la puissance LH
if ~isempty(wt.Actuators_Heating_LHCD_phase_c1)
    npar1 = 1.83 + wt.Actuators_Heating_LHCD_phase_c1 ./ 200;
else
    npar1 = zeros(size(wt.time));
end
if ~isempty(wt.Actuators_Heating_LHCD_phase_c2)
    npar2 = 2.03 + (wt.Actuators_Heating_LHCD_phase_c2 + 90) ./ 310;
else
    npar2 = zeros(size(wt.time));
end
wt.plhtot = zeros(size(wt.time));
wt.npari  = zeros(size(wt.time));
if ~isempty(wt.Actuators_Heating_LHCD_power_c1)
  wt.plhtot = wt.plhtot + wt.Actuators_Heating_LHCD_power_c1;
  npar1 = npar1 .*  wt.Actuators_Heating_LHCD_power_c1;
end
if ~isempty(wt.Actuators_Heating_LHCD_power_c2)
  wt.plhtot = wt.plhtot + wt.Actuators_Heating_LHCD_power_c2;
  npar2 = npar2 .*  wt.Actuators_Heating_LHCD_power_c2;
end
npari  = trapz(wt.time,npar1 + npar2) ./ max(1e-6,trapz(wt.time,wt.plhtot));

% calcul de la puissance FCI
wt.picrhtot = zeros(size(wt.time));
if ~isempty(wt.Actuators_Heating_ICRH_power_c1)
	wt.picrhtot = wt.picrhtot + wt.Actuators_Heating_ICRH_power_c1;
end
if ~isempty(wt.Actuators_Heating_ICRH_power_c2)
	wt.picrhtot = wt.picrhtot + wt.Actuators_Heating_ICRH_power_c2;
end
if ~isempty(wt.Actuators_Heating_ICRH_power_c3)
	wt.picrhtot = wt.picrhtot + wt.Actuators_Heating_ICRH_power_c3;
end

% others data
if shot == 0
    wt.shot = tsdernier_shot + 1;
else
    wt.shot = shot;
end
% gaz
% to be done
warning('gas type not read');
gaz = 2;

% parametre
z0dinput.option = map_option_west_preshot(z0dinput.option);
%
z0dinput.cons.temps    = wt.time;
z0dinput.cons.ip       = max(1e3,wt.Plasma_Ip);
z0dinput.cons.flux     = zeros(size(temps));
z0dinput.geo.a         = wt.amin;
z0dinput.geo.R         = wt.Plasma_Rgeo;
z0dinput.geo.z0        = wt.Plasma_Zgeo;
z0dinput.geo.K         = (wt.Kl + wt.Kh) ./ 2;
z0dinput.geo.d         = (wt.dl + wt.dh) ./ 2;
if any(wt.gas_puff ~= 0)
    z0dinput.cons.nbar     = wt.nbar .* 1e19 +  2 .* sqrt(-1) .* wt.gas_puff .* phys.pam3;
else
    z0dinput.cons.nbar     = wt.nbar .* 1e19;
end
z0dinput.cons.picrh    = wt.picrhtot;
z0dinput.cons.plh      = wt.plhtot;
z0dinput.cons.pnbi     = zeros(size(temps)); 
z0dinput.cons.pecrh    = zeros(size(temps)); 
z0dinput.cons.hmore    = ones(size(temps));
z0dinput.cons.zeff    = (1.2 + (gaz == 4)) .* ones(size(temps));
z0dinput.option.li     = 0.5;
z0dinput.option.mino   = 'H';
	
% donnee calculee dans le zerod
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];
	
	
% les parametres
ffci = tsmat(shot,'DFCI;PILOTAGE;FREQUENC');
if ~isempty(ffci)
	z0dinput.option.freq = max(1,mean(ffci));
else
	z0dinput.option.freq = 55.5;
end

z0dinput.option.npar0 = 2;
if ~isempty(npari)
	z0dinput.option.npar0 = npari;	
end

z0dinput.machine     = 'WEST';
z0dinput.shot        = wt.shot;
z0dinput.option.first_wall = '';
z0dinput.option.available_flux = 9.8  - (-7.82);

% LCFS generation
exp0d = z0dinput.exp0d;
for k=1:length(wt.time)
    
    sepa_option.rxup      = wt.dh(k);     % upper triangularity (minor radius unit)
    sepa_option.zxup      = wt.Kh(k);    % upper altitude X point (minor radius unit)
    if ~isempty(wt.Plasma_dXUp) && any(isfinite(wt.Plasma_dXUp(k)))
        sepa_option.apup      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
        sepa_option.amup      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
    else
        sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
        sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
    end
    sepa_option.ra        = wt.Plasma_Rgeo(k);       % major radius R0 (m) [6.2]
    sepa_option.za        = wt.Plasma_Zgeo(k);       % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = wt.amin(k);         % minor radius (m) [2]
    sepa_option.rxdo      = wt.dl(k);     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = wt.Kl(k);       % lower altitude X point (minor radius unit)
    if ~isempty(wt.Plasma_dXLow) && any(isfinite(wt.Plasma_dXLow))
        sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
        sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
    else
        sepa_option.apdo      = 0;   % lower separatrix angle (R,X)  (LFS, degrees)
        sepa_option.amdo      = 0;   % lower separatrix angle (-R,X)  (HFS, degrees)
    end
    sepa_option.delta     = 0.73;      % magnetic field at R0
    sepa_option.b0        = 11;
    sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
    sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
    
    z0d1t = zerod_get1t(z0dinput,k);
    rep = z0separatrix(z0d1t,sepa_option,0);
    z0dinput.geo.K(k) = rep.geo.K;
    z0dinput.geo.d(k) = rep.geo.d;
    z0dinput.geo.R(k) = rep.geo.R;
    z0dinput.geo.a(k) = rep.geo.a;
    z0dinput.geo.z0(k) = rep.geo.z0;
    exp0d.Rsepa(k,:) = rep.exp0d.Rsepa(1,:);
    exp0d.Zsepa(k,:) = rep.exp0d.Zsepa(1,:);
end
z0dinput.exp0d = exp0d;

% toroidal magnetic field
if shot == 0
    sdn = datestr(now - 7e-4,'dd-mm-yyyy HH:MM:SS');
    [sdn,hn1] = strtok(sdn,' ');
    hn1 = deblank(hn1);
    sdn = datestr(now,'dd-mm-yyyy HH:MM:SS');
    [sdn,hn2] = strtok(sdn,' ');
    hn2 = deblank(hn2);
    hn1(hn1<=' ') = [];
    hn2(hn2<=' ') = [];
    sdn(sdn<' ')  = [];
else
    [sdnx,hnx] = tsdate(shot);
    sdn = char(abs(sdnx));
    hn = char(abs(hnx));
    sdn(sdn <= ' ') = [];
    sdn(sdn >= '@') = [];
    sdn = strrep(sdn,'/','-');
    hn(hn <= ' ') = [];
    hn(hn >= '@') = [];
    [hh,r] = strtok(hn,':');
    [mm,ss] =  strtok(r(2:end),':');
    mm = str2num(mm);
    if mm < 59
      mm = mm + 1;
      hn1 = hn;
      hn2 = sprintf('%s:%d:%s',hh,mm,ss(2:end));
    else
      mm = mm - 1;
      hn2 = hn;
      hn1 = sprintf('%s:%d:%s',hh,mm,ss(2:end));    
    end
end
[itor,titor] = tsbase('gitoro',sdn,hn1,hn2);
if isempty(itor)
  warning('ITOR can''t be read: B0 must be set manually');
  if refshot > 0
      [itor,titor]           = tsbase(refshot,'gmag_itor');
  else
      [itor,titor]           = tsbase(tsdernier_choc,'gmag_itor'); 
      if isempty(itor)
	 itor = 1247;
      end
  end
end
itor = mean(max(itor,[],2));
rb0                    = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;

% securite densite
%nsat = min(1e20 .* z0dinput.cons.ip ./ 1e6 ./ pi ./ max(z0dinput.geo.a) .^2,0.06e20 .* (z0dinput.cons.ip ./ 1e6) .* z0dinput.geo.R .* sqrt(gaz) ./ z0dinput.geo.a .^ (5/2));
%z0dinput.cons.nbar(z0dinput.cons.nbar < 3e17) = (2/3) * nsat(z0dinput.cons.nbar < 3e17); 
if all(real(z0dinput.cons.nbar) < 3e17)
  [vp,sp,sext,peri] = zgeo0(z0dinput.geo);
  z0dinput.cons.nbar = 5e17 +sqrt(-1) .* imag(z0dinput.cons.nbar);
  if all(imag(z0dinput.cons.nbar) == 0)
	ind_gp = max(2,find(vp == max(vp),1));
	gp = z0dinput.option.p_prefill ./ 1.3806503e-23 ./ z0dinput.option.temp_vac .* z0dinput.option.VV_volume ./  ...
	     (z0dinput.cons.temps(ind_gp) - z0dinput.cons.temps(1));
	gp = gp * ones(size(z0dinput.cons.temps));
	gp(ind_gp:end) = 0;
	z0dinput.cons.nbar = real(z0dinput.cons.nbar) + sqrt(-1) .* gp;
  end
  z0dinput.option.natural = 2;
  warning('No reference for density: natural density enforced and prefill pressure used to compute density waveform reference');
end
z0dinput.cons.nbar(real(z0dinput.cons.nbar) < 3e17) = 3e17 + sqrt(-1) .* imag(z0dinput.cons.nbar(real(z0dinput.cons.nbar) < 3e17)); 
	
% choc de test 
if refshot > 0 
	void = zerod_init_west(12,refshot,((z0dinput.option.gaz == 4) + 1),z0dinput.cons.temps);
	if isempty(void)
	  disp('No data in reference shot');
	  return
	end
	if isfield(void.exp0d,'Rsepa')
		void.exp0d = rmfield(void.exp0d,'Rsepa');
		void.exp0d = rmfield(void.exp0d,'Zsepa');
	end
	z0dinput.exp0d = void.exp0d;
	z0dinput.cons.zeff = max(1.1,interp10d(void.cons.temps,void.cons.zeff,z0dinput.cons.temps,'nearest'));
	ipbad = find(void.cons.ip < 3e5);
	if ~isempty(ipbad)
		[n,z]  = hist(void.cons.zeff(find(void.cons.ip > 3e5)));
		zeff0  = z(find(n==max(n),1));
		z0dinput.cons.zeff(ipbad) = zeff0;
	end
	z0dinput.option.cmin = void.option.cmin ;
	if isfield(z0dinput.exp0d,'XDURt') & isfield(z0dinput.exp0d,'XDURx') & isfield(z0dinput.exp0d,'XDURv')
		indok            = find(all(z0dinput.exp0d.XDURv>=0,2));
		if length(indok) >3
			z0dinput.option.dlh       = 0;
		end
	end
	
	% etalonage consigne sur reference cons = etalon(exp,cons_exp,cons)
	if shot ~= refshot
		etal = zerod_init(13,refshot,[],z0dinput.cons.temps);
	else 
		etal = z0dinput;
	end
	fprintf('correction nbar = ')
	% ajout de la densite naturelle
	z0dinput.cons.nbar  =  etalon(real(void.cons.nbar),real(etal.cons.nbar),real(z0dinput.cons.nbar)) + sqrt(-1) .* imag(z0dinput.cons.nbar);
	z0dinput.cons.nbar  = max(min(real(void.cons.nbar)),real(z0dinput.cons.nbar))+ sqrt(-1) .* imag(z0dinput.cons.nbar);
	fprintf('correction gas puff = ')
    z0dinput.cons.nbar  =  real(z0dinput.cons.nbar) + sqrt(-1) .* etalon(imag(void.cons.nbar),imag(etal.cons.nbar),imag(z0dinput.cons.nbar));

	fprintf('correction plh = ')
	z0dinput.cons.plh   =  etalon(void.cons.plh,etal.cons.plh,z0dinput.cons.plh);
	fprintf('correction pecrh = ')
	z0dinput.cons.pecrh =  etalon(void.cons.pecrh,etal.cons.pecrh,z0dinput.cons.pecrh);
	fprintf('correction picrh = ')
	z0dinput.cons.picrh =  etalon(void.cons.picrh,etal.cons.picrh,z0dinput.cons.picrh); 
	fprintf('correction ip = ')
	z0dinput.cons.ip    =  etalon(void.cons.ip,etal.cons.ip,z0dinput.cons.ip);
	
	% on recopie la consigne de idn
	% provisoirement on recopie le choc de reference
	if any(void.cons.pnbi > 10e3) 
	    %z0dinput.cons.pnbi = interp10d(void.cons.temps,void.cons.pnbi,z0dinput.cons.temps,'nearest');
	    z0dinput.option.rtang = 1.47;
	    z0dinput.option.angle_nbi = 90;
	    z0dinput.option.einj = void.option.einj; 
	    z0dinput.cons.ftnbi = ones(size(z0dinput.cons.temps)) .* mean(void.cons.ftnbi);    
    end
    
    
    % geometry 
    if all(abs(z0dinput.geo.K -1) < 0.05) &&  all(abs(z0dinput.geo.d) < 0.05)
       % geometry is no available from PCS
       fform = (wt.IXh + wt.IXh) .* z0dinput.cons.ip ./ 1e6;
       fform = fform ./ max(fform);
       indok = find((void.geo.K > 0.8) & (void.geo.K < 2) & (void.geo.d < 1) & (void.geo.d > - 0.5));
       if length(indok) > 2
                Kref  = trapz(temps(indok), void.cons.ip(indok) .* void.geo.K(indok)) ./ max(eps,trapz(temps(indok), void.cons.ip(indok)));
                dref  = trapz(temps(indok), void.cons.ip(indok) .* void.geo.d(indok)) ./ max(eps,trapz(temps(indok), void.cons.ip(indok)));
                Rref  = trapz(temps(indok), void.cons.ip(indok) .* void.geo.R(indok)) ./ max(eps,trapz(temps(indok), void.cons.ip(indok)));
                aref  = trapz(temps(indok), void.cons.ip(indok) .* void.geo.a(indok)) ./ max(eps,trapz(temps(indok), void.cons.ip(indok)));
                z0ref  = trapz(temps(indok), void.cons.ip(indok) .* void.geo.z0(indok)) ./ max(eps,trapz(temps(indok), void.cons.ip(indok)));
                z0dinput.geo.K = 1 + (Kref - 1) .* fform;
                z0dinput.geo.d = dref .* fform;
                z0dinput.geo.a = aref .* fform + z0dinput.geo.a .* (1 - fform);
                rb0 =  z0dinput.geo.R .* z0dinput.geo.b0;
                z0dinput.geo.R = Rref .* fform + z0dinput.geo.R .* (1 - fform);
                z0dinput.geo.z0 = z0ref .* fform + z0dinput.geo.z0 .* (1 - fform);
                z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;
       end
    end
end






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

% etalonnage des consigne sure le choc de reference
function cons = etalon(mesure_exp,cons_exp,cons)

% cas 0
if all(cons == 0)
	fprintf('1 * consigne + 0   (consigne nulle)\n');
	return
end

% 1 recherche des point valides
delta = abs(mesure_exp - cons_exp) ./ max(cons);
indok = find((delta <= 0.5) & (cons_exp > eps) & (mesure_exp > eps));
if length(indok) < 7
	fprintf('1 * consigne + 0 (donnees insuffisantes)\n');
	return
end

% etalonnage lineaire
mm = mean(cons_exp(indok));
ss = std(cons_exp(indok));
if ss ==0
	ss = 1;
end
xx = (cons_exp(indok) - mm) ./ ss;
warning off
pp = polyfit(xx,mesure_exp(indok),1);
warning on
if any(~isfinite(pp))
	fprintf('1 * consigne + 0 (donnees insuffisantes)\n');
	return
end
if pp(1) < 0
	fprintf('1 * consigne + 0 (pente negative)\n');
	return
end
% nouvelle consigne corrigee
xx = (cons - mm - min(cons)) ./ ss;
cons = max(0,polyval(pp,xx) .* (cons > 0)) + min(cons);

fprintf('%g * consigne + %g\n',pp(1) ./ ss , pp(2) - pp(1) .* mm ./ ss)



function time = time_waveforms(wf,time)

if nargin < 2
    time = [];
end

if isstruct(wf)
    if isfield(wf,'time') && isfield(wf,'value')
        time = union(time,wf.time);
    else
        noms = fieldnames(wf);
        for k=1:length(noms)
            time = time_waveforms(wf.(noms{k}),time);
        end
    end
else
    disp('What is this data ?');
    keyboard
end


function wt = data_waveforms(wf,t,root)

if nargin < 3
    root  = '';
    rootb = '';
else
    root = strrep(root,'waveform','');
    root = strrep(root,'env_waveform','env');
    root = strrep(root,'WEST_PCS','');
    root = strrep(root,'val','');
    root = strrep(root,'dp','');
    if isempty(root)
        %nothing
    elseif root(1) == '_'
        root = root(2:end);
    end
    rootb = root;
    root = sprintf('%s_',root);
end
wt = [];

if isstruct(wf)
    if isfield(wf,'time') && isfield(wf,'value')
        if ~isempty(wf.value)
            time = wf.time;
            value = wf.value;
            [time,idt]   = sort(time);
            value        = value(idt);
            [time,ia,ic] = unique(time);
            value        = value(ia);
            wt = interp1(time,value,t,'linear',0);
        else
            wt = NaN * ones(size(t));
        end
    else
        noms = fieldnames(wf);
        for k=1:length(noms)
            out = data_waveforms(wf.(noms{k}),t,sprintf('%s%s',root,noms{k}));
            if isstruct(out)
                noms_out = fieldnames(out);
                if length(noms_out)  == 1
                    if isfield(wt,noms{k})
                        wt.(rootb) = out.(noms_out{1});                        
                    else
                        switch noms{k}
                            case {'waveform','dp','val'}
                                wt.(rootb) = out.(noms_out{1});
                            otherwise
                                wt.(sprintf('%s%s',root,noms{k})) = out.(noms_out{1});
                        end
                    end
                else
                    for l=1:length(noms_out)
                        wt.(noms_out{l}) = out.(noms_out{l});                            
                    end
                end
            else
                if isfield(wt,noms{k})
                    wt.(rootb) = out;
                else
                    switch noms{k}
                        case {'waveform','dp','val'}
                            wt.(rootb) = out;
                        otherwise
                            wt.(sprintf('%s%s',root,noms{k})) = out;
                    end
                end
            end
        end
    end
else
    disp('What is this data ?');
    keyboard
end





