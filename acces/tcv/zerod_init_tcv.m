%============================================================
% acces au donnees de TCV
%============================================================
function z0dinput = zerod_init_tcv(mode_exp,shot,gaz,temps,z0dinput,tunnel)

if nargin < 6
  tunnel = 0;
end
langue                 =  'anglais';
% cas donnees TCV
z0dinput.exp0d         = [];
% 1 - numero du choc
if isempty(gaz) || isempty(shot) || ((nargin > 5) && isempty(tunnel))

    try 
	numshot = fix(evalin('base','param.from.shot.num'));
    catch
	numshot = 55555;
    end 

    prompt={'shot number :','charge of main gas:','tunnel port:'};
    def={sprintf('%d',numshot),'1',tunnel};
    dlgTitle='Access to TCV data';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
            z0dinput = [];
	    return
    end
    shot  = str2num(answer{1});
    gaz   = str2num(answer{2});
    tunnel =  str2num(answer{3});
 
end


% connection to TCV data base
cr = NaN;
try
  if tunnel > 0
      % bypass error appening on some system where the first try return an error 
      try
	  cr = mdsconnect(sprintf('localhost:%d',tunnel));
      catch
	  cr = mdsconnect(sprintf('localhost:%d',tunnel));      
      end
      if isempty(cr)
	cr =0;			
      end
      
  else
      % bypass error appening on some system where the first try return an error 
      try
	  cr = mdsconnect('tcvdata');
      catch
	  cr = mdsconnect('tcvdata');      
      end
      if isempty(cr)   
	cr =0;			
      end
  end
catch
     try 
	port = detect_tunnel;
	if ~isempty(port)
            % bypass error appening on some system where the first try return an error 
	    try
		cr = mdsconnect(sprintf('localhost:%d',port));
            catch
		cr = mdsconnect(sprintf('localhost:%d',port));            
            end
	    if isempty(cr)
	      cr =0;			
	    end	   
	else
	    disp('you can try to redirect MDS+ access if you have an account on EPFL server:  ssh -L <port_number>:tcvdata:8000 <user_name>@lac911.epfl.ch or ssh -L <port_number>:tcvdata:8000 <user_name>@lac912.epfl.ch (with port number > 5500)');
	end
     catch
      cr = -7;
     end
end
if cr < 0
  disp('ubable to connect to TCV data base');
  disp('you can try to redirect MDS+ access if you have an account on EPFL server:  ssh -L <port_number>:tcvdata:8000 <user_name>@lac911.epfl.ch or ssh -L <port_number>:tcvdata:8000 <user_name>@lac912.epfl.ch (with port number > 5500)');
  z0dinput = [];
  return
else
  rep = mdsopen('tcv_shot',shot);
  if isempty(rep) || (rep ~= shot)
    fprintf('unable to access shot %d\n',shot);
    z0dinput = [];
    return
  else
    fprintf('Reading data for TCV shot %d\n',rep);
  end 
  
end

% lecture des donnees
ip = mdsvalue(eval(sprintf(['''\\','magnetics::iplasma:trapeze',''''])));
if ~isnumeric(ip)
  ip = mdsvalue(eval(sprintf(['''\\','results::I_p',''''])));
  if ~isnumeric(ip)
	  disp('No plasma')
	  z0dinput = [];
	  return
  end

  tip = mdsvalue(eval(sprintf(['''dim_of(\\','results::I_p',')'''])));
else
  tip = mdsvalue(eval(sprintf(['''dim_of(\\','magnetics::iplasma:trapeze',')'''])));
end
tip(~isfinite(ip)) = 0;
ip(~isfinite(ip)) = 0;
sign_ip = sign(mean(sign(ip)));
ip = abs(ip);
if isempty(ip)
	disp('No plasma')
        z0dinput = [];
	return
end

if gaz < 0
  tt = tip(tip >=0);
else
  tt = tip(tip >=0 & ip > 1e4);
end
if isempty(tt)
	disp('No plasma current');
        z0dinput = [];
    	return
end

if ~isempty(temps)
	% on prend le vecteur donnee en entree
elseif gaz < 0
    gaz = abs(gaz);
    temps = tt;
    temps(find(diff(temps)<=0)) =[];
else
    temps = tt;
    temps(find(diff(temps)<=0)) =[];
    if length(temps) > 1001
	temps = linspace(min(tt),max(tt),1001)';
    end	
end
%
% probleme tip
%	
indtip = find(diff(tip) <= 0);
if ~isempty(indtip)
    tip(indtip) = [];
    ip(indtip) = [];
end
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = max(1,interp10d(tip,ip,temps,'nearest'));
ip   = z0dinput.cons.ip;

% flux au bord
% must be checked !
vl = mdsvalue(eval(sprintf(['''\\','magnetics::vloop[*,"001"]',''''])));
if isnumeric(vl)
    vl(~isfinite(vl)) = 0;
    tvl = mdsvalue(eval(sprintf(['''dim_of(\\','magnetics::vloop[*,"001"]',')'''])));
    fluxbord = cumtrapz(tvl,vl);
    %fluxbord = mdsvalue(eval(sprintf(['''\\','magnetics::flux[*,"001"]',''''])));
    %tfluxbord = mdsvalue(eval(sprintf(['''dim_of(\\','magnetics::flux[$1,"001"]',')'''])));
    z0dinput.cons.flux    = - interp10d(tvl,fluxbord,temps,'nearest') ./ 2 ./ pi;
    vloop = interp10d(tvl,vl,temps,'nearest');
else
    disp('No Vloop');
    vloop = NaN * temps;
end
% LCFS
RLCFS = mdsvalue(eval(sprintf(['''\\','r_contour',''''])));
if ~isnumeric(RLCFS)
	disp('No LCFS')
        z0dinput = [];
	return
end
ZLCFS = mdsvalue(eval(sprintf(['''\\','z_contour',''''])));
tLCFS = mdsvalue(eval(sprintf(['''dim_of(\\','results::r_axis',')'''])));

% cast to double
RLCFS = double(RLCFS);
ZLCFS = double(ZLCFS);
tLCFS = double(tLCFS);

if ~isempty(RLCFS)
        indbad = [];
        % put the LCFS on same number of point
        for k = 1:size(RLCFS,2);
	    if any(~isfinite(RLCFS(:,k))) || any(~isfinite(ZLCFS(:,k)))
		Rext = RLCFS(:,k);
		Zext = ZLCFS(:,k);
		indok = find(isfinite(Rext) & isfinite(Zext));
		if length(indok) < 5
		    indbad(end + 1) = k;
		else
		    nbp = length(Rext) + 1;
		    r = Rext(indok);
		    r = r(:);
		    z = Zext(indok);
		    z = z(:);		    
		    KH = sort(unique(convhull(r,z)));
		    if (length(KH) ~= length(r))
		    	    r = r(KH);
	                    z = z(KH);
		            r = r(:);
		            z = z(:);
		    end
		    if length(r) < 5
			indbad(end + 1) = k;
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
			  indbad(end+1) = k;
			else
			  RLCFS(:,k) = R(1:end-1);
			  ZLCFS(:,k) = Z(1:end-1);
			end
			%figure(21);clf;plot(Rext,Zext,'b',R,Z,'.r');drawnow
			%pause(0.1);
		    end
		end
	    end
        end
        if ~isempty(indbad)
	      RLCFS(:,indbad) = [];
	      ZLCFS(:,indbad) = [];
	      tLCFS(:,indbad) = [];
        end
        if length(tLCFS) < 2
              disp('No valid LCFS')
	      z0dinput = [];
	      return
	end
        
	Rext = interp1(tLCFS,RLCFS',temps,'nearest','extrap');
	Zext = interp1(tLCFS,ZLCFS',temps,'nearest','extrap');		
	Rext(:,end)    = Rext(:,1);
	Zext(:,end)    = Zext(:,1);
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
	z0dinput.exp0d.Rsepa = Rext;
	z0dinput.exp0d.Zsepa = Zext - za * ones(1,size(Zext,2));
	z0dinput.geo.a       = a;
	z0dinput.geo.R       = ra;
	z0dinput.geo.z0      = za;
	z0dinput.geo.K       = k;
	z0dinput.geo.d       = d;

end


% Btor
RB0 = mdsvalue(eval(sprintf(['''\\','magnetics::rbphi',''''])));
if ~isnumeric(RB0)
	disp('No B_tor')
        z0dinput = [];
	return
end
sign_b0 = sign(mean(sign(RB0)));
RB0 = abs(RB0);
tRB0 = mdsvalue(eval(sprintf(['''dim_of(\\','magnetics::rbphi',')'''])));
rb0                    = interp10d(tRB0,RB0,temps,'nearest');
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;


% density
nbar = mdsvalue('\results::fir:n_average');
if isnumeric(nbar)
    tnbar = mdsvalue(eval(sprintf(['''dim_of(\\','results::fir:n_average',')'''])));
else
    nbar = [];
end
ne_fit   = mdsvalue(eval(sprintf(['''\\','results::proffit.avg_time:neft',''''])));
tne_fit  = mdsvalue(eval(sprintf(['''dim_of(\\','results::proffit.avg_time:neft',')'''])));
if isnumeric(ne_fit)
    ne_fit    = interp10d(tne_fit,ne_fit,temps,'nearest');
    for k=1:length(temps)
      ind_ok = find(isfinite(ne_fit(k,:)) & (ne_fit(k,:) >= 1e17));
      if ~isempty(ind_ok)
	  nbar_alt(k)     = mean(ne_fit(k,ind_ok),2);
      else
	  nbar_alt(k)     = NaN;
      end
    end
    ne_fit(~isfinite(ne_fit)) = 1e13;
else
    nbar_alt = [];
end
ne_thom = mdsvalue('\results::ne_thomson:foo');
if isnumeric(ne_thom)
    t_thom = mdsvalue('\results::thomson:times');
    ne_thom    = interp10d(t_thom,ne_thom,temps,'nearest');
    for k=1:length(temps)
      ind_ok = find(isfinite(ne_thom(k,:)) & (ne_thom(k,:) >= 1e17));
      if ~isempty(ind_ok)
	  nbar_thom(k)     = mean(ne_thom(k,ind_ok),2);
      else
	  nbar_thom(k)     = NaN;
      end
    end
    ne_thom(~isfinite(ne_thom)) = 1e13;
else 
  nbar_thom =[];
end
if ~isempty(nbar)
    z0dinput.cons.nbar     = interp10d(tnbar,nbar,temps,'nearest');
    %figure;plot(temps,z0dinput.cons.nbar,'r',temps,nbar_alt,'b',temps,nbar_thom,'g');
    indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
    if ~isempty(indbad)
	if ~isempty(nbar_thom)
	    z0dinput.cons.nbar(indbad)     = nbar_thom(indbad);
	end    
    end
%      indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
%      if ~isempty(indbad)
%  	if ~isempty(nbar_alt)
%  	    z0dinput.cons.nbar(indbad)     = nbar_alt(indbad);
%  	end    
%      end
    indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
    if ~isempty(indbad)
	    z0dinput.cons.nbar(indbad)     = max(1e17,0.3e14 .*  ip(indbad) ./ pi ./ z0dinput.geo.a(indbad) .^ 2);
    end
elseif ~isempty(nbar_thom)
    z0dinput.cons.nbar     = nbar_thom;
%      indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
%      if ~isempty(indbad)
%  	if ~isempty(nbar_alt)
%  	    z0dinput.cons.nbar(indbad)     = nbar_alt(indbad);
%  	end    
%      end
    indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
    if ~isempty(indbad)
	    z0dinput.cons.nbar(indbad)     = max(1e17,0.4e14 .*  ip(indbad) ./ pi ./ z0dinput.geo.a(indbad) .^ 2);
    end
elseif ~isempty(nbar_alt)
    z0dinput.cons.nbar     = nbar_alt;
    indbad = find(~isfinite(z0dinput.cons.nbar) | (z0dinput.cons.nbar < 1e17));
    if ~isempty(indbad)
	    z0dinput.cons.nbar(indbad)     = max(1e17,0.4e14 .*  ip(indbad) ./ pi ./ z0dinput.geo.a(indbad) .^ 2);
    end
else
    disp('No density')
    z0dinput = [];
    return
end

% profiles
te_fit   = mdsvalue(eval(sprintf(['''\\','results::proffit.avg_time:teft',''''])));
tte_fit  = mdsvalue(eval(sprintf(['''dim_of(\\','results::proffit.avg_time:teft',')'''])));
ne_fit   = mdsvalue(eval(sprintf(['''\\','results::proffit.avg_time:neft',''''])));
tne_fit  = mdsvalue(eval(sprintf(['''dim_of(\\','results::proffit.avg_time:neft',')'''])));
if isnumeric(ne_fit)
  ne_fit    = interp10d(tne_fit,ne_fit,temps,'nearest');
  ane = max(ne_fit,[],2) ./ mean(ne_fit,2);
  ane = ane -1;
else
  ane = 0.5 * ones(size(temps));
end
if isnumeric(te_fit)
  te_fit    = interp10d(tte_fit,te_fit,temps,'nearest');
end
% confinement
z0dinput.cons.hmore    = ones(size(temps));

% zeff
zeff  = mdsvalue('\results::ibs:z_eff');
if ~isnumeric(zeff)
    zeff  = mdsvalue('\results::zx:foo');
    if isnumeric(zeff)
      tzeff = mdsvalue(eval(sprintf(['''dim_of(\\','results::zx:foo',')'''])));
    end
else
  tzeff = mdsvalue(eval(sprintf(['''dim_of(\\','results::ibs:z_eff',')'''])));
end

if isnumeric(zeff)
  zeff  = interp10d(tzeff,zeff,temps,'nearest');
  for k = 1:size(zeff,1)
      if any(isfinite(zeff(k,:)))
          zeff_ok = zeff(k,isfinite(zeff(k,:)));
          if length(zeff_ok) > 1
	    zeffm(k) = mean(zeff_ok);
	  else
	    zeffm(k) = zeff_ok;	  
	  end
      elseif k > 1
	  zeffm(k) = zeffm(k-1);
      else
	  zeffm(k) = 2.3;
      end
  end
else
  disp('No Zeff measurement - using 2.3')
  zeffm = 2.3 * ones(size(temps));
end
if gaz == 2
	z0dinput.cons.zeff    = max(2,min(16,zeffm));
else
	z0dinput.cons.zeff    = max(1,min(7,zeffm));
end

% power
z0dinput.cons.plh      = zeros(size(temps));
z0dinput.cons.pnbi     = zeros(size(temps)); 
z0dinput.cons.pecrh    = zeros(size(temps)); 
z0dinput.cons.picrh    = zeros(size(temps)); 


% les parametres
z0dinput.option = map_option_tcv(z0dinput.option);
z0dinput.option.signe = sign_ip * sign_b0;


% ECRH
pecrh  = mdsvalue('\results::toray.input:P_GYRO');
if isnumeric(pecrh)
    tpecrh = mdsvalue(eval(sprintf(['''dim_of(\\','results::toray.input:P_GYRO',')'''])));
    pecrh  = interp10d(tpecrh,pecrh,temps,'nearest');
    pecrh_x2 = pecrh(:,1:6);
    pecrh_x3 = pecrh(:,7:9);
    % 10 is sum of contribution
    for k=1:length(temps)
	px2 = pecrh_x2(k,:);
	if any(isfinite(px2))
	    z0dinput.cons.pecrh(k) = sum(px2(isfinite(px2))) * 1e3;
	end
	px3 = pecrh_x3(k,:);
	if any(isfinite(px3))
	    z0dinput.cons.plh(k) = sum(px2(isfinite(px3))) * 1e3;
	end
    end
    
    % try to determine position of deposition
    rho_toray = mdsvalue('\TCV_SHOT::TOP.RESULTS.TORAY.INPUT.RHO_TORAY');
    t_toray = mdsvalue('dim_of(\TCV_SHOT::TOP.RESULTS.TORAY.INPUT.RHO_TORAY,1)');
    toray_pdens_tdi = mdsvalue('\TCV_SHOT::TOP.RESULTS.TORAY.OUTPUT_X:PDENS');
    if isnumeric(toray_pdens_tdi)
	% compute x2 and x3 sum profile
	px3 = zeros(size(rho_toray));
	px2 = zeros(size(rho_toray));
	for k = 1:size(rho_toray,2)
	  plx2 = squeeze(toray_pdens_tdi(:,1:6,k));
	  if any(all(isfinite(plx2),1))
	      plx2 = plx2(:,all(isfinite(plx2),1));
	      px2(:,k) = sum(plx2,2);
	  end
	  plx3 = squeeze(toray_pdens_tdi(:,7:9,k));
	  if any(all(isfinite(plx3),1))
	      plx3 = plx3(:,all(isfinite(plx3),1));
	      px3(:,k) = sum(plx3,2);
	  end
	end
	x = linspace(0,1,size(rho_toray,1))';
	rhox2_center = trapz(x,rho_toray .* px2,1) ./ trapz(x,max(eps,px2),1);
	width_x2     = sqrt(trapz(x,(rho_toray -  ones(size(rho_toray,1),1) * rhox2_center) .^ 2 .* px2,1) ./ trapz(x,max(eps,px2),1));
%      for k = 1:size(rho_toray,2)
%        figure(21);
%        clf
%        width_ecrh = width_x2(k)*sqrt(2);    
%        width_lh   = width_x2(k)*sqrt(2);  
%        %
%        deccd    = width_ecrh / sqrt(2);
%        jeccd    = exp(-(rho_toray(:,k) - rhox2_center(k)).^ 2 ./2 ./ deccd .^ 2);  
%        %
%        jlh    = exp(-(rho_toray(:,k) - rhox2_center(k)).^ 2 ./ width_lh .^ 2);  
%        plot(rho_toray(:,k),px2(:,k)./ max(px2(:,k)),'b',rho_toray(:,k),jlh,'r',rho_toray(:,k),jeccd,'g');
%        drawnow
%        pause(0.1);
%      end
	rhox3_center = trapz(x,rho_toray .* px3,1) ./ trapz(x,max(eps,px3),1);
	width_x3     = sqrt(trapz(x,(rho_toray -  ones(size(rho_toray,1),1) * rhox3_center) .^ 2 .* px3,1) ./ trapz(x,max(eps,px3),1));
	
	
	% interpolation on time
	rhox2_center  = interp1(t_toray,rhox2_center,temps,'nearest','extrap');
	rhox3_center  = interp1(t_toray,rhox3_center,temps,'nearest','extrap');
	width_x2      = interp1(t_toray,width_x2,temps,'nearest','extrap');
	width_x3      = interp1(t_toray,width_x3,temps,'nearest','extrap');

	% power averaged value
	z0dinput.cons.xece         = rhox2_center;
	z0dinput.option.width_ecrh = sqrt(2) .* trapz(temps,z0dinput.cons.pecrh .* width_x2) ./ max(eps,trapz(temps,z0dinput.cons.pecrh));
	z0dinput.option.dlh        = sqrt(2) .* trapz(temps,z0dinput.cons.plh .* width_x3) ./ max(eps,trapz(temps,z0dinput.cons.plh));
	z0dinput.option.xlh        = trapz(temps,z0dinput.cons.plh .* rhox3_center) ./ max(eps,trapz(temps,z0dinput.cons.plh));
    end
end

%NBI

%mdsvalue('\ATLAS::NBH.DATA.MAIN_ADC:DATA');
NBH_in_TCV = false;
if (shot >= 51641)
    lcs_mode   = mdsvalue('data(\VSYSTEM::TCV_NBHDB_I["NBI:LCS_MODE"])');
    % NBH in TCV equiv lcs_mode = 9
    NBH_in_TCV = (lcs_mode == 9);
else
    % Nodes used in previous block only exist outside of Vista for shots after 51641
    if any(shot == [51411 51412 51415 51420 ...                         % 27.JAN.2016
            51458 51459 51460 51461 51463 51465 51470 51472 ...         % 29.JAN.2016
            51628 51629 51631 51632 51633 ...                           % 09.FEB.2016
            51639 51640 ... 51641                                       % 10.FEB.2016
            ]),
        NBH_in_TCV = true;
    end
end
if NBH_in_TCV
    nbh_data      = mdsvalue('\ATLAS::NBH.DATA.MAIN_ADC:DATA');
    if isnumeric(nbh_data)
	neutral_power = max(0,nbh_data(:,37) .* 1e6);
	inj_energy    = nbh_data(:,33) .* 1e3;
	nbh_time      = mdsvalue('dim_of(\ATLAS::NBH.DATA.MAIN_ADC:DATA,0)');
	z0dinput.option.einj = trapz(nbh_time,inj_energy .* neutral_power) ./ max(eps,trapz(nbh_time,neutral_power));
	z0dinput.cons.pnbi   = interp10d(nbh_time,neutral_power,temps,'nearest');
	z0dinput.cons.pnbi(~isfinite(z0dinput.cons.pnbi)) = 0;
    end
end

% NBI diag
% not yet


% other references
z0dinput.cons.ftnbi = zeros(size(temps));
z0dinput.cons.iso  = zeros(size(temps));



% autres donnees
li = mdsvalue('\results::l_i');
if isnumeric(li)
    tli =  mdsvalue(eval(sprintf(['''dim_of(\\','results::l_i',')'''])));
    li    = interp10d(tli,li,temps,'nearest');
else
    li = ones(size(temps));
end

% donnee calculee dans le zerod
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];	
	
% only this two option are available for TCV
if gaz == 1
	z0dinput.option.gaz = 2;
else
	z0dinput.option.gaz = 4;
end
if isfinite(vl(1))
  z0dinput.option.breakdown = - vl(1);
end
if isfinite(li(1))
  z0dinput.option.li = li(1);
elseif isfinite(li(2))
  z0dinput.option.li = li(2);
else
  z0dinput.option.li = 0.5;
end
if ~isempty(isfinite(ane))
    z0dinput.option.vane = mean(ane(isfinite(ane)));
end
% donnees experimentales
z0dinput.exp0d.temps  = temps;
pohm = mdsvalue('\results::oh_power');
if isnumeric(pohm)
  tpohm =  mdsvalue(eval(sprintf(['''dim_of(\\','results::oh_power',')'''])));
  pohm    = interp10d(tpohm,pohm,temps,'nearest');
  z0dinput.exp0d.pohm   = pohm;
else
  z0dinput.exp0d.pohm   = abs(vloop .* ip);    
end
z0dinput.exp0d.vloop    = vloop;
z0dinput.exp0d.pin   = real(z0dinput.cons.pnbi) + imag(z0dinput.cons.pnbi) + z0dinput.cons.plh + z0dinput.cons.pecrh + z0dinput.exp0d.pohm;
z0dinput.exp0d.pw    = z0dinput.exp0d.pin;		      
if ~isempty(zeffm)
	z0dinput.exp0d.zeff  = z0dinput.cons.zeff;
end
z0dinput.exp0d.ane  = ane;
if isnumeric(te_fit)
  vx = ones(length(temps),1) * linspace(0,1,size(ne_fit,2));
  z0dinput.exp0d.nem  = trapz(vx(1,:),ne_fit .* vx,2) ./ trapz(vx(1,:),vx,2);
  z0dinput.exp0d.ne0  = max(ne_fit,[],2);
  z0dinput.exp0d.nea  = min(ne_fit,[],2);
end
wtot = mdsvalue('\results::total_energy');
if isnumeric(wtot)
    twtot =  mdsvalue(eval(sprintf(['''dim_of(\\','results::total_energy',')'''])));
    wtot    = interp10d(twtot,wtot,temps,'nearest');
    z0dinput.exp0d.w     = wtot;
end
%z0dinput.exp0d.taue  = 
if isnumeric(te_fit)
    z0dinput.exp0d.tem   = trapz(vx(1,:),te_fit .* vx,2) ./ trapz(vx(1,:),vx,2);
    z0dinput.exp0d.te0   = max(te_fit,[],2);
    z0dinput.exp0d.tea   = min(te_fit,[],2);
end
%z0dinput.exp0d.wrad  =  
betap = mdsvalue('\results::beta_pol');
if isnumeric(betap)
    tbetap =  mdsvalue(eval(sprintf(['''dim_of(\\','results::beta_pol',')'''])));
    betap    = interp10d(tbetap,betap,temps,'nearest');
    z0dinput.exp0d.betap = betap;
end
qa = mdsvalue('\results::q_psi');
if isnumeric(qa)
  tqa =  mdsvalue(eval(sprintf(['''dim_of(\\','results::q_psi',')'''])));
  qa    = interp10d(tqa,qa,temps,'nearest');
  z0dinput.exp0d.qa    = qa;
end
z0dinput.exp0d.ip    = z0dinput.cons.ip;
prad    = mdsvalue('\results::bolo:prad:total');
if isnumeric(prad)
    tprad   =  mdsvalue(eval(sprintf(['''dim_of(\\','results::bolo:prad:total',')'''])));
    prad    = interp10d(tprad,prad,temps,'nearest');
    z0dinput.exp0d.prad  = prad;
else
    z0dinput.exp0d.prad  = NaN * temps;
end
z0dinput.exp0d.ploss = z0dinput.exp0d.pin - z0dinput.exp0d.prad / 3;

z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
z0dinput.exp0d.li    = li;
z0dinput.exp0d.picrh = z0dinput.cons.picrh;
z0dinput.exp0d.plh   = z0dinput.cons.plh;
z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;
z0dinput.exp0d.edgeflux    = z0dinput.cons.flux;
z0dinput.machine     = 'TCV';
z0dinput.shot        = shot;
z0dinput.option.first_wall = '';     
z0dinput.option.available_flux =   NaN;

mdsclose;
mdsdisconnect;


function option = map_option_tcv(option)

option.gaz = 2;
option.neasser = 1;
option.Recycling = 0.7;
option.natural = 1;
option.ane = 11;
option.fn0a = 1;
option.fn0a_div = 0.1;
%
option.scaling = 12;
option.dilution = 1;
option.tau_limitation = 'On';
option.l2hscaling = 3;
option.pl2h_mass_charge = 1;
option.modeh = 1;
option.hysteresis = 0;
option.configuration = 2;
option.l2hslope =  0.5;
option.usepped_scl = 1;
option.taurotmul = 0;
option.fintrinsic = 0.2;
option.xiioxie = -4.5;
option.kishape = 0;
option.xieorkie = 0;
option.omega_shape = 0;
option.fstiff = 1;
option.ploss_exp = 'max_power';
option.xiioxie_ped = 0;
option.hmore_pped  = 2;
option.fpl2h_lim   = 2;
option.ki_expo     = 2;
option.plhthr      = 'P_LCFS';
option.grad_ped    = 3;
option.ode_pped    = 1;
option.adiabatic   = 1;

%
option.qdds = 1;
option.kidds = 3;
option.sitb = 2;
option.itb_sensitivity = 1;
option.itb_slope_max = 2;
option.smhd = 100;
option.tmhd = 0;
%
option.runaway = 5;
option.modeboot = 2;
%
option.li = 1;
option.breakdown = - 10;
option.berror = 0; % no breakdown simulation by default
option.L_eddy = 0;
option.R_eddy = 0;
option.C_eddy = 0;
option.B_eddy = 0;
option.I_eddy = 0;
option.p_prefill = 0;
option.VV_volume = 0;
%
option.zeff = 0;
option.faccu = 0;
option.heat_acc = 0;
option.fne_acc = 0;
option.zmax = 8;
option.zimp = 6;
option.rimp = 0.3;
option.density_model ='minconv';
%
option.frad = 1;
option.matthews = 1;
option.fte_edge = 1;
option.gaunt = 1;
option.noncoronal = 0;
option.z_prad = 'Stangeby';
%
option.sol_lscale = 0;
option.eioniz     = 0;
option.fnesol     = 0;
option.sol_model  = 'scaling';
option.lcx = 1;
option.fcond = -1;
option.fmom = 0;
option.lambda_scale = 3;
option.sol_rad = 'decoupled';
option.fzmax_div = -1;
option.W_effect = 0;
option.cw_factor = 0;
option.cw_offset = 5e-5;
option.factor_scale = 1;
option.fpower = 0.6000;
option.fR_target = 1;
option.mach_corr = 1;

%
option.angle_ece = 90;
option.synergie  = 1;
option.sens      = 1;
option.eccdmul   = 1;
%
option.angle_nbi = 90;
option.rtang     = 0.65;
option.zext      = 0.2;
option.einj      = 25e3;
option.nbicdmul  = 1;
option.nb_nbi    = 1;
option.e_shielding = 'Honda-Sauter';
option.drs1        = 0;
option.dzs1        = 0;
%
option.angle_nbi2 = 0;
option.rtang2 = 2.85;
option.zext2  = 0.1;
option.einj2  = 85000;
option.nbicdmul2 = 1;
option.drs2   = 0;
option.dzs2   = 0;
%
option.lhmode = 5;  % used as ECCD system for breakdown
option.etalh  = 1;  % perpendicular
option.wlh = 0;
option.xlh = 0.3;
option.dlh = 0.2;
option.angle_ece2 = 90;
option.npar0 = 2;
option.npar_neg = -4;

% used as third injector
option.fwcd = 0;
option.mino = 'H';
option.cmin = 0.15;
option.nphi = 30;
option.freq = 55.5;
option.icrh_width = 0.7;
option.fact_mino  = 0;
option.orbit_width  = 1;
option.icrh_model = 'PION_fit-Stix';
%
option.equi_ppar = 3;
option.signe = 1;
option.cronos_regul= 2;
option.available_flux =  NaN ; %WB
option.machine = 'TCV';
option.evolution = 0;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UAL writing control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


option.init_output_cpo = 0;
option.restart= '';
option.coreprof= 1;
option.coretransp= 1;
option.coresource_lhcd= 1;
option.coresource_eccd= 1;
option.coresource_icrh= 1;
option.coresource_nbicd= 1;
option.coresource_fusion= 1;
option.coreneutrals= 1;
option.coresource_radiation= 1;
option.coresource_cyclotron= 1;
option.neoclassic= 1;
option.coresource_full= 1;
option.equilibrium= 1;
option.grid_equi= 0;
option.scenario_occurrence= '';
option.coreprof_occurrence= '';
option.coretransp_occurrence= '';
option.coreneutrals_occurrence= '';
option.neoclassic_occurrence= '';
option.equilibrium_occurrence= '';
option.coresources_occurrence= '';

% default value for ITM
option.COCOS  = 13;




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


function port = detect_tunnel

if ispc
  port = [];
  return
end

[s,t] = unix('lsof -i | grep epfl');
if s ~= 0
  port = [];
  return
end
[s,t] = unix('netstat -tpln | grep ssh | grep 127.0.0.1');
if s ~= 0
  port = [];
else
  tt = tseparec(t);
  info = t(end,:);
  pos = findstr(info,'127.0.0.1:');
  if isempty(pos)
    port = [];
  else
      [s,r] = strtok(info(pos:end),':');
      data  = strtok(r(2:end),' ');
      port  = str2num(data); 
      fprintf('Port redirection detected on %d\n',port);
  end
end
