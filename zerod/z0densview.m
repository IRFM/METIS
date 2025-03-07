% visulaisation 3D des profils 0D
function z0densview

post = evalin('base','post');

% core profile for IMAS profiles at WEST
cp =[];

%  switch post.z0dinput.machine
%  case 'JET'
%      try
%      if nargin('mdsdisconnect') > 0
%  	  mdsdisconnect( 'mdsplus.jet.efda.org');
%      else 
%  	  mdsdisconnect;
%      end
%      end
%  end


zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo = post.z0dinput.geo;
ts    = zs.temps;

profli = post.profil0d;
xli    = profli.xli;
t      = profli.temps;


h = findobj(0,'type','figure','tag','z0densview');
if isempty(h)
       h=figure('tag','z0densview');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')


subplot(1,3,1)
switch post.z0dinput.machine
case 'TS'
   zplotprof(gca,t,xli,profli.nep./1e19,'color','r');
   leg= {'simulation'};
   occ = 8;
   ne = [];
   while isempty(ne) & (occ >= 0)
   	[ne,tne,rne,cne] = tsbase(post.z0dinput.shot + occ /10,'gprofnefit');
	occ = occ -1;
   end
   if ~isempty(ne)
 	zplotprof(gca,tne,rne,ne./1e19,'color','k','linestyle','none','marker','o');  
	leg{end+1} = 'data fit';
   end
   [gne,tgne,rgne,cgne] = tsbase(post.z0dinput.shot,'gne');
   if ~isempty(gne)
 	zplotprof(gca,tgne,rgne,gne./1e19,'color',[0,0.7,0],'linestyle','none','marker','+');     
	leg{end+1} = 'GNE';
   end
   [neth,tth,yth,cth] = tsbase(post.z0dinput.shot,'gnethom');
   if ~isempty(neth)
	vt = ones(size(tth));
   	rlas = tsbase(post.z0dinput.shot,'rthomlas');
	if isempty(rlas)
		rlas = 2.37;
	end
	xxth = cat(2,-profli.xli(end:-1:2),profli.xli);   
        amat = interp1(zs.temps,geo.a,tth,'nearest');
        z0 = interp1(zs.temps,geo.z0,tth,'nearest');
        rmat = interp1(profli.temps,profli.Raxe,tth,'nearest');
	rmat = cat(2,rmat(:,end:-1:2),rmat);
	zth  = (amat * xxth) .* sqrt(1 - ((rlas - rmat) ./  (amat * max(eps,xxth))) .^ 2) + z0 * ones(size(xxth));
	zvoid = 1e38 + vt * xxth .* 1e38;
	zth(:,21) = 0;
	indi      = find(imag(zth));
	zth(indi) = zvoid(indi);
	xth  = abs(tsplinet(zth,vt *xxth,vt * yth));
	neth(neth<=0) =NaN;
	xth(xth>=0.98) =NaN;
	
 	zplotprof(gca,tth,xth,neth./1e19,'color',[0,0.7,0],'linestyle','none','marker','x');  
        leg{end+1} = 'N_e_T_h_o_m';
  
   end
   legend(leg)
   xlabel('x (normalized radius)')

case 'WEST'

    % two point model information
    twopts = twopts_dc(post.z0dinput.shot);

    tigni = tsbase(post.z0dinput.shot,'rignitron');
    if isempty(tigni)
      tigni = 0;
    end
    leg = {};
    % read profile in coreprofile
    cp = imas_west_get(post.z0dinput.shot, 'core_profiles',0, 0,'imas_public','west');
    if isfield(cp,'profiles_1d') && ~isempty(cp.profiles_1d{1}.electrons.density)
	rhoprof = NaN * ones(length(cp.profiles_1d),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	psinprof = NaN * ones(length(cp.profiles_1d),length(cp.profiles_1d{1}.grid.psi));
	neprof = NaN * ones(length(cp.profiles_1d),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	for k=1:length(cp.profiles_1d)
	    tprof(k)  = cp.profiles_1d{k}.time - tigni;
	    rhoprof(k,:) = cp.profiles_1d{k}.grid.rho_tor_norm;
	    psi = cp.profiles_1d{k}.grid.psi;
	    psinprof(k,:) =  (psi - psi(1)) ./(psi(end) - psi(1)); 
	    neprof(k,:) =   cp.profiles_1d{k}.electrons.density;
	    neprof_upper(k,:) =   cp.profiles_1d{k}.electrons.density + ...
			cp.profiles_1d{k}.electrons.density_error_upper;
	    neprof_lower(k,:) =   cp.profiles_1d{k}.electrons.density - ...
			cp.profiles_1d{k}.electrons.density_error_lower;
	end
	zplotprof(gca,tprof,rhoprof,neprof./1e19, ...
	  'color','b','linestyle','-','marker','none');
	  leg{end+1} = 'N_e Fast NF';

	
    end
    % read coreprofile from NICE
    nice_on = 0;
    profil_nice = imas_west_get(post.z0dinput.shot, 'core_profiles',0, 1,'imas_public','west');
    if isfield(profil_nice,'profiles_1d')
	tprof_nice = profil_nice.time - tigni;
	rhoprof_nice = [];
	psinprof_nice = [];
	neprof_nice = [];
	for k=1:length(profil_nice.profiles_1d)
	    
	    if ~isempty(profil_nice.profiles_1d{k}.grid.rho_tor_norm)
		if isempty(rhoprof_nice)
		      rhoprof_nice = NaN * ones(length(profil_nice.profiles_1d),1001);
		end
		rhoprof_nice(k,1:length(profil_nice.profiles_1d{k}.grid.rho_tor_norm)) = profil_nice.profiles_1d{k}.grid.rho_tor_norm;
	    end
	    if ~isempty(profil_nice.profiles_1d{k}.grid.psi)
	      if isempty(psinprof_nice)
		  psinprof_nice = NaN * ones(length(profil_nice.profiles_1d),1001);
	      end
	      psi_nice = profil_nice.profiles_1d{k}.grid.psi;
	      psinprof_nice(k,1:length(psi_nice)) =  (psi_nice - psi_nice(1)) ./ (psi_nice(end) - psi_nice(1));
	    end
	    if ~isempty(profil_nice.profiles_1d{k}.electrons.density)
		if isempty(neprof_nice)
		    neprof_nice = NaN * ones(length(profil_nice.profiles_1d),1001);
		    nice_on = 1; 
		end
		neprof_nice(k,1:length(profil_nice.profiles_1d{k}.electrons.density)) =   profil_nice.profiles_1d{k}.electrons.density;
	    end
	end
	%
	% 
	% resampling on a fixed grid
	if nice_on
	    nb = 101;
	    for k=1:length(tprof_nice)
		x = linspace(0,1,nb);
		indf = find(isfinite(rhoprof_nice(k,:)));
		if isempty(indf)
		    neprof_nice_rho(k,:) = NaN * x;
        else
            if sum(rhoprof_nice(k,indf) == 0) > 1
                ind0 = max(find(rhoprof_nice(k,indf) == 0));
                indf = indf(ind0:end);
            end
            if length(indf) > 1
                neprof_nice_rho(k,:) = interp1(rhoprof_nice(k,indf),neprof_nice(k,indf),x,'linear',0);
            else
                neprof_nice_rho(k,:) = NaN * x;               
            end
		end
		indf = find(isfinite(psinprof_nice(k,:)));
		if isempty(indf)
		    neprof_nice_psin(k,:) = NaN * x;
		else          
		    neprof_nice_psin(k,:) = interp1(psinprof_nice(k,indf),neprof_nice(k,indf),x,'linear',0);
	      end
	  end
	  indok = find(tprof_nice >= 0);
	  zplotprof(gca,tprof_nice(indok),ones(size(tprof_nice(indok))) * x,neprof_nice_rho(indok,:)./1e19, ...
		'color','c','linestyle','-','marker','none');
	  leg{end+1} = 'N_e Nice';

	end
    end

  %psin = profli.psi - profli.psi(:,1) * ones(size(profli.xli));
  %psin = psin ./ ( (psin(:,end) - psin(:,1)) * ones(size(profli.xli)));
  rtn = profli.rmx ./ (profli.rmx(:,end) * ones(1,size(profli.rmx,2)));
  zplotprof(gca,t,rtn,profli.nep./1e19,'color','r');
  leg{end+1} = 'N_e METIS';
  
  if isfield(cp,'profiles_1d') && ~isempty(cp.profiles_1d{1}.electrons.density)
	zplotprof(gca,tprof,rhoprof,neprof_upper ./ 1e19 +  ...
	sqrt(-1) .* neprof_lower ./1e19, ...
	'facecolor',[0.7,0.7,0.85]);
	  leg{end+1} = 'error on N_e Fast NF';
	
  end
   
  % two points model result 
  if ~isempty(twopts)
      zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),twopts.ne_lcfs * cat(2,1,1) .* 1e-19,'color','k','linestyle','none','marker','*');
      leg{end+1} = 'n_{e,LCFS} from two points model';
      zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),cat(2,twopts.ne_lcfs_o,twopts.ne_lcfs_i) .* 1e-19,'color','k','linestyle','-','marker','none');
      leg{end+1} = 'inner/outer n_{e,LCFS} from two points model span';  
  end

  
  legend(leg,'Location','best')
  xlabel('\rho_{tor,norm}')
  %ylabel('10^{19} m^{-3}');
  

case 'JET'
  leg = {};
  liste        = {};
  liste{end+1} = 'ppf/@shot/HRTS/NE';
  liste{end+1} = 'ppf/@shot/HRTS/DNE';
  liste{end+1} = 'ppf/@shot/HRTS/RMID';     
  da_ne          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da_ne          = da_ne.ppf.HRTS;
  if isempty(da_ne.NE.t) || ischar(da_ne.NE.data) 
	  liste        = {};
	  liste{end+1} = 'ppf/@shot/HRTS/NE?uid=chain1';
	  liste{end+1} = 'ppf/@shot/HRTS/DNE?uid=chain1';
	  liste{end+1} = 'ppf/@shot/HRTS/RMID?uid=chain1';     
	  da_ne           = cgcgetjet(post.z0dinput.shot,liste,'','');
	  da_ne           = da_ne.ppf.HRTS;
  end
  if ~isempty(da_ne.NE.t) && ~isempty(da_ne.DNE.data)
        thrts = da_ne.NE.t * ones(size(da_ne.NE.x(:)'));
  elseif ~isempty(da_ne.NE.t) && ~isempty(da_ne.NE.data)
        thrts = da_ne.NE.t;
  end
  if ~isempty(da_ne.RMID.t)
	  rhrts = da_ne.RMID.data;
  elseif ~isempty(da_ne.NE.t) && ~isempty(da_ne.DNE.data)
	  rhrts = ones(size(da_ne.NE.t(:))) * da_ne.NE.x(:)';
  else 
      rhrts = [];
  end
  if ~isempty(da_ne.DNE.t) &&  ~isempty(da_ne.DNE.data)
    nehrts=da_ne.NE.data;
    ne_plus  = da_ne.NE.data + da_ne.DNE.data;
    ne_moins = da_ne.NE.data - da_ne.DNE.data;
    ne_max = 3 * ceil(max(profli.nep(:)./1e19))*1e19;
    indbad = find((ne_moins < 0) | (ne_plus > ne_max));
    ne_plus(indbad) = NaN;
    ne_moins(indbad) = NaN;
    nehrts_e = make_data4errorbar(ne_moins,ne_plus);
    rhrts_e  = make_data4errorbar(rhrts);
    thrts_e  = make_data4errorbar(thrts);
    try
	  zplotprof(gca,thrts,rhrts,nehrts./1e19, ...
	  'color','b','linestyle','none','marker','o');
	  zplotprof(gca,thrts_e,rhrts_e,nehrts_e./1e19, ...
	  'color','b','linestyle','-','marker','none');         
	  leg{end+1} = 'N_e HRTS';
	  leg{end+1} = 'error N_e HRTS';

    end
  else
      nehrts=da_ne.NE.data;
      try
	  zplotprof(gca,thrts,rhrts,nehrts./1e19, ...
	  'color','b','linestyle','none','marker','o');
	  leg{end+1} = 'N_e HRTS';
      end
  end
  liste        = {};
  liste{end+1} = 'ppf/@shot/LIDR/NE';
  da           = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf.LIDR;
  try
 	 zplotprof(gca,da.NE.t,da.NE.x,da.NE.data./1e19, ...
            'color','m','linestyle','none','marker','o');
         leg{end+1} = 'N_e LIDAR';
  end
  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.nep./1e19,'color','r');
  leg{end+1} = 'N_e METIS';
  legend(leg,'Location','best')
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.nep./1e19,'color','r');
  xlabel('R (m)')
otherwise
        zplotprof(gca,t,xli,profli.nep./1e19,'color','r');
        xlabel('r / a')
	
end
ylabel('10^1^9 m^-^3')

subplot(1,3,2)
switch post.z0dinput.machine
case 'TS'
   zplotprof(gca,t,xli,profli.tep./1e3,'color','r');
   leg = {'T_e'};
   occ = 8;
   te = [];
   while isempty(te) & (occ >= 0)
   	[te,tte,rte,cte] = tsbase(post.z0dinput.shot + occ /10,'gproftefit');
	occ = occ -1;
   end
   if ~isempty(te)
 	zplotprof(gca,tte,rte,te,'color','k','linestyle','none','marker','o');  
        leg{end+1} = 'T_e_f_i_t';
   end
   [tece,ttece] = tsbase(post.z0dinput.shot,'gshtenv');
   if isempty(tece)
    	[tece,ttece] = tsbase(post.z0dinput.shot,'gshte');
   end
   if ~isempty(tece)
   	[rtece,ttece] = tsbase(post.z0dinput.shot,'gshr');
   else
        [rtece,ttece] = tsbase(post.z0dinput.shot,'gtrshr');
    	[tece,ttece] = tsbase(post.z0dinput.shot,'gtrshte'); 
        tece = tece ./1e3;
        rtece = rtece ./1e3;
   end
   if ~isempty(tece)
   	ttece = ttece(:,1);
	vt = ones(size(ttece));
        amat = interp1(zs.temps,geo.a,ttece,'nearest');
        rmat = interp1(profli.temps,profli.Raxe,ttece,'nearest');
	xxece = vt * cat(2,-profli.xli(end:-1:2),profli.xli);   
	Rxece = cat(2,rmat(:,end:-1:2) - amat * profli.xli(end:-1:2),rmat + amat * profli.xli);
	xece  = abs(tsplinet(Rxece,xxece,rtece));
	tece(tece <=0) = NaN;
	zplotprof(gca,ttece,xece,tece,'color',[0 0.7 0],'linestyle','none','marker','+');  
        leg{end+1} = 'T_e_c_e';
	
   end 
   
   [teth,tth,yth,cth] = tsbase(post.z0dinput.shot,'gtethom');
   if ~isempty(teth)
	vt = ones(size(tth));
   	rlas = tsbase(post.z0dinput.shot,'rthomlas');
	if isempty(rlas)
		rlas = 2.37;
	end
	xxth = cat(2,-profli.xli(end:-1:2),profli.xli);   
        amat = interp1(zs.temps,geo.a,tth,'nearest');
        z0 = interp1(zs.temps,geo.z0,tth,'nearest');
        rmat = interp1(profli.temps,profli.Raxe,tth,'nearest');
	rmat = cat(2,rmat(:,end:-1:2),rmat);
	zth  = (amat * xxth) .* sqrt(1 - ((rlas - rmat) ./  (amat * max(eps,xxth))) .^ 2) + z0 * ones(size(xxth));
	zvoid = 1e38 + vt * xxth .* 1e38;
	zth(:,21) = 0;
	indi      = find(imag(zth));
	zth(indi) = zvoid(indi);
	xth  = abs(tsplinet(zth,vt *xxth,vt * yth));
	teth(teth<=0) = NaN;
	xth(xth>=0.98) = NaN;
	
 	zplotprof(gca,tth,xth,teth,'color',[0 0.7 0],'linestyle','none','marker','x');  
        leg{end+1} = 'T_T_h_o_m';
   end
   
   [teb,tb] = tsbase(post.z0dinput.shot,'stebrag');
   if ~isempty(teb)
   	teb(teb >7| teb < 0) = NaN;
	
   	zplotprof(gca,tb,ones(size(tb)) *cat(2,0,0.3) ,cat(2,teb,teb),'color','m');
        leg{end+1} = 'T_e_B_r_a_g_g';
   	
   end
   
   legend(leg);
   xlabel('x (normalized radius)')
   
case 'WEST'
   leg = {};
    tigni = tsbase(post.z0dinput.shot,'rignitron');
    if isempty(tigni)
      tigni = 0;
    end
    try
	ece_ids = imas_west_get(post.z0dinput.shot,'ece');
    catch
	ece_ids = [];
    end
    if ~isempty(ece_ids)
      time_ece = ece_ids.time - tigni;
      te_ece = NaN * ones(length(time_ece),length(ece_ids.channel));
      dte_p  = NaN * ones(length(time_ece),length(ece_ids.channel));
      dte_m  = NaN * ones(length(time_ece),length(ece_ids.channel));
      r_ece = NaN * ones(length(time_ece),length(ece_ids.channel));
      rho_ece     = NaN * ones(length(time_ece),length(ece_ids.channel));
      opd_ece     = NaN * ones(length(time_ece),length(ece_ids.channel));
      for k=1:length(ece_ids.channel)
	  if ~isempty(ece_ids.channel{k}.t_e.data)
	      te_ece(:,k) = ece_ids.channel{k}.t_e.data(:);  
	      dte_p(:,k)  = ece_ids.channel{k}.t_e.data(:) + ece_ids.channel{k}.t_e.data_error_upper(:);
	      dte_m(:,k)  = ece_ids.channel{k}.t_e.data(:) - ece_ids.channel{k}.t_e.data_error_lower(:);
	  end
	  if  ~isempty(ece_ids.channel{k}.position.r.data)
	      r_ece(:,k) = ece_ids.channel{k}.position.r.data(:);  
	  end
	  if  ~isempty(ece_ids.channel{k}.position.rho_tor_norm.data)
	      rho_ece(:,k) = ece_ids.channel{k}.position.rho_tor_norm.data(:);
	  end
	  if  ~isempty(ece_ids.channel{k}.optical_depth.data)
	      opd_ece(:,k) = ece_ids.channel{k}.optical_depth.data(:);
	  end
      end
      % compute rho_cut and ind_cut
      fdrhodr = NaN * ones(size(time_ece));
      for k=1:length(time_ece)
	      % whith R span LFH and HFS
	      ind_ok = find(isfinite(r_ece(k,:)));
	      if ~isempty(ind_ok)
		  fdrhodr(k) = (max(rho_ece(k,ind_ok)) - min(rho_ece(k,ind_ok))) ./ ...
		      (max(r_ece(k,ind_ok)) - min(r_ece(k,ind_ok))) ./ 2;
	      end
      end

      % compute locla width of emission with optical depth
      fte_opd   = te_ece / 1e3;
      delta_opd = max(0,r_ece ./ opd_ece ./ fte_opd ./ (2*sqrt(2*log(2))));
      % in rho unity, error in localisation can't be greater than
      % plasma size
      delta_opd = min(1,(fdrhodr * ones(1,size(delta_opd,2))) .* delta_opd);

      te_ece(~isfinite(r_ece)) = 0;
      dte_p(~isfinite(r_ece)) = 0;
      dte_m(~isfinite(r_ece)) = 0;
      r_ece(~isfinite(r_ece)) = -1;
      rho_ece(~isfinite(r_ece)) = -1;
      opd_ece(~isfinite(r_ece)) = NaN;
      te_ece(abs(te_ece) > 1e37) = 0;
      dte_p(abs(dte_p) > 1e37) = 0;
      dte_m(abs(dte_m) > 1e37) = 0;
      r_ece(abs(r_ece) > 1e37) = -1;
      rho_ece(abs(rho_ece) > 1e37) = -1;
      opd_ece(abs(opd_ece) > 1e37) = NaN;
      ind_ok = find((time_ece >= 0) & (time_ece <= (zs.temps(end) + 2)));
      zplotprof(gca,time_ece(ind_ok),rho_ece(ind_ok,:),te_ece(ind_ok,:)./1e3, ...
	      'color','b','linestyle','none','marker','o');
      leg{end+1} = 'T_e ECE ';
    
      te_plus  = abs(dte_p(ind_ok,:))./1e3;
      te_moins = abs(dte_m(ind_ok,:))./1e3;
      tehrts_e = make_data4errorbar(te_moins,te_plus);
      rhrts    = rho_ece(ind_ok,:);
      rhrts_e  = make_data4errorbar(rhrts);
      thrts    = time_ece(ind_ok) * ones(1,size(te_plus,2));
      thrts_e  = make_data4errorbar(thrts);

      zplotprof(gca,thrts_e,rhrts_e,tehrts_e, ...
	  'color','b','linestyle','-','marker','none');
    	  leg{end+1} = 'error T_e ECE';
    	  
      tehrts_e   = make_data4errorbar(te_ece(ind_ok,:) / 1e3);
      rhrts_m    = rho_ece(ind_ok,:) - delta_opd(ind_ok,:) / 2;
      rhrts_p    = rho_ece(ind_ok,:) + delta_opd(ind_ok,:) /2;
      rhrts_e  = make_data4errorbar(rhrts_m,rhrts_p);
      thrts    = time_ece(ind_ok) * ones(1,size(te_ece,2));
      thrts_e  = make_data4errorbar(thrts);
      zplotprof(gca,thrts_e,rhrts_e,tehrts_e, ...
	  'color','c','linestyle','-','marker','none');
	   leg{end+1} = 'error T_e ECE (position ~ optical depth)';
	   
	   
   
    	  

    end
    
    
    % read profile data
    if isempty(cp)
      try
	    cp = imas_west_get(post.z0dinput.shot,'core_profiles');
      catch    
	  cp = [];
      end
   end
    % swaptime
    if isempty(cp)
	tprof         = [];
	teprof        = [];
	dteprof_upper = [];
	dteprof_lower = [];
	rho_prof      = [];   
    elseif ~isempty(cp.time)
	tprof = cp.time;
    else
	tprof = NaN * ones(length(cp.profiles_1d),1);
	for k=1:length(cp.profiles_1d);
	    tprof(k) = cp.profiles_1d{k}.time;
	end
    end
    if ~isempty(cp) && (length(cp.profiles_1d{1}.electrons.temperature) > 1)
	tprof         = tprof - tigni;
	teprof        = NaN * ones(length(tprof),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	dteprof_upper = NaN * ones(length(tprof),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	dteprof_lower = NaN * ones(length(tprof),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	rho_prof      = NaN * ones(length(tprof),length(cp.profiles_1d{1}.grid.rho_tor_norm));
	
	for k=1:length(cp.profiles_1d)
	    teprof(k,:)        = cp.profiles_1d{k}.electrons.temperature;
	    dteprof_upper(k,:) = cp.profiles_1d{k}.electrons.temperature + cp.profiles_1d{k}.electrons.temperature_error_upper;
	    dteprof_lower(k,:) = cp.profiles_1d{k}.electrons.temperature - cp.profiles_1d{k}.electrons.temperature_error_lower;
	    rho_prof(k,:)      = cp.profiles_1d{k}.grid.rho_tor_norm;
	end
	
	zplotprof(gca,tprof,rho_prof, teprof ./1e3,'color','k','linestyle','-','marker','none');
        leg{end+1} = 'fit T_e ECE ';

    end

    
    
    
    
    %amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
    %rli   = profli.Raxe + amat * profli.xli;
    %zplotprof(gca,t,rli,profli.tep./1e3,'color','r');
    rtn = profli.rmx ./ (max(profli.rmx,[],2) * ones(1,size(profli.rmx,2)));
    zplotprof(gca,t,rtn,profli.tep./1e3,'color','r');
    leg{end+1} = 'T_e METIS';
    
    
    if  ~isempty(cp) && (length(cp.profiles_1d{1}.electrons.temperature) > 1)
	  zplotprof(gca,tprof,rho_prof,dteprof_upper ./ 1e3 +  ...
	  sqrt(-1) .* dteprof_lower ./1e3, ...
	  'facecolor',[0.7,0.7,0.7]);
	leg{end+1} = 'error fit T_e ECE';
    end

    % two points model result 
    if ~isempty(twopts)
      zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),twopts.te_lcfs * cat(2,1,1) .* 1e-3,'color','k','linestyle','none','marker','*');
      leg{end+1} = 'n_{e,LCFS} from two points model';
      zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),cat(2,twopts.te_lcfs_o,twopts.te_lcfs_i) .* 1e-3,'color','k','linestyle','-','marker','none');
      leg{end+1} = 'inner/outer n_{e,LCFS} from two points model span';  
    end

   
    
    legend(leg,'Location','best')
    %rli   = profli.Raxe - amat * profli.xli;
    %zplotprof(gca,t,rli,profli.tep./1e3,'color','r');
    %xlabel('R (m)')
    xlabel('\rho_{tor,norm}')
    set(gca,'xlim',[0 1]);
   
case 'JET'
  leg  = {};
  liste        = {};
  liste{end+1} = 'ppf/@shot/HRTS/TE';
  liste{end+1} = 'ppf/@shot/HRTS/DTE';
  liste{end+1} = 'ppf/@shot/HRTS/RMID';     
  da_te          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da_te          = da_te.ppf.HRTS;
  if isempty(da_te.TE.t) || ischar(da_te.TE.data) 
	  liste        = {};
	  liste{end+1} = 'ppf/@shot/HRTS/TE?uid=chain1';
	  liste{end+1} = 'ppf/@shot/HRTS/DTE?uid=chain1';
	  liste{end+1} = 'ppf/@shot/HRTS/RMID?uid=chain1';     
	  da_te           = cgcgetjet(post.z0dinput.shot,liste,'','');
	  da_te           = da_te.ppf.HRTS;
  end
  if isempty(da_te.TE.t) || ischar(da_te.TE.data) 
	  liste        = {};
	  liste{end+1} = 'ppf/@shot/LIDR/TE';
	  liste{end+1} = 'ppf/@shot/LIDR/DTE';
	  liste{end+1} = 'ppf/@shot/LIDR/RMID';   
	  da_te           = cgcgetjet(post.z0dinput.shot,liste,'','');
	  da_te           = da_te.ppf.LIDR;
  end
  if ~isempty(da_te.TE.t)  &&  ~isempty(da_te.DTE.data)
        thrts = da_te.TE.t * ones(size(da_te.TE.x(:)'));
  elseif ~isempty(da_te.TE.t) && ~isempty(da_te.TE.data)
    thrts = da_te.TE.t;
  end
  if ~isempty(da_te.RMID.t)
	  rhrts = da_te.RMID.data;
  elseif ~isempty(da_te.TE.t)  &&  ~isempty(da_te.TE.x)
	  rhrts = ones(size(da_te.TE.t(:))) * da_te.TE.x(:)';
  else
      rhrts = [];
  end
  if ~isempty(da_te.DTE.t) &&  ~isempty(da_te.DTE.data)
    tehrts=da_te.TE.data;
    te_plus  = da_te.TE.data + da_te.DTE.data;
    te_moins = da_te.TE.data - da_te.DTE.data;
    te_max = 3 * ceil(max(profli.tep(:)./1e3))*1e3;
    indbad = find((te_moins < 0) | (te_plus > te_max));
    te_plus(indbad) = NaN;
    te_moins(indbad) = NaN;
    tehrts_e = make_data4errorbar(te_moins,te_plus);
    rhrts_e  = make_data4errorbar(rhrts);
    thrts_e  = make_data4errorbar(thrts);
    try
	  zplotprof(gca,thrts,rhrts,tehrts./1e3, ...
	  'color','b','linestyle','none','marker','o');
	  leg{end+1} = 'error T_e HRTS';
	  zplotprof(gca,thrts_e,rhrts_e,tehrts_e./1e3, ...
	  'color','b','linestyle','-','marker','none');
	  leg{end+1} = 'T_e HRTS';
    end
  else
      tehrts=da_te.TE.data;
      try
	  zplotprof(gca,thrts,rhrts,tehrts./1e3, ...
	  'color','b','linestyle','none','marker','o');
	  leg{end+1} = 'T_e HRTS';
      end
  end
  liste        = {};
  liste{end+1} = 'ppf/@shot/LIDR/TE';
  da           = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf.LIDR;
  try 
  	zplotprof(gca,da.TE.t,da.TE.x,da.TE.data./1e3, ...
            'color','m','linestyle','none','marker','o');
        leg{end+1} = 'T_e LIDAR';
  end
  
  % ECE
  liste        = {};
  liste{end+1} = 'ppf/@shot/KK3/TPRF';
  liste{end+1} = 'ppf/@shot/KK3/CPRF'; 
  da           = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf.KK3;
  if ~isempty(da.TPRF.data) && ~isempty(da.CPRF.data)
    tesh         = da.TPRF.data;
    ttesh        = da.TPRF.t;
    rcsh         = da.TPRF.x;
    Rsh          = da.CPRF.data;
    try 
	  zplotprof(gca,ttesh,Rsh,tesh./1e3, ...
	      'color','g','linestyle','none','marker','o');
	  leg{end+1} = 'T_e ECE (KK3)';
    end
  end
  
  
  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.tep./1e3,'color','r');
  leg{end+1} = 'T_e METIS';
  legend(leg,'Location','best')
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.tep./1e3,'color','r');
  xlabel('R (m)')

otherwise
        zplotprof(gca,t,xli,profli.tep./1e3,'color','r');
	legend('T_e')	
        xlabel('r / a')
end
ylabel('keV');
title(sprintf('Metis : %s@%d/Profiles: comparison to experiment', ...
             post.z0dinput.machine,post.z0dinput.shot));

subplot(1,3,3)
switch post.z0dinput.machine
case 'TS'
   zplotprof(gca,t,xli,profli.tip./1e3,'color','r');
   leg = {'T_i'};   
   [tib,tb] = tsbase(post.z0dinput.shot,'stibrag');
   if ~isempty(tib)
   	tib(tib >7 | tib < 0) = NaN;
	
   	zplotprof(gca,tb,ones(size(tb)) *cat(2,0,0.3) ,cat(2,tib,tib),'color','c');
        leg{end+1} = 'T_i_B_r_a_g_g';
   	
   end
   
   legend(leg);
   xlabel('x (normalized radius)')
   
case 'JET'

  leg = {};
  liste        = {};
  liste{end+1}   = 'ppf/@shot/CXG6/TI';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXG6/RCOR';      % R of measurment
  liste{end+1}   = 'ppf/@shot/CXG6/TILO';    %    eV  minimum  Ti 	
  liste{end+1}   = 'ppf/@shot/CXG6/TIHI';    %	eV maximum	 Ti 	
  liste{end+1}   = 'ppf/@shot/CXG6/TICR';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXG6/TIRH';      % profil de Ti
  liste = switchuid(liste,post.z0dinput.shot);
  da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da_ti          = da_ti.ppf.CXG6;
  if isempty(da_ti.TI.data)
      liste        = {};
      liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
      liste{end+1}   = 'ppf/@shot/CXFM/RCOR';      % R of measurment
      liste{end+1}   = 'ppf/@shot/CXFM/TILO';    %    eV  minimum  Ti 	
      liste{end+1}   = 'ppf/@shot/CXFM/TIHI';    %	eV maximum	 Ti 	
      liste{end+1}   = 'ppf/@shot/CXFM/TICR';      % profil de Ti
      liste{end+1}   = 'ppf/@shot/CXFM/TIRH';      % profil de Ti
      liste = switchuid(liste,post.z0dinput.shot);
      da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
      da_ti          = da_ti.ppf.CXFM;
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CXSM/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXSM/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CXSM/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CXSM/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CXSM/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXSM/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CXSM;  
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CXGM/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXGM/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CXGM/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CXGM/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CXGM/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXGM/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CXGM;  
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CXD6/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXD6/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CXD6/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CXD6/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CXD6/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXD6/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CXD6;  
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CXRS/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXRS/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CXRS/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CXRS/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CXRS/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXRS/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CXRS;  
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CX7P/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CX7P/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CX7P/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CX7P/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CX7P/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CX7P/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CX7P;  
  end
  if isempty(da_ti.TI.data)
    liste        = {};
    liste{end+1}   = 'ppf/@shot/CXSE/TI';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXSE/RCOR';      % R of measurment
    liste{end+1}   = 'ppf/@shot/CXSE/TILO';    %    eV  minimum  Ti 	
    liste{end+1}   = 'ppf/@shot/CXSE/TIHI';    %	eV maximum	 Ti 	
    liste{end+1}   = 'ppf/@shot/CXSE/TICR';      % profil de Ti
    liste{end+1}   = 'ppf/@shot/CXSE/TIRH';      % profil de Ti
    liste = switchuid(liste,post.z0dinput.shot);
    da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
    da_ti          = da_ti.ppf.CXSE;  
  end
   
   

  
  
  if ~isempty(da_ti.TI.t)
    fprintf('Wraning: first CX measurement @ t = %g s\n',min(da_ti.TI.t)); 
  end
  if ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TILO.data)
        thrts = da_ti.TI.t * ones(size(da_ti.TI.x(:)'));
  elseif ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TI.data)
	thrts = da_ti.TI.t;
  end
  if ~isempty(da_ti.RCOR.t)
	  rhrts = da_ti.RCOR.data;
  elseif ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TILO.data)
	  rhrts = ones(size(da_ti.TI.t)) * da_ti.TI.x(:)';
  else
	  rhrts = [];
  end
  if ~isempty(da_ti.TILO.t) && ~isempty(da_ti.TILO.data)
	  tihrts=da_ti.TI.data;
	  ti_plus  = da_ti.TIHI.data;
	  ti_moins = da_ti.TILO.data;
	  ti_max = 3 * ceil(max(profli.tip(:)./1e3))*1e3;
	  indbad = find((ti_moins < 0) | (ti_plus > ti_max));
	  ti_plus(indbad) = NaN;
	  ti_moins(indbad) = NaN;
	  tihrts_e  = make_data4errorbar(ti_moins,ti_plus);
	  rhrts_e   = make_data4errorbar(rhrts);
	  thrts_e   = make_data4errorbar(thrts);
	  try
		  zplotprof(gca,thrts,rhrts,tihrts./1e3, ...
		  'color','b','linestyle','none','marker','o');
		  leg{end+1} = sprintf('CX (first @ %g s)',min(da_ti.TI.t));
		  zplotprof(gca,thrts_e,rhrts_e,tihrts_e./1e3, ...
		  'color','b','linestyle','-','marker','none');
		  leg{end+1} ='error CX';
	  end
  else
	  try
		  zplotprof(gca,thrts,rhrts,da_ti.TI.data ./1e3,'color','b','marker','o','linestyle','none');
		  leg{end+1} = sprintf('CX (first @ %g s)',min(da_ti.TI.t));
	  end 
 
  end
  try
	  zplotprof(gca,thrts,rhrts,da_ti.TICR.data ./1e3,'color','m','marker','+','linestyle','none');
	  leg{end+1} ='CX corrected';
  end 
  

  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.tip./1e3,'color','r');	
  leg{end+1} = 'T_i METIS';
  legend(leg,'Location','best')
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.tip./1e3,'color','r');	
  xlabel('R (m)')
  
case 'WEST'
        zplotprof(gca,t,profli.rmx./ (max(profli.rmx,[],2) * ones(1,size(profli.rmx,2))),profli.tip./1e3,'color','b');
	legend('T_i METIS')	
        xlabel('\rho_{tor,norm}')
 
otherwise
        zplotprof(gca,t,xli,profli.tip./1e3,'color','b');
	legend('T_i METIS')	
        xlabel('r / a')
end
ylabel('keV');


% pour TS plot avec reflex
switch post.z0dinput.machine
    case 'TS'
        h = findobj(0,'type','figure','tag','z0densview2');
        if isempty(h)
            h=figure('tag','z0densview2');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
        colormap('hot')
        
        rli  = cat(2,profli.Raxe(:,end:-1:2) - profli.rmx(:,end:-1:2),profli.Raxe + profli.rmx);
        nep  = cat(2,profli.nep(:,end:-1:2),profli.nep)./1e19;
        zplotprof(gca,t,rli,nep,'color','r');
        leg = {'simulation'};
        
        [ner,tner,rner,cner] = tsbase(post.z0dinput.shot,'grefnex');
        [rner,rtner,rrner,rcner] = tsbase(post.z0dinput.shot,'grefrx');
        if ~isempty(ner) & ~isempty(rner)
	
		zplotprof(gca,tner,rner,ner./1e19,'color','b','linestyle','none','marker','o'); 
		
%  		% recherche du decalage
%  		rla = interp1(t,rli(:,end),tner,'nearest');
%  		nea_ = interp1(t,profli.nep(:,end),tner,'nearest');
%  		dd   = abs(ner - nea_ * ones(1,size(ner,2)));
%  		mask = (dd == min(dd,[],2) * ones(1,size(ner,2)));
%  		delta = sum((rner - rla * ones(1,size(ner,2))) .* mask,2) ./ max(1,sum(mask,2));
		
%  		zplotprof(gca,tner,rner-delta * ones(1,size(ner,2)) ,ner./1e19,'color','c','linestyle','none','marker','o'); 
		
		leg{end+1} = 'Reflex';
%  		leg{end+1} = 'Reflex with shift';
		

   	end
  	[ner,tner,rner,cner] = tsbase(post.z0dinput.shot,'gneflucint');
   	[rner,rtner,rrner,rcner] = tsbase(post.z0dinput.shot,'grflucint');
   	if ~isempty(ner) & ~isempty(rner)
		zplotprof(gca,tner,rner,ner./1e19,'color','m','linestyle','none','marker','+'); 
			
		leg{end+1} = 'Refluc (gne)';
	end	
  	[ner,tner,rner,cner] = tsbase(post.z0dinput.shot,'gneflucref');
   	[rner,rtner,rrner,rcner] = tsbase(post.z0dinput.shot,'grflucref');
   	if ~isempty(ner) & ~isempty(rner)
		zplotprof(gca,tner,rner,ner./1e19,'color','c','linestyle','none','marker','x'); 
			
		leg{end+1} = 'Refluc (ref)';
	end	
  	[ner,tner,rner,cner] = tsbase(post.z0dinput.shot,'gnebordint');
   	[rner,rtner,rrner,rcner] = tsbase(post.z0dinput.shot,'grbordint');
   	if ~isempty(ner) & ~isempty(rner)
		zplotprof(gca,tner,rner,ner./1e19,'color','k','linestyle','none','marker','*'); 
			
		leg{end+1} = 'Refluc (gne)';
	end	
	ylabel('10^1^9 m^3');
	xlabel('x (normalized radius)')
	legend(leg);
		
	title(sprintf('Metis : %s@%d/reflectometry  profile', ...
		post.z0dinput.machine,post.z0dinput.shot));
    
    case 'WEST'
    
       % try to read reflectometers data
       [reflec,tref,neref,rho_ref,psi_ref,psin_ref] = get_west_reflectometer_data(post.z0dinput.shot,0,0);
       [reflec1,tref1,neref1,rho_ref1,psi_ref1,psin_ref1] = get_west_reflectometer_data(post.z0dinput.shot,0,1);
       [reflec2,tref2,neref2,rho_ref2,psi_ref2,psin_ref2] = get_west_reflectometer_data(post.z0dinput.shot,0,2);
       
       if isempty(tref) && isempty(tref1) && isempty(tref2)
           return
       end
       
       tplusdip = ftplusdip(post.z0dinput.shot);

       
        h = findobj(0,'type','figure','tag','z0densview2');
        if isempty(h)
            h=figure('tag','z0densview2');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
        colormap('hot')
 
       leg = {};
        
       psin = (profli.psi - profli.psi(:,1) * ones(1,size(profli.psi,2))) ./ ((profli.psi(:,end) - profli.psi(:,1)) * ones(1,size(profli.psi,2)));
       zplotprof(gca,t,psin,profli.nep./1e19,'color','r');
         leg = {'METIS (simulation)'};
   
       
       if ~isempty(tref)
	  indok = find((tref >= 0) &(tref <= tplusdip));
	  zplotprof(gca,tref(indok),psin_ref(indok,:), neref(indok,:)/1e19,'color','b','linestyle','-','marker','.');
	  leg{end+1} = 'n_e from reflectometer DREFRAP';
       end
       if ~isempty(tref1)
	  indok = find((tref1 >= 0) &(tref1 <= tplusdip));
	  zplotprof(gca,tref1(indok),psin_ref1(indok,:), neref1(indok,:)/1e19,'color','c','linestyle','-','marker','.');
	  leg{end+1} = 'n_e from reflectometer DREFLUC initialised with DREFRAP';
       end
       if ~isempty(tref2)
	  indok = find((tref2 >= 0) &(tref2 <= tplusdip));
	  zplotprof(gca,tref2(indok),psin_ref2(indok,:), neref2(indok,:)/1e19,'color','m','linestyle','-','marker','.');
	  leg{end+1} = 'n_e from reflectometer DREFLUC initialised with interferometer';
       end
      
       
%         if isfield(profile_out,'NF_inversion')  && isfield(profile_out.NF_inversion,'density') && ~isempty(profile_out.NF_inversion.density)
%             indok = find((profile_out.NF_inversion.time >= 0) &(profile_out.NF_inversion.time <= tplusdip));
%             zplotprof(gca,profile_out.NF_inversion.time(indok),profile_out.NF_inversion.psin(indok,:), profile_out.NF_inversion.density(indok,:)/1e19,'color','c','linestyle','-','marker','none');
%             leg{end+1} = 'n_e from interferometer identified by fastNF';
%         end
%         if isfield(profile_out,'NICE_inversion') && isfield(profile_out.NICE_inversion,'density') && ~isempty(profile_out.NICE_inversion.density)
%           indok = find((profile_out.NICE_inversion.time >= 0) &(profile_out.NICE_inversion.time <= tplusdip));
%           zplotprof(gca,profile_out.NICE_inversion.time(indok),profile_out.NICE_inversion.psin(indok,:),profile_out.NICE_inversion.density(indok,:)/1e19,'color','m','linestyle','-','marker','none');
%           leg{end+1} = 'n_e from interferometer idnetified by NICE';
%         end
%         
        % two points model result 
	if ~isempty(twopts)
	    zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),twopts.ne_lcfs * cat(2,1,1) .* 1e-19,'color','k','linestyle','none','marker','*');
	    leg{end+1} = 'n_{e,LCFS} from two points model';
	    zplotprof(gca,twopts.time,cat(2,ones(size(twopts.time)),ones(size(twopts.time))),cat(2,twopts.ne_lcfs_o,twopts.ne_lcfs_i) .* 1e-19,'color','k','linestyle','-','marker','none');
	    leg{end+1} = 'inner/outer n_{e,LCFS} from two points model span';  
	end

      
 	ylabel('10^1^9 m^3');
	xlabel('normalized poloidal flux')
	legend(leg);
  
 	title(sprintf('Metis : %s@%d/reflectometry  profile', ...
		post.z0dinput.machine,post.z0dinput.shot));
    
end


switch post.z0dinput.machine
case 'TS'
	try
		z0dtsticxs;
	end
end


function liste = switchuid(liste,shot)

switch shot
case 86614
  for k=1:length(liste)
    liste{k} = sprintf('%s?uid=cgiroud+seq=977',liste{k});
  end
end

