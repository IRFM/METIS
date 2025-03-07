% script pour le test du model complet de SOL
if isfield(post.z0dinput.option,'Sn_fraction') && (post.z0dinput.option.Sn_fraction > 0)
    error('JET do not contain any tin (Sn) in plasma composition (option.Sn_fraction should be 0)');
end


zerod 	 = post.zerod;
profil0d = post.profil0d;
cons     = post.z0dinput.cons;
geo      = post.z0dinput.geo;
exp0d    = post.z0dinput.exp0d;
z0plotsc;
drawnow	
[tend_meas,void] = ginput(1);
drawnow
nbt     = 101;
ne_max     = log10(max(post.zerod.ne0));
nlcfs_list = logspace(17,ne_max,nbt)';
pmax    = max(post.zerod.pin) ./ min(1,post.z0dinput.option.fpower);
% offset pour tenir compte du rayonnement 
p_offset = max(0,max(post.zerod.plim - post.zerod.pin)) ./ min(1,post.z0dinput.option.fpower);
pin_list = linspace(1,pmax + p_offset,nbt+2);



temps = profil0d.temps;
ind_sol = find(temps >= tend_meas,1)-1;
vt      = ones(nbt,1);

% extraction et replication
times = cons.temps;
% modification des donnees
noms = fieldnames(cons);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(cons,nomc);
	valn  = interp1(times,val,tend_meas,'linear');
	indbad      = find(any(~isfinite(valn),2));
	if ~isempty(indbad)
		valn(indbad,:) = ones(length(indbad),1) * val(end,:);
	end
	cons = setfield(cons,nomc,vt*valn);
end
cons.temps = linspace(1,nbt,nbt)';

noms = fieldnames(geo);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(geo,nomc);
	if ~isempty(val)
		valn  = interp1(times,val,tend_meas,'linear');
		indbad      = find(any(~isfinite(valn),2));
		if ~isempty(indbad)
			valn(indbad,:) = ones(length(indbad),1) * val(end,:);
		end
		geo = setfield(geo,nomc,vt*valn);
	end
end

noms = fieldnames(zerod);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(zerod,nomc);
	if length(val) == length(times)
		val  = interp1(times,val,tend_meas,'nearest');
		zerod = setfield(zerod,nomc,vt*val);
	end
end

noms = fieldnames(exp0d);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(exp0d,nomc);
	if length(val) == length(times)
                warning off
		val  = interp1(times,val,tend_meas,'nearest');
		warning on
		exp0d = setfield(exp0d,nomc,vt*val);
	end
end

% modification des donnees
noms = fieldnames(profil0d);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(profil0d,nomc);
	if all(size(val) > 1)
		warning off
		valn  = interp1(temps,val,tend_meas,'linear');
		warning on
		indbad      = find(any(~isfinite(valn),2));
		if ~isempty(indbad)
			valn(indbad,:) = ones(length(indbad),1) * val(end,:);
		end
		profil0d = setfield(profil0d,nomc,vt*valn);
	end
end
profil0d.temps = cons.temps;

%rampe de densite
ve           = ones(size(profil0d.xli));
nebord_mem   = zerod.nebord;
% reglage ici
zerod.nebord = nlcfs_list;
%
zerod.nibord = zerod.nibord .* zerod.nebord  ./ nebord_mem;
nep_mem      = profil0d.nep;
profil0d.nep = profil0d.nep - profil0d.nep(:,end) * ve + zerod.nebord * ve;
profil0d.nip = profil0d.nep ./ nep_mem  .* profil0d.nip;
profil0d.n1p = profil0d.nep ./ nep_mem  .* profil0d.n1p;
profil0d.nhep = profil0d.nep ./ nep_mem  .* profil0d.nhep;
profil0d.nzp = profil0d.nep ./ nep_mem  .* profil0d.nzp;
profil0d.nwp = profil0d.nep ./ nep_mem  .* profil0d.nwp;
vol          = trapz(profil0d.xli,profil0d.vpr,2);
zerod.nem    = trapz(profil0d.xli,profil0d.vpr .* profil0d.nep,2) ./ vol;
zerod.nbar   = trapz(profil0d.xli,profil0d.nep,2);
zerod.nhem   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nhep,2) ./ vol;
zerod.nim   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nip,2) ./ vol;
m1m_mem     = zerod.n1m;
zerod.n1m   = trapz(profil0d.xli,profil0d.vpr .* profil0d.n1p,2) ./ vol;
zerod.nDm   = zerod.nDm .* zerod.n1m ./ m1m_mem;
zerod.nTm   = zerod.nTm .* zerod.n1m ./ m1m_mem;
zerod.nimpm   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nzp,2) ./ vol;
zerod.nwm   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nwp,2) ./ vol;
zerod.neped = profil0d.nep(:,end-1);
zerod.niped = profil0d.nip(:,end-1);
zerod.ne0 = profil0d.nep(:,1);
zerod.ni0 = profil0d.nip(:,1);
option = z0dinput.option;
option.sol_model ='2_points';


% reservation des sorties
tebord_tab = NaN * ones(nbt,nbt+2);
nelim_tab  = NaN * ones(nbt,nbt+2);
telim_tab  = NaN * ones(nbt,nbt+2);
plim_tab   = NaN * ones(nbt,nbt+2);
		
for k=1:nbt+2
	zerod.pin(:) = pin_list(k);
	[tebord,nelim,telim,qpl,err,nb,indbad,fmom,delta_qpl_rad,delta_qpl_neutral,qpl_in_loc,pl, ...
          	zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond] = ...
          	z0convergence_2points_dic(option,cons,geo,zerod,profil0d);
	tebord_tab(:,k) = tebord;
	nelim_tab(:,k)  = nelim;
	telim_tab(:,k)  = telim;
	plim_tab(:,k)   = pl;
end 
indbad = find(all(plim_tab <= 1,1));
tebord_tab(:,indbad) = [];
nelim_tab(:,indbad)  = [];
telim_tab(:,indbad)  = [];
plim_tab(:,indbad)   = [];
pin_list(indbad)     = [];
figure;
subplot(2,2,1);
mesh(nlcfs_list,pin_list,nelim_tab');
subplot(2,2,2);
mesh(nlcfs_list,pin_list,telim_tab');
subplot(2,2,3);
mesh(nlcfs_list,pin_list,tebord_tab');
subplot(2,2,4);
mesh(nlcfs_list,pin_list,plim_tab');
drawnow

% lecture des donnees des sondes de Langmuir
% target data : Langmuir probes outer in iner bord
liste        = {};
liste{end+1}   = 'ppf/@shot/Y4PO/NES?uid=KY4D';    % Ti0
liste{end+1}   = 'ppf/@shot/Y4PO/TES?uid=KY4D';      % profil de Ti
liste{end+1}   = 'ppf/@shot/Y4PI/NES?uid=KY4D';    % Ti0
liste{end+1}   = 'ppf/@shot/Y4PI/TES?uid=KY4D';      % profil de Ti
target         = cgcgetjet(post.z0dinput.shot,liste,'','');
% outbord data
telang = sgolayfilt(target.ppf.Y4PO.TES.data,1,7);

% etalonage telang from C. Guillemaut et al NF 54 (2014) 093012
%pp = polyfit([18.03 10.4],[19.32 4.7],1) 
%telang = min(telang,max(0.1,polyval(pp,telang)));

nelang = sgolayfilt(target.ppf.Y4PO.NES.data,1,7);
% creation des tableaux de parametres
pin_tab = ones(nbt,1) * pin_list;
nlcfs_tab = nlcfs_list * ones(size(pin_list));
% interpolation
nebord_lang_ident    = griddata(nelim_tab,telim_tab,nlcfs_tab,nelang,telang,'linear');
%ind_NaN = find(~isfinite(nebord_lang_ident));
%nebord_lang_ident(ind_NaN)   = griddata(nelim_tab,telim_tab,nlcfs_tab,nelang(ind_NaN),telang(ind_NaN),'nearest');

pout_lcfs_lang_ident = griddata(nelim_tab,telim_tab,plim_tab,nelang,telang,'linear');
%ind_NaN = find(~isfinite(pout_lcfs_lang_ident));
%pout_lcfs_lang_ident(ind_NaN)   = griddata(nelim_tab,telim_tab,plim_tab,nelang(ind_NaN),telang(ind_NaN),'nearest');

pin_lang_ident       = griddata(nelim_tab,telim_tab,pin_tab,nelang,telang,'linear');
%ind_NaN = find(~isfinite(pin_lang_ident));
%pin_lang_ident(ind_NaN)   = griddata(nelim_tab,telim_tab,pin_tab,nelang(ind_NaN),telang(ind_NaN),'nearest');

tebord_lang_ident       = griddata(nelim_tab,telim_tab,tebord_tab,nelang,telang,'linear');
%ind_NaN = find(~isfinite(pin_lang_ident));
%pin_lang_ident(ind_NaN)   = griddata(nelim_tab,telim_tab,pin_tab,nelang(ind_NaN),telang(ind_NaN),'nearest');




h = findobj(0,'type','figure','tag','JET_lang1');
if isempty(h)
       h=figure('tag','JET_lang1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'toolbar','figure')
k    = 6;
subplot(k,1,1)
plot(target.ppf.Y4PO.TES.t,sgolayfilt(target.ppf.Y4PO.TES.data,1,7));
hold on
plot(post.zerod.temps,post.zerod.telim,'k');
ylabel('T_{e,target} (eV)');
title('color = Langmuir/identification, black = METIS');

subplot(k,1,2)
plot(target.ppf.Y4PO.NES.t,sgolayfilt(target.ppf.Y4PO.NES.data,1,7) ./ 1e19);
hold on
plot(post.zerod.temps,post.zerod.nelim./ 1e19,'k');
ylabel('N_{e,target} (1e19 m^{-3})');
set(gca,'ylim',[0,max(target.ppf.Y4PO.NES.data(:))/1e19]);

subplot(k,1,3)
hold on
plot(post.zerod.temps,post.zerod.nebord./ 1e19,'k');
ylabel('N_{e,LCFS} (1e19 m^{-3})');
plot(target.ppf.Y4PO.NES.t,nebord_lang_ident./ 1e19);

subplot(k,1,4)
hold on
plot(post.zerod.temps,post.zerod.tebord ,'k');
ylabel('T_{e,LCFS} (eV)');
plot(target.ppf.Y4PO.NES.t,tebord_lang_ident);


subplot(k,1,5)
hold on
plot(post.zerod.temps,post.zerod.plim ./ 1e6,'k');
ylabel('P_{out,LCFS} (MW)');
plot(target.ppf.Y4PO.NES.t,pout_lcfs_lang_ident ./ 1e6);

subplot(k,1,6)
hold on
plot(post.zerod.temps,post.zerod.pin ./ 1e6,'k');
ylabel('P_{in} (MW)');
plot(target.ppf.Y4PO.NES.t,pin_lang_ident ./ 1e6);


xlabel('time (s)');
%z0loglin(gca);
joint_axes(h,k);
drawnow
edition2


