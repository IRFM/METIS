% script pour le test du model complet de SOL
zerod_loc 	 = post.zerod;
profil0d = post.profil0d;
cons     = post.z0dinput.cons;
geo      = post.z0dinput.geo;
exp0d    = post.z0dinput.exp0d;
z0plotsc;
drawnow	
[tend_meas,void] = ginput(1);
drawnow
prompt={'minimum','maximum','number of points'};
name='Range of density : input power of 10 given the minimum and the maximum density at the LCFS';
numlines=1;
defaultanswer={'18','21','301'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
        nlcfs_list = logspace(max(13,min(str2num(answer{1}),str2num(answer{2}))), ...
                              min(30,max(str2num(answer{1}),str2num(answer{2}))), ...
                              min(3001,max(11,str2num(answer{3}))))';
	nbt     = length(nlcfs_list);
else
	nbt     = 301;
	nlcfs_list = logspace(18,21,nbt)';
end


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

noms = fieldnames(zerod_loc);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(zerod_loc,nomc);
	if length(val) == length(times)
		val  = interp1(times,val,tend_meas,'nearest');
		zerod_loc = setfield(zerod_loc,nomc,vt*val);
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
nebord_mem   = zerod_loc.nebord;
% reglage ici
zerod_loc.nebord = nlcfs_list;
%
zerod_loc.nibord = zerod_loc.nibord .* zerod_loc.nebord  ./ nebord_mem;
nep_mem      = profil0d.nep;
profil0d.nep = profil0d.nep - profil0d.nep(:,end) * ve + zerod_loc.nebord * ve;
profil0d.nip = profil0d.nep ./ nep_mem  .* profil0d.nip;
profil0d.n1p = profil0d.nep ./ nep_mem  .* profil0d.n1p;
profil0d.nhep = profil0d.nep ./ nep_mem  .* profil0d.nhep;
profil0d.nzp = profil0d.nep ./ nep_mem  .* profil0d.nzp;
profil0d.nwp = profil0d.nep ./ nep_mem  .* profil0d.nwp;
vol          = trapz(profil0d.xli,profil0d.vpr,2);
zerod_loc.nem    = trapz(profil0d.xli,profil0d.vpr .* profil0d.nep,2) ./ vol;
zerod_loc.nbar   = trapz(profil0d.xli,profil0d.nep,2);
zerod_loc.nhem   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nhep,2) ./ vol;
zerod_loc.nim   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nip,2) ./ vol;
m1m_mem     = zerod_loc.n1m;
zerod_loc.n1m   = trapz(profil0d.xli,profil0d.vpr .* profil0d.n1p,2) ./ vol;
zerod_loc.nDm   = zerod_loc.nDm .* zerod_loc.n1m ./ m1m_mem;
zerod_loc.nTm   = zerod_loc.nTm .* zerod_loc.n1m ./ m1m_mem;
zerod_loc.nimpm   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nzp,2) ./ vol;
zerod_loc.nwm   = trapz(profil0d.xli,profil0d.vpr .* profil0d.nwp,2) ./ vol;
zerod_loc.neped = profil0d.nep(:,end-1);
zerod_loc.niped = profil0d.nip(:,end-1);
zerod_loc.ne0 = profil0d.nep(:,1);
zerod_loc.ni0 = profil0d.nip(:,1);
option = z0dinput.option;
option.sol_model ='2_points';

[tebord,nelim,telim,qpl,err,nb,indbad,fmom,delta_qpl_rad,delta_qpl_neutral,qpl_in_loc,pl, ...
          zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond] = ...
          z0convergence_2points_dic(option,cons,geo,zerod_loc,profil0d);

h = findobj(0,'type','figure','tag','z02p_test');
if isempty(h)
       h=figure('tag','z02p_test');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'toolbar','figure')
k    = 6;
subplot(k,1,1)
semilogy(zerod_loc.nebord/1e19,zerod_loc.nebord/1e19,zerod_loc.nebord/1e19,nelim/1e19);
ylabel('10^{19} m^{-3}');
legend('n_{e,LCFS}','n_{e,t}')
z0loglin(gca);
title(sprintf('METIS : %s@%d/ Two Points model (ramp of density, use GUI set of parameters)',post.z0dinput.machine,post.z0dinput.shot));
subplot(k,1,2)
semilogy(zerod_loc.nebord/1e19,tebord,zerod_loc.nebord/1e19,telim);
ylabel('eV');
legend('T_{e,LCFS}','T_{e,t}')
z0loglin(gca);
subplot(k,1,3)
semilogy(zerod_loc.nebord/1e19,qpl./1e6,zerod_loc.nebord/1e19,qpl_in_loc./1e6,zerod_loc.nebord/1e19,delta_qpl_rad./1e6,zerod_loc.nebord/1e19,delta_qpl_neutral./1e6);
ylabel('MW/m^2');
legend('q_{//,t}','q_{//,LCFS}','\delta(q_{//,rad,div })','\delta(q_{//,neutral,div })');
z0loglin(gca);
subplot(k,1,4)
semilogy(zerod_loc.nebord/1e19,fmom,zerod_loc.nebord/1e19,mach_target,zerod_loc.nebord/1e19,fcond);
legend('f_{mom}','M_t','f_{cond} * f_{pe}');
z0loglin(gca);
subplot(k,1,5)
semilogy(zerod_loc.nebord/1e19,zeff_div,zerod_loc.nebord/1e19,gamma);
legend('Z_{eff} div','\gamma');
z0loglin(gca);
subplot(k,1,6)
semilogy(zerod_loc.nebord/1e19,prad_loc./1e6,zerod_loc.nebord/1e19,pradsol_loc./1e6);
legend('Prad core (line)','Prad SOL (line + brem)');
ylabel('MW');
xlabel('LCFS density (10^{19} m^{-3})');
z0loglin(gca);
joint_axes(h,k);
drawnow
edition2



h = findobj(0,'type','figure','tag','z02p_test_extract');
if isempty(h)
       h=figure('tag','z02p_test_extract');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'toolbar','figure')
k    = 2;
subplot(k,1,1)
semilogy(zerod_loc.nebord/1e19,zerod_loc.nebord/1e19,zerod_loc.nebord/1e19,nelim/1e19);
ylabel('10^{19} m^{-3}');
legend('n_{e,LCFS}','n_{e,t}')
z0loglin(gca);
title(sprintf('METIS : %s@%d/ Two Points model (ramp of density, use GUI set of parameters)',post.z0dinput.machine,post.z0dinput.shot));
subplot(k,1,2)
semilogy(zerod_loc.nebord/1e19,tebord,zerod_loc.nebord/1e19,telim);
ylabel('eV');
legend('T_{e,LCFS}','T_{e,t}')

xlabel('LCFS density (10^{19} m^{-3})');
z0loglin(gca);
joint_axes(h,k);
drawnow
edition2

