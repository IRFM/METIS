[info,ref_model] = imas_call_metis_lh_model;
option = info.valeur;
imas_ref_in = ref_model;
%imas_ref_in.shot = 54719;
imas_ref_in.shot = 55564;
imas_ref_in.run  = 22;
imas_ref_in.tokamak  = 'west';
imas_ref_in.occurrence.equilibrium  = 0;
imas_ref_out = imas_ref_in;
imas_ref_out.run = 21;
option.gaz      = 2;
option.etalh    = 0.81;
option.wlh      = 0.58;
option.npar_neg = -4;
option.fupshift = 1;
option.npar0    = 1.9;
option.freqlh    = 3.7;
core_sources    = imas_call_metis_lh_model(imas_ref_in,imas_ref_out,option);

% search for existing sources
index_sources = [];
for k=1:length(core_sources.source)
    if ~isempty(core_sources.source{k}.identifier.index)
        index_sources(end+1) = core_sources.source{k}.identifier.index;
    end
end
if isempty(index_sources)
   error('No lH source')
else
   index_lh = find(index_sources == 4);
   index_lh_prev = find(index_sources == 99999);
end
  
lhs = core_sources.source{index_lh};
lhs_prev = core_sources.source{index_lh_prev};


time = core_sources.time;
nbt  = length(time);
nbx  = length(lhs.profiles_1d{1}.grid.rho_tor_norm);
%
rho      = NaN * ones(nbt,nbx);
rho_pol  = NaN * ones(nbt,nbx);
psi      = NaN * ones(nbt,nbx);
plh      = NaN * ones(nbt,nbx);
jlh      = NaN * ones(nbt,nbx);
rot      = NaN * ones(nbt,nbx);
ilh      = NaN * ones(nbt,1);
ptot     = NaN * ones(nbt,1);
tor      = NaN * ones(nbt,1);
%
rho_prev      = NaN * ones(nbt,nbx);
rho_pol_prev  = NaN * ones(nbt,nbx);
psi_prev      = NaN * ones(nbt,nbx);
plh_prev      = NaN * ones(nbt,nbx);
jlh_prev      = NaN * ones(nbt,nbx);
rot_prev      = NaN * ones(nbt,nbx);
ilh_prev      = NaN * ones(nbt,1);
ptot_prev     = NaN * ones(nbt,1);
tor_prev      = NaN * ones(nbt,1);

for k=1:length(time)
    rho(k,:)      = lhs.profiles_1d{k}.grid.rho_tor_norm;
    rho_pol(k,:)  = lhs.profiles_1d{k}.grid.rho_pol_norm;
    psi(k,:)      = lhs.profiles_1d{k}.grid.psi;
    plh(k,:)      = lhs.profiles_1d{k}.electrons.energy;
    jlh(k,:)      = lhs.profiles_1d{k}.j_parallel;
    rot(k,:)      = lhs.profiles_1d{k}.momentum_tor;
    ilh(k)        = lhs.global_quantities{k}.current_parallel;
    ptot(k)       = lhs.global_quantities{k}.power;
    tor(k)        = lhs.global_quantities{k}.torque_tor;

    rho_prev(k,:)      = lhs_prev.profiles_1d{k}.grid.rho_tor_norm;
    psi_prev(k,:)      = lhs_prev.profiles_1d{k}.grid.psi;
    if ~isempty(lhs_prev.profiles_1d{k}.grid.rho_pol_norm)
        rho_pol_prev(k,:)  = lhs_prev.profiles_1d{k}.grid.rho_pol_norm;
    else
        rho_pol_prev(k,:)  = sqrt((psi_prev(k,:) - psi_prev(k,1)) ./ (psi_prev(k,1) - psi_prev(k,end)));
    end
    plh_prev(k,:)      = lhs_prev.profiles_1d{k}.electrons.energy;
    jlh_prev(k,:)      = lhs_prev.profiles_1d{k}.j_parallel;
    rot_prev(k,:)      = lhs_prev.profiles_1d{k}.momentum_tor;
    ilh_prev(k)        = lhs_prev.global_quantities{k}.current_parallel;
    ptot_prev(k)       = lhs_prev.global_quantities{k}.power;
    tor_prev(k)        = lhs_prev.global_quantities{k}.torque_tor;
 
end

figure
subplot(2,2,1)
plot(rho',psi','r',rho_prev',psi_prev','.b')
title(sprintf('%s@%d #%d',imas_ref_in.tokamak,imas_ref_in.shot ,imas_ref_in.run));
xlabel('rho_tor_norm')
ylabel('psi (Wb)')
subplot(2,2,2)
plot(rho',plh','r',rho_prev',plh_prev','.b')
xlabel('rho_{tor,norm}')
ylabel('plh (W/m^{-3})')
title('red = IMAS actor & blue = IMAS database');
subplot(2,2,3)
plot(rho',jlh','r',rho_prev',jlh_prev','.b')
xlabel('rho_{tor,norm}')
ylabel('jlh (A/m^{-2})')
subplot(2,2,4)
plot(rho',rot','r',rho_prev',rot_prev','.b')
xlabel('rho_{tor,norm}')
ylabel('momentum  (kg m^{-1} s^{-2})')

figure
subplot(3,1,1)
plot(time,ptot/1e6,'r',time,ptot_prev/1e6,'.b');
title(sprintf('%s@%d #%d',imas_ref_in.tokamak,imas_ref_in.shot ,imas_ref_in.run));
%xlabel('time (s)');
ylabel('P_{LH} (MW)');
legend('IMAS actor','IMAS database')

subplot(3,1,2)
plot(time,ilh/1e3,'r',time,ilh_prev/1e3,'.b');
%xlabel('time (s)');
ylabel('I_{LH} (kA)');

subplot(3,1,3)
plot(time,tor,'r',time,tor_prev,'.b');
xlabel('time (s)');
ylabel('Torque_{LH} (kg m ^2 s^{-2})');









