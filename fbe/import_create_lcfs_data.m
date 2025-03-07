%Import CREATE LCFS data
[fname,pathtofile] = uigetfile('*.mat','Select a CREATE LCFS data file');
if isnumeric(fname)
    % action canceled
    return
end
filename = fullfile(pathtofile,fname);
ref = load(filename);
nbp_max = 0;
clear fbe_lcfs 
if isfield(ref,'LCFS')
    for k=1:length(ref.LCFS)
        fbe_lcfs{k}.time = ref.LCFS(k).time;
        fbe_lcfs{k}.ip   = ref.LCFS(k).ip .* 1e6;
        fbe_lcfs{k}.R    = ref.LCFS(k).r_lim(:);
        fbe_lcfs{k}.Z    = ref.LCFS(k).z_lim(:);
        nbp_max = max(nbp_max,length(ref.LCFS(k).r_lim));
    end
elseif isfield(ref,'SH')
    for k=1:length(ref.SH)
        fbe_lcfs{k}.time = ref.SH(k).time;
        fbe_lcfs{k}.ip   = ref.SH(k).Ip;
        fbe_lcfs{k}.R    = ref.SH(k).rc;
        fbe_lcfs{k}.Z    = ref.SH(k).zc;
        nbp_max = max(nbp_max,length(ref.SH(k).rc));
    end
else
    for k=1:length(ref.Shape_sequence)
        fbe_lcfs{k}.time = ref.Shape_sequence(k).time;
        fbe_lcfs{k}.ip   = ref.Shape_sequence(k).Ip;
        fbe_lcfs{k}.R    = ref.Shape_sequence(k).rc;
        fbe_lcfs{k}.Z    = ref.Shape_sequence(k).zc;
        nbp_max = max(nbp_max,length(ref.Shape_sequence(k).rc));
    end
end
km = menu('Choose a time rsampling method:','kept time slices from CREATE','kept time slices from METIS','resampling on plasma current','no resmapling');
list = {'input','time','ip','none'};
z0dinput = z0separatrix_fromfbetometis(z0dinput,fbe_lcfs,list{km},[],1);

% make matrix for figure
time_fig = NaN * ones(length(fbe_lcfs),1);
R_fig    = NaN * ones(length(fbe_lcfs),nbp_max);
Z_fig    = NaN * ones(length(fbe_lcfs),nbp_max);
for k=1:length(fbe_lcfs)
    time_fig(k) = fbe_lcfs{k}.time;
    ip_fig(k)   = fbe_lcfs{k}.ip;
    R_fig(k,1:length(fbe_lcfs{k}.R))    = fbe_lcfs{k}.R;
    Z_fig(k,1:length(fbe_lcfs{k}.Z))    = fbe_lcfs{k}.Z;
end
indbad = find(R_fig <= 0);
if ~isempty(indbad)
    R_fig(indbad) = NaN;
    Z_fig(indbad) = NaN;
end

% figure for control
hold on
zplotprof(gca,time_fig,R_fig,Z_fig,'color','k','marker','.','linestyle','none');
legend('METIS moments','METIS (R,Z)','CREATE (R,Z)');
if ~isempty(strfind(upper(z0dinput.machine),'JT-60SA'))
    load('jt60sawall.mat');
    plot(wall.data(:,1),wall.data(:,2),'k');
end
axis('square')
axis('equal')
xlabel('R (m)');
ylabel('Z (m)');
title(sprintf('importation of LCFS from CREATE (file = %s) in METIS with resampling method ''%s'' for %s',fname,list{km},z0dinput.machine));

% figure current
h = findobj(0,'type','figure','tag','sepacur');
if isempty(h)
h=figure('tag','sepacur');
else
figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
plot(z0dinput.cons.temps,z0dinput.cons.ip./1e6,'r',time_fig,ip_fig./1e6,'.k');
xlabel('Time (s)');
ylabel('Ip (MA)');
title(sprintf('importation of LCFS from CREATE (file = %s) in METIS with resampling method ''%s'' for %s',fname,list{km},z0dinput.machine));

