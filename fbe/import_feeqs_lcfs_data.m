%Import CREATE LCFS data
clear fbe_lcfs
[fname,pathtofile] = uigetfile('*.mat','Select a FEEQS.M LCFS data file');
filename = fullfile(pathtofile,fname);
load(filename);
nbp_max = 0;
for k=1:length(fbe_lcfs)
    nbp_max = max(nbp_max,length(fbe_lcfs{k}.R));
end
km = menu('Choose a time rsampling method:','kept time slices from FEEQS.M','kept time slices from METIS','resampling on plasma current','no resmapling');
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
title(sprintf('importation of LCFS from FEEQS.M (file = %s) in METIS with resampling method ''%s'' for %s',fname,list{km},z0dinput.machine));

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
title(sprintf('importation of LCFS from FEEQS.M (file = %s) in METIS with resampling method ''%s'' for %s',fname,list{km},z0dinput.machine));
