


% charge fichier CXS prive
shot = input('CXS shot number ? ');
% ancienne version par fichier
%fname = ['load /home/sccp/gcbd/fenzi/CXRS/pourThierry/TICXS_',int2str(shot)];
%eval(fname)

% avec tsbase
[Ticxs,tcxs,pos]=tsbase(shot,'GTICXS'); % shot+0.1 pour occurrence 1
[erTicxs]=tsbase(shot,'GERTICXS');


disp(['CXS data for time t = ',num2str(tcxs'),' s'])
if length(tcxs > 1)
   icxs = input('Which CXS time slice do you want to see (1,2,3, ...) ? ');
   %pos = pos(:,icxs);   % pour fichiers
   %Ticxs = Ticxs(:,icxs);
   %erTicxs = erTicxs(:,icxs);
   pos = pos(icxs,:)';    % pour TSBASE
   Ticxs = Ticxs(icxs,:)';
   erTicxs = erTicxs(icxs,:)';
   
   tcxs = tcxs(icxs);
end

time = input('CRONOS time for comparison ? (s) ');

icronos = iround(data.gene.temps,time);
rhoCXS_LFS = interp1(squeeze(double(data.equi.R(icronos,:,1))),squeeze(double(data.equi.rhoRZ(icronos,:)))./data.equi.rhomax(icronos),pos(pos>double(data.equi.R(icronos,1,1))));  % low field side
rhoCXS_HFS = interp1(double(squeeze(data.equi.R(icronos,:,floor(size(data.equi.R,3)/2)))),squeeze(double(data.equi.rhoRZ(icronos,:)))./data.equi.rhomax(icronos),pos(pos<double(data.equi.R(icronos,1,1))));  % high field side
rhoCXS = [rhoCXS_HFS' rhoCXS_LFS'];


figure
plot(param.gene.x,data.prof.ti(icronos,:),'r')
hold on
errorbar(rhoCXS,Ticxs,erTicxs)
legend('Cronos','CXS')
xlabel('\rho')
ylabel('Ti (eV)')
title(['Shot ',int2str(shot),' , t = ',num2str(tcxs),' s'])
