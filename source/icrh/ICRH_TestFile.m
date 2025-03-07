load METIS_ICRH_WEST1
% load METIS_ICRH_LH_WEST1
% 0d 
% pabs_ME = post.zerod.picrh_th; 
% Pe_ME = post.zerod.pel_icrh;
% Pi_ME = post.zerod.pion_icrh;

%%% post.z0dinput.option.gaz = 2 correspond to D Majority
% % post.z0dinput.option.mino = H  for minority

tune.ichcd.active = 1;
tune.ichcd.uindices = 2;
% ICRH_config
ipa.active = tune.ichcd.active;
ipa.uindices = tune.ichcd.uindices;
ipa.check = false;

%% METIS Test Find all in  post.z0dinput.option
ipa.eps_Xi = 0.03;
ipa.upshift = 1;
ipa.p = 1;
ipa.f = post.z0dinput.option.freq*1e6; % Hz 5.5e10; 5.5e7;
ipa.fz = 0.4; 0.5;
ipa.rate_ni =  post.z0dinput.option.cmin; % nm/nM
ipa.n = post.z0dinput.option.nphi;
% ipa.widthtau = 1;2; 1.22; %*delR == post.z0dinput.option.icrh_width
ipa.P_lossTau = 0;  
ipa.Nu = 1000;
ipa.Nz = 100;
% In model.atom
ipa.Zm = 1;% H in METIS; 
ipa.ZM = 1;% Deuteron 1 % charge of Majority with electron
ipa.Am = 1;% number of mass of minority
ipa.AM = 2;% number of mass of Majority

ichcd_params = ipa;

g = 0;

model.geom.rhogauss = post.profil0d.xli';
rho = model.geom.rhogauss;
Pe = zeros(size(post.profil0d.qjli'));
PM = Pe;
Pm = Pe; 
pabs = Pe;
Wpert = Pe;
Wpar = Pe;
R = zeros(1,length(Pe));
next = [];

for it = 1:296 % from it = 85
U = post.z0dinput.cons.picrh(it);

if U >  1e2 % for icrh
    
    q = post.profil0d.qjli(it,:)';
    te = post.profil0d.tep(it,:)';
    ne = post.profil0d.nep(it,:)';
    ti = post.profil0d.tip(it,:)';
    ni = post.profil0d.nip(it,:)';
    
%     figure(1);
%     subplot(311)
%     plot(rho,te,rho,ti)
%     legend ('te','ti')
%     subplot(312)
%     plot(rho,ne,rho,ni)
%     legend ('ne','ni')
%     subplot(313)
%     plot(rho,q)
%     legend ('q')
    
    model.equi.R0 = post.z0dinput.geo.R(it);
    model.equi.epsilon = post.z0dinput.geo.a(it)/model.equi.R0;
    model.equi.B0 = post.z0dinput.geo.b0(it);
    model.equi.kappa = post.z0dinput.geo.K(it);
    ichcd_params.picrh = (post.profil0d.picrh(it,:))';
    ichcd_params.Ploss =  post.z0dinput.cons.picrh(it) - post.zerod.picrh(it);
    last = next;
    [Pe(:,it),PM(:,it),pabs(:,it),R(it),Wpert(:,it),Wpar(:,it),next] = ichcd_colli_test(q,te,ne,ti,ni,U,model,ichcd_params,g,last);
%     plot([Pe(:,it), pabs(:,it),(post.profil0d.picrh_ion(it,:))']) % ~sum/5/ipa.widthtau
%     plot([pabs(:,it),(post.profil0d.picrh(it,:))'])
%    figure(2);
%    plot([(post.profil0d.picrh_ion(it,:))',(post.profil0d.picrh(it,:))'])
%     legend('Pe','pabs','picrh_{METIS}')
%     hold on

end
end
%% 
% % Wpar_METIS == post.zerod.esup_icrh;
figure(9);plot(max(Wpar)); hold on;
plot(max(Wpert/2),'c'); 
plot(post.zerod.esup_icrh,'r');
plot(post.zerod.picrh_th/1e2,'k');
legend('Wpar','Wpert/2','esup_{icrh}METIS','picrh/100')
xlabel('t')
% compare to METIS
% load METIS_ICRH_WEST1
% trapz(rho,2*pi*R*pabs)%%== post.zerod.picrh_th(it);
% trapz(rho,2*pi*R*PM) %%== post.zerod.pion_icrh(it);
% trapz(rho,2*pi*R*Pe) == post.zerod.pel_icrh(it);
