function z0plot_shot_digest(post,tprofi,ts_list)


% selection of the time slice
if (nargin == 1)
    evalin('base','z0plot_reference_post;');
    subplot(3,1,1);
    title('choose 1 time slice for profiles visualization')
    drawnow
    [tprofi,void] = ginput(1);
    subplot(3,1,1);
    title('choose 3 time slices for safety factor evolution visualization')
    drawnow
    [ts_list,void] = ginput(3);
    close(gcf);
    drawnow
end
ts_list = sort(ts_list);


% choose between mode reactor and non nuclear experiment
disp(' ')
switch post.z0dinput.option.gaz
    case 3
        % DEMO or ITER
        z0plot_shot_digest_dt(post,tprofi,ts_list)
        
    case 5
        % D-He3
        z0plot_shot_digest_dhe3(post,tprofi,ts_list)
        
    case 11
        % p-B11
        z0plot_shot_digest_pb11(post,tprofi,ts_list)
        
    otherwise
        % other tokamak
        z0plot_shot_digest_other(post,tprofi,ts_list)
end

function z0plot_shot_digest_other(post,tprofi,ts_list)

disp('===================================================')
disp('Non nuclear tokamak')
%combines wwmetis and figMetis

%tmax=120; tprofi=90;
tmax=max(post.profil0d.temps);
tmin_round = fix(min(post.profil0d.temps));


nfig=0;lls='-';posfig=0;

a=post.z0dinput.geo.a;
R0=post.z0dinput.geo.R;
eps=a./R0;
k=post.z0dinput.geo.K;
delta=post.z0dinput.geo.d;
B=post.z0dinput.geo.b0;
%betapol=post.zerod.betap;
betapol=post.zerod.betaptot;
rho=post.profil0d.rmx./(max(post.profil0d.rmx,[],2)*ones(1,21));
t=post.z0dinput.cons.temps;
ttt=post.profil0d.temps;
nbar=post.z0dinput.cons.nbar;
tau=post.zerod.taue;
tauh=post.zerod.tauh;
tauhe=post.zerod.tauhe;
wth=post.zerod.wth;
V=post.zerod.vp;
S=post.zerod.sp;
P_tot=post.zerod.pin;  %0-d quantity
P_loss=post.zerod.ploss;  %0-d quantity
P_fus=post.zerod.pfus;  %0-d quantity, puissance alpha
P_nbi=post.zerod.pnbi;  %0-d quantity
P_ic=post.zerod.picrh;  %0-d quantity
P_ec=post.zerod.pecrh;  %0-d quantity
P_lh=post.zerod.plh;  %0-d quantity
P_aux=P_tot-P_fus;  %0-d quantity
Prad=post.zerod.prad;  %0-d quantity
Pradsol=post.zerod.pradsol;  %0-d quantity
Pbrem=post.zerod.pbrem;  %0-d quantity
Pcyclo=post.zerod.pcyclo;  %0-d quantity
Pradtot=Prad+Pbrem+Pcyclo;
Pradcore=Prad+Pbrem;
Q=5.0286*P_fus./(P_tot-P_fus);
Zeff=post.zerod.zeff;

q=post.profil0d.qjli;
j=post.profil0d.jli;
jbs=post.profil0d.jboot + post.profil0d.jfus;
jpar=post.profil0d.jeff;
Bpol=post.profil0d.bpol;
epar=post.profil0d.epar;
Te=post.profil0d.tep/1000;
Ti=post.profil0d.tip/1000;
xie=post.profil0d.xie;
xii=post.profil0d.xii;
ne=post.profil0d.nep;
ni=post.profil0d.nip;
Pfce=post.profil0d.pecrh;
Plhp=post.profil0d.plh;
Pnbi_i=post.profil0d.pnbi_ion;
Pnbi=post.profil0d.pnbi;
Pnbi_e=Pnbi-Pnbi_i;
Pfci_i=post.profil0d.picrh_ion;
Pfci=post.profil0d.picrh + post.profil0d.pfweh;
Pfci_e=Pfci-Pfci_i;
Pfus_i=post.profil0d.pfus_ion;
Pfus=post.profil0d.pfus;
Pfus_e=Pfus-Pfus_i;
P_i=post.profil0d.source_ion;
P_e=post.profil0d.source_el;
jfce=post.profil0d.jeccd;
jlh=post.profil0d.jlh;
jnbi=post.profil0d.jnbicd;
jfci=post.profil0d.jfwcd;
Pohm=post.profil0d.pohm;
q_ei=post.profil0d.qei;
Psyn=post.profil0d.pcyclo;
P_rad=post.profil0d.prad;
P_brem=post.profil0d.pbrem;

Ip=post.zerod.ip;
Ipar=post.zerod.ipar;
li=post.zerod.li;
Vmes=post.zerod.vmes;
Vloop=post.zerod.vloop;
Ioh=post.zerod.iohm;
Ilh=post.zerod.ilh;
Ifce=post.zerod.ieccd;
Iboot=post.zerod.iboot;
Ibs=post.zerod.iboot;
Ini=post.zerod.ini;
Inbi=post.zerod.inbicd;
psi=post.profil0d.psi;
ptot=post.profil0d.ptot; %total pressure

nHe=post.profil0d.nhep;
betaN=post.zerod.betan*100;
beta=1e6*Ip.*betaN./(a.*B);
Ptot=Pnbi+Pfci+Pfus+Pohm+Pfce+Plhp;
f_b_s=Iboot./Ipar;
f_n_i=Ini./Ipar;
nG=(Ip./1e6)/pi./a.^2;
f_G=(nbar/1e20)./nG;
ka=V./(2*pi^2*R0.*a.^2);
M=post.z0dinput.option.gaz;
tauIPB=0.0562*(Ip/1e6).^0.93.*B.^0.15.*(P_loss/1e6).^(-0.69).*(nbar/1e19).^0.41 ...
    .*M^0.19.*R0.^1.97.*eps.^0.58.*ka.^0.78;
H=tau./tauh;
tauIPB=tauh;

if exist('post.profil0d.Rsepa')==1
    Rse=post.profil0d.Rsepa;
    Zse=post.profil0d.Zsepa;
else
end
q95=post.zerod.q95;
nem=post.zerod.nem;
nhem=post.zerod.nhem;
nimpm=post.zerod.nimpm;
negr=post.zerod.negr;
tem=post.zerod.tem/1e3;
tite=post.zerod.tite;
Pnbi1=real(post.z0dinput.cons.pnbi)/1e6;
Pnbi2=imag(post.z0dinput.cons.pnbi)/1e6;
Pnbi1_abs=real(post.zerod.pnbi)/1e6;
Pnbi2_abs=imag(post.zerod.pnbi)/1e6;
Pecrh=post.z0dinput.cons.pecrh/1e6;
Picrh=post.z0dinput.cons.picrh/1e6;
Plh=post.z0dinput.cons.plh/1e6;
tauj=post.zerod.tauj;
Sn=post.zerod.ndd;
Teped=post.zerod.teped/1e3;
neped=post.zerod.neped/1e19;
rimp=post.z0dinput.option.rimp;

eval('nne=ne;');eval('pptot=ptot;');eval('jj=j;');eval('bbeta=beta;');

% GLOBAL QUANTITIES


dpos=-20; posfig=0; pos1=256; pos2=247; pos3=364; pos4=425;
%psize1=0.7; psize10=0.8; psize2=1.1; psize3=1.3; psize4=1.1; psize5=1.05; psize6=0.97;
psize1=1; psize10=1; psize2=1; psize3=1; psize4=1; psize5=1; psize6=1;

tit=['t = ',num2str(tprofi),' s'];
elon=k;

nfig=nfig+1;
figuren(nfig); set(gcf,'name','Ip,ne0,Q,P vs t','Position',[pos1+posfig pos2+posfig pos3*psize3 pos4*psize3]);
subplot(211); set(gca,'fontsize',fix(14*psize3));
plot(t,Ip/1e6,ttt,nne(:,1)/1e19,'linewidth',2,'linestyle',lls)
xlabel('time (s)')
axis([tmin_round  tmax 0 Inf])
legend('Ip (MA)','ne(0) (10^{19} m^{-3})','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize3));
%legend('Pnbi','Pec','Location','Best')
%plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,'linewidth',2,'linestyle',lls)
%legend('Pnnbi','Ppnbi','Pec','Palpha')
%plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
if post.z0dinput.option.lhmode == 5
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,(P_ec + P_lh)/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Pic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,(P_ec + P_lh)/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Location','Best')
    end
elseif any(P_lh > 1e3)
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_lh/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Plh','Pic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_lh/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Plh','Location','Best')
    end
else
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','P_ic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Location','Best')
    end
end
axis([tmin_round  tmax 0 Inf])
xlabel('time (s)')
ylabel('(MW)')

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','fbs,fni,fG,betaN,li,H,Zeff vs t','Position',[pos1+posfig pos2+posfig pos3*psize2 pos4*psize2]);
subplot(211); set(gca,'fontsize',fix(14*psize2));
plot(t,f_b_s,t,f_n_i,t,f_G,'linewidth',2,'linestyle',lls)
xlabel('time (s)')
axis([tmin_round  tmax 0 1.2])
legend('f_b_s','f_n_i','f_G','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize2));
HH=smooth(H,4);
plot(t,betaN,t,4*li,t,HH,t,Zeff,'linewidth',2,'linestyle',lls)
hold on;plot([0 tmax],[1 1],'k--');hold off;
axis([tmin_round  tmax 0 5.5])
legend('betaN','4 li','H','Zeff','Location','Best')
xlabel('time (s)')

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','gCD,Vl,Prad,Pbrem,Psyn vs t','Position',[pos1+posfig pos2+posfig pos3*psize2 pos4*psize2]);
subplot(211); set(gca,'fontsize',fix(14*psize2));
gCD=(Ini-Ibs)./(real(P_nbi)+imag(P_nbi)+P_lh+P_ec).*nbar.*R0/1e20;
plot(t,gCD,t,smooth(Vloop,4),'linewidth',2,'linestyle',lls)
hold on;plot([0 tmax],[0 0],'k--');hold off;
xlabel('time (s)')
axis([tmin_round  tmax -0.05 Inf])
legend('gammaCD','Vloop (V)','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize2));
plot(t,Pradtot/1e6,t,Pbrem/1e6,t,Pcyclo/1e6,'linewidth',2,'linestyle',lls)
legend('Pradtot','Pbrem','Psyn','Location','Best')
axis([tmin_round  tmax 0 Inf])
xlabel('time (s)')
ylabel('(MW)')

% PROFILES

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','Te,Ti,ne profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
k=iround(ttt,tprofi);
h=plot(rho(k,:),Te(k,:),rho(k,:),Ti(k,:),rho(k,:),nne(k,:)/1e19); set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
legend('Te (keV)','Ti (keV)','ne (10^{19} m^{-3})','Location','Best')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','chi_e,chi_i profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
h=plot(rho(k,:),xie(k,:),rho(k,:),xii(k,:)); set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
legend('chie (m^2/s)','chii (m^2/s)','Location','NorthWest')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','P profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize4 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
l=(1:21);
%h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
%legend('Palpha','Pnbi','Pec','Location','Best')
if any(Plhp(k,l) > 1e3)
    if any(Pfci(k,l) > 1e3)
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6+Plhp(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Pnbi','Pec','Pic','Location','Best');
        else
            h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Plhp(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Pnbi','Pec','Plh','Pic','Location','Best');
        end
    else
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6+Plhp(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Pnbi','Pec','Location','Best');
        else
            h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Plhp(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Pnbi','Pec','Plh','Location','Best');
        end
    end
else
    if any(Pfci(k,l) > 1e3)
        h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('Pnbi','Pec','Pic','Location','Best')
    else
        h=plot(rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('Pnbi','Pec','Location','Best')
    end
end
%h=plot(rho(k,l),Pnbi_i(k,l)/1e6,rho(k,l),Pnbi_e(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
%legend('Pnbi-i','Pnbi-e','Pec')
axis([0 1 0 Inf])
xlabelrho
ylabel('(MW/m^3)')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','j profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize4 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
l=(1:21);
if any(jlh(k,l) > 1e3)
    if any(jfci(k,l) > 1e3)
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),(jfce(k,l)+ jlh(k,l))/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jfwcd','Location','Best')
        else
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jlh(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jlh','jfwcd','Location','Best')
        end
    else
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),(jfce(k,l)+ jlh(k,l))/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','Location','Best')
        else
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jlh(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jlh','Location','Best')
        end
    end
else
    if any(jfci(k,l) > 1e3)
        h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('j','jnbi','jec','jbs','jfwcd','Location','Best')
    else
        h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('j','jnbi','jec','jbs','Location','Best')
    end
end
%h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
ylabel('(MA/m^2)')
%legend('j','jnbi','jec','jbs','Location','Best')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','q profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize5 pos4*psize1*psize6]);
set(gca,'FontSize',fix(18*psize10))
kk=iround(ttt,ts_list);
h=plot(rho(kk,:)',q(kk,:)',[0 1],[1 1],'k--');set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
ylabel('q')
legend(sprintf('t = %0.2f s',ts_list(1)),sprintf('t = %0.2f s',ts_list(2)),sprintf('t = %0.2f s',ts_list(3)),'Location','Best')

qlim = ceil(min(max(post.zerod.qa),2.1 .* trapz(post.zerod.temps,post.zerod.ip .* post.zerod.qa) ./ max(1,trapz(post.zerod.temps,post.zerod.ip))));
posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','q vs t','Position',[pos1+posfig pos2+posfig pos3*psize1*psize5 pos4*psize1*psize6]);
set(gca,'FontSize',fix(18*psize10))
h=plot(ttt,q(:,1),ttt,q(:,11),ttt,q(:,20),[0 tmax],[1 1],'k--');set(h,'linewidth',2,'linestyle',lls);
axis([tmin_round tmax 0 qlim])
xlabel('time (s)')
ylabel('q(rho)')
legend('q(0)','q(0.5)','q(0.95)','Location','Best')

%Poloidal current density (units to be checked)

d2=post.profil0d.df2dpsi(k,l);
F=post.profil0d.fdia(k,l);
dFdpsi=0.5*d2./F;
Bp=post.profil0d.bpol(k,l);
jpol=dFdpsi.*Bp/(4*pi*1e-7);
%figuren(99)
%plot(rho(k,l),jpol)

%  TABLE

t0= tprofi;

ii=iround(t,t0);
iipr=iround(ttt,t0);

if exist('Rse')==1
    
    R=Rse(iipr,:);
    Z=Zse(iipr,:);
    
    imax=max(find(Z==max(Z)));
    imin=max(find(Z==min(Z)));
    Zmax=Z(imax); Zmin=Z(imin);
    
    jmax=max(find(R==max(R)));
    jmin=max(find(R==min(R)));
    Rmax=R(jmax); Rmin=R(jmin);
    
    RZmax=R(imax); RZmin=R(imin);
    
    RR0=(Rmax+Rmin)/2;
    aa=(Rmax-Rmin)/2;
    k=(Zmax-Zmin)/(Rmax-Rmin);
    dup=(RR0-RZmax)/aa;
    dlo=(RR0-RZmin)/aa;
    delta=(dup+dlo)/2;
else
    k=elon(ii);
    delta=delta(ii);
    
end


time=t(ii);
Ip0=Ip(ii)/1e6;
B0=B(ii);
q950=q95(ii);
a0=a(ii);
R00=R0(ii);
A=R00/a0;
%k=z0dinput.geo.K(ii);
%delta=z0dinput.geo.d(ii);
Sf=q950*Ip0/(a0*B0);
V0=V(ii);
S0=S(ii);
betapol0=betapol(ii);
betaN0=betaN(ii);
beta0=Ip0.*betaN0./(a0.*B0);
nbar0=nbar(ii)/1e19;
nem0=nem(ii)/1e19;
nhem0=nhem(ii)/1e19;
nimpm0=nimpm(ii)/1e19;
ne0=nne(iipr,1)/1e19;
negr0=negr(ii)/1e19;
fgr=nbar0/negr0;
te0=Te(iipr,1);
tem0=tem(ii);
ti0=Ti(iipr,1);
tite0=tite(ii);
tim=tite0*tem0;
wth0=wth(ii)/1e6;
Pnbi10=Pnbi1(ii);
Pnbi20=Pnbi2(ii);
Pecrh0=Pecrh(ii);
tau0=tau(ii);
tauh0=tauh(ii);
H0=tau0/tauh0;
Vloop0=Vloop(ii);
Iboot0=Iboot(ii)/1e6;
Ipar0=Ipar(ii)/1e6;
fbs=Iboot0/Ipar0;
Ini0=Ini(ii)/1e6;
fni=Ini0/Ipar0;
tauj0=tauj(ii);
tauR=tauj0/14.682;
Sn0=Sn(ii);

ka=V0/(2*pi^2*R00*a0^2);
Ploss=P_loss(ii)/1e6;
Pradcore0=Pradcore(ii)/1e6;
Prad0=Prad(ii)/1e6;
Pradsol0=Pradsol(ii)/1e6;
Pbrem0=Pbrem(ii)/1e6;
Pcyclo0=Pcyclo(ii)/1e6;
Pnbi1_abs0=Pnbi1_abs(ii);
Pnbi2_abs0=Pnbi2_abs(ii);
Picrh0=Picrh(ii);
Plh0=Plh(ii);
Zeff0=Zeff(ii);
Palpha=P_fus(ii)/1e6;
Padd=Pnbi10+Pnbi20+Pecrh0+Palpha;
Pfus=5.0286*Palpha;
peakTe=te0/tem0;
peakTi=ti0/tim;
peakn=ne0/nem0;
fHe=nhem0/nem0*100;

flightimp=nimpm0/nem0*100;
fheavyimp=rimp*nimpm0/nem0*100;

M=2;
eps=1/A;
tauIPB=0.0562*Ip0.^0.93.*B0.^0.15.*Padd.^(-0.69).*nbar0.^0.41 ...
    .*M^0.19.*R00.^1.97.*eps.^0.58.*ka.^0.78;

tau_padd=H0*tauIPB;
tau_He=tauhe(ii);

% McDonald NF2007

Teped0=Teped(ii);
neped0=neped(ii);
Wped3=0.024*Ip0^1.64*R00^1.03*Padd^0.56*nem0^(-0.18)*eps^(-0.39);
Teped_sc=208*Wped3/(0.92*V0*neped0);

%tab=sprintf(\t)
Pradcyclo=Pradcore0+Pcyclo0;
if post.z0dinput.option.lhmode == 5
    aa=strvcat('time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC1(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','PEC2(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ');
else
    aa=strvcat('time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','Plh(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ');
end
bb=num2str([time ,Ip0 ,B0 ,q950 ,R00 ,a0 ,A ,k ,delta ,Sf ,V0 ,S0 ,beta0 ,betapol0 ,betaN0 ,nbar0 ,nem0 ,ne0 ,negr0 ,fgr ,tem0 ,te0 ,tim ,ti0 ,wth0 ,Padd ,Pnbi10 ,Pnbi20,Pecrh0 ,tau0 ,H0 ,Vloop0 ,fbs ,fni ,tauR ,Sn0 ,0 ,ka ,Ploss ,tau_padd ,tau_He ,Pnbi1_abs0 ,Pnbi2_abs0 ,Picrh0 ,Plh0 ,neped0 ,Teped0 ,Teped_sc ,Zeff0 ,Palpha ,Pfus ,peakTe ,peakTi ,peakn ,fHe ,fheavyimp ,flightimp ,Pradcore0 ,Pradsol0 ,Prad0 ,Pbrem0 ,Pcyclo0 ,Pradcyclo]',4);
table=[aa bb];
disp(table)

% make csv file
if post.z0dinput.option.lhmode == 5
    aa={'time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC1(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','PEC2(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) '};
else
    aa={'time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','Plh(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) '};
end
bb=[time ,Ip0 ,B0 ,q950 ,R00 ,a0 ,A ,k ,delta ,Sf ,V0 ,S0 ,beta0 ,betapol0 ,betaN0 ,nbar0 ,nem0 ,ne0 ,negr0 ,fgr ,tem0 ,te0 ,tim ,ti0 ,wth0 ,Padd ,Pnbi10 ,Pnbi20,Pecrh0 ,tau0 ,H0 ,Vloop0 ,fbs ,fni ,tauR ,Sn0 ,0 ,ka ,Ploss ,tau_padd ,tau_He ,Pnbi1_abs0 ,Pnbi2_abs0 ,Picrh0 ,Plh0 ,neped0 ,Teped0 ,Teped_sc ,Zeff0 ,Palpha ,Pfus ,peakTe ,peakTi ,peakn ,fHe ,fheavyimp ,flightimp ,Pradcore0 ,Pradsol0 ,Prad0 ,Pbrem0 ,Pcyclo0 ,Pradcyclo];
filename = sprintf('METIS_simulation_digest_table_TOK_%s@shot_%d.txt',post.z0dinput.machine,post.z0dinput.shot);
fid = fopen(filename,'w');
if fid > 0
    for k=1:length(bb)
        fprintf(fid,'%-25s\t%g\n',deblank(aa{k}),bb(k));
    end
    fclose(fid);
    disp('===================================================')
    fprintf('File %s as been created\n',filename);
    disp('===================================================')
end



            
function z0plot_shot_digest_dt(post,tprofi,ts_list)

% DEMO or ITER
disp('===================================================')
disp('Reactor tokamak');
%combines wwmetis and figMetis

tmax=max(post.profil0d.temps); %tprofi=0.95*tmax;
tmin_round = fix(min(post.profil0d.temps));
nfig=0;lls='-'; posfig=0;
% if (nargin == 1)
%     % selection of the time slice
%     evalin('base','z0plot_reference_post;');
%     subplot(3,1,1);
%     title('choose 1 time slice for profiles visualization')
%     drawnow
%     [tprofi,void] = ginput(1);
%     subplot(3,1,1);
%     title('choose 3 time slices for safety factor evolution visualization')
%     drawnow
%     [ts_list,void] = ginput(3);
%     close(gcf);
%     drawnow
% end
% ts_list = sort(ts_list);


dpos=-20; pos1=256; pos2=247; pos3=364; pos4=425;
%psize1=0.7; psize10=0.8; psize2=1.1; psize3=1.3; psize4=1.1; psize5=1.05; psize6=0.97;
psize1=1; psize10=1; psize2=1; psize3=1; psize4=1; psize5=1; psize6=1;

a=post.z0dinput.geo.a;
R0=post.z0dinput.geo.R;
eps=a./R0;
k=post.z0dinput.geo.K;
delta=post.z0dinput.geo.d;
B=post.z0dinput.geo.b0;
betapol=post.zerod.betaptot;
%betapol=post.zerod.betap;
rho=post.profil0d.rmx./(max(post.profil0d.rmx,[],2)*ones(1,21));
t=post.z0dinput.cons.temps;
ttt=post.profil0d.temps;
nbar=post.z0dinput.cons.nbar;
tau=post.zerod.taue;
tauh=post.zerod.tauh;
tauhe=post.zerod.tauhe;
wth=post.zerod.wth;
V=post.zerod.vp;
S=post.zerod.sp;
P_tot=post.zerod.pin;  %0-d quantity
P_loss=post.zerod.ploss;  %0-d quantity
P_l2h=post.zerod.plossl2h;  %0-d quantity
P_lhthr=post.zerod.plhthr;  %0-d quantity
P_lim=post.zerod.plim;  %0-d quantity
P_fus=post.zerod.pfus;  %0-d quantity, puissance alpha
P_nbi=post.zerod.pnbi;  %0-d quantity
P_ic=post.zerod.picrh;  %0-d quantity
P_ec=post.zerod.pecrh;  %0-d quantity
P_lh=post.zerod.plh;  %0-d quantity
P_aux=P_tot-P_fus;  %0-d quantity
Prad=post.zerod.prad;  %0-d quantity
Pradsol=post.zerod.pradsol;  %0-d quantity
Pbrem=post.zerod.pbrem;  %0-d quantity
Pcyclo=post.zerod.pcyclo;  %0-d quantity
Pradtot=Prad+Pbrem+Pcyclo;
Pradcore=Prad+Pbrem;
Q=5.0286*P_fus./(P_tot-P_fus);
Zeff=post.zerod.zeff;

q=post.profil0d.qjli;
j=post.profil0d.jli;
jbs=post.profil0d.jboot + post.profil0d.jfus;
jpar=post.profil0d.jeff;
Bpol=post.profil0d.bpol;
epar=post.profil0d.epar;
Te=post.profil0d.tep/1000;
Ti=post.profil0d.tip/1000;
xie=post.profil0d.xie;
xii=post.profil0d.xii;
ne=post.profil0d.nep;
ni=post.profil0d.nip;
Pfce=post.profil0d.pecrh;
Plhp=post.profil0d.plh;
Pnbi_i=post.profil0d.pnbi_ion;
Pnbi=post.profil0d.pnbi;
Pnbi_e=Pnbi-Pnbi_i;
Pfci_i=post.profil0d.picrh_ion;
Pfci=post.profil0d.picrh + post.profil0d.pfweh;
Pfci_e=Pfci-Pfci_i;
Pfus_i=post.profil0d.pfus_ion;
Pfus=post.profil0d.pfus;
Pfus_e=Pfus-Pfus_i;
P_i=post.profil0d.source_ion;
P_e=post.profil0d.source_el;
jfce=post.profil0d.jeccd;
jlh=post.profil0d.jlh;
jnbi=post.profil0d.jnbicd;
jfci=post.profil0d.jfwcd;
Pohm=post.profil0d.pohm;
q_ei=post.profil0d.qei;
Psyn=post.profil0d.pcyclo;
P_rad=post.profil0d.prad;
P_brem=post.profil0d.pbrem;

Ip=post.zerod.ip;
Ipar=post.zerod.ipar;
li=post.zerod.li;
Vmes=post.zerod.vmes;
Vloop=post.zerod.vloop;
Ioh=post.zerod.iohm;
Ilh=post.zerod.ilh;
Ifce=post.zerod.ieccd;
Iboot=post.zerod.iboot;
Ibs=post.zerod.iboot;
Ini=post.zerod.ini;
Inbi=post.zerod.inbicd;
psi=post.profil0d.psi;
ptot=post.profil0d.ptot; %total pressure

nHe=post.profil0d.nhep;
betaN=post.zerod.betan*100;
beta=1e6*Ip.*betaN./(a.*B);
Ptot=Pnbi+Pfci+Pfus+Pohm+Pfce+Plhp;
f_b_s=Iboot./Ipar;
f_n_i=Ini./Ipar;
nG=(Ip./1e6)/pi./a.^2;
f_G=(nbar/1e20)./nG;
ka=V./(2*pi^2*R0.*a.^2);
M=2.5;
tauIPB=0.0562*(Ip/1e6).^0.93.*B.^0.15.*(P_loss/1e6).^(-0.69).*(nbar/1e19).^0.41 ...
    .*M^0.19.*R0.^1.97.*eps.^0.58.*ka.^0.78;
H=tau./tauh;
tauIPB=tauh;

if exist('post.profil0d.Rsepa')==1
    Rse=post.profil0d.Rsepa;
    Zse=post.profil0d.Zsepa;
else
end
q95=post.zerod.q95;
nem=post.zerod.nem;
nhem=post.zerod.nhem;
nimpm=post.zerod.nimpm;
negr=post.zerod.negr;
tem=post.zerod.tem/1e3;
tite=post.zerod.tite;
Pnbi1=real(post.z0dinput.cons.pnbi)/1e6;
Pnbi2=imag(post.z0dinput.cons.pnbi)/1e6;
Pnbi1_abs=real(post.zerod.pnbi)/1e6;
Pnbi2_abs=imag(post.zerod.pnbi)/1e6;
Pecrh=post.z0dinput.cons.pecrh/1e6;
Picrh=post.z0dinput.cons.picrh/1e6;
Plh=post.z0dinput.cons.plh/1e6;
tauj=post.zerod.tauj;
Sn=post.zerod.ndd;
Teped=post.zerod.teped/1e3;
neped=post.zerod.neped/1e19;
rimp=post.z0dinput.option.rimp;

eval('nne=ne;');eval('pptot=ptot;');eval('jj=j;');eval('bbeta=beta;');

% GLOBAL QUANTITIES


tit=['t = ',num2str(tprofi),' s'];
elon=k;

nfig=nfig+1;
figuren(nfig); set(gcf,'name','Ip,ne0,Q,P vs t','Position',[pos1+posfig pos2+posfig pos3*psize3 pos4*psize3]);
subplot(211); set(gca,'fontsize',fix(14*psize3));
plot(t,Ip/1e6,ttt,nne(:,1)/1e19,t,Q,'linewidth',2,'linestyle',lls)
xlabel('time (s)')
axis([tmin_round  tmax 0 Inf])
legend('Ip (MA)','ne(0) (10^{19} m^{-3})','Q','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize3));
if post.z0dinput.option.lhmode == 5
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,(P_ec + P_lh)/1e6,t,P_fus/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Pic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,(P_ec + P_lh)/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Location','Best')
    end
elseif any(P_lh > 1e3)
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_fus/1e6,t,P_lh/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Plh','Pic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_fus/1e6,t,P_lh/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Plh','Location','Best')
    end
else
    if any(P_ic > 1e3)
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_fus/1e6,t,P_ic/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Pic','Location','Best')
    else
        plot(t,(real(P_nbi)+imag(P_nbi))/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        %legend('Pnnbi','Ppnbi','Pec','Palpha')
        %plot(t,P_nbi/1e6,t,P_ec/1e6,t,P_fus/1e6,'linewidth',2,'linestyle',lls)
        legend('Pnbi','Pec','Palpha','Location','Best')
    end
end
axis([tmin_round  tmax 0 Inf])
xlabel('time (s)')
ylabel('(MW)')

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','fbs,fni,fG,betaN,li,H,Zeff vs t','Position',[pos1+posfig pos2+posfig pos3*psize2 pos4*psize2]);
subplot(211); set(gca,'fontsize',fix(14*psize2));
plot(t,f_b_s,t,f_n_i,t,f_G,'linewidth',2,'linestyle',lls)
hold on;plot([0 tmax],[1 1],'k--');hold off;
xlabel('time (s)')
axis([tmin_round  tmax 0 1.5])
legend('f_b_s','f_n_i','f_G','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize2));
HH=smooth(H,4);
plot(t,betaN,t,4*li,t,HH,t,Zeff,'linewidth',2,'linestyle',lls)
hold on;plot([0 tmax],[1 1],'k--');hold off;
axis([tmin_round  tmax 0 5.5])
legend('betaN','4 li','H','Zeff','Location','Best')
xlabel('time (s)')

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','gCD,Vl,Prad,Pbrem,Psyn vs t','Position',[pos1+posfig pos2+posfig pos3*psize2 pos4*psize2]);
subplot(211); set(gca,'fontsize',fix(14*psize2));
gCD=(Ini-Ibs)./(real(P_nbi)+imag(P_nbi)+P_lh+P_ec).*nbar.*R0/1e20;
plot(t,gCD,t,smooth(Vloop,4),'linewidth',2,'linestyle',lls)
hold on;plot([0 tmax],[0 0],'k--');hold off;
xlabel('time (s)')
axis([tmin_round  tmax -0.05 Inf])
legend('gammaCD','Vloop (V)','Location','Best')
subplot(212); set(gca,'fontsize',fix(14*psize2));
plot(t,Pradtot/1e6,t,Pbrem/1e6,t,Pcyclo/1e6,'linewidth',2,'linestyle',lls)
legend('Pradtot','Pbrem','Psyn','Location','Best')
axis([tmin_round  tmax 0 Inf])
xlabel('time (s)')
ylabel('(MW)')

% PROFILES

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','Te,Ti,ne profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
k=iround(ttt,tprofi);
h=plot(rho(k,:),Te(k,:),rho(k,:),Ti(k,:),rho(k,:),nne(k,:)/1e19); set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
legend('Te (keV)','Ti (keV)','ne (10^{19} m^{-3})','Location','Best')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','chi_e,chi_i profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
h=plot(rho(k,:),xie(k,:),rho(k,:),xii(k,:)); set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
legend('chie (m^2/s)','chii (m^2/s)','Location','NorthWest')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','P profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize4 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
l=(1:21);
if any(Plhp(k,l) > 1e3)
    if any(Pfci(k,l) > 1e3)
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6+Plhp(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Palpha','Pnbi','Pec','Pic','Location','Best');
        else
            h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Plhp(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Palpha','Pnbi','Pec','Plh','Pic','Location','Best');
        end
    else
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6+Plhp(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Palpha','Pnbi','Pec','Location','Best');
        else
            h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Plhp(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('Palpha','Pnbi','Pec','Plh','Location','Best');
        end
    end
else
    if any(Pfci(k,l) > 1e3)
        h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6,rho(k,l),Pfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('Palpha','Pnbi','Pec','Pic','Location','Best')
    else
        h=plot(rho(k,l),Pfus(k,l)/1e6,rho(k,l),Pnbi(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('Palpha','Pnbi','Pec','Location','Best')
    end
end
%h=plot(rho(k,l),Pnbi_i(k,l)/1e6,rho(k,l),Pnbi_e(k,l)/1e6,rho(k,l),Pfce(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
%legend('Pnbi-i','Pnbi-e','Pec')
axis([0 1 0 Inf])
xlabelrho
ylabel('(MW/m^3)')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','j profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize4 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10))
%k=iround(t,tprofi);
l=(1:21);
if any(jlh(k,l) > 1e3)
    if any(jfci(k,l) > 1e3)
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),(jfce(k,l)+ jlh(k,l))/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jfwcd','Location','Best')
        else
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jlh(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jlh','jfwcd','Location','Best')
        end
    else
        if post.z0dinput.option.lhmode == 5
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),(jfce(k,l)+ jlh(k,l))/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','Location','Best')
        else
            h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jlh(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
            legend('j','jnbi','jec','jbs','jlh','Location','Best')
        end
    end
else
    if any(jfci(k,l) > 1e3)
        h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6,rho(k,l),jfci(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('j','jnbi','jec','jbs','jfwcd','Location','Best')
    else
        h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
        legend('j','jnbi','jec','jbs','Location','Best')
    end
end
%h=plot(rho(k,l),jj(k,l)/1e6,rho(k,l),jnbi(k,l)/1e6,rho(k,l),jfce(k,l)/1e6,rho(k,l),jbs(k,l)/1e6);set(h,'linewidth',2,'linestyle',lls);
%legend('j','jnbi','jec','jbs','Location','Best')
axis([0 1 0 Inf])
xlabelrho
ylabel('(MA/m^2)')
title(tit);

posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','q profiles','Position',[pos1+posfig pos2+posfig pos3*psize1*psize5 pos4*psize1*psize6]);
set(gca,'FontSize',fix(18*psize10))
kk=iround(ttt,ts_list);
%kk=iround(ttt,[500 1000 2000]);
%kk=iround(ttt,[1000 2000 2900]);
h=plot(rho(kk,:)',q(kk,:)',[0 1],[1 1],'k--');set(h,'linewidth',2,'linestyle',lls);
axis([0 1 0 Inf])
xlabelrho
ylabel('q')
%legend('t = 1000 s','t = 2000 s','t = 3000 s')
%legend('t = 500 s','t = 1000 s','t = 2000 s')
legend(sprintf('t = %0.2f s',ts_list(1)),sprintf('t = %0.2f s',ts_list(2)),sprintf('t = %0.2f s',ts_list(3)),'Location','Best')
%legend('t = 500 s','t = 2000 s','t = 10000 s','Location','Best')

qlim = ceil(min(max(post.zerod.qa),2.1 .* trapz(post.zerod.temps,post.zerod.ip .* post.zerod.qa) ./ max(1,trapz(post.zerod.temps,post.zerod.ip))));
posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','q vs t','Position',[pos1+posfig pos2+posfig pos3*psize1*psize5 pos4*psize1*psize6]);
set(gca,'FontSize',fix(18*psize10))
h=plot(ttt,q(:,1),ttt,q(:,11),ttt,q(:,20),[0 tmax],[1 1],'k--');set(h,'linewidth',2,'linestyle',lls);
axis([tmin_round tmax 0 qlim])
xlabel('time (s)')
ylabel('q(rho)')
legend('q(0)','q(0.5)','q(0.95)','Location','Best')

%  TABLE

t0= tprofi;

ii=iround(t,t0);
iipr=iround(ttt,t0);

if exist('Rse')==1
    
    R=Rse(iipr,:);
    Z=Zse(iipr,:);
    
    imax=max(find(Z==max(Z)));
    imin=max(find(Z==min(Z)));
    Zmax=Z(imax); Zmin=Z(imin);
    
    jmax=max(find(R==max(R)));
    jmin=max(find(R==min(R)));
    Rmax=R(jmax); Rmin=R(jmin);
    
    RZmax=R(imax); RZmin=R(imin);
    
    RR0=(Rmax+Rmin)/2;
    aa=(Rmax-Rmin)/2;
    k=(Zmax-Zmin)/(Rmax-Rmin);
    dup=(RR0-RZmax)/aa;
    dlo=(RR0-RZmin)/aa;
    delta=(dup+dlo)/2;
else
    k=elon(ii);
    delta=delta(ii);
    
end


time=t(ii);
Ip0=Ip(ii)/1e6;
B0=B(ii);
q950=q95(ii);
a0=a(ii);
R00=R0(ii);
A=R00/a0;
%k=z0dinput.geo.K(ii);
%delta=z0dinput.geo.d(ii);
Sf=q950*Ip0/(a0*B0);
V0=V(ii);
S0=S(ii);
betapol0=betapol(ii);
betaN0=betaN(ii);
beta0=Ip0.*betaN0./(a0.*B0);
nbar0=nbar(ii)/1e19;
nem0=nem(ii)/1e19;
nhem0=nhem(ii)/1e19;
nimpm0=nimpm(ii)/1e19;
ne0=nne(iipr,1)/1e19;
negr0=negr(ii)/1e19;
fgr=nbar0/negr0;
te0=Te(iipr,1);
tem0=tem(ii);
ti0=Ti(iipr,1);
tite0=tite(ii);
tim=tite0*tem0;
wth0=wth(ii)/1e6;
Pnbi10=Pnbi1(ii);
Pnbi20=Pnbi2(ii);
Pecrh0=Pecrh(ii);
tau0=tau(ii);
tauh0=tauh(ii);
H0=tau0/tauh0;
Vloop0=Vloop(ii);
Iboot0=Iboot(ii)/1e6;
Ipar0=Ipar(ii)/1e6;
fbs=Iboot0/Ipar0;
Ini0=Ini(ii)/1e6;
fni=Ini0/Ipar0;
tauj0=tauj(ii);
tauR=tauj0/14.682;
Sn0=Sn(ii);
gCD0=gCD(ii);
ka=V0/(2*pi^2*R00*a0^2);
Ploss=P_loss(ii)/1e6;
Pl2h=P_l2h(ii)/1e6;
Plhthr=P_lhthr(ii)/1e6;
Plim=P_lim(ii)/1e6;
Pradcore0=Pradcore(ii)/1e6;
Prad0=Prad(ii)/1e6;
Pradsol0=Pradsol(ii)/1e6;
Pbrem0=Pbrem(ii)/1e6;
Pcyclo0=Pcyclo(ii)/1e6;
Pnbi1_abs0=Pnbi1_abs(ii);
Pnbi2_abs0=Pnbi2_abs(ii);
Picrh0=Picrh(ii);
Plh0=Plh(ii);
Zeff0=Zeff(ii);
Palpha=P_fus(ii)/1e6;
Padd=Pnbi10+Pnbi20+Pecrh0+Picrh0+Plh0+Palpha;
Pfus0=5.0286*Palpha;
Q0=Q(ii);
peakTe=te0/tem0;
peakTi=ti0/tim;
peakn=ne0/nem0;
fHe=nhem0/nem0*100;

flightimp=nimpm0/nem0*100;
fheavyimp=rimp*nimpm0/nem0*100;

M=2.5;
eps=1/A;

tauIPB=0.0562*Ip0.^0.93.*B0.^0.15.*Padd.^(-0.69).*nbar0.^0.41 ...
    .*M^0.19.*R00.^1.97.*eps.^0.58.*ka.^0.78;
%tauIPB=0.0562*Ip0.^0.93.*B0.^0.15.*Ploss.^(-0.69).*nbar0.^0.41 ...
%           .*tmin_roundM^0.19.*R00.^1.97.*eps.^0.58.*ka.^0.78
tau_padd=H0*tauIPB;
tau_He=tauhe(ii);

% McDonald NF2007
Teped0=Teped(ii);
neped0=neped(ii);
Wped3=0.024*Ip0^1.64*R00^1.03*Padd^0.56*nem0^(-0.18)*eps^(-0.39);
Teped_sc=208*Wped3/(0.92*V0*neped0);

%tab=sprintf(\t)
Pradcyclo=Pradcore0+Pcyclo0;
%Pdiv0=Padd-Pradcyclo;
Pdiv0=Plim;
radfrac0=Pradcyclo/Padd;
radfrac1=(Pcyclo0+Pbrem0)/Padd;

%l2hfrac0=1-Pl2h/Padd;
l2hsploss0=Ploss/Pl2h;
PsPthr=Plhthr/Pl2h;


if post.z0dinput.option.lhmode == 5
    aa=strvcat('time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC1(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','PEC2(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ','Pdiv(MW) ','Pl2h(MW) ','core rad. fraction ','Ploss/Pl2h ','P/Pthr ','gCD(10^20 A/W/m??) ','Q ');
else
    aa=strvcat('time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','Plh(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ','Pdiv(MW) ','Pl2h(MW) ','core rad. fraction ','Ploss/Pl2h ','P/Pthr ','gCD(10^20 A/W/m??) ','Q ');
end
bb=num2str([time ,Ip0 ,B0 ,q950 ,R00 ,a0 ,A ,k ,delta ,Sf ,V0 ,S0 ,beta0 ,betapol0 ,betaN0 ,nbar0 ,nem0 ,ne0 ,negr0 ,fgr ,tem0 ,te0 ,tim ,ti0 ,wth0 ,Padd ,Pnbi10 ,Pnbi20,Pecrh0 ,tau0 ,H0 ,Vloop0 ,fbs ,fni ,tauR ,Sn0 ,0 ,ka ,Ploss ,tau_padd ,tau_He ,Pnbi1_abs0 ,Pnbi2_abs0 ,Picrh0 ,Plh0 ,neped0 ,Teped0 ,Teped_sc ,Zeff0 ,Palpha ,Pfus0 ,peakTe ,peakTi ,peakn ,fHe ,fheavyimp ,flightimp ,Pradcore0 ,Pradsol0 ,Prad0 ,Pbrem0 ,Pcyclo0 ,Pradcyclo ,Pdiv0 ,Pl2h ,radfrac0 ,l2hsploss0 ,PsPthr ,gCD0 ,Q0]',4);
table=[aa bb];
disp(table)
if post.z0dinput.option.lhmode == 5
    aa={'time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC1(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','PEC2(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ','Pdiv(MW) ','Pl2h(MW) ','core rad. fraction ','Ploss/Pl2h ','P/Pthr ','gCD(10^20 A/W/m??) ','Q '};
else
    aa={'time(s) ','Ip(MA) ','B(T) ','q95 ','R(m) ','a(m) ','A ','k ','delta (av) ','Sf ','V(m^3) ','S(m^2) ','beta ','betap ','betaN ','nbar(10^19 m^-3) ','nem(10^19 m^-3) ','ne0(10^19 m^-3) ','nG(10^19 m^-3) ','fG ','Tem(keV) ','Te0(keV) ','Tim(keV) ','Ti0(keV) ','Wth(MJ) ','Padd(MW) ','PNBI1(MW) ','PNBI2(MW) ','PEC(MW) ','taue(s) ','H ','Vloop(V) ','fbs ','fni ','tauR(s) ','Sn(n/s) ','  ','ka ','Ploss(MW) ','tau_padd(s) ','tau_He(s) ','Pnbi1_abs(MW) ','Pnbi2_abs(MW) ','Picrh(MW) ','Plh(MW) ','neped(10^19 m^-3) ','Teped(keV) ','Teped_sc(keV) ','Zeff ','Palpha (MW) ','Pfus (MW) ','peak(Te) ','peak(Ti) ','peak(n) ','fHe(%) ','f_heavy_imp(%) ','f_light_imp(%) ','Prad_core(MW) ','Prad_SOL(MW) ','Pline(MW) ','Pbrem(MW) ','Pcyclo(MW) ','Pline+brem+cyclo(MW) ','Pdiv(MW) ','Pl2h(MW) ','core rad. fraction ','Ploss/Pl2h ','P/Pthr ','gCD(10^20 A/W/m??) ','Q '};
end
bb=[time ,Ip0 ,B0 ,q950 ,R00 ,a0 ,A ,k ,delta ,Sf ,V0 ,S0 ,beta0 ,betapol0 ,betaN0 ,nbar0 ,nem0 ,ne0 ,negr0 ,fgr ,tem0 ,te0 ,tim ,ti0 ,wth0 ,Padd ,Pnbi10 ,Pnbi20,Pecrh0 ,tau0 ,H0 ,Vloop0 ,fbs ,fni ,tauR ,Sn0 ,0 ,ka ,Ploss ,tau_padd ,tau_He ,Pnbi1_abs0 ,Pnbi2_abs0 ,Picrh0 ,Plh0 ,neped0 ,Teped0 ,Teped_sc ,Zeff0 ,Palpha ,Pfus0 ,peakTe ,peakTi ,peakn ,fHe ,fheavyimp ,flightimp ,Pradcore0 ,Pradsol0 ,Prad0 ,Pbrem0 ,Pcyclo0 ,Pradcyclo ,Pdiv0 ,Pl2h ,radfrac0 ,l2hsploss0 ,PsPthr ,gCD0 ,Q0];
filename = sprintf('METIS_simulation_digest_table_TOK_%s@shot_%d.txt',post.z0dinput.machine,post.z0dinput.shot);
fid = fopen(filename,'w');
if fid > 0
    for k=1:length(bb)
        fprintf(fid,'%-25s\t%g\n',deblank(aa{k}),bb(k));
    end
    fclose(fid);
    disp('===================================================')
    fprintf('File %s as been created\n',filename);
    disp('===================================================')
end


% recover figure handle
function h = figuren(n)

tag = sprintf('z0plot_shot_digest_%d',n);
h = findobj(0,'type','figure','tag',tag);
if isempty(h)
    h=figure('tag',tag,'color',[1 1 1]);
else
    figure(h);
    set(h,'color',[1 1 1]);
end
clf
setcoloroldcolororder(h);

function setcoloroldcolororder(hfig)

if nargin == 0
    hfig = gcf;
end

co = [
    0         0    1.0000
    0    0.5000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
    ];



set(hfig,'defaultAxesColorOrder',co);



function xlabelrho

if verLessThan('matlab','8.6')
    xlabel('r','fontname','symbol');
else
    xlabel('\rho');
end
