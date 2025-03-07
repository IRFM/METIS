function x=batchcronos(data,tstart,tend,comment,tsamp)
% plots the cronos batchs
% use: >>batchcronos('mat name', start time(opt), end time(opt), comment(80 char., string),marker time(opt)) 
if nargin==1,tstart=0;tend=30;comment=' ',tsamp=tstart;end 
if nargin==2,tend=30;comment=' ',tsamp=tstart;end 
if nargin==3,comment=' ',tsamp=tstart;end 
if nargin==4,tsamp=tstart;end 
 
clf
 
hbox=.12; 
dt=tend-tstart; 
if dt<300,dti=30;end
if dt<60,dti=5;end 
if dt<20,dti=1;end 
if dt<2,dti=.1;end 
if dt<.2,dti=.01;end 
if dt<.02,dti=.001;end 
%global scenario parameters 
figure(1) 
clf
set(1,'renderer','paint','doublebuffer','on') 
h=axes('position',[0 0 1 1]); 
set(h,'visible','off') 
t00=text(0.3,0.98,['CRONOS - #' nommat]);set(t00,'fontsize',14);set(t00,'fontweight','bold')
t01=text(0.1,0.05,comment);set(t01,'fontsize',12);set(t01,'fontweight','bold') 


eval(['load ' nommat '_resultat'])


tps=data.gene.temps;

%toroidal field
Btor1=data.geo.b0;Btor=eval(Btor1);
%ICRH power 
Picrh=data.gene.paddfci;Picrh=Picrh/1e6;
%LHCD power 
Plh=data.gene.paddhyb;Plh=Plh/1e6;
%IDN
Pidn=data.gene.paddidn;Pidn=Pidn/1e6;
%FCE
Pfce=data.gene.paddfce;Pfce=Pfce/1e6;
%Ploss
Ploss=data.gene.ploss;Ploss=Ploss/1e6;
ploss_rec=data.gene.ploss+data.gene.dwdiadt-gradient(data.gene.wdia,tps(2)-tps(1));
%plasma current 
Ipla=data.gene.ip;Ipla=Ipla/1e6;

h=axes('position',[0.05 .94-hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
p1=plot(tps,Picrh,'b',tps,Plh,'r',tps,Pfce,'g',tps,ploss_rec/1e6,'-.k',[tsamp tsamp],[0 20],'--r'); 
set(p1,'linewidth',2); 
axis([tstart tend 0 20])

title([' Bt=' sprintf('%0.2g',Btor(1)) 'T']) 
t14=text(tstart+1,17,'Ploss(MW)');set(t14,'color','black')
t13=text(tstart+1,11,'Picrh(MW)');set(t13,'color','blue')
t12=text(tstart+1,8,'Plhcd(MW)');set(t12,'color','red') 
t10=text(tstart+1,5,'Pecrh(MW)');set(t10,'color','green') 
t11=text(tstart+1,14,'Pnbi(MW)');set(t11,'color','cyan')
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[5:5:15]) 
set(h,'yticklabel',[5:5:15]) 
 
%ohmic power 
Pohm=data.gene.paddohm;Pohm=Pohm/1e6; 
%Radiated power 
Prad=data.gene.prad;Prad=Prad/1e6;
 
 
h=axes('position',[0.05 .94-2*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
p2=plot(tps,Pohm,'g',tps,Prad,'b',[tsamp tsamp],[0 10],'--r'); 
set(p2,'linewidth',2); 
axis([tstart tend 0 10]) 
t22=text(tstart+1,5,'Pohm(MW)');set(t22,'color','green') 
t23=text(tstart+1,2,'Prad(MW)');set(t23,'color','blue') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[2:2:8]) 
set(h,'yticklabel',[2:2:8]) 
 
%volume averaged density 
nvol=data.gene.nemoy;
nvol=nvol/1e19;
ne=data.prof.ne;ne=ne/1e19;
 
h=axes('position',[0.05 .94-3*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
p3=plot(tps,ne(:,1),'r',tps,nvol,'g',[tsamp tsamp],[0 10],'--r'); 
set(p3,'linewidth',2); 
axis([tstart tend 0 10]) 
t31=text(tstart+1,6.5,'ne0(e19m-3)');set(t31,'color','red') 
t32=text(tstart+1,4.5,'<n>(e19m-3)');set(t32,'color','green') 
xlabel('t(s)') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[1:9]) 
set(h,'yticklabel',[1:9]) 


%major radius 
R=data.equi.raxe;
%minor radius 
a=data.equi.a;
 
h=axes('position',[0.05 .94-4*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
 
p4=plot(tps,R(:,1)-6,'r',tps,a(:,101)-2,'g',[tsamp tsamp],[-.5 .5],'--r'); 
set(p4,'linewidth',2); 
axis([tstart tend -.5 .5]) 
t41=text(tstart+1,.1,'Rmag-6(m)');set(t41,'color','red') 
t42=text(tstart+1,-.1,'a-2(m)');set(t42,'color','green') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[-.5:.2:.5]) 
set(h,'yticklabel',[-.5:.2:.5]) 

 
%diamagnetic energy 
wdia=data.gene.wdia;wdia=wdia/1e6; 
%neutrons 
Rnt=data.source.totale.neutron.dt;

h=axes('position',[0.55 .94-hbox 0.4 hbox]);
set(h,'drawmode','fast') 
set(h,'box','on') 
 
p5=plot(tps,wdia,'r',tps,log10(Rnt(:,1)),'g',[tsamp tsamp],[0 16],'--r'); 
set(p5,'linewidth',2); 
axis([tstart tend 0 16]) 
t51=text(tstart+1,6,'Wdia(MJ)');set(t51,'color','red') 
t52=text(tstart+1,9,'log(Rnt(0))');set(t52,'color','green')
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:3:15]) 
set(h,'yticklabel',[0:3:15]) 

%betap 
betap=data.gene.betap; 
%betan 
betan=(4*betap.*Ipla)./(a(:,101).*Btor);
 
h=axes('position',[0.55 .94-2*hbox 0.4 hbox]); 
set(h,'drawmode','fast') 
set(h,'box','on') 
 
p6=plot(tps,betap,'r',tps,betan,'g',[tsamp tsamp],[0 1.3],'--r'); 
set(p6,'linewidth',2); 
axis([tstart tend 0 1.3]) 
t61=text(tstart+1,0.8,'betap');set(t61,'color','red') 
t62=text(tstart+1,1,'betaN');set(t62,'color','green') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0.2:.2:1.3]) 
set(h,'yticklabel',[0.2:.2:1.3]) 
 

%Hfactors 
M=2

%scaling ITER L-mode - tau L thermique - NF39(1999)2206
tauL=0.023*(eval(data.geo.r0)).^1.89*M^0.2.*(data.cons.ip/1e6).^0.96.*(eval(data.geo.b0)).^0.03.*(data.equi.e(:,101)).^0.64.*...
(eval(data.geo.a)).^-0.06.*(data.gene.nbar./2e19./eval(data.geo.a)).^0.4.*(ploss_rec/1e6).^-0.73;
tauL(find(isnan(tauL)==1))=0;
tauL(find(iscomplex(tauL)==1))=0;

%scaling ITER H-mode - IPB98(y,2) - NF39(1999)2208
tauH=0.0562*(eval(data.geo.r0)).^1.97*M^0.19.*(data.cons.ip/1e6).^0.93.*(eval(data.geo.b0)).^0.15.*(data.equi.e(:,101)).^0.78.*...
(eval(data.geo.a)./eval(data.geo.r0)).^0.58.*(data.gene.nbar./2e19./eval(data.geo.a)).^0.41.*(ploss_rec/1e6).^-0.69;
tauH(find(isnan(tauH)==1))=0;
tauH(find(iscomplex(tauH)==1))=0;

%scaling ITER H-mode - ITER89-P - NF37(1997)1303
tauH2=0.038*(eval(data.geo.r0)).^1.5*M^0.5.*(data.cons.ip/1e6).^0.85.*(eval(data.geo.b0)).^0.2.*(data.equi.e(:,101)).^0.5.*...
(eval(data.geo.a)./eval(data.geo.r0)).^0.3.*(data.gene.nbar./2e19./eval(data.geo.a)).^0.1.*(ploss_rec/1e6).^-0.5;
tauH2(find(isnan(tauH2)==1))=0;
tauH2(find(iscomplex(tauH2)==1))=0;

 
h=axes('position',[0.55 .94-3*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
set(h,'box','on') 
 
p7=plot(tps,data.gene.taue./tauL,'r',tps,data.gene.taue./tauH,'g',tps,data.gene.taue./tauH2,'b',[tsamp tsamp],[0 3],'--r'); 
set(p7,'linewidth',2); 
axis([tstart tend 0 3]) 
t71=text(tstart+1,1.5,'H-L');set(t71,'color','red') 
t72=text(tstart+1,1,'H98th');set(t72,'color','green') 
t73=text(tstart+1,0.5,'H89-P');set(t72,'color','blue') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[.5:.5:2.5]) 
set(h,'yticklabel',[.5:.5:2.5]) 


%Vloop 
vloop=data.gene.vloop; 
h=axes('position',[0.55 .94-4*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
 
p8=plot(tps,vloop,'r',[tsamp tsamp],[0 2],'--r'); 
set(p8,'linewidth',2); 
axis([tstart tend 0 2]) 
t81=text(tstart+1,0.5,'Vloop');set(t81,'color','red') 
grid 
 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0.1:0.4:1.9]) 
set(h,'yticklabel',[.1:.4:1.9]) 
 

%Te 
te=data.prof.te;te=te/1e3;  

h=axes('position',[0.05 .94-5*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
 
p9=plot(tps,te(:,1),'r',[tsamp tsamp],[0 10],'--r'); 
set(p9,'linewidth',2); 
axis([tstart tend 0 10]) 
t91=text(tstart+1,7,'Te0(keV)');set(t91,'color','red')   
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[2:2:8]) 
set(h,'yticklabel',[2:2:8]) 
 
 
%Ti 
ti=data.prof.ti;ti=ti/1e3;  

h=axes('position',[0.05 .94-6*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
 
pa=plot(tps,ti(:,1),'+r',[tsamp tsamp],[0 10],'--r'); 
set(pa,'linewidth',2); 
axis([tstart tend 0 10]) 
ta1=text(tstart+1,5,'Ti0(keV)');set(ta1,'color','red')  
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[2:2:8]) 
set(h,'yticklabel',[2:2:8]) 

 
% 
 
h=axes('position',[0.05 .94-7*hbox 0.4 hbox]);
set(h,'drawmode','fast') 

pb=plot([tsamp tsamp],[0 1e5],'--r'); 
%set(pb,'linewidth',2); 
axis([tstart tend 0 1e5]) 
tb1=text(tstart+1,8e4,' ');set(tb1,'color','red') 
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'ytick',[.2:.2:0.8]*1e5) 
set(h,'yticklabel',[.2:.2:0.8]*1e5) 
xlabel('t(s)') 

 
%bootstrap current

iboot=data.gene.iboot/1e6;
icd=data.gene.icd/1e6;
 
h=axes('position',[0.55 .94-5*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
pc=plot(tps,Ipla,'r',tps,iboot,'g',tps,icd,'b',[tsamp tsamp],[0 10],'--r'); 
set(pc,'linewidth',2); 
axis([tstart tend 0 10]) 
tc1=text(tstart+1,9,'Ip(MA)');set(tc1,'color','red'); 
tc2=text(tstart+1,5,'Iboot(MA)');set(tc2,'color','green');
tc3=text(tstart+1,2,'ICD(MA)');set(tc3,'color','blue'); 
 
grid 

set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[1:9]) 
set(h,'yticklabel',[1:9]) 
 

 
%q95 
q=data.prof.q; q95=q(:,95);q0=q(:,1);
%li 
li=data.equi.li;
iq=size(q);
qmin=zeros(iq(1));
for ii=1:iq(1)
qmin(ii)=min(q(ii,2:iq(2)));
end
 
 
h=axes('position',[0.55 .94-6*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
pd=plot(tps,li,'c',tps,q0,'r',tps,q95,'g',tps,qmin,'b',[tsamp tsamp],[0 6],'--r'); 
set(pd,'linewidth',2); 
axis([tstart tend 0 6]) 
td2=text(tstart+1,5,'q95');set(td2,'color','green'); 
td3=text(tstart+1,3,'q0');set(td3,'color','red'); 
td4=text(tstart+1,1,'li');set(td4,'color','cyan');
td5=text(tstart+1,2,'qmin');set(td5,'color','blue');
 
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[1:5]) 
set(h,'yticklabel',[1:5]) 

 
%Zeff 
zeff=data.impur.zeff; 
 
h=axes('position',[0.55 .94-7*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
pe=plot(tps,zeff(:,1),'r',[tsamp tsamp],[0 3],'--r'); 
set(pe,'linewidth',2); 
axis([tstart tend 0 3]) 
te1=text(tstart+1,2,'Zeff(0)');set(te1,'color','red'); 
 
 
grid 
set(h,'xtick',[floor(tstart):dti:floor(tend)]) 
%set(h,'xticklabel',[]) 
set(h,'ytick',[1:2]) 
set(h,'yticklabel',[1:2]) 
xlabel('t(s)')
