function zbatchcronos
   
   
data = evalin('base','data');
param = evalin('base','param');
   
tstart = param.gene.tdeb;
tend   = param.gene.t;
tsamp  = tstart;
comment= param.from.machine;
 
h = findobj(0,'type','figure','tag','batchcronos');
if isempty(h)
       h=figure('tag','batchcronos');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

hbox=.12; 
dt=tend-tstart; 
dti = round(dt/5*100)/100;
%if dt<300,dti=30;end
%if dt<60,dti=5;end 
%if dt<20,dti=1;end 
%if dt<2,dti=.2;end 
%if dt<.2,dti=.01;end 
%if dt<.02,dti=.001;end 
%global scenario parameters 

 
h=axes('position',[0 0 1 1]); 
set(h,'visible','off') 
nom = param.gene.file;
ind = find(nom=='_')
nom(ind)='-';
lim = max(find(nom=='/'));
t00=text(0.3,0.98,['CRONOS : ' nom((lim+1):end)]);
set(t00,'fontsize',14);
set(t00,'fontweight','bold')
t01=text(0.1,0.05,param.from.creation.com);
set(t01,'fontsize',12);
set(t01,'fontweight','bold') 




tps=data.gene.temps;

%toroidal field
Btor=data.geo.b0;
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
Pmax = ceil(max([max(Picrh) max(Plh) max(Pfce) max(ploss_rec/1e6)])*1.1);
if Pmax > 20
   Pens = [max(Picrh) max(Plh) max(Pfce) max(ploss_rec/1e6)];
   ind = find(Pens > 20);
   Pens(ind) = [];
   Pmax = ceil(max(Pens)*1.1);
end
if Pmax < 1
  Pmax=5;
end
p1=plot(tps,Picrh,'b',tps,Plh,'r',tps,Pfce,'g',tps,ploss_rec/1e6,'-.k',[tsamp tsamp],[0 Pmax],'--r'); 
set(p1,'linewidth',2); 

axis([tstart tend 0 Pmax])

title([' Bt=' sprintf('%1.3g',Btor(1)) 'T, machine :',param.from.machine]) 
t14=text(tstart+1,Pmax*0.9,'Ploss(MW)');set(t14,'color','black')
t13=text(tstart+1,Pmax*0.7,'Picrh(MW)');set(t13,'color','blue')
t12=text(tstart+1,Pmax*0.5,'Plhcd(MW)');set(t12,'color','red') 
t10=text(tstart+1,Pmax*0.3,'Pecrh(MW)');set(t10,'color','green') 
t11=text(tstart+1,Pmax*0.1,'Pnbi(MW)');set(t11,'color','cyan')
grid 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
int = max(fix(Pmax/4),0.5);
set(h,'ytick',[0:int:Pmax]) 
set(h,'yticklabel',[0:int:Pmax]) 
 
%ohmic power 
Pohm=data.gene.paddohm;Pohm=Pohm/1e6; 
%Radiated power 
Prad=data.gene.prad;Prad=Prad/1e6;
Pmax = ceil(max(max(Pohm),max(Prad))*1.1); 
if Pmax > 20
   Pens = [max(Pohm) max(Prad)];
   ind = find(Pens > 20);
   Pens(ind) = [];
   Pmax = ceil(max(Pens)*1.1);
end 
h=axes('position',[0.05 .94-2*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
p2=plot(tps,Pohm,'r',tps,Prad,'b',[tsamp tsamp],[0 Pmax],'--r'); 
set(p2,'linewidth',2); 
axis([tstart tend 0 Pmax]) 
t22=text(tstart+1,Pmax*0.9,'Pohm(MW)');set(t22,'color','red') 
t23=text(tstart+1,Pmax/2,'Prad(MW)');set(t23,'color','blue') 
grid 
 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
int = Pmax/4;
set(h,'ytick',[0:int:Pmax]) 
set(h,'yticklabel',[0:int:Pmax]) 
 
%volume averaged density 
nvol=data.gene.nemoy;
nvol=nvol/1e19;
ne=data.prof.ne;ne=ne/1e19;
 
h=axes('position',[0.05 .94-3*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
dmax=ceil(max(max(ne(:,1)),max(nvol))*1.1);
p3=plot(tps,ne(:,1),'r',tps,nvol,'b',[tsamp tsamp],[0 dmax],'--r'); 
set(p3,'linewidth',2); 
axis([tstart tend 0 dmax]) 
t31=text(tstart+1,dmax*0.8,'ne0(e19m-3)');set(t31,'color','red') 
t32=text(tstart+1,dmax*0.6,'<n>(e19m-3)');set(t32,'color','blue') 
xlabel('t(s)') 
grid 
int = 0:(dmax/5):dmax 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',int) 
set(h,'yticklabel',int) 


%major radius 
R=data.equi.raxe;
%minor radius 
a=data.equi.a;
 
h=axes('position',[0.05 .94-4*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
Rmean = fix(mean(data.geo.r0)*1000)/1000;
amean = fix(mean(data.equi.rhomax)*1000)/1000; 
p4=plot(tps,R(:,1)-Rmean,'r',tps,a(:,end)-amean,'g',[tsamp tsamp],[-.5 .5],'--r'); 
set(p4,'linewidth',2); 
v2=ceil(max(max(R(:,1)-Rmean),max(a(:,101)-amean))*10)/10;
if v2==0

  v2 = 0.1;

end
if length(R) > 10
  v1=ceil(min(min(R(10:(end-1),1)-Rmean),min(a(10:(end-1),101)-amean))*10)/10;
else
  v1=ceil(min(min(R(1:(end-1),1)-Rmean),min(a(1:(end-1),101)-amean))*10)/10;
end
if v1==0

  v1 = -0.1;

end
if v1 >= v2
  v2 = 2*v1;
end
axis([tstart tend v1 v2]) 
t41=text(tstart+1,v2*0.8,['R0-',num2str(Rmean,3),' (m)']);set(t41,'color','red') 
t42=text(tstart+1,v2*0.6,['a-',num2str(amean,3),' (m)']);set(t42,'color','green') 
grid 
int=fix((v2-v1)/4*100)/100;
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[v1:int:v2]) 
set(h,'yticklabel',[v1:int:v2]) 

 
%diamagnetic energy 
wdia=data.gene.wdia;wdia=wdia/1e6; 
%neutrons 
Rnt=data.source.totale.neutron.dt;
Rnd=data.source.totale.neutron.dd;

h=axes('position',[0.55 .94-hbox 0.4 hbox]);
set(h,'drawmode','fast') 
set(h,'box','on') 
echwd = max(wdia(2:(end-1)))*1.2; 
p5=plot(tps,wdia,'r',[tsamp tsamp],[0 echwd],'--r'); 
set(p5,'linewidth',2); 
axis([tstart tend 0 echwd]) 
t51=text(tstart+1,max(wdia),'Wdia(MJ)');set(t51,'color','red')
if isempty(param.fonction.coefb) 
  title(['modele : ',param.fonction.coefa])
else
  title(['modele : ',param.fonction.coefa,' + ',param.fonction.coefb])
end
grid 
 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:(echwd/4):echwd]) 
set(h,'yticklabel',round(100*[0:(echwd/4):echwd])/100) 

%betap 
betap=data.gene.betap; 
%betan 
betan=(4*betap.*Ipla)./(a(:,101).*Btor);
 
h=axes('position',[0.55 .94-2*hbox 0.4 hbox]); 
set(h,'drawmode','fast') 
set(h,'box','on') 
bem = ceil(max(max(betap(2:(end-1))),max(betan(2:(end-1))))*1.2); 
 
p6=plot(tps,betap,'r',tps,betan,'g',[tsamp tsamp],[0 bem],'--r'); 
set(p6,'linewidth',2);
axis([tstart tend 0 bem]) 
t61=text(tstart+1,bem*0.9,'betap');set(t61,'color','red') 
t62=text(tstart+1,bem*0.7,'betaN');set(t62,'color','green') 
grid 
 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:(bem/4):bem]) 
set(h,'yticklabel',[0:(bem/4):bem]) 
 

%Hfactors 
M=2

%scaling ITER L-mode - tau L thermique - NF39(1999)2206
tauL=0.023*(data.geo.r0).^1.89*M^0.2.*(data.cons.ip/1e6).^0.96.*(data.geo.b0).^0.03.*(data.equi.e(:,101)).^0.64.*...
(data.geo.a).^-0.06.*(data.gene.nbar./2e19./data.geo.a).^0.4.*(ploss_rec/1e6).^-0.73;
tauL(find(isnan(tauL)==1))=0;
tauL(find(iscomplex(tauL)==1))=0;

%scaling ITER H-mode - IPB98(y,2) - NF39(1999)2208
tauH=0.0562*(data.geo.r0).^1.97*M^0.19.*(data.cons.ip/1e6).^0.93.*(data.geo.b0).^0.15.*(data.equi.e(:,101)).^0.78.*...
(data.geo.a./data.geo.r0).^0.58.*(data.gene.nbar./2e19./data.geo.a).^0.41.*(ploss_rec/1e6).^-0.69;
tauH(find(isnan(tauH)==1))=0;
tauH(find(iscomplex(tauH)==1))=0;

%scaling ITER H-mode - ITER89-P - NF37(1997)1303
tauH2=0.038*(data.geo.r0).^1.5*M^0.5.*(data.cons.ip/1e6).^0.85.*(data.geo.b0).^0.2.*(data.equi.e(:,101)).^0.5.*...
(data.geo.a./data.geo.r0).^0.3.*(data.gene.nbar./2e19./data.geo.a).^0.1.*(ploss_rec/1e6).^-0.5;
tauH2(find(isnan(tauH2)==1))=0;
tauH2(find(iscomplex(tauH2)==1))=0;
[a,b]=butter(2,0.1);
tauH2 = filtfilt(a,b,real(tauH2));
tauH = filtfilt(a,b,real(tauH));
tauL = filtfilt(a,b,real(tauL));
taue = filtfilt(a,b,real(data.gene.taue));

 
h=axes('position',[0.55 .94-3*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
set(h,'box','on') 
 
p7=plot(tps,taue./tauL,'r',tps,taue./tauH,'k',tps,taue./tauH2,'b',[tsamp tsamp],[0 3],'--r'); 
set(p7,'linewidth',2); 
axis([tstart tend 0 3]) 
t71=text(tstart+1,2.5,'H-L');set(t71,'color','red') 
t72=text(tstart+1,2,'H98th');set(t72,'color','black') 
t73=text(tstart+1,1.5,'H89-P');set(t72,'color','blue') 
grid 
 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
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
 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:0.4:1.9]) 
set(h,'yticklabel',[0:.4:1.9]) 
 

%Te 
te=data.prof.te;te=te/1e3;  

h=axes('position',[0.05 .94-5*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
Tem = min(ceil(max(te(2:(end-1),1))*1.2),15);
ind50 = min(find(param.gene.x>0.5));
p9=plot(tps,te(:,1),'r',tps,te(:,ind50),'b',[tsamp tsamp],[0 Tem],'--r'); 
set(p9,'linewidth',2); 
axis([tstart tend 0 Tem]) 
t91=text(tstart+1,Tem*0.4,'Te0(keV)');set(t91,'color','red')   
t92=text(tstart+1,Tem*0.2,'Te(x=0.5)');set(t92,'color','blue')   
grid 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:(Tem/4):Tem]) 
set(h,'yticklabel',[0:(Tem/4):Tem]) 
 
 
%Ti 
ti=data.prof.ti;ti=ti/1e3;  

h=axes('position',[0.05 .94-6*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
 
Tim = min(ceil(max(ti(2:(end-1),1))*1.2),10);
pa=plot(tps,ti(:,1),'+r',[tsamp tsamp],[0 Tim],'--r'); 
set(pa,'linewidth',2); 
axis([tstart tend 0 Tim]) 
ta1=text(tstart+1,Tim*0.2,'Ti0(keV)');set(ta1,'color','red')  
grid 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[0:(Tim/4):Tim]) 
set(h,'yticklabel',[0:(Tim/4):Tim]) 

 
% 
 
h=axes('position',[0.05 .94-7*hbox 0.4 hbox]);
set(h,'drawmode','fast') 
pb=plot(tps,log10(Rnt(:,1)),'r',tps,log10(Rnd(:,1)),'b'); 
set(pb,'linewidth',2); 
valdt = log10(Rnt(2:(end-1),1));
valdn = log10(Rnd(2:(end-1),1));
valdt(~isfinite(valdt)) = 0;
valdn(~isfinite(valdn)) = 0;
fluxM = ceil(max(max(valdt)*1.2,max(valdn)*1.2));
if fluxM < 8 
  fluxM = 12;
end 
fluxm = 8;
ecflux = fluxM-fluxm;
axis([tstart tend fluxm fluxM]) 
pos=ceil(max(max(valdt)*1.2,max(valdn)*1.2))/1.1;
tb1=text(tstart+1,pos,'log(Rnt(0))');set(tb1,'color','red')
tb2=text(tstart+1,pos/1.1,'log(Rnd(0))');set(tb2,'color','blue')

grid 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'ytick',[fluxm:(ecflux/4):fluxM]) 
set(h,'yticklabel',round(100*[fluxm:(ecflux/4):fluxM])/100) 
xlabel('t(s)') 

 
%bootstrap current

iboot=data.gene.iboot/1e6;
icd=data.gene.icd/1e6;
ini=data.gene.ini/1e6;
 
h=axes('position',[0.55 .94-5*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
pc=plot(tps,Ipla,'r',tps,iboot,'g',tps,icd,'b',tps,ini,'k',[tsamp tsamp],[0 max(Ipla)*1.1],'--r'); 
set(pc,'linewidth',2); 
axis([tstart tend 0 max(Ipla)*1.2]) 
tc1=text(tstart+1,max(Ipla),'Ip(MA)');set(tc1,'color','red'); 
tc2=text(tstart+1,max(Ipla)/1.5,'Iboot(MA)');set(tc2,'color','green');
tc3=text(tstart+1,max(Ipla)/2.25,'ICD(MA)');set(tc3,'color','blue'); 
tc4=text(tstart+1,max(Ipla)/4,'INI(MA)');set(tc3,'color','black'); 
 
grid 

set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
int = ceil(max(Ipla)*1.2)/5
set(h,'ytick',[0:int:ceil(max(Ipla)*1.2)]) 
set(h,'yticklabel',[0:int:ceil(max(Ipla)*1.2)]) 
 

 
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
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
set(h,'xticklabel',[]) 
set(h,'ytick',[1:5]) 
set(h,'yticklabel',[1:5]) 

 
%Zeff 
zeff=data.impur.zeff; 
zefm=max(ceil(max(zeff(:,1))*1.5),2);
h=axes('position',[0.55 .94-7*hbox 0.4 hbox]);
set(h,'drawmode','fast'); 
set(h,'box','on') 
 
pe=plot(tps,zeff(:,1),'r',[tsamp tsamp],[0 zefm],'--r'); 
set(pe,'linewidth',2); 
axis([tstart tend 0 zefm]) 
te1=text(tstart+1,zefm*0.8,'Zeff(0)');set(te1,'color','red'); 
 
 
grid 
set(h,'xtick',[floor(tstart):dti:ceil(tend)]) 
%set(h,'xticklabel',[]) 
set(h,'ytick',[1:zefm]) 
set(h,'yticklabel',[1:zefm]) 
xlabel('t(s)')
