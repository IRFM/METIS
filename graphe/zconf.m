temps=data.gene.temps;
x = param.gene.x;
Ipf = data.gene.ip/1e6;
bt = data.geo.b0;
nbarf=data.gene.nbar/1e19;
[taue,ttaue]=tsbase(fix(param.from.shot.num),'staudia');
inder = find(diff(taue) <= 0);
taue(inder)=[];
ttaue(inder)=[];
taue = interp1(ttaue,taue,data.gene.temps);
%
%  plth = puissance totale absorbee (on enleve la puissance perdue dans le ripple)
%
plth=data.gene.ploss/1e6;
plth=(data.gene.paddfci+data.gene.paddhyb*0.95+data.gene.paddohm)/1e6;
prad = data.gene.prad/3/1e6;
dwdt = gradient(temps,data.gene.wdia);
tauth = data.gene.wth./plth/1e6;
ntauth = data.gene.wth./plth/1e6;
R0=data.geo.r0;
af=data.geo.a;
%
% plth1 = puissance inject�e
%
plth1 = sum(abs(data.cons.fci'))'+sum(abs(data.cons.hyb'))'+data.gene.paddohm;
plth1=plth1/1e6;
if strcmp(deblank(param.cons.fci.mode{1}),'FWEH')
    % pour plaier au probleme d'absorb, VB, 24 aout 2005
    plth = plth1;
    tauth = data.gene.wth./plth/1e6;
end
    
tauth1 = data.gene.wth./plth1/1e6;
ntauth1 = data.gene.wth./plth1/1e6;
temps=data.gene.temps;
tauL=zeros(size(temps));
tauL1=tauL;
clear totsup totel totion
[a,b]=butter(2,0.05);    
wdia=data.gene.wdia;                
nwdia=filtfilt(a,b,wdia);
dwdt=gradient(filtfilt(a,b,wdia),temps)/1e6;
temps1=temps(param.gene.kmin:param.gene.k);
for k=param.gene.kmin:param.gene.k
  Pel   = data.prof.pe(k,:);
  Pion  = data.prof.pion(k,:);
  Psup  = data.source.fci.psupra(k,:);
  totsup(k) =  zintvol(Psup,x,data.equi.vpr(k,:),data.equi.rhomax(k))/1e6;
  totion(k) =  zintvol(Pion,x,data.equi.vpr(k,:),data.equi.rhomax(k))/1e6;
  totel(k)  =  zintvol(Pel,x,data.equi.vpr(k,:),data.equi.rhomax(k))/1e6;
  if ~isnan(dwdt(k))
      pcalc = plth(k)-prad(k)-dwdt(k);
      pcalc1 = plth1(k)-prad(k)-dwdt(k);
  else
      pcalc=plth(k)-prad(k);
      pcalc1=plth1(k)-prad(k);
  end
  tauL(k) = H97(Ipf(k),bt(k),nbarf(k),pcalc,R0(k),1,af(k),2);
  tauL1(k)=H97(Ipf(k),bt(k),nbarf(k),pcalc1,R0(k),1,af(k),2);

end
oui=0;
if oui
  if ~exist('tauf','var')
    [tauf,ttauf]=tsbase(33951,'staudia');
  end
%
% tauth1 : taudia mesur� � partir de Wdia et ploss
%
  tauth1 = interp1(ttauf,tauf,data.gene.temps);
%
% tauth : taudia corrig� des pertes ripple
%
  tauth = tauth1 .* plth1 ./ plth;
end
h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
ver = (plth < 15 & plth > 1) | tauL<1;
ver(ver==0)=-1;
sup = abs(1-data.gene.wth./data.gene.wdia);
sup(data.gene.paddfci<1e5)=0;
sup(sup>0.5)=0;
rap=tauth./tauL;
rap1=tauth1./tauL1;
ver(rap>2)=-1;
ver(rap<0.1)=-1;
subplot(3,1,1)
plot(temps,rap.*ver,temps,(tauth.*(1-sup)./tauL).*ver,'--')
axis([param.gene.tdeb param.gene.t 0 2])
grid
title(['shot ',int2str(param.from.shot.num),' Pabs = Pinj-Prip'])
ylabel('H96')
legend('no sup.','sup.')
hold on
subplot(3,1,2)
plot(temps1,totel,'r',temps1,totion,'b.-',temps1,totsup,'k--')
legend('Pe','Pion','Psup')
ylabel('Pa')
axis([param.gene.tdeb param.gene.t 0 inf])
hold off
grid

title('pressure')


subplot(3,1,3)
frrip=(sum(abs(data.cons.fci'))'-data.gene.paddfci)./sum(abs(data.cons.fci'))'*100;
frrip(frrip>50) = 0;
frrip(frrip<0) = 0;
plot(temps,frrip,temps,sup*100)
legend('rippl','sup')
ylabel('%')
axis([param.gene.tdeb param.gene.t 0 40])
grid
xlabel('time (s)')
