
Ipf = data.gene.ip/1e6;
bt = data.geo.b0;
nbarf=data.gene.nbar/1e19;
%
%  plth = puissance totale absorbée (on enleve la puissance perdue dans le ripple) 
%
plth=data.gene.ploss/1e6;
plth=(data.gene.paddfci+data.gene.paddhyb*0.95+data.gene.paddohm)/1e6;
tauth = data.gene.wth./plth/1e6;
ntauth = data.gene.wth./plth/1e6;
R0=data.geo.r0;
af=data.geo.a;
%
% plth1 = puissance injectée
%
plth1 = sum(abs(data.cons.fci'))'+sum(abs(data.cons.hyb'))'+data.gene.paddohm;
plth1=plth1/1e6;
tauth1 = data.gene.wth./plth1/1e6;
ntauth1 = data.gene.wth./plth1/1e6;
temps=data.gene.temps;
for k=param.gene.kmin:param.gene.k

  tauL(k) = H97(Ipf(k),bt(k),nbarf(k),plth(k),R0(k),1,af(k),1);
  tauL1(k)=H97(Ipf(k),bt(k),nbarf(k),plth1(k),R0(k),1,af(k),1);

end
oui=1;
if oui
  if ~exist('tauf','var')
    [tauf,ttauf]=tsbase(33951,'staudia');
  end
%
% tauth1 : taudia mesuré à partir de Wdia et ploss
%
  tauth1 = interp1(ttauf,tauf,data.gene.temps);
%
% tauth : taudia corrigé des pertes ripple
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
ver = (plth < 15 & plth > 1) | tauL'<1;
sup = abs(data.gene.wdia-data.gene.wth)./data.gene.wth;
sup(data.gene.paddfci<1e5)=0;
sup(sup>0.5)=0;
rap=tauth./tauL';
rap1=tauth1./tauL1';

ind = find(ver~=0 & rap<2 & rap>0.1);

subplot(3,1,1)
plot(temps(ind),rap(ind),temps(ind),(tauth(ind).*(1-sup(ind))./tauL(ind)'),'--',...
      temps(ind),ntauth(ind)./tauL(ind)','+')
axis([param.gene.tdeb param.gene.t 1 2])
grid
title(['shot ',int2str(param.from.shot.num),' Pabs = Pinj-Prip'])
ylabel('H96')
legend('exp no sup.','exp sup.','cronos')
hold on
subplot(3,1,2)
plot(temps(ind),rap1(ind),temps(ind),(tauth1(ind).*(1-sup(ind))./tauL1(ind)'),'--',...
      temps(ind),ntauth1(ind)./tauL1(ind)','+')
axis([param.gene.tdeb param.gene.t 1 2])
hold off
grid
xlabel('time (s)')
ylabel('H96')
title('Pabs=Pinj')
subplot(3,2,5)

plot(temps,sup*100)
ylabel('frsup%')
axis([param.gene.tdeb param.gene.t 0 40])


subplot(3,2,6)
frrip=(sum(abs(data.cons.fci'))'-data.gene.paddfci)./sum(abs(data.cons.fci'))'*100;
frrip(frrip>50) = 0;
frrip(frrip<0) = 0;
plot(temps,frrip)
axis([param.gene.tdeb param.gene.t 0 40])

ylabel('frrip %')

