%
% ce script trace le bilan de puissance 
%
h = findobj(0,'type','figure','tag','power-b');
if isempty(h)
       h=figure('tag','power-b');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1],'defaultlinemarkersize',4)
%
% nombre de subplot
k = 6;
% comptage subplot
l =1;
% le temps
t = data.gene.temps;
%
% pour tot,fus ohm
%
subplot(k,1,l); l = l + 1;
idn  = sum(real(data.cons.hyb),2);
plot(t,data.gene.paddohm/1e6,'b',t,data.gene.paddtot/1e6,'r',t,data.gene.paddfus/1e6,'k');
text('units','normalized','position',[1.01 0.75],'string','Ohm','color','b')
text('units','normalized','position',[1.01 0.5],'string','total','color','r')
text('units','normalized','position',[1.01 0.25],'string','Fus','color','k')
title(sprintf('%s, #%d, Injected power checking',param.from.machine,param.from.shot.num))  
ylabel('MW') 
set(gca,'xticklabel','');
%
% pour fci
%
subplot(k,1,l); l = l + 1;
fci  = sum(abs(data.cons.fci),2);
plot(t,fci/1e6,'b',t,data.gene.paddfci/1e6,'or');
text('units','normalized','position',[1.01 0.25],'string','ICRH','color','r')
text('units','normalized','position',[1.01 0.75],'string','ICRH_c_o_n_s','color','b')
ylabel('MW') 
set(gca,'xticklabel','');

%
% pour fce
%
subplot(k,1,l); l = l + 1;
fce  = sum(abs(data.cons.fce),2);
plot(t,fce/1e6,'b',t,data.gene.paddfce/1e6,'or');
text('units','normalized','position',[1.01 0.25],'string','ECRH','color','r')
text('units','normalized','position',[1.01 0.75],'string','ECRH_c_o_n_s','color','b')
ylabel('MW') 
set(gca,'xticklabel','');

%
% pour hybride
%
subplot(k,1,l); l = l + 1;
hyb  = sum(abs(data.cons.hyb),2);
plot(t,hyb/1e6,'b',t,data.gene.paddhyb/1e6,'or');
text('units','normalized','position',[1.01 0.25],'string','LH','color','r')
text('units','normalized','position',[1.01 0.75],'string','LH_c_o_n_s','color','b')
ylabel('MW') 
set(gca,'xticklabel','');

%
% pour idn
%
subplot(k,1,l); l = l + 1;
if strcmp(param.from.machine,'JET') 
  user     = getenv('HOME');
  if strcmp(user,'/usr/drfc/cgc')
    chemin = strcat(user,'/cgc_data/jet/data/');
  else
    chemin =  strcat(user,'zineb/data/jet/');
  end
  chemin   = strcat(chemin,int2str(param.from.shot.num),'/transp',int2str(param.from.shot.num),'.mat');
  if exist(chemin,'file')
    transp = load(chemin);
	 if ~isempty(transp.qbitransp)
	   nbitransp = transp.qbitransp + transp.qbetransp;
		nbitemps  = transp.ttransp;
	 else
	   nbitransp = [];
		nbitemps  = [];
	 end
  else
	 nbitransp = [];
    nbitemps  = [];
  end
else
  nbitransp = [];
  nbitemps  = [];
end
idn  = sum(real(data.cons.idn),2);
plot(t,idn/1e6,'b',t,data.gene.paddidn/1e6,'or',nbitemps,nbitransp,'m--');
text('units','normalized','position',[1.01 0.25],'string','NBI','color','r')
text('units','normalized','position',[1.01 0.5],'string','NBI_c_o_n_s','color','b')
if ~isempty(nbitransp)
  text('units','normalized','position',[1.01 0.75],'string','NBI_t_r_a_n_s','color','m')
end
ylabel('MW') 
set(gca,'xticklabel','');

%
% pour ext
%
subplot(k,1,l); l = l + 1;
ext   = sum(abs(data.cons.ext),2);
exti  = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .*  ...
                                (data.source.ext.el + data.source.ext.ion),2);
plot(t,ext/1e6,'b',t,exti/1e6,'or');
text('units','normalized','position',[1.01 0.25],'string','External','color','r')
text('units','normalized','position',[1.01 0.75],'string','Ext_c_o_n_s','color','b')
ylabel('MW') 

% apres le dernier plot
xlabel('time (s)')


