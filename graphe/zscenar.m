%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
h = findobj(0,'type','figure','tag','scenario');
if isempty(h)
       h=figure('tag','scenario');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
ind = (param.gene.kmin+1):param.gene.k;
indt = find(param.gene.file=='_');
param.gene.file(indt) = '-';


 
wrlw = zrlw(data.gene.ip(ind)/1e6,data.geo.r0(ind),data.equi.rhomax(ind),data.gene.nemoy(ind)/1e19,...
            data.geo.b0(ind),data.gene.zeffm(ind),smooth(data.gene.ploss(ind),5)/1e6,data.geo.e1(ind))*1e6; 

Hrlw = data.gene.we(ind)./wrlw;


subplot(5,1,1)
if mean(abs(data.cons.fce(:))) > 0
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,1))/1e6,'b')
  hold on
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,2))/1e6,'r+')
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,3))/1e6,'k--')
  ylabel('P_F_C_E (MW)')
  hold off
  legend([' ant 1, phi =',num2str(param.cons.fce.angle_tor(1)),' °'],...
         [' ant 2, phi=',num2str(param.cons.fce.angle_tor(2)),' °'],...
       [' ant 3, phi=',num2str(param.cons.fce.angle_tor(2)),' °'])
  title([' simu : ',param.gene.file(22:end)])
end
if mean(abs(data.cons.fci(:))) > 0
  plot(data.gene.temps(ind),sum(abs(data.cons.fci(ind,:)'))/1e6,'b')
  ylabel('P_F_C_i (MW)')
  hold off
  legend([' ant. freq=',num2str(param.cons.fci.frequence(1)),' MHz'])
  betapmax = max(data.gene.betap(ind));
  frboot   = max(smooth(data.gene.iboot(ind),1)./smooth(data.gene.ip(ind),1))*100;
  title([' choc FWEH, densité type 18805, bp =',num2str(betapmax,2),' fr_b_o_o_t=',num2str(frboot,2),' %, H_r_l_w= ',num2str(max(Hrlw))])
end

subplot(5,1,2)
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,'b')
ylabel('P_L_H (MW)')
subplot(5,1,3)
plot(data.gene.temps(ind),data.gene.nbar(ind)/1e19,'b',...
     data.gene.temps(ind),data.gene.temoy(ind)/1e3,'r',...
     data.gene.temps(ind),Hrlw,'k--')
legend('nbar 10^1^9 m^-^3','<Te> keV','H_r_l_w')
subplot(5,1,4)
plot(data.gene.temps(ind),data.gene.ip(ind)/1e6,'b',...
     data.gene.temps(ind),data.gene.ini(ind)/1e6,'r',...
     data.gene.temps(ind),data.gene.iboot(ind)/1e6,'g')
     
legend('Ip','Ini','Iboot')
ylabel('MA')
ind1=15;
ind2=30;
ind3=40;
subplot(5,3,13)

plot(param.gene.x,data.prof.q(iround(data.gene.temps,ind1),:),'b',...
     param.gene.x,data.prof.q(iround(data.gene.temps,ind2),:),'r',...
     param.gene.x,data.prof.q(iround(data.gene.temps,ind3),:),'k')

ylabel('q')
subplot(5,3,14)

plot(param.gene.x,data.prof.te(iround(data.gene.temps,ind1),:)/1e3,'b',...
     param.gene.x,data.prof.te(iround(data.gene.temps,ind2),:)/1e3,'r',...
     param.gene.x,data.prof.te(iround(data.gene.temps,ind3),:)/1e3,'k')

ylabel('Te (keV)')
subplot(5,3,15)

plot(param.gene.x,data.prof.ne(iround(data.gene.temps,ind1),:)/1e19,'b',...
     param.gene.x,data.prof.ne(iround(data.gene.temps,ind2),:)/1e19,'r',...
     param.gene.x,data.prof.ne(iround(data.gene.temps,ind3),:)/1e19,'k')
legend(['t=',num2str(ind1,2),' s'],['t=',num2str(ind2,2),' s'],['t=',num2str(ind3,2),' s']);
ylabel('ne (10^1^9 m^-^3')

h = findobj(0,'type','figure','tag','courant');
if isempty(h)
       h=figure('tag','courant');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

plotprof(gca,data.gene.temps(ind),param.gene.x,data.source.hyb.j(ind,:),'color','r')
plotprof(gca,data.gene.temps(ind),param.gene.x,data.source.fce.j(ind,:),'color','b')
plotprof(gca,data.gene.temps(ind),param.gene.x,data.source.totale.j(ind,:),'color','k')
plotprof(gca,data.gene.temps(ind),param.gene.x,data.prof.jmoy(ind,:),'color','m','linestyle','o')
ylabel('A/m^2')
legend('LH','FCE','NI','Jmoy')

