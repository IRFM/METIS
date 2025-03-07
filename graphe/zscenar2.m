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


 
wrlw  = zrlw(data.gene.ip(ind)/1e6,data.geo.r0(ind),data.equi.rhomax(ind),data.gene.nemoy(ind)/1e19,...
            data.geo.b0(ind),data.gene.zeffm(ind),smooth(data.gene.ploss(ind),5)/1e6,data.geo.e1(ind))*1e6; 

Hrlw  = data.gene.we(ind)./wrlw;
temps = data.gene.temps(ind);
Pfce1 = abs(data.cons.fce(ind,1))/1e6;
Pfce2 = abs(data.cons.fce(ind,2))/1e6;
Pfce3 = abs(data.cons.fce(ind,3))/1e6;
Pfci  = abs(sum(data.cons.fci(ind,:)')')/1e6;
Plh   = data.gene.paddhyb(ind)/1e6;
nbar  = data.gene.nbar(ind)/1e19;
Temoy = data.gene.temoy(ind)/1e3;
Temax = data.prof.te(ind,1)/1e3;
p    = data.gene.ip(ind)/1e6;
Ini   = data.gene.ini(ind)/1e6;
Iboot = data.gene.iboot(ind)/1e6;
subplot(5,1,1)
if mean(Pfce1+Pfce2+Pfce3) > 0
  plot(temps,Pfce1,'b')
  hold on
  plot(temps,Pfce2,'r+')
  plot(temps,Pfce3,'k--')
  ylabel('P_F_C_E (MW)')
  hold off
  legend([' ant 1, phi =',num2str(param.cons.fce.angle_tor(1)),' °'],...
         [' ant 2, phi=',num2str(param.cons.fce.angle_tor(2)),' °'],...
       [' ant 3, phi=',num2str(param.cons.fce.angle_tor(2)),' °'])
  title([' simu : ',param.gene.file(22:end)])
end
if mean(Pfci) > 0
  plot(temps,Pfci,'b')
  ylabel('P_F_C_i (MW)')
  hold off
  legend([' ant. freq=',num2str(param.cons.fci.frequence(1)),' MHz'])
  betapmax = max(data.gene.betap(ind));
  frboot   = max(smooth(data.gene.iboot(ind),1)./smooth(data.gene.ip(ind),1))*100;
  title([' choc FWEH, densité type 18805, bp =',num2str(betapmax,2),' fr_b_o_o_t=',num2str(frboot,2),' %, H_r_l_w= ',num2str(max(Hrlw))])
end
if mean(Pfce1+Pfce2+Pfce3) > 0 & mean(Pfci) > 0
  plot(temps,(Pfce1+Pfce2+Pfce3)*10,'r')
  hold on 
  plot(temps,Pfci,'b')
  ylabel('P_i_n_j (MW)')
  legend([' P*10, ant. ECRH, phi =',num2str(param.cons.fce.angle_tor(1)),' °'],...
         [' ant. ICRH, freq=',num2str(param.cons.fci.frequence(1)),' MHz'])
  hold off
end
subplot(5,1,2)
plot(temps,Plh,'b')
ylabel('P_L_H (MW)')
subplot(5,1,3)
plot(temps,nbar,'b',...
     temps,Temax,'r',...
     temps,Hrlw,'k--')
legend('nbar 10^1^9 m^-^3','Te(0) keV','H_r_l_w')
subplot(5,1,4)
plot(temps,Ip,'b',...
     temps,Ini,'r',...
     temps,Iboot,'g')
     
legend('Ip','Ini','Iboot')
ylabel('MA')
long  = length(temps);

ind1  = long/5;
ind2  = long/2;
ind3  = long/1.2;
q     = data.equi.q(ind,:);
Te    = data.prof.te(ind,:)/1e3;
ne    = data.prof.ne(ind,:)/1e19;
Jmoy  = data.prof.jmoy(ind,:)/1e6;
Jmoy1 = Jmoy(ind1,:)';
Jmoy2 = Jmoy(ind2,:)';
Jmoy3 = Jmoy(ind3,:)';

Jlh   = data.source.hyb.j(ind,:)/1e6;
Jlh1  = Jlh(ind1,:)';
Jlh2  = Jlh(ind2,:)';
Jlh3  = Jlh(ind3,:)';

Jfce  = data.source.fce.j(ind,:)/1e6;
Jfce1 = Jfce(ind1,:)';
Jfce2 = Jfce(ind2,:)';
Jfce3 = Jfce(ind3,:)';


q1    = q(ind1,:)';
q2    = q(ind2,:)';
q3    = q(ind3,:)';
Te1   = Te(ind1,:)';
Te2   = Te(ind2,:)';
Te3   = Te(ind3,:)';
ne1   = ne(ind1,:)';
ne2   = ne(ind2,:)';
ne3   = ne(ind3,:)';
x     = param.gene.x';
tempsana = [temps(ind1) temps(ind2) temps(ind3)]';

saveigor = 0;
saveigor = inputd('sauvegarde pour fichier igor [1 -> oui]',saveigor);
if saveigor == 1
  nom    = input('nom de la sauvegarde [defaut : cronos_igor]','s');
  if isempty(nom)
    nom  = 'cronos_igor';
  end
  eval(['save ',nom,' temps tempsana Plh Pfce1 Pfce2 Pfce3 Pfci Ip Ini Iboot Hrlw nbar Temoy x Jmoy1 Jmoy2 Jmoy3 Jlh1 Jlh2 Jlh3 Jfce1 Jfce2 Jfce3 q1 q2 q3 Te1 Te2 Te3 ne1 ne2 ne3'])

end
subplot(5,3,13)

plot(x,q1,'b',...
     x,q2,'r',...
     x,q3,'k')

ylabel('q')
subplot(5,3,14)

plot(x,Te1,'b',...
     x,Te2,'r',...
     x,Te3,'k')

ylabel('Te (keV)')
subplot(5,3,15)

plot(x,ne1,'b',...
     x,ne2,'r',...
     x,ne3,'k')
legend(['t=',num2str(temps(ind1),2),' s'],['t=',num2str(temps(ind2),2),' s'],['t=',num2str(temps(ind3),2),' s']);
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

