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
set(h,'defaultaxesfontsize',10,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])
ind = (param.gene.kmin+1):param.gene.k;
indt = find(param.gene.file=='_');
param.gene.file(indt) = '-';


 
wrlw = zrlw(data.gene.ip(ind)/1e6,data.geo.r0(ind),data.equi.rhomax(ind),data.gene.nemoy(ind)/1e19,...
            data.geo.b0(ind),data.gene.zeffm(ind),smooth(data.gene.ploss(ind),5)/1e6,data.geo.e1(ind))*1e6; 

Hrlw = data.gene.we(ind)./wrlw;
frboot = data.gene.iboot(ind)./data.gene.ip(ind)*100;


if mean(abs(data.cons.fce(:))) > 0
  subplot(6,1,1)
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,1))/1e6,'b')
  hold on
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,2))/1e6,'r+')
  plot(data.gene.temps(ind),abs(data.cons.fce(ind,3))/1e6,'k--')
  ylabel('P_F_C_E (MW)')
  hold off
  legend([' ant 1, phi =',num2str(param.cons.fce.angle_tor(1)),' °'],...
         [' ant 2, phi=',num2str(param.cons.fce.angle_tor(2)),' °'],...
       [' ant 3, phi=',num2str(param.cons.fce.angle_tor(2)),' °'])
  title([' simu : ',param.gene.file(22:end),' Ip=',num2str(max(data.gene.ip)/1e6,2),' MA'])
  grid
end
if mean(abs(data.cons.fci(:))) > 0
  subplot(6,1,2)
  plot(data.gene.temps(ind),sum(abs(data.cons.fci(ind,:)'))/1e6,'b')
  ylabel('P_F_C_i (MW)')
  hold off
  legend([' ant. freq=',num2str(param.cons.fci.frequence(1)),' MHz'])
  betapmax = max(data.gene.betap(ind));
  grid
end

subplot(6,1,3)
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,'b')
ylabel('P_L_H (MW)')
grid
subplot(6,1,4)
plot(data.gene.temps(ind),data.gene.nbar(ind)/1e19,'b',...
     data.gene.temps(ind),data.prof.te(ind,1)/1e3,'r')
legend('nbar 10^1^9 m^-^3','Te(0) keV')
grid
subplot(6,1,5)
plot(data.gene.temps(ind),Hrlw,'b-')
ylabel('H_r_l_w')
grid
subplot(6,1,6)
plot(data.gene.temps(ind),frboot,'r-')
ylabel('%')
legend('fr_b_o_o_t')
grid
