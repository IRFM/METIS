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
indkk = param.gene.kmin:param.gene.k
t = data.gene.temps(indkk);
netot = data.gene.netot(indkk);
dnetotdt = cat(1,0,diff(netot)./diff(t));
neout = data.prof.flux.sortant.ge(indkk,end);
nein  = data.gene.sne_bord(indkk) + data.gene.sne_idn(indkk);
neinb = data.gene.sne_bord(indkk);
neinm = data.bord.fluxmur_c(indkk)+ data.bord.fluxmur_f(indkk);
bilane = data.gene.bilan_e(indkk);
neinint  = cumtrapz(t,nein,1);
neoutint  = cumtrapz(t,neout,1);

h = findobj(0,'type','figure','tag','bilanmat');
if isempty(h)
       h=figure('tag','bilanmat');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',0.5,'color',[1 1 1])

subplot(5,1,1)
plot(t,netot)
ylabel('Netot (electrons)')

st =sprintf('choc %s@%8.1f , t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine,...
             param.from.shot.num,...   
             param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);

subplot(5,1,2)
plot(t,dnetotdt,'or',t,neout,'+b',t,nein,'xg',t,bilane,'k')
legend('dNe_t_o_t/dt','Ne_o_u_t','Ne_i_n','bilan_e');
ylabel('electrons/s')

subplot(5,1,3)
plot(t,(dnetotdt+neout-nein) ./ netot)
title(num2str(max(abs(dnetotdt+neout-nein) ./ netot)))
ylabel('bilan ponctuel')
max(abs(dnetotdt+neout-nein) ./ netot)

subplot(5,1,4)
plot(t,neinb - neinm)
legend('int(source)- flux du mur');
ylabel('electrons/s')

subplot(5,1,5)
plot(t,neinb - neinm)
plot(t,(netot+neoutint-neinint)./netot-1)
ylabel('bilan integral')
xlabel('temps (s)')

