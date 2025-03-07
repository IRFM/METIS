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

h = findobj(0,'type','figure','tag','asserqpol');
if isempty(h)
       h=figure('tag','asserqpol');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
  x = param.gene.x;
  parametre = param.cons.asser;
  indq = fix(interp1(param.gene.x,1:length(param.gene.x),parametre.qpos,'nearest'));
  tempsvrai  = data.gene.temps;
  tempsvrai(data.mode.asser<2) = 0;
  tempsvrai(tempsvrai > 0) = tempsvrai(tempsvrai > 0) - min(tempsvrai(tempsvrai > 0));

  eq        = ones(size(tempsvrai))*parametre.qref - data.prof.q(:,indq);

  inte      = cumtrapz(tempsvrai,eq);
  phyb      = parametre.pmin + sum( (ones(size(tempsvrai))*parametre.K) .* (parametre.gp .* (eq + parametre.gi .* inte)),2);


  subplot(2,2,1)
  plotprof(gca,data.gene.temps,x,data.prof.q,'color','r');
  hold on
  plot(param.cons.asser.qpos,param.cons.asser.qref,'ro');
  hold off
  subplot(2,2,2)
  plot(data.gene.temps,phyb,'color','r');
  
  subplot(2,5,6)
  
  plot(data.gene.temps,param.cons.asser.qref(1)*ones(size(data.gene.temps)),'bo')
  hold on
  plot(data.gene.temps,data.prof.q(:,min(find(x>param.cons.asser.qpos(1)))),'b')
  hold off  
  
  subplot(2,5,7)
  
  plot(data.gene.temps,param.cons.asser.qref(2)*ones(size(data.gene.temps)),'bo')
  hold on
  plot(data.gene.temps,data.prof.q(:,min(find(x>param.cons.asser.qpos(2)))),'b')
  hold off
  
  subplot(2,5,8)
  
  plot(data.gene.temps,param.cons.asser.qref(3)*ones(size(data.gene.temps)),'bo')
  hold on
  plot(data.gene.temps,data.prof.q(:,min(find(x>param.cons.asser.qpos(3)))),'b')
  hold off  
  
  subplot(2,5,9)
  
  plot(data.gene.temps,param.cons.asser.qref(4)*ones(size(data.gene.temps)),'bo')
  hold on
  plot(data.gene.temps,data.prof.q(:,min(find(x>param.cons.asser.qpos(4)))),'b')
  hold off

  subplot(2,5,10)
  
  plot(data.gene.temps,param.cons.asser.qref(5)*ones(size(data.gene.temps)),'bo')
  hold on
  plot(data.gene.temps,data.prof.q(:,min(find(x>param.cons.asser.qpos(5)))),'b')
  hold off
  
