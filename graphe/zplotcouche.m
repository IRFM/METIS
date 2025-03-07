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
% fonction ecrite par V. Basiuk , poste 61-26  
% version  2.2  du  25/03/2004  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  

h = findobj(0,'type','figure','tag','couche FCI');
if isempty(h)
       h=figure('tag','couche FCI');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
freq  = param.cons.fci.frequence;
for k=1:param.nombre.fci
  list_opt{k} = sprintf('ant. %d (%g MHz)',k,freq(k));
end
opt_def = 1;
callback='plotcouche3(data,param,temps_zslider,option_zslider)';
if param.gene.k == param.gene.kmin

  fin = param.gene.kmin+1;

else

  fin = param.gene.k;

end
zslider(h,data.gene.temps(param.gene.kmin:fin),list_opt,opt_def,callback)
