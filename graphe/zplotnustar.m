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
% version  3.0  du  9/11/2005
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  

h = findobj(0,'type','figure','tag','nustar');
if isempty(h)
       h=figure('tag','nustar');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])


list_opt = [];
opt_def = 1;
if ~exist('jeux1')
 jeux1=[];
end
chargedata
callback='nurhostar(data,param,jeux1,temps_zslider)';
if param.gene.k == param.gene.kmin

  fin = param.gene.kmin+1;

else

  fin = param.gene.k;

end
zslider(h,data.gene.temps(param.gene.kmin:fin),[],[],callback)
