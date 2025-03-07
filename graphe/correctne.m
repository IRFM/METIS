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
for k=1:length(data.gene.temps)

  [pa,pq]=piquage(param.gene.x,data.prof.ni(k,:)',0,100);
  newni(k,:) = (data.prof.ni(k,1)-data.prof.ni(k,end)) * ...
               (1-param.gene.x.^(2*pa)).^pq + max(data.prof.ni(k,end),data.prof.ni(k,1)/20);
end
