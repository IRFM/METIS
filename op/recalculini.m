%  RECALCULINI  courte description  
%------------------------------------------------------------------------------- 
% fichier :  recalculini.m  ->  recalculini ,  zintsurf 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   data=recalculini(data,param) 
%  
% entrees :  
%  data  = 
%  param = 
%  
% sorties :  
%   data = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function data=recalculini(data,param)
for k=1:length(data.gene.temps)
  data.source.totale.j(k,1)=data.source.totale.j(k,2);
  data.source.ext.j(k,1)=0;
  
  data.gene.ini(k)       = zintsurf(data.source.totale.j(k,:),param.gene.x,data.equi.spr(k,:),data.equi.rhomax(k,:));

end
function s=zintsurf(e,x,sp,rhomax)   

    s = rhomax .* trapz(x,sp .* e,2);
  
