%  NESTAT  courte description  
%------------------------------------------------------------------------------- 
% fichier :  nestat.m  ->  nestat 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   dndx = nestat(x,ne,void,s,d,v,xv,rm,vpr,grho2) 
%  
% entrees :  
%  x     = 
%  ne    = 
%  void  = 
%  s     = 
%  d     = 
%  v     = 
%  xv    = 
%  rm    = 
%  vpr   = 
%  grho2 = 
%  
% sorties :  
%   dndx  = 
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
function dndx = nestat(x,ne,void,s,d,v,xv,rm,vpr,grho2)

f  = rm .* v./d;
g  = rm .^ 2 .* cumtrapz(xv,vpr .* s) ./ max(eps,vpr) ./grho2 ./ d;
  
ind  = fix(interp1(xv,1:length(xv),x,'nereast'));
if isempty(ind)
   if x < min(xv)
      ind = 1;
   else
      ind = length(xv);
   end
end
dndx = - ne .* f(ind) - g(ind);
