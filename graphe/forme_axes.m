%  FORME_AXES  courte description  
%------------------------------------------------------------------------------- 
% fichier :  forme_axes.m  ->  forme_axes 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   forme_axes(ha) 
%  
% entrees :  
%  ha = 
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
function forme_axes(ha)

lim = get(ha,'xlim');
xmin = 1e300;
xmax = -1e300;
for k = 1:length(lim)
   lc   = lim{k}; 
   if all(isfinite(lc))
      xmin = min(xmin,lc(1));
      xmax = max(xmax,lc(2));
   end
end
set(ha,'xlim',[xmin,xmax]);
set(ha(1:(end-1)),'xticklabel','');

for k = 1:(length(ha)-1)
      lab = get(ha(k),'yticklabel');
      if ~isempty(lab)
         lab(1,:) = char(32*ones(size(lab(1,:))));
         set(ha(k),'yticklabel',lab,'ytickmode','manual');
      end
end
      
