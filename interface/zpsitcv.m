%  ZPSITCV  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zpsitcv.m  ->  zpsitcv 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [volloc,psiloc]=zpsitcv(rpsi,zpsi,psi) 
%  
% entrees :  
%  rpsi = 
%  zpsi = 
%  psi  = 
%  
% sorties :  
%  volloc = 
%  psiloc = 
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
function [volloc,psiloc]=zpsitcv(rpsi,zpsi,psi)

[n,m,lon]=size(psi);

volloc = zeros(lon,41);

for k=1:lon
 psiloc(k,:)=linspace(eps,max(max(psi(:,:,k))),41);
 disp(['reste :',int2str(lon-k)])
 cc=contour(rpsi,zpsi,squeeze(psi(:,:,k))',psiloc(k,:));
 for l=1:40
   n=cc(2,1);
   r = cc(1,2:(n+1));
   z= cc(2,2:(n+1));
   cc(:,1:(n+1))=[];
   ang = unwrap(angle((r-mean(r))...
            +sqrt(-1)*(z-mean(z))));
   ds = z.*pdederive(ang,r,2,2,2,1);
   volloc(k,l) = abs(trapz(ang,2*pi*r.*ds,2));
  end

 end