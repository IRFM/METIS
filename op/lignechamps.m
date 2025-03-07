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


angc = linspace(0,20,41)*pi/180;
RRana = RRref(:,1);
clear bt2 br2 bz2
for k1=1:length(RRana)
  for l=1:length(angc)
    [bt2(k1,l),br2(k1,l),bz2(k1,l)]= bcronoscall1p(RRana(k1),0,angc(l),itor);
  end
end
[Rt1,Pt1] = meshgrid(RRana,angc);
Rt1       = Rt1';
val       = cumsum(br2./bt2.*Rt1.*mean(diff(angc))) + Rt1;


drip      = val(:,10)-val(:,3);

angc = linspace(0,20,41)*pi/180;
clear bt3 br3 bz3
RRana = RRref(1:5:size(RRref,1),:);
ZZana = ZZref(1:5:size(RRref,1),:);
[n,m] = size(RRana);
for i=1:n
  for j=1:m
    for k1=1:length(angc)
      [bt3(i,j,k1),br3(i,j,k1),bz3(i,j,k1)]= bcronoscall1p(RRana(i,j),ZZana(i,j),angc(k1),itor);
      Rt1(i,j,k1) = RRana(i,j);
    end
  end
end
val       = cumsum(br3./bt3.*Rt1.*mean(diff(angc))) + Rt1;
drip      = squeeze(val(:,:,10))-squeeze(val(:,:,3));
