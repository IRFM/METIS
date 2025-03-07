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
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
t=input('temps ? (s) ');
k = min(find(t<= data.gene.temps));
dk = zget1t(data,k);
xmin = input('xmin ? ');
xmax = input('xmax ? ');
ind = find((param.gene.x>=xmin)&(param.gene.x<=xmax));
x  = param.gene.x(ind);
figure('color',[1 1 1]);
subplot(3,3,1)
plot(x,dk.prof.pe(ind),'r',x,dk.prof.pion(ind),'b');
ylabel('Pe et Pion')

subplot(3,3,2)
plot(x,dk.prof.ped1(ind),'r',x,dk.prof.piond1(ind),'b');
ylabel('d/dx(Pe) et d/dx(Pion)')

subplot(3,3,3)
plot(x,dk.prof.ped2(ind),'r',x,dk.prof.piond2(ind),'b');
ylabel('d^2/dx^2(Pe) et d^2/dx^2(Pion)')

subplot(3,3,4)
plot(x,dk.prof.ne(ind),'r',x,dk.prof.ne(ind).* dk.prof.ae(ind),'b');
ylabel('Ne et Ni')

subplot(3,3,5)
plot(x,dk.prof.ned1(ind)./dk.prof.ne(ind),'r', ...
     x,(dk.prof.ned1(ind).*dk.prof.ae(ind)+dk.prof.ne(ind).*dk.prof.aed1(ind)) ./ ...
     (dk.prof.ne(ind).* dk.prof.ae(ind)),'b');
ylabel('d/dx(Ne)./Ne et d/dx(Ni)./Ni')

subplot(3,3,6)
plot(x,dk.prof.ned2(ind)./ dk.prof.ne(ind),'r', ...
     x,(dk.prof.ned2(ind).*dk.prof.ae(ind)+dk.prof.ne(ind).*dk.prof.aed2(ind) + ...
     2.* dk.prof.ned1(ind) .* dk.prof.aed1(ind)) ./ ...
     (dk.prof.ne(ind).* dk.prof.ae(ind)),'b');
ylabel('d^2/dx^2(Ne) ./Ne et d^2/dx^2(Ni)./Ni')

subplot(3,3,7)
plot(x,dk.prof.ned2(ind)./ dk.prof.ne(ind) .^ 2,'r', ...
     x,(dk.prof.ned2(ind).*dk.prof.ae(ind)+dk.prof.ne(ind).*dk.prof.aed2(ind) + ...
     2.* dk.prof.ned1(ind) .* dk.prof.aed1(ind)) ./ ...
     (dk.prof.ne(ind).* dk.prof.ae(ind)) .^ 2,'b');
ylabel('d^2/dx^2(Ne) ./Ne^2 et d^2/dx^2(Ni)./Ni^2')

subplot(3,3,8);
xie        = (dk.coef.ee+dk.neo.coef.ee)./ dk.prof.ne;
xii        = (dk.coef.ii+dk.neo.coef.ii)./ dk.prof.ne ./ dk.prof.ae;
plot(x,xie(ind),'r',x,xii(ind));
ylabel('Xie et Xii')

subplot(3,3,9);
xied1      = rpdederive(param.gene.x,xie,0,2,2,1);
xiid1      = rpdederive(param.gene.x,xii,1,2,2,2);
plot(x,xied1(ind),'r',x,xiid1(ind));
ylabel('d/dx(Xie) et d/dx(Xii)')


