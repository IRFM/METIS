function plotHarm(parammish,EV,IEV,EW,IEW,ivar);
%
% plotHarm(parammish,EV,IEV,EW,IEW);
%
% ivar = 2 ; vitesse perpendiculaire au surface de flux
% ivar = 1 ; densite
% ivar = 3 ; vitesse poloidal
% ivar = 4 ; vitesse parallele au champ magnetique
% ivar = 5 ; temperature
% ivar = 6 ; potentiel vecteur, composant perpendiculaire
% ivar = 7 ; potentiel vecteur, composant poloidal
% ivar = 8 ; potentiel vecteur, composant toroidal
  
colormap(jet);
cmap = colormap; szmap = size(cmap);

manz = 5;
rfour(1:manz) = parammish(2):(parammish(2)+manz-1);
ng            = parammish(1);
cs            = linspace(0,1,ng);
neq           = 8;
nbg=2*neq*manz;
zma = zeros(2,nbg);
for i=1:ng
  zma = EV(:,i);
  zmai = IEV(:,i);
  ev(1,:,i) = zma;
  ev(2,:,i) = zmai;
end
ev(1,:,1) = 0.;

eigvec = zeros(2*neq,2,manz,ng);
for i=1:neq
   i1 = 2*(i-1)*manz + 1;
   i2 = 2*(i-1)*manz + 2;
   m1 = i1 + 2*manz - 1;
%--------------------------------------------- v1,V2 ------------
  eigvec(2*i-1,:,1:manz,:) = ev(:,i1:2:m1,:);
%--------------------------------------------- d(v1)/ds, v2 midnode ------
  eigvec(2*i,:,1:manz,:) = ev(:,i2:2:m1,:);
end
phase = 0;
hwm = zeros(manz,3);
ns=2;
hold off;
for m=1:manz
  icol = 1;
  if ( manz > 1 )
    icol = fix((m-1)/(manz-1) * (szmap(1)-1)) + 1;
  end
  vec(ns:ng) = (squeeze(eigvec(2*ivar-1,1,m,ns:ng)) * cos(phase) ...
             - squeeze(eigvec(2*ivar-1,2,m,ns:ng)) * sin(phase));
  plot(cs(ns:ng),vec(ns:ng),'-o','Color',cmap(icol,:),'MarkerSize',4);
  mlab = num2str(rfour(m));
  [ym,iym] = max(abs(vec(ns:ng)));
  text(cs(ns+iym-1),vec(ns+iym-1),mlab,'FontSize',12,'Color',cmap(icol,:));
  hold on;
end
xlabel('s');
switch ivar
case 1
  ylabel('N');
case 2
  ylabel('Vper');
case 3
  ylabel('Vpol');
case 4
  ylabel('V//');
case 5
  ylabel('T');
case 6
  ylabel('Aper');
case 7
  ylabel('Apol');
case 8
  ylabel('Ator');
end

zoom on;
str = [' eigenvalue = (',num2str(EW),',',num2str(IEW),')'];
title(str);
hold off;

