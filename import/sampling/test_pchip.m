% script de test de pchip
figure(21);
for k = 1:1e6
  nin  = max(2,fix(11*rand(1)));
  nout = max(1,fix(11*rand(1)));
  x    = sort(randn(1,nin));
  y    = randn(1,nin);
  xi   = sort(randn(1,nout));
  yi = pchip(x,y,xi);
  %plot(x,y,'or',xi,yi,'b');
  %drawnow
end
