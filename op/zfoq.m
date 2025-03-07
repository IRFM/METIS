% calcul de la reduction du transport sur les surface rationnelle a partir du flux helicoidal
function frs = zfoq(x,q,psi,mmax)

% si mmax n'est pas specifie
if nargin <4
   mmax = 30;
elseif isempty(mmax)
   mmax = 30;
end

% derivee de psi
psid1  = pdederive(x,psi,0,1,2,1);

% calcul de l'integrant
% boucle sur m et n
count = 0;
inte  = ones(size(x));
for m = 1:mmax
   for n = 1:mmax
      if (m./n < max(q)) &((m./n) >=1)
         poid  =  m .* n;
         inte  = inte .* (1 + tanh(abs(q - m./n) .*poid)) ./ 2;
         count = count +1;
      end
   end
end

frs   = max(0.1,min(inte,1));
