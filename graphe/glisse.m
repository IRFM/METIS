%
%	Lisse.m
%               Yl = Glisse(Yb,Alpha) 
%
%	Lissage glissant d'un vecteur Yb avec coefficient alpha
%	Si matrice : lissage suivant la plus grande dimension
%		

function Yl = Glisse(Yb,Alpha);

Yt = Yb; ss = size(Yb); Beta = (1-Alpha);
if ss(1) < ss(2), Yt = Yt'; end
Yl = Yt*0; Yl(1,:) = Yt(1,:);

for n = 2:length(Yt),
  Yl(n,:) = Yl(n-1,:)*Beta+Yt(n,:)*Alpha;
end

if ss(1) < ss(2), Yl = Yl'; end