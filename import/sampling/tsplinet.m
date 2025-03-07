% TSPLINET		interpolation par spline cubic sur l'espace pour une matrice (temps,epasce)
%-----------------------------------------------------------------------------------------------
%
% fonction tsplinet.m -> utilise un MexFile (source tsplinet.c)
%
%	cette fonction interpole, sur l'espace, par spline cubic une matrice du 
%	type matrice(temps,espace). cette fonction permet connaissant x et y , 
%	tel que y(t)=F(t,x(t)) de deduire yi en donnant xi, tel que 
%	yi(t)=F(t,xi(t)), a partir d'une interpolation par spline cubic.
%
% syntaxe :
%
%		yi=tsplinet(x,y,xi)
%
% entree :
%
%		x = matrice des coordonnes d'espace 
%			x(indice des temps, indice d'espace)
%
%		y = valeurs de la fonction F pour les memes temps
%			que x, evalues en x (matrice de memes dimensions
%			que x, y(indice de temps,indice d'espace))
%
%		xi = matrice des nouvelles coordonees d'espace
%			 xi(indice de temps, indice d'espace). 
%			 xi doit avoir la meme premiere dimension que x 
%
% sortie :
%
%		yi = valeurs de la fonction F pour les position xi  
%
% remarque :
%
% Le contenu de x doit etre range par ordre croissant !
%
% fonction ecrit par J-F Artaud , poste 46-78
% version 1, derniere mise a jour le 10/01/95
%-----------------------------------------------------------------------------------------------
%
%
function yi=tsplinet(x,y,xi)

nbt = size(y,1);
nbx = size(x,2);
nbi = size(xi,2);

yi = NaN * ones(nbt,nbi);

for k = 1:nbt
	if size(x,1) == nbt
		kx =k;
	else
		kx = 1
	end
	if size(xi,1) == nbt
		ki =k;
	else
		ki = 1
    end
    if any(~isfinite(y(k,:)))
        var = y(k,:);
        var(~isfinite(var)) = 0;
        yi(k,:) = spline(x(kx,:),var,xi(ki,:));
        var     = interp1(x(kx,:),y(k,:),xi(ki,:),'nearest','extrap');
        yi(k,~isfinite(var)) = var(~isfinite(var));
    else
        yi(k,:) = spline(x(kx,:),y(k,:),xi(ki,:));
    end

end







