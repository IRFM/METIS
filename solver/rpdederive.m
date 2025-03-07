% RPDEDERIVE calcule la derivee d'une grandeur y (reel) par rapport a x (reel)
%----------------------------------------------------------------
% fichier rpdederive.m -> rpdederive
%
% Fonction Matlab 5 : attention compilee avec mcc -> ne marche que sur des reel
%
% Cette fonction calcule la derivee de y (reel) par rapport a x (reel) sur 3 points
% en fixant des conditions de prolongation aux bornes de l'intervalle. 
% La derivee est effectuee selon la dimension d des matrices.
% 
% voir aussi pdederive.m
%
% syntaxe :
% 
%    dydx = rpdederive(x,y,c0,c1,d,o); 
%
% entree :
%        
%   x  = matrice reelle  des X  de memes dimensions que y ou verifiant
%        length(x)=K. (les X, le long de la dimension de derivation
%        doivent etre equidistants pour obtenir une prolongation au 
%        bord de pra continuite precise)   
% 
%   y  = matrice  reelle des Y de dimension dim = size (y), 
%        dim etant un vecteur aynt au moins d elements.
%        et dim(d)=K
%   
%   c0 = matrice donnant les conditions de prolongation k=1,
%        k etant l'indice de x et y sur la dimension d. c0 
%        peut etre un scalaire (condition uniforme pour tous
%        les points) ou une matrice de memes dimension que x, 
%        sauf pour la dimsion d :
%            dim0 =size(c0), dim0(d)=1.
%        Cette matrice code pour des conditions de prolongation
%        variables selon les points.Les valeur de c0 possible sont :
%          * pour la derive premiere :
%                 0 -> derivee nulle en k=1
%                 1 -> derivee est calculee sur deux points  en k=1 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3 en k=1 
%          * pour la derive seconde :
%                 0 -> derivee nulle en k=1
%                 1 -> derivee 1ere nulle en k=K 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3  en k=1      
%        
%   c1 = matrice donnant les conditions de prolongation k=K,
%        k etant l'indice de x et y sur la dimension d. c1
%        peut etre un scalaire (condition uniforme pour tous
%        les points) ou une matrice de memes dimension que x, 
%        sauf pour la dimsion d :
%            dim1 =size(c1), dim1(d)=1.
%        Cette matrice code pour des conditions de prolongation
%        variables selon les points. Les valeur de c1 possible sont :
%          * pour la derive premiere :
%                 0 -> derivee nulle en k=K
%                 1 -> derivee est calculee sur deux points  en k=K 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3 en k=K 
%          * pour la derive seconde :
%                 0 -> derivee 2ieme nulle en k=K
%                 1 -> derivee 1ere nulle en k=K 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3  en k=K      
%
%   d = dimension d'espace sur laquelle porte la derivee
%
%   o = ordre de la derivee (1 ou 2)
%   
% sortie : 
%
%   dydx = derivee de y par rapport a x (memes dimensions que x)
%
% exemple :
% 
%   x=linspace(0,pi/2);
%   y=cos(x);
%   dydx=rpdederive(x,y,0,1,2);
%   plot(x,-sin(x),'o',x,dydx,'+');
%
% ordre de compilation : mcc -V1.2 -r rpdederive
%                        ( chemin vers le compilateur ->addpath /usr/local/matlab5/toolbox/compiler)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.5, du 14/06/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%

function dydx=derive(x,y,c0,c1,d,o)

% pas de verification des dimensions
if nargin == 5
	o=1;
end
% nb dimension des matrices
ts = size(y);
n = length(ts);

% mode 1 (x est un vecteur)
if length(x) == prod(size(x))
	mode = 1;
	x    = x(:);
else
	mode = 0;
end

% permutation si on ne derive pas selon la dimension 1 
if n == 2 & d == 2
	if mode == 0
		x = x.';
	end
	y  = y.';
	c0 = c0.';
	c1 = c1.';

elseif d>1
	if mode == 0
		x = shiftdim(x,d-1);
	end
	y = shiftdim(y,d-1);
	c0 = shiftdim(c0,d-1);
	c1 = shiftdim(c1,d-1);

end

% dimension d'espace
K  = size(x,1);
k0 = 1:(K-1);
k1 = 2:K;

% calcule de dx
% si x est un vecteur
if mode == 1
			% prolongement au bord
			x1 = 2 .* x(K) - x(K-1);
			% prolongement au centre
			x0 = 2 .* x(1) - x(2);
			% concatenation des matrices
			dx= (cat(1,x(k1),x1) - cat(1,x0,x(k0))) ./ 2;

			% complement pour la division
			if n == 1
					comp = 1;
			elseif n == 2
					comp = ones(1,size(y,2));
			elseif n == 3
					comp = ones(1,size(y,2),size(y,3));

			elseif n == 4
					comp = ones(1,size(y,2),size(y,3),size(y,4));

			else
					error('Pas encore implante - a vous de l''ecrire')

			end

else
	% selon le nombre de dimension (1 a 4)
	if n == 1
			% prolongement au bord
			x1 = 2 .* x(K) - x(K-1);
			% prolongement au centre
			x0 = 2 .* x(1) - x(2);
			% concatenation des matrices
			dx= (cat(1,x(k1),x1) - cat(1,x0,x(k0))) ./ 2;

	elseif n == 2
			% prolongement au bord
			x1 = 2 .* x(K,:) - x(K-1,:);
			% prolongement au centre
			x0 = 2 .* x(1,:) - x(2,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:),x1) - cat(1,x0,x(k0,:))) ./ 2;

	elseif n == 3
			% prolongement au bord
			x1 = 2 .* x(K,:,:) - x(K-1,:,:);
			% prolongement au centre
			x0 = 2 .* x(1,:,:) - x(2,:,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:,:),x1) - cat(1,x0,x(k0,:,:))) ./ 2;
			
	elseif n == 4
			% prolongement au bord
			x1 = 2 .* x(K,:,:,:) - x(K-1,:,:,:);
			% prolongement au centre
			x0 = 2 .* x(1,:,:,:) - x(2,:,:,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:,:,:),x1) - cat(1,x0,x(k0,:,:,:))) ./ 2;
			
	else 
			error('Pas encore implante - a vouds de l''ecrire')
			
	end	
	
end

% calcule de dy 
% en fonction de l'ordre de la derivee
if o ==1
	% selon le nombre de dimension (1 a 4)
	if n == 1
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K) - y(K-1)) + ...
	     	     (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1) - y(2)) +  ...
	     	     (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
			% concatenation des matrices
			dy=cat(1,y(k1),y1)-cat(1,y0,y(k0));
			
	elseif n == 2
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:) - y(K-1,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:) - y(2,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:),y1)-cat(1,y0,y(k0,:));
			
	elseif n == 3
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:,:) - y(K-1,:,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:,:) - y(2,:,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:,:),y1)-cat(1,y0,y(k0,:,:));
			
	elseif n == 4
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:,:,:) - y(K-1,:,:,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:,:,:) - y(2,:,:,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:,:,:),y1)-cat(1,y0,y(k0,:,:,:));
			
	else 
			error('Pas encore implante - a vous de l''ecrire')
			
	end	
else
	% selon le nombre de dimension (1 a 4)
	if n == 1
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1) + ...
	     	     (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2) +  ...
	     	     (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
			% concatenation des matrices
			dy=cat(1,y(k1),y1) - 2 .* y + cat(1,y0,y(k0));
			
	elseif n == 2
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:) + ...
	     	     (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:),y1) - 2 .* y + cat(1,y0,y(k0,:));
			
	elseif n == 3
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:,:) + ...
	     	     (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:));
			
	elseif n == 4
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:,:,:) + ...
	     	     (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:,:,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
			% concatenation des matrices
			dy=cat(1,y(k1,:,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:,:));
			
	else 
			error('Pas encore implante - a vous de l''ecrire')
			
	end	
end


% calcul de la derivee
if mode ==1
	dx = dx * comp;
end

if o ==1
	dydx = dy ./ (2 .* dx) ;
else
	dydx = dy ./ (dx .^ 2);
end

% mise a zeros selon les demandes des bornes
% selon le nombre de dimension (1 a 4)
if n == 1
		% 0 au bord
		dydx(K) =(c1~=0) .* dydx(K);
		% 0 au centre
		dydx(1) = (c0~=0) .* dydx(1);
		
elseif n == 2
		% 0 au bord
		dydx(K,:) =(c1~=0) .* dydx(K,:);
		% 0 au centre
		dydx(1,:) = (c0~=0) .* dydx(1,:);
		
elseif n == 3
		% 0 au bord
		dydx(K,:,:) =(c1~=0) .* dydx(K,:,:);
		% 0 au centre
		dydx(1,:,:) = (c0~=0) .* dydx(1,:,:);
		
elseif n == 4
		% 0 au bord
		dydx(K,:,:,:) =(c1~=0) .* dydx(K,:,:,:);
		% 0 au centre
		dydx(1,:,:,:) = (c0~=0) .* dydx(1,:,:,:);
		
else 
		error('Pas encore implante - a vous de l''ecrire')
		
end	


% permutation inverse si on ne derive pas selon la dimension 1 
if n == 2 & d ==2
	dydx  = dydx.';

elseif d>1
	dydx=shiftdim(dydx,n-(d-1));
end
