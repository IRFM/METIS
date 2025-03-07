% PDEDERIVE calcule la derivee d'une grandeur y par rapport a x
%----------------------------------------------------------------
% fichier pdederive.m -> pdederive
%
% Fonction Matlab 5 :
%
% Cette fonction calcule la derivee de y par rapport a x sur 3 points
% en fixant des conditions de prolongation aux bornes de l'intervalle. 
% La derivee est effectuee selon la dimension d des matrices.
%
% syntaxe :
% 
%    dydx = pdederive(x,y,c0,c1,d,o,v0,v1); 
%
% entree :
%        
%   x  = matrice  des X  de memes dimensions que y ou verifiant
%        length(x)=K. (les X, le long de la dimension de derivation
%        doivent etre equidistants pour obtenir une prolongation au 
%        bord de pra continuite precise)   
% 
%   y  = matrice des Y de dimension dim = size (y), 
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
%                 4 -> d?rivee donnee en k=1 
%          * pour la derive seconde :
%                 0 -> derivee nulle en k=1
%                 1 -> derivee 1ere nulle en k=K 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3  en k=1      
%                 4 -> d?rivee donnee en k=1 
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
%                 4 -> d?rivee donnee en k=K 
%          * pour la derive seconde :
%                 0 -> derivee 2ieme nulle en k=K
%                 1 -> derivee 1ere nulle en k=K 
%                 2 -> prolongation par continuite de la derivee 
%                      d'ordre 3  en k=K      
%                 4 -> d?rivee donnee en k=K 
%
%   d = dimension d'espace sur laquelle porte la derivee
%
%   o = ordre de la derivee (1 ou 2)
%
%   v0 = valeur de la d?riv?e en k=1
%
%   v1 = valeur de la d?riv?e en k=K
%   
% sortie : 
%
%   dydx = derivee de y par rapport a x (memes dimensions que x)
%
% exemple :
% 
%   x=linspace(0,pi/2);
%   y=cos(x);
%   dydx=pdederive(x,y,0,1,2);
%   plot(x,-sin(x),'o',x,dydx,'+');
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.4, du 14/05/2001.
% 
% 
% liste des modifications : 
%
% 04/08/99 -> option 1 pour la derivee 2ieme
% 14/05/2001 -> acceleration dans le cas de 2 dimensions
%--------------------------------------------------------------
%

function dydx=derive(x,y,c0,c1,d,o,v0,v1)

% pas de verification des dimensions
if nargin < 6
	o=1;
end
if nargin < 7
        sy = size(y);
        sy(d) = 1;
	v0 = zeros(sy);
end
if nargin < 8
        sy = size(y);
        sy(d) = 1;
	v1 = zeros(sy);
end
% nb dimension des matrices
ts = size(y);
n = length(ts);

% mode 1 (x est un vecteur)
if length(x) == prod(size(x))
	mode = 1;
	x     = x(:);
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
	v0 = v0.';
	v1 = v1.';

elseif d>1
	if mode == 0
		x = shiftdim(x,d-1);
	end
	y = shiftdim(y,d-1);
	c0 = shiftdim(c0,d-1);
	c1 = shiftdim(c1,d-1);
	v0 = shiftdim(v0,d-1);
	v1 = shiftdim(v1,d-1);
	
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
			switch n
				
				case 1
					comp = 1;
				case 2
					comp = ones(1,size(y,2));
				case 3 
					comp = ones(1,size(y,2),size(y,3));
					
				case 4 
					comp = ones(1,size(y,2),size(y,3),size(y,4));
					
				otherwise 
					error('Pas encore implante - a vouds de l''ecrire')
			
			end
				
else
	% selon le nombre de dimension (1 a 4)
	switch n
		
	     case 1
			% prolongement au bord
			x1 = 2 .* x(K) - x(K-1);
			% prolongement au centre
			x0 = 2 .* x(1) - x(2);
			% concatenation des matrices
			dx= (cat(1,x(k1),x1) - cat(1,x0,x(k0))) ./ 2;
			
		case 2
			% prolongement au bord
			x1 = 2 .* x(K,:) - x(K-1,:);
			% prolongement au centre
			x0 = 2 .* x(1,:) - x(2,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:),x1) - cat(1,x0,x(k0,:))) ./ 2;
			
		case 3
			% prolongement au bord
			x1 = 2 .* x(K,:,:) - x(K-1,:,:);
			% prolongement au centre
			x0 = 2 .* x(1,:,:) - x(2,:,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:,:),x1) - cat(1,x0,x(k0,:,:))) ./ 2;
			
		case 4
			% prolongement au bord
			x1 = 2 .* x(K,:,:,:) - x(K-1,:,:,:);
			% prolongement au centre
			x0 = 2 .* x(1,:,:,:) - x(2,:,:,:);
			% concatenation des matrices
			dx= (cat(1,x(k1,:,:,:),x1) - cat(1,x0,x(k0,:,:,:))) ./ 2;
			
		otherwise 
			error('Pas encore implante - a vouds de l''ecrire')
			
	end	
	
end

% calcule de dy 
% en fonction de l'ordre de la derivee
if o == 1
	% selon le nombre de dimension (1 a 4)
	switch n
		
	   case 1
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K) - y(K-1)) + ...
	     	     (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
            	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .*  (2 .* (x(end) - x(end - 1)) .* v1 + y(end - 1));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1) - y(2)) +  ...
	     	     (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
            	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (y(2) + 2 .* (x(1) - x(2)) .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1),y1)-cat(1,y0,y(k0));
			
		case 2
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:) - y(K-1,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
            	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .*  (2 .* (x(end,:) - x(end - 1,:)) .* v1 + y(end - 1,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:) - y(2,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
            	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (y(2,:) + 2 .* (x(1,:) - x(2,:)) .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1,:),y1)-cat(1,y0,y(k0,:));
			
		case 3
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:,:) - y(K-1,:,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
            	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .*  (2 .* (x(end,:) - x(end - 1,:)) .* v1 + y(end - 1,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:,:) - y(2,:,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
            	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (y(2,:,:) + 2 .* (x(1,:,:) - x(2,:,:)) .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1,:,:),y1)-cat(1,y0,y(k0,:,:));
			
		case 4
			% prolongement au bord
	     	y1 = (c1~=2) .* (2 .* y(K,:,:,:) - y(K-1,:,:,:)) + ...
	     	     (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
            	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .*  (2 .* (x(end,:,:,:) - x(end - 1,:,:,:)) .* v1 + y(end - 1,:,:,:));
			% prolongement au centre
	     	y0 = (c0~=2)  .* (2 .* y(1,:,:,:) - y(2,:,:,:)) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
            	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (y(2,:,:,:) + 2 .* (x(1,:,:,:) - x(2,:,:,:)) .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1,:,:,:),y1)-cat(1,y0,y(k0,:,:,:));
			
		otherwise 
			error('Pas encore implante - a vous de l''ecrire')
			
	end	
else
	% selon le nombre de dimension (1 a 4)
	switch n

        case 1
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1) + ...
	     	     (c1==2) .* (4 .* y(K) - 6 .* y(K-1) + 4 .* y(K-2) - y(K-3));
 	  	dxend = (x(end) - x(end -1)) .^ 2;
           	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .* (2 .* y(end) - y(end - 1)  + dxend .* v1);
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2) +  ...
	     	     (c0==2) .* (4 .* y(1) - 6 .* y(2) + 4 .* y(3) - y(4));
	        dx21 = (x(2) - x(1)) .^ 2;
           	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (2 .* y(1) - y(2)  + dx21 .* v0);
		% concatenation des matrices
		dy=cat(1,y(k1),y1) - 2 .* y + cat(1,y0,y(k0));
			
		case 2
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:) + ...
	     	     (c1==2) .* (4 .* y(K,:) - 6 .* y(K-1,:) + 4 .* y(K-2,:) - y(K-3,:));
 	  	dxend = (x(end,:) - x(end -1,:)) .^ 2;
           	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .* (2 .* y(end,:) - y(end - 1,:)  + dxend .* v1);
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:) - 6 .* y(2,:) + 4 .* y(3,:) - y(4,:));
	        dx21 = (x(2,:) - x(1,:)) .^ 2;
           	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (2 .* y(1,:) - y(2,:)  + dx21 .* v0);
		
		% concatenation des matrices
		dy=cat(1,y(k1,:),y1) - 2 .* y + cat(1,y0,y(k0,:));
			
		case 3
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:,:) + ...
	     	     (c1==2) .* (4 .* y(K,:,:) - 6 .* y(K-1,:,:) + 4 .* y(K-2,:,:) - y(K-3,:,:));
 	  	dxend = (x(end,:,:) - x(end -1,:,:)) .^ 2;
           	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .* (2 .* y(end,:,:) - y(end - 1,:,:)  + dxend .* v1);
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:) - 6 .* y(2,:,:) + 4 .* y(3,:,:) - y(4,:,:));
	        dx21 = (x(2,:,:) - x(1,:,:)) .^ 2;
           	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (2 .* y(1,:,:) - y(2,:,:)  + dx21 .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:));
			
		case 4
			% prolongement au bord
	     	y1 = (c1~=2) .* y(K-1,:,:,:) + ...
	     	     (c1==2) .* (4 .* y(K,:,:,:) - 6 .* y(K-1,:,:,:) + 4 .* y(K-2,:,:,:) - y(K-3,:,:,:));
 	  	dxend = (x(end,:,:,:) - x(end -1,:,:,:)) .^ 2;
           	y1 = (c1 ~= 4) .* y1 + (c1 == 4) .* (2 .* y(end,:,:,:) - y(end - 1,:,:,:)  + dxend .* v1);
			% prolongement au centre
	     	y0 = (c0~=2)  .* y(2,:,:,:) +  ...
	     	     (c0==2) .* (4 .* y(1,:,:,:) - 6 .* y(2,:,:,:) + 4 .* y(3,:,:,:) - y(4,:,:,:));
	        dx21 = (x(2,:,:,:) - x(1,:,:,:)) .^ 2;
           	y0 = (c0 ~= 4) .* y0 + (c0 == 4) .* (2 .* y(1,:,:,:) - y(2,:,:,:)  + dx21 .* v0);
			% concatenation des matrices
			dy=cat(1,y(k1,:,:,:),y1) - 2 .* y + cat(1,y0,y(k0,:,:,:));
			
		otherwise 
			error('Pas encore implante - a vous de l''ecrire')
			
	end	
end


% calcul de la derivee
if mode ==1
	try
		dx = dx(:,ones(size(comp)));
	catch
		dx = dx * comp;
        end
end

if o ==1
	dydx = dy ./ (2 .* dx) ;
else
	dydx = dy ./ (dx .^ 2);
end

% mise a zeros selon les demandes des bornes
% selon le nombre de dimension (1 a 4)
switch n
	
     case 1
		% 0 au bord
		dydx(K) =(c1~=0) .* dydx(K);
		% 0 au centre
		dydx(1) = (c0~=0) .* dydx(1);
		
	case 2
		% 0 au bord
		dydx(K,:) =(c1~=0) .* dydx(K,:);
		% 0 au centre
		dydx(1,:) = (c0~=0) .* dydx(1,:);
		
	case 3
		% 0 au bord
		dydx(K,:,:) =(c1~=0) .* dydx(K,:,:);
		% 0 au centre
		dydx(1,:,:) = (c0~=0) .* dydx(1,:,:);
		
	case 4
		% 0 au bord
		dydx(K,:,:,:) =(c1~=0) .* dydx(K,:,:,:);
		% 0 au centre
		dydx(1,:,:,:) = (c0~=0) .* dydx(1,:,:,:);
		
	otherwise 
		error('Pas encore implante - a vous de l''ecrire')
		
end	


% permutation inverse si on ne derive pas selon la dimension 1 
if n == 2 & d ==2
	dydx  = dydx.';

elseif d>1
	dydx=shiftdim(dydx,n-(d-1));
end
