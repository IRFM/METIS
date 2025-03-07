% ZEXTRAITNOEUD extrait les noeuds d'une consigne ou d'un profil
%---------------------------------------------------------------------
% fichier zextraitnoeud.m ->  zextraitnoeud
%
%
% fonction Matlab 5 :
%
% Cette fonction extrait les noeuds d'une consigne ou d'un profil. 
% Elle est utilisee par l'editeur de consignes
%  
% syntaxe  :
%  
%     [xx,yy]= zextraitnoeud(x,y,nb,mode);
%    
% entree :
%
%     x     = vecteur d'entree (temps, rho)
%     y     = donnees (x) [vecteur]
%     nb    = nombre de points en sortie (si mode =1)
%     mode  = 0   -> mode lineaire (consigne)
%             1   -> mode spline   (profil)
%             2   -> comme 0 mais retour a zero forcer en fin d'intervalle
%
% sorties :
% 
%     xx =  abscisses des noeuds [vecteur]
%     yy =  valeur de la consigne ou du profils aux noeuds yy(x) [vecteur]
%     
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.9, du 29/08/2002.
% 
% 
% liste des modifications :
%
% * 21/03/2002 -> test d'un nouvel algorithme 
% * 14/05/2002 -> petit bug detection alignement point final
% * 03/06/2002 -> bug si nb topr grand
% * 29/08/2002 -> ajout de mode = 2 
%
%--------------------------------------------------------------
%
function  [xx,yy]= zextraitnoeud(x,y,nb,mode)

xmemin = x;
ymemin = y;

if mode  > 1
	xx = linspace(min(x),max(x),nb);
	yy = pchip(x,y,xx);
elseif mode == 1
	xx = linspace(min(x),max(x),nb);
	yy = spline(x,y,xx);

elseif mode == 0
	% orientation
	if size(x,1) == 1;
	  or = 0;
	else
	  or = 1;
	  x  = x';
	  y  = y';
	end

	% nombre de points en entree
	nbtot  = length(y);
	
	% prefiltrage
	%y = medfilt1(y,3);
	%y = sgolayfilt(y,1,3);
	
	% echelle
	scale  = max(y) - min(y);
	if (std(y) ./ max(eps,abs(mean(y))))  < (1./nbtot./2)
		xx = x([1,end]);
		yy = y([1,end]);
		if or == 1
	 		 xx  = xx(:);
	  		 yy  = yy(:);
		else
	 		 xx  = xx(:)';
	  		 yy  = yy(:)';
		end
		return
	end

	% detection des points alignes
	indalign = NaN;
	ym = cat(2,Inf,y(1:(end-1)));
	yp = cat(2,y(2:end),Inf);
	dy =  max(abs(y-ym),abs(y-yp));
	indalign = find((dy./scale) < (1./nbtot./2));
	if ~isempty(indalign)
			x(indalign) = [];
			y(indalign) = [];
	end
	nb = min(nb,length(x)-2);

	% boucle de recherche des point pertinents
	xx = x([1,end]);
	yy = y([1,end]);

	fin   = 0;
	nbc   = 3 .* nb;
	while (nbc >0) & (fin ==0)
		yt = pchip(xx,yy,x);
		%figure(61);clf;plot(x,y,xx,yy,'o',x,yt,'+');drawnow
		%d  = abs(medfilt1(y,3) - yt) ./ scale;
		d  = abs(y - yt) ./ scale;
		d([1,end]) = 0;
		e  = trapz(x,abs(yt -y),2) ./ scale ./ (max(x) - min(x));
		if max(e) < (1./nb./10)
			fin = 1;
		end
		indadd = min(find(d == max(d)));
		if any(xx == x(indadd))
			x(indadd) = [];
			y(indadd) = [];
		else
			xx     = cat(2,xx,x(indadd));
			yy     = cat(2,yy,y(indadd));
			[xx,ind] = sort(xx);
			yy       = yy(ind);
		end
		nbc    = nbc - 1;
	end


	if or == 1
		xx  = xx(:);
		yy  = yy(:);
	else
		xx  = xx(:)';
		yy  = yy(:)';
	end

else
	% orientation
	if size(x,1) == 1;
	  or = 0;
	else
	  or = 1;
	  x  = x';
	  y  = y';
	end

	dy    = pdederive(x,y,2,2,2,1);
	yp    = cumtrapz(x,abs(dy),2);
	if mode == 2
		seuil = 2 .* (max(y) - min(0,min(y))) /nb;
	else
		seuil = 2 .* (max(y) - min(y)) /nb;
	end
	if seuil == 0
	   xx(1)  = x(1);
	   yy(1)  = y(1);
		xx(2)  = x(end);
	   yy(2)  = y(end);
	else
	    yp    = yp ./ seuil;
       %nb    = ceil(max(yp))+1;
		 xx     = x;
		 yy     = y;
		 k      = 1;
		 indc   = 0;
		 l      = 1;
	    while ((l-1) <= max(yp))
	   	 if k == 1
	   		 xx(1)  = x(1);
	   		 yy(1)  = y(1);
	   		 zz(1)  = yp(1);
				 k      = k + 1;
	   	 else
	   		 d 	 = abs(yp -(l-1));
	   		 ind   = find(d == min(d));
				 if length(ind) == 1
				    indc  = ind;
				 else
					 indc = min(ind);
				 end
	   	    %
				 % detection d'alignement
				 %
			    if k > 2
				   x1 = xx(k - 2);
				   x2 = xx(k - 1);
				   y1 = yy(k - 2);
				   y2 = yy(k - 1);
					a  = (y2 - y1) ./ (x2 -x1);
					b  = y1 - a .* x1;
					ys = a .* x(indc) + b;
					if ((abs(ys - y(indc))) ./ (max(y) - min(y))) < 1e-2;
	   		      xx(k-1) = x(indc);
	   		      yy(k-1) = y(indc);
	   		      zz(k-1) = yp(indc);
					elseif 	 xx(k - 1) ~= x(indc);
	   		      xx(k) = x(indc);
	   		      yy(k) = y(indc);
	   		      zz(k) = yp(indc);
					   k     = k + 1;
					end
				 elseif 	 xx(k - 1) ~= x(indc);
	   		   xx(k) = x(indc);
	   		   yy(k) = y(indc);
	   		   zz(k) = yp(indc);
					k     = k + 1;
				 end
			 end
			 l  = l + 1;
	    end
	   %
	   % le point final
	   %
	   indc = length(y);
	   %
	   % detection d'alignement
	   %
	   if k > 2
	      x1 = xx(k - 2);
		   x2 = xx(k - 1);
		   y1 = yy(k - 2);
		   y2 = yy(k - 1);
		   a  = (y2 - y1) ./ (x2 -x1);
		   b  = y1 - a .* x1;
		   ys = a .* x(indc) + b;
		   %if (abs(ys - y(indc))) ./ seuil < 0.1
			if ((abs(ys - y(indc))) ./ (max(y) - min(y))) < 1e-2;
	   	   xx(k-1) = x(indc);
	   	   yy(k-1) = y(indc);
	   	   zz(k-1) = yp(indc);
		   elseif 	 xx(k - 1) ~= x(indc);
	   	   xx(k) = x(indc);
	   	   yy(k) = y(indc);
	   	   zz(k) = yp(indc);
			   k     = k + 1;
		   end
	   else
	     xx(k) = x(indc);
	     yy(k) = y(indc);
	     zz(k) = yp(indc);
		  k     = k + 1;
	   end
      %
	   % on garde que les points utils
	   %
	   xx     = xx(1:(k-1));
	   yy     = yy(1:(k-1));
	   zz     = zz(1:(k-1));

	end


	indxx  = find(diff(xx) <=0);
	if ~isempty(indxx)
	   xx(indxx) =[];
	   yy(indxx) =[];
	   zz(indxx) =[];
	end
	if or == 1
	  xx  = xx';
	  yy  = yy';

	end
end

