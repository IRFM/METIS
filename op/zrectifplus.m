% corrige les profils negatif au bord
function [cor,crx] = zrectifplus(x,old,new)

crx   = 0;	
cor   = new;
if any(new <= 0)
  % au bord seulement
  ind = find(new <= 0);
  if x(min(ind)) < 0.93
      % pas de correction
		crx = 1;
		return
  end
  if any(new <= new(end))
     ind = find(new > new(end));
	  if ~isempty(ind)
	      cor =interp1(cat(2,x(ind),x(end)),cat(2,new(ind),new(end)),x,'linear');
	  else
	     zverbose('Probleme grave -> valeur au bord trop grande ou profil negatif\n');
		  crx = 2;
		  cor = new; 
	  end
   end	  
end
    
mask  = (cor <= 0);
if any(mask)
   cor = cor .* (~mask) + old .* mask; 
end
