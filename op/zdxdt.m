% ZDXDT calul la derivee temporelle (difference simple)
%----------------------------------------------------------
% fichier zdxdt.m ->  zdxdt
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la derivee premiere temporelle :
% 
%  dx/dt(k) = (x(k+1) -x(k)) / (t(k+1) -t(k));
%  
% La derivee est prolongee lineairement a la fin de l'intervalle
%
% syntaxe  :
%  
%       xp = zdxdt(x,t)
%       
% entrees :
%
%       t = vecteur temps {size(t) = [N,1]}
%       x = vecteur ou matrice {size(x) = [N,M]}                       
% sorties :
% 
%       xp = derivee temporelle de x par rapport a t {size(xp) = [N,M]}  
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 12/10/1999.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function xp = zdxdt(x,t)

% dimension des entrees
sx=size(x);
nd = length(sx);
t=t(:);

% prolongation lineaire
t(end+1)  =  2.* t(end) -t(end-1);

% selon les dimensions
switch nd
   case 2
   	x(end+1,:) = 2.*x(end,:) - x(end-1,:);
   	ve = ones(1,size(x,2));
   	xp = diff(x,1,1) ./ diff(t*ve,1,1);
   	
	otherwise 
			error('Pas encore implante - a vous de l''ecrire')
			
end	
	
