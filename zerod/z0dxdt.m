% Z0DXDT calul la derivee temporelle (difference simple)
%----------------------------------------------------------
% fichier z0dxdt.m ->  z0dxdt
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
%       xp = z0dxdt(x,t)
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
function xp = z0dxdt(x,t)

% dimension des entrees
sx  = size(x);
nd  = length(sx);
t   = t(:);
ve = ones(1,size(x,2));
xp = diff(x,1,1) ./ diff(t*ve,1,1);
if length(t) < 3
    xp  = cat(1,xp(1),xp(1));
elseif  length(t) < 5   
    xp = cat(1,pchip(t(2:end)',xp',t(1))',xp);  
else
    %xp_ = xp;
    if ndims(xp) < 3
        xp = cat(1,pchip(t(2:5)',xp(1:4,:)',t(1))',xp); 
    elseif ndims(xp) == 3
        xp = cat(1,pchip(t(2:5)',xp(1:4,:,:)',t(1))',xp); 
    elseif ndims(xp) == 4
        xp = cat(1,pchip(t(2:5)',xp(1:4,:,:,:)',t(1))',xp); 
    else
        error('not yet implemanted');
    end     
    %xp_ = cat(1,pchip(t(2:end)',xp_',t(1))',xp_);  
    %std(xp - xp_)    
end
