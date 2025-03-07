% ZINOUT determine les points interieurs et exterieurs au plasma
%---------------------------------------------------------------
% fichier zinout.m ->  zinout
%
%
% fonction Matlab 5 :
%
% Cette fonction determine les points interieurs et exterieurs au plasma.
% La frontiere appartient au plasma.
%  
% syntaxe  :
%  
%     mask = zinout(R,Z,r,z,stric);
%    
% entree :
%
%     R, Z =  coordonnees de la derniere surface magnetique (vecteur)
%     r, z =  coordonneees des points a tester [M,N]
%     strict = si = 0 la frontiere est contenu dans la surface;
%             si = 1 que les points strictement dans la surface
%             si = -1 que les points sur le contour
%             si = -2 que les points du contour
%
% sorties :
% 
%     mask =  mask  des points interieurs (1 = interieur, 0 = exterieur) [M,N]
%
% fonction ecrite par J-F Artaud , poste 62-15
% version 4.1 , du 06/10/2008.
% 
% 
%--------------------------------------------------------------
%
function mask_out = zinout(R,Z,r,z,strict)

if nargin < 5
	strict = 1;
elseif isempty(strict)
	strict = 1;
end


% mise en ligne
R    = R(:);
Z    = Z(:);

% securite ITM
u = R + sqrt(-1) * Z;
indbad = find(diff(u) == 0);
while(~isempty(indbad))
  u(indbad) = [];
  indbad = find(diff(u) == 0);
end
R = real(u);
Z = imag(u);

if (R(1) ~= R(end)) | (Z(1) ~= Z(end))
	R(end+1) = R(1);
	Z(end+1) = Z(1);
end

% inpolygon is faster
[IN ON] = inpolygon(r, z, R, Z);
switch strict
  case 0
    mask_out = IN | ON;
  case 1
    mask_out = IN;
  case -1
    mask_out = ON;
  case -2
    [void,indice] = intersect(r + sqrt(-1) * z , R + sqrt(-1) * Z);
    mask_out = zeros(size(r));
    mask_out(indice) = 1;
  otherwise
    error('input strict value not supported')
end
