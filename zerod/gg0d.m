% coefficient de geometrie
function [phi,dphidx,vpr,grho2r2,r2i,ri,C2,C3,grho,grho2,phiplasma,badjac] = ...
           gg0d(x,Raxe,a,kx,d,Vp,Sp,Fin,Rsepa,Zsepa,rho_in,expo)

% securite sur les valeurs d'entree
Raxe  = max(0.01,min(1e3,Raxe));
a     = max(0.01,min(0.9 .* Raxe(:,end-1),a));
kx    = max(0.5,min(10,kx));
d     = max(0,min(0.999,d));
% delta = Raxe - Raxe(:,end) * ve;

% rb0
rb0 =  Fin(:,end);

%  % derivee
%  dkx    = pdederive(x(1,:),kx,0,2,2,1);
%  ddx    = pdederive(x(1,:),t,0,2,2,1);
%  dsh    = pdederive(x(1,:),delta,0,2,2,1);
%  drhodx = pdederive(x(1,:),rho_in,2,2,2,1);

% x est la coordonnee de Lao
% extention vers l'exterieur
xin    = x;
x      = cat(2,x(1:end - 1),x(end) .* 0.99,x(end),x(end) .* 1.01);
ve     = ones(size(x));
vt     = ones(size(Raxe,1),1);
rho_in = pchip(xin,rho_in,x); 
kx     = pchip(xin,kx,x); 
Raxe   = pchip(xin,Raxe,x); 
t      = asin(pchip(xin,d,x)); 
Fin    = cat(2,Fin(:,1:end - 1),Fin(:,end),Fin(:,end),Fin(:,end)); 



% vecteur utils
amin  = a  * ve;
x      = vt * x;
Zc    = 0 .* vt;
R0    = Raxe(:,end - 1) * ve;

% initialidation :
vpr      = zeros(size(Raxe));
grho2r2  = zeros(size(Raxe));
r2i      = zeros(size(Raxe));
ri       = zeros(size(Raxe));
C2       = zeros(size(Raxe));
C3       = zeros(size(Raxe));
C4       = zeros(size(Raxe));
Cgrho    = zeros(size(Raxe));
Cgrho2   = zeros(size(Raxe));
dphidx   = zeros(size(Raxe));
dphidxw  = zeros(size(Raxe));
badjac   = zeros(size(Raxe));
% dphidx_alt   = zeros(size(Raxe));
% dphidxw_alt  = zeros(size(Raxe));

% angle theta
th   = linspace(0,2.*pi,35);
vh   = ones(size(th));
sh   = sin(th);
ch   = cos(th);
%sina = vt * sh;
%cosa = vt * ch;

if ~isempty(Rsepa) & ~isempty(Zsepa)
	% separatrice donnees par les moments
	Rmom   = R0(:,end-1) * vh + (a * vh) .* cos(vt * th + t(:,end-1) * sin(th));
	Zmom   = (a .* kx(:,end-1)) * sin(th);	
end


% flag
lmem = [];
rhogp = [];
warning_0 = false;
% boucle sur les surface
for l = 1:size(x,2)

	if l > 1
		Rm = R;
		Zm = Z;
		R  = Rp;
		Z  = Zp;
		dSphim = dSphi;
		dSphi  = dSphip;
		dSphimw = dSphiw;
		dSphiw  = dSphipw;
		rhogm = rhog;
		rhog  = rhogp;
	else
		R = Raxe(:,1) * vh;
		Z = 0 .* (vt *vh);
		dSphi   = Z;	
		dSphiw   = Z;	
		rhog    = Z;		
	end
	if l < size(x,2)
        if size(Raxe,1) == 1
            Rp = Raxe(:,l+1) * vh + (amin(:,l+1) * vh) .* (x(:,l+1) * vh) .*  ...
	   	   	     cos( vt * th + (t(:,l+1) * vh) .* sin(vt * th));
		    Zp = (kx(:,l+1) * vh) .*(amin(:,l+1) * vh) .* (x(:,l+1) * vh) .* sin(vt * th);
        else
            Rp = squeeze(Raxe(:,l+1,ones(size(vh)))) + squeeze(amin(:,l+1,ones(size(vh)))) .* squeeze(x(:,l+1,ones(size(vh)))) .*  ...
                    cos(th(ones(size(vt)),:) + squeeze(t(:,l+1,ones(size(vh)))) .* sin(th(ones(size(vt)),:)));
            Zp = squeeze(kx(:,l+1,ones(size(vh)))) .* squeeze(amin(:,l+1,ones(size(vh)))) .* squeeze(x(:,l+1,ones(size(vh)))) .* sin(th(ones(size(vt)),:));
        end
		% calcul de la correction de forme
        Rp_no = Rp;
        Zp_no = Zp;        
		if ~isempty(Rsepa) & ~isempty(Zsepa) & (expo >0) 		
			[Rp,Zp]   = z0morph(Rp,Zp,Rsepa,Zsepa,Rmom,Zmom,x(1,l+1),expo,l);
		end
		[Rp,Zp,rhogp,warning_0] = chang(Raxe(:,1),Zc,Rp,Zp,th,rhogp);
		
        if any(~isfinite(Rp(:))) || any(~isfinite(Zp(:))) ||any(imag(Rp(:))) || any(imag(Zp(:)))
                warning('gg0d: miss defined LCFS');
                [Rp_no,Zp_no,rhogp_no] = chang(Raxe(:,1),Zc,Rp_no,Zp_no,th);  
                indbad_morph = find(any(~isfinite(Rp),2) | any(~isfinite(Zp),2) | any(imag(Rp),2) | any(imag(Zp),2));
                Rp(indbad_morph,:) = Rp_no(indbad_morph,:);
                Zp(indbad_morph,:) = Zp_no(indbad_morph,:);
                rhogp(indbad_morph,:) = rhogp_no(indbad_morph,:);
        end
        if size(Raxe,1) == 1
                dSphip   = rhogp .* squeeze(Fin(:,l+1,ones(size(vh))))' ./ Rp;
        else
                dSphip   = rhogp .* squeeze(Fin(:,l+1,ones(size(vh)))) ./ Rp;
        end
        inter    = Fin(:,l+1) - rb0;			
		dSphipw  = rhogp .* inter(:,ones(size(vh))) ./ Rp;			
 	end
	if l > 1 & l < size(x,2)
                drho_in = rho_in(:,l+1) - rho_in(:,l-1);
		drho_in = drho_in(:,ones(size(vh)));
		dRdx = (Rp - Rm) ./ drho_in;
		dZdx = (Zp - Zm) ./ drho_in;
		dRdth = pdederive(th,R,2,2,2,1);
		dZdth = pdederive(th,Z,2,2,2,1);
		% periodisation
		dRdth(:,1)  = (dRdth(:,1) + dRdth(:,end)) ./ 2;
		dRdth(:,end)  =dRdth(:,1);
		dZdth(:,1)  = (dZdth(:,1) + dZdth(:,end)) ./ 2;
		dZdth(:,end)  =dZdth(:,1);
		
		
		% jacobien
		RZjac      = (dRdx .* dZdth - dZdx .* dRdth);
		badjac(:,l) = double(any(RZjac < eps,2));
                sign_RZjac = sign(RZjac);
                RZjac      = abs(RZjac);
		% grad rho 
     		dxdR   =  dZdth  ./ max(eps,RZjac) .* sign_RZjac;
		dxdZ   = -dRdth  ./ max(eps,RZjac) .* sign_RZjac;
 		
                gradrho = sqrt(dxdR .^ 2 + dxdZ .^ 2);
		
		% vpr 
		factinte      = 2 .* pi .* R .* RZjac;

		vpr(:,l)      = trapz(th,factinte,2);
		C2(:,l)       = trapz(th,gradrho .^ 2 ./ R .^ 2 .* factinte,2);
		% securite pour haut beta
		if ~isempty(lmem)
			if any(C2(:,l) < C2(:,lmem))
				indtm = find(C2(:,l) < C2(:,lmem));
				C2(indtm,l) = 1.1 .* C2(indtm,lmem);
			end
		end
		lmem     = l;
		C3(:,l)       = trapz(th,factinte ./ R .^ 2,2);
		C4(:,l)       = trapz(th,factinte ./ R,2);		
		Cgrho(:,l)    = trapz(th,gradrho .* factinte,2); 
		Cgrho2 (:,l)  = trapz(th,gradrho .^ 2 .* factinte,2); 
						
%  		dphidx(:,l+1) = trapz(th,abs(dSphip + dSphi) ./ 2  .* (rhogp - rhog) ./ ((x(:,l+1) - x(:,l)) * vh) ,2); 
%  		if l == 2
%  			dphidx(:,l) = trapz(th,abs(dSphi  + dSphim)./ 2  .* (rhog - rhogm) ./ ((x(:,l) - x(:,l-1)) * vh) ,2); 
%  		end
	end
	
	if l == 1
		dphidx(:,l) = 0;
		dphidxw(:,l) = 0;
	elseif l < size(x,2)
  		dphidx(:,l) = trapz(th,(abs(dSphi  + dSphim)./ 2  .* (rhog - rhogm) +  ...
		                        abs(dSphip + dSphi) ./ 2  .* (rhogp - rhog)) ./  ...
					((x(:,l+1) - x(:,l-1)) * vh) ,2); 
 		dphidxw(:,l) = trapz(th,(abs(dSphiw  + dSphimw)./ 2  .* (rhog - rhogm) +  ...
		                        abs(dSphipw + dSphiw) ./ 2  .* (rhogp - rhog)) ./  ...
					((x(:,l+1) - x(:,l-1)) * vh) ,2); 

		%dvdx_alt(:,l)     = trapz(th, pi .* (R + Rm) .* (rhog .^ 2 - rhogm .^ 2),2) ./ (x(:,l) - x(:,l-1)) ./ 2;
		%dsdx_alt(:,l)     = trapz(th, rhog .^ 2 - rhogm .^ 2,2) ./ (x(:,l) - x(:,l-1)) ./ 2;
%                 if all((x(:,l) - x(:,l-1)) > 0)
% 		    dphidx_alt(:,l)   = trapz(th,((Fin(:,l) * vh) ./ R + (Fin(:,l-1) * vh) ./ Rm ) .* (rhog .^ 2 - rhogm .^ 2),2) ./ (x(:,l) - x(:,l-1)) ./ 4;
% 		    dphidxw_alt(:,l)  = trapz(th,((Fin(:,l) * vh  - rb0 * vh) ./ R + (Fin(:,l-1) * vh - rb0 * vh) ./ Rm) .* (rhog .^ 2 - rhogm .^ 2),2) ./ (x(:,l) - x(:,l-1)) ./ 4;
%                 end
		%volume(:,l)   = trapz(th,pi .* deuxd.R(k,:) .* deuxd.rhog(k,:) .^ 2,2);
		%surface(:,l)  = trapz(th,deuxd.rhog(k,:) .^ 2 / 2,2);
               
	else
		dphidx(:,l) = NaN;
		dphidxw(:,l) = NaN;
	end	
 	%save(sprintf('contexte_gg0d_sat_%d',l));

%  	keyboard
end
if warning_0
    fprintf('0');
end

% valeur centrale
vpr(:,1) = 0;
C2(:,1) = 0;
C3(:,1) = 0;
C4(:,1) = 0;
Cgrho(:,1) = 0;
Cgrho2(:,1) = 0;
dphidx(:,1) = 0;
dphidxw(:,1) = 0;

% donnees finale
ll       = cat(2,1:(size(x,2) - 4),(size(x,2) - 1));
vpr      = pchip(x(1,ll),vpr(:,ll),xin);
C2       = pchip(x(1,ll),C2(:,ll),xin);
C3       = pchip(x(1,ll),C3(:,ll),xin);
C4       = pchip(x(1,ll),C4(:,ll),xin);
Cgrho    = pchip(x(1,ll),Cgrho(:,ll),xin);
Cgrho2   = pchip(x(1,ll),Cgrho2(:,ll),xin);
dphidx   = pchip(x(1,ll),dphidx(:,ll),xin);
dphidxw  = pchip(x(1,ll),dphidxw(:,ll),xin);

% rmx nouvelle definition
phi        = cumtrapz(xin,dphidx,2);
phiplasma  = trapz(xin,dphidxw,2);

% calcul des donnees derivees
grho2r2      = C2 ./ max(eps,vpr);
r2i          = C3 ./ max(eps,vpr);
ri           = C4 ./ max(eps,vpr);
grho         = Cgrho ./ max(eps,vpr);
grho2        = Cgrho2 ./ max(eps,vpr);

% continuite au centre
r2i(:,1)     = r2i(:,2);
ri(:,1)      = ri(:,2);
% valeur limite de grho2
ll           = 5:size(xin,2);
grho2r2      = pchip(cat(2,-1,xin(ll)),cat(2,r2i(:,2),grho2r2(:,ll)),xin);
grho2        = pchip(cat(2,-1,xin(ll)),cat(2,vt,grho2(:,ll)),xin);
grho         = pchip(cat(2,-1,xin(ll)),cat(2,vt,grho(:,ll)),xin);

% securite ultime
ll        = 1:size(xin,2);
r2i       = max(0.25 ./ R0(:,ll) .^ 2 ,min(r2i,4 ./ R0(:,ll) .^ 2));
grho2r2   = max(0.25 ./ R0(:,ll) .^ 2 ,min(grho2r2,25 ./ R0(:,ll) .^ 2));
ri        = max(0.25 ./ R0(:,ll) ,min(ri,4 ./ R0(:,ll)));
vpr       = max(0,min(vpr,8 .* pi .^ 2 .* R0(:,ll) .^ 2));
C2        = max(0,C2);
C3        = max(0,C3);
grho      = max(0.5,min(5,grho));
grho2     = max(0.25,min(25,grho2));

%save contexte_gg0d_loc

function [R,Z,rhoo,warning_0] = chang(Rc,Zc,R,Z,th,rhoo_mem)

%output
warning_0 = false;

%R_in = R;
%Z_in = Z;
if nargin  == 5
  rhoo_mem = [];
end

% donnees 2D (grille droite)
indbad = find((Rc <= min(R,[],2)) |(Rc >= max(R,[],2)));
if ~isempty(indbad)
	Rc(indbad) = (min(R(indbad,:),[],2) + max(R(indbad,:),[],2)) ./ 2;
end 

cth  = cos(th);
sth  = sin(th);
vh   = ones(size(th));
vt   = ones(size(R,1),1);
cx   = (R - Rc(:,ones(size(vh)))) + sqrt(-1) .* (Z-Zc(:,ones(size(vh))));
%thx  = unwrap(angle(cx),[],2);
thx = angle(cx);
thx = thx .* (thx>=0) + (2 .* pi + thx) .* (thx <0);
rhox = abs(cx);
thx(thx<0) = thx(thx<0) + 2 .* pi;
[thx,indx] = sort(thx,2);
for k = 1:size(indx,1)
	rhox(k,:)       = rhox(k,indx(k,:));
end
rhox = cat(2,rhox,rhox(:,2:end-1),rhox);
thx = cat(2,thx -2.*pi,thx(:,2:end-1),thx+2.*pi);
% correction of non strictely monotonic coordinate
mask0 = logical(diff(thx,1,2) <= 0);
mask0 = cat(2,zeros(size(thx,1),1),mask0);
if any(mask0(:))
    warning_0 = true;
    mask0(:,end-1) = mask0(:,end-1) | mask0(:,end);
    mask0(:,end)   = ~mask0(:,end-1) & mask0(:,end);
    ind0 = find(any(mask0,2));
    indx = 1:size(thx,2);
    for k=1:length(ind0)
        indc = ind0(k);
        if sum(~mask0(indc,:)) >= 2
            %disp(mask0(indc,:))
            %disp(thx(indc,:))
            %disp(interp1(indx(~mask0(indc,:)),thx(indc,~mask0(indc,:)),indx,'linear','extrap'));
            thx(indc,:) = interp1(indx(~mask0(indc,:)),thx(indc,~mask0(indc,:)),indx,'linear','extrap');          
        else
            thx(indc,:) = linspace(-2.*pi,4.*pi,size(thx,2));
        end
    end
end

%rhoo   = tsplinet(thx,rhox,vt*th);
rhoo   = tsplinet(thx,rhox,th(ones(size(vt)),:));
if ~isempty(rhoo_mem)
  %figure(21);hist(abs(max(rhoo(:),rhoo_mem(:) .* 1.001) - rhoo(:)));drawnow
  rhoo   = max(rhoo,rhoo_mem .* 1.001);
end
R      = Rc(:,ones(size(vh))) + rhoo .* cth(ones(size(vt)),:);
Z      = Zc(:,ones(size(vh))) + rhoo .* sth(ones(size(vt)),:);


