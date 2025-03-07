function [rayon,surface,L] = lcordeyp(Rp,ap,epsip,hp,Rm,Ri,ai,fi,tv,kv,mv,nv0,nc);
%
%	Calculation of the transfer matrix for a given plasma configuration
%	between line integrated emission and local emission. Each matrix coefficient T(i,j)
%	is corresponding to the length of the chord i in the plasma layer j.
%
%	Input:
%
%		- Rp: plasma major radius (mm)  [1,1]
%		- ap: plasma minor radius (mm) [1,1]
%		- epsip: plasma ellipticity (at three positions: 0, 1/2 and 1) [3,1]
%		- hp: plasma vertical shift (mm) [1,1]
%		- Rm: magnetic axis (mm) or Shafranov shift profile [1,1] or [n,2]
%		- Ri: chamber major radius (mm)  [1,1]
%		- ai: chamber minor radius (mm) [1,1]
%		- fi: chamber ellipticity, elongation and triangularity  [1,3]
%		- tv: angular position of the intersection between detector axis direction 
%			  and the chamber (degrees) [1,m] 
% 		- kv: horizontal rotation (degrees) [1,m] 
% 		- mv: vertical rotation (degrees) [1,m] 
%		- nv0: number of chords viewing the plasma in a poloidal plane [1,m]
%		- nc: number of flux surfaces [1,1]
%		  (default = 100)
%
%	Output:
%
%		- rayon: radius of each layer (mm) [1,nc]
%		- surface: surface of each layer (mm^2) [1,nc]
%		- L: sparse transfer matrix (mm) [nv0,nc]
%
%	Warning: If length(Rm) > 1, then Rm = [rshfit(mm),shift(mm)]. Its a way to use 
%			 the experimental Shafranov shift profile in the calculations.
%
%by Y.PEYSSON CEA-DRFC 19/09/1991 <peysson@fedv09.cad.cea.fr>
%revised for MatLab 4.1 (25/08/1994) 
%
if nargin < 12,
	infoyp(2,'Wrong number of input arguments for lcordeyp');
	return;
end
if nargin == 12,
	nc = 100;
end
%
conv = pi/180.0;
rayon = [ap/nc:ap/nc:ap];
xepsip = fepsipyp(epsip,length(rayon));
surfaceT = xepsip.*pi.*rayon.^2;
if length(Rm) > 1,
	delta = spline(Rm(:,1),Rm(:,2),rayon);
else
	delta0 = Rm-Rp
	delta = delta0*(1-rayon.^2/ap^2);
end

xi1 = Ri+ai*(cos(conv*tv)+fi(2)*cos(2*conv*tv));
zi1 = ai*(fi(1)*sin(conv*tv)+fi(3)*sin(2*conv*tv));
%
dx = Rp+delta;dz = hp;
H = 0.0*ones(nv0,nc);
%H = sparse(nv0,nc);
L = H;
sL = size(L);
for i=1:nv0,
 	for j=1:nc,
		mhu = sin(conv*(90+mv(i)))^2+cos(conv*(90+mv(i)))^2/xepsip(j)^2;
  		if mv(i)==0,
   			H(i,j) = 2*max(rayon(j)^2-(zi1(i)-dz)^2/xepsip(j)^2,0)^0.5;
  		else
   			bb = -(xi1(i)-dx(j))*tan(conv*(90+mv(i)))-zi1(i)*tan(conv*(90+mv(i)))^2-dz^2/xepsip(j)^2;
   			aa = tan(conv*(90+mv(i)))^2+1.0/xepsip(j)^2;
   			cc = xi1(i)^2+dx(j)^2+zi1(i)^2*tan(conv*(90+mv(i)))^2+2*xi1(i)*zi1(i)*tan(conv*(90+mv(i)));
   			cc = cc-2*dx(j)*xi1(i)-2*dx(j)*zi1(i)*tan(conv*(90+mv(i)))+dz^2/xepsip(j)^2-rayon(j)^2;
			H(i,j) = 2*abs(cos(conv*(90+mv(i))))*max(bb^2-aa*cc,0)^0.5/mhu;
 		end
 	end
end
%
for j=nc:-1:2,
	L(:,j) = H(:,j)-H(:,j-1);
 	surface(j) = surfaceT(j)-surfaceT(j-1);
end
%
L(:,1) = H(:,1);
surface(1) = surfaceT(1);





