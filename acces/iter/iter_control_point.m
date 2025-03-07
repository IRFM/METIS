% point de control pour la separatrice de ITER
function [e,rs,zs] = iter_control_point(e) 

if nargin == 0
	e.Rx      =  5.0811;
	e.Zx      = -3.3421;
	e.R_leg_r =  5.562;
	e.Z_leg_r = -4.436;
	e.R_leg_l =  4.0541;
	e.Z_leg_l = -3.7763;
	e.R_rmax  =  8.1995;
	e.Z_rmax  =  0.4211;
	e.R_rmin  =  4.2005;
	e.Z_rmin  =  1.0033;
	e.R_zmax  =  5.4213;
	e.Z_zmax  =  4.0476; 
	e.xp_onoff = 1;
end

%
rc = (e.R_rmax + e.R_rmin) ./ 2;
ac = (e.R_rmax - e.R_rmin) ./ 2;
zc = e.Z_rmax;
ku = (e.Z_zmax - zc) ./ ac;
du = (rc - e.R_zmax) ./ ac;
kd = (e.Zx - zc) ./ ac;
dd = (rc - e.Rx) ./ ac;
% partie haute de la separatrice
u    = linspace(pi-pi/12,pi/2 + pi/16,11);
t    = asin(du);
Rup  = rc + ac .* cos(u + t * sin(u));
Zup  = zc + ac .* ku * sin(u);
u    = linspace(pi/2 - pi/16,pi/32,15);
t    = asin(du);
Ruq  = rc + ac .* cos(u + t * sin(u));
Zuq  = zc + ac .* ku * sin(u);
u    = linspace(pi-pi/16,pi+pi/32,5);
t    = asin(dd);
Rl  = rc + ac .* cos(u + t * sin(u));
Zl  = zc + ac .* kd * sin(u);
u    = linspace(pi/32,0.25 .* pi,9);
t    = asin(dd);
Rr  = rc + ac .* cos(u + t * sin(u));
Zr  = zc + ac .* kd * sin(u);

% liste pour rechantionnage
Rct  = [e.R_leg_r,e.Rx,e.R_rmin,e.R_zmax,e.R_rmax,e.Rx,e.R_leg_l];
Zct  = [e.Z_leg_r,e.Zx,e.Z_rmin,e.Z_zmax,e.Z_rmax,e.Zx,e.Z_leg_l];
R    = [e.R_leg_r,e.Rx,Rl,e.R_rmin,Rup,e.R_zmax,Ruq,e.R_rmax,Rr,e.Rx,e.R_leg_l];
Z    = [e.Z_leg_r,e.Zx,Zl,e.Z_rmin,Zup,e.Z_zmax,Zuq,e.Z_rmax,Zr,e.Zx,e.Z_leg_l];

C  = (R -rc) + sqrt(-1) .* (Z -zc);
S  = cumsum(abs(C));
s  = linspace(S(1),S(end),201);
r = pchip(S,R,s);
z = pchip(S,Z,s);
%
if e.xp_onoff ~= 0
	ss  = linspace(S(2),S(end - 1),201);
	rs = pchip(S,R,s);
	zs = pchip(S,Z,s);
else
	u    = linspace(0,2 .* pi,201);
	t    = asin(abs(du));
	rs  = rc + ac .* cos(u + t * sin(u));
	zs  = zc + ac .* abs(ku) * sin(u);
	t    = asin(abs(dd));
	R2  = rc + ac .* cos(u + t * sin(u));
	Z2  = zc + ac .* abs(kd) * sin(u);
        indl = find(zs < zc);
	rs(indl) = R2(indl);
	zs(indl) = Z2(indl);	
end
%

if nargout ~= 0
	return
end

[rl,zl] = iterlimiteur;
[rref,zref] = ref_rz;
% 
figure
plot(Rct,Zct,'or',r,z,'c',rs,zs,'b',rl,zl,'m',Rup,Zup,'og',Ruq,Zuq,'og',Rl,Zl,'og',Rr,Zr,'og',rc,zc,'+r', ...
     rref,zref,'.k');
axis('equal')


function [rref,zref] = ref_rz


ref = [
6.0094 -5.4566
5.9133 -5.2071
5.8119 -4.9570
5.7104 -4.7219
5.5798 -4.4308
5.4460 -4.1396
5.3113 -3.8486
5.0811 -3.3421
4.9261 -2.9753
4.8127 -2.6841
4.7097 -2.3930
4.6179 -2.1018
4.5424 -1.8314
4.4669 -1.5196
4.3879 -1.1315
4.3245 -0.7434
4.2750 -0.3552
4.2386 0.0329
4.2191 0.3241
4.2066 0.6152
4.2005 1.0033
4.2044 1.2944
4.2160 1.5855
4.2365 1.8767
4.2664 2.1677
4.3080 2.4589
4.3639 2.7500
4.4448 3.0646
4.5380 3.3322
4.6829 3.6233
4.8354 3.8218
4.9330 3.9069
5.0306 3.9681
5.1283 4.0078
5.2259 4.0335
5.3236 4.0460
5.4213 4.0476
5.5188 4.0403
5.6165 4.0252
5.7142 4.0035
5.9095 3.9441
6.1976 3.8174
6.3726 3.7204
6.5930 3.5782
6.7864 3.4293
6.9836 3.2548
7.1789 3.0534
7.3511 2.8470
7.4718 2.6844
7.6181 2.4589
7.7783 2.1677
7.9098 1.8767
8.0155 1.5855
8.0969 1.2944
8.1554 0.9973
8.1801 0.8092
8.1952 0.6152
8.1995 0.4211
8.1930 0.2270
8.1753 0.0329
8.1464 -0.1611
8.1059 -0.3552
8.0222 -0.6463
7.9094 -0.9374
7.7648 -1.2294
7.5855 -1.5196
7.3647 -1.8108
7.1789 -2.0175
6.9836 -2.2066
6.7629 -2.3930
6.5015 -2.5870
6.3001 -2.7209
6.0316 -2.8782
5.7142 -3.0486
5.4569 -3.1693
5.0811 -3.3421
4.5777 -3.5574
4.3444 -3.6544
4.0541 -3.7763
3.6636 -3.9426
];

rref = ref(:,1);
zref = ref(:,2);