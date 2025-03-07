% test 2 points
geo = post.z0dinput.geo;
option = post.z0dinput.option;
zs  = post.zerod;
% fraction perdue en volume dans la sol par rayonnement:
fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,zs.ploss)));
%
% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
%lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 5 .* zs.qa) .^ 2);   % le 5 pour tenir compte du point X
lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 7 .* zs.qa) .^ 2);   % le 7 pour tenir compte du point X (valeur etalonnee sur A.S Kukushkin et al NF 43 p 716-723)
switch option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lcpol;
case 3
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lclim;
otherwise
	lc  = lcx;
end
% puissance conduite a la separatrice
pl        = max(1,zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz);

if option.sol_lscale  == 0
    dsol        = geo.a ./ 100;
elseif option.sol_lscale  > 0
    dsol        = geo.a .* option.sol_lscale;
else
    dsol        = - geo.R .* option.sol_lscale;
end
vpr_tor = interp1(post.profil0d.temps,post.profil0d.vpr_tor(:,end),zs.temps,'pchip','extrap');
qpl = pl .* (1 - fesol) ./ (vpr_tor .* dsol);
% reference : P.C. Stangeby NF 51 (2011) 063001
[tebord_x,telim_x,nelim_x] = z0twopoints(qpl,zs.nebord,lc,zs.zeff,zs.meff,option.fpe,option.fei);
