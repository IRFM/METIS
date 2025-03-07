% script de run du rampup de ITER
% valeur pour la creation du run vide
b0    = 5.3;
R     = 6.2;
a     = 2;
K95     = 1.7;
d95     = 0.33;
ip    = 15e6;
nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
plh   = 0;
picrh = 0;
pecrh = 0;
pnbi  = 0;
zeff  = 1.23;
li    = 0.85;
hmore = 1.0;
iso   = 1;
ftnbi = 0.5;
xece  = 0.6;


% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);
d  = d95 ./ (0.95 .^ 2);

% limite de stabilite verticale
k95lim = 2.22 - 0.17 .* R ./ a;


% base temps 
temps    = cat(2,2:1:100,105:5:600,601:2:660)';
	
% creation du jeu de donnees
z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi, ...
           zeff,xece,hmore,iso,ftnbi,0, ...
           'gaz',3,'frhe0',0,'tauhemul',5,'ane',0,'vane',1, ...
	   'scaling',1,'l2hscalin',0,'modeh',1,'l2hmul',5,'fpped',0.5,'fprad',1/3, ...
	   'xiioxie',0,'kishape',3,'qdds',0.95,'kidds',3, ...
	   'vloop',0,'runaway',1,'modeboot',1,'vref',0,'li',li, ...
	   'zeff',0,'zmax',10,'zimp',4,'rimp',0.06,'matthews',0,'frad',0.1,'rw',0.7, ...
	   'angle_ece',90,'synergie',0, ...
	   'sens',1,'angle_nbi',60,'einj',1e6,'rtang',R, ...
	   'lhmode',3,'etalh',0.75,'npar0',2,'freqlh',3.7,'wlh',0.545,'xlh',0,'dlh',0, ...
	   'fwcd',0,'mino','T','cmin',1,'nphi',25,'freq',72, ...
	   'sitb',0,'tae',0,'smhd',0,'tmhd',0,'rip',0,'configuration',2);

% creation des consignes  temporelles
tv     = [1.6	4.61	7.82	11.38	15.24	19.52	24.17	29.37	35.25	42.12	49.26	56.21 70    100    600    660];
ipv    = [0.40	1.50	2.50	3.50	4.50	5.50	6.50	7.50	8.50	9.50	10.50	11.50 13    15     15     2];
Rv     = [7.480	6.889	6.656	6.524	6.443	6.373	6.324	6.230	6.227	6.219	6.219	6.215 6.215 6.215  6.215  6.8];
av     = [0.800	1.393	1.621	1.751	1.831	1.906	1.93	1.950	1.953	1.96	1.966	1.972 1.972 1.972  1.972  1];
K95v   = [1.005	1.008	1.113	1.230	1.343	1.438	1.527	1.631	1.641	1.648	1.662	1.672 1.7   1.7    1.7    1];
d95v   = [0.026	0.058	0.098	0.137	0.183	0.228	0.284	0.30	0.31	0.322	0.326	0.327 0.327 0.327  0.327  0.327];
pnbiv  = [0	0	0	0	0	0	0	0	0	0	0	0     15    28     28     28];
picrhv = [0	0	0	0	0	0	8	17	17	17	17	17    17    17     17     17];
plhv   = [0	1       3	5	5	5	5	5	5	5	5	5     5     5      5      5];
%nv     = [0.05	0.1	0.2	0.3	0.4	0.5	0.4	0.4	0.5	%0.5	0.5	0.6   0.7   1      1      0.2];
nv       =ipv ./ 13;
isov   = [0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.5	1	1     1     1      1      0.6];

Kv  = 2 .* (K95v - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);
dv  = d95v ./ (0.95 .^ 2);

z0dinput.cons.ip = pchip(tv,ipv,temps) .* 1e6;
z0dinput.cons.picrh = pchip(tv,picrhv,temps) .* 1e6;
z0dinput.cons.pnbi = pchip(tv,pnbiv,temps) .* 1e6;
z0dinput.cons.plh = pchip(tv,plhv,temps) .* 1e6;
z0dinput.cons.pecrh = 0.*pchip(tv,plhv,temps) .* 1e6;
z0dinput.geo.a   = pchip(tv,av,temps);
z0dinput.geo.R   = pchip(tv,Rv,temps);
z0dinput.geo.b0  = b0.*R./z0dinput.geo.R;
z0dinput.geo.K   = pchip(tv,Kv,temps);
z0dinput.geo.d   = pchip(tv,dv,temps);
z0dinput.cons.nbar = pchip(tv,nv,temps).*1e20;
z0dinput.cons.iso = pchip(tv,isov,temps);

z0dinput.cons.zeff  = 1.7 + 2.3 .*  (min(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) .^ 2.6;





% utilisation des moments
if isfield(z0dinput.exp0d,'Rsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end

% appel a METIS
[zs,infovoid,profli] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;
op0d =post.z0dinput.option;
z0plotsc
