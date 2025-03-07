% fonction de calcul de l'estimation des sources de neutrons DD 
function profli = z0neutron_profil(zs,profli)

x = profli.xli;
temps = profli.temps;
ux  = 1 - x .^ 2;
ve  = ones(size(x));
vt  = ones(size(temps));
spr = profli.spr;
vpr = profli.vpr;
nep = profli.nep;
tep = profli.tep;
tip = profli.tip;
nip = profli.nip;
n1p = profli.n1p;
nD  = profli.n1p .* ((zs.nDm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
nT  = profli.n1p .* ((zs.nTm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
pnbi = real(profli.pnbi) + imag(profli.pnbi);
zeffp = profli.zeff;
pnbi_th = trapz(x,vpr .* pnbi,2);

profli.nD = nD;
profli.nT = nT;

if isfield(zs,'temps')
	ndd          = pchip(zs.temps,zs.ndd,temps);
	ndd_th       = pchip(zs.temps,zs.ndd_th,temps);
	ndd_nbi_th   = pchip(zs.temps,zs.ndd_nbi_th,temps);
	ndd_nbi_nbi  = pchip(zs.temps,zs.ndd_nbi_nbi,temps);
else
	ndd = zs.ndd;
	ndd_th = zs.ndd_th;
	ndd_nbi_th = zs.ndd_nbi_th;
	ndd_nbi_nbi = zs.ndd_nbi_nbi;
end


% fusion DD du plasma thermique (correction ,doit etre petit, l'enrichissement du milieu en T, H et He3 n'est pas prise en compte)
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(max(tip(:),13.6));
profli.sn0_th       = 0.5 .* nD .^ 2 .* reshape(dd_n.sv,size(tip)) .* (tip >= dd_n.timin) .*  (tip <= dd_n.timax);
profli.sn0_th       = profli.sn0_th .* ((ndd_th ./ max(1,trapz(x,profli.sn0_th .* vpr,2))) * ve);

% source de neutron beam-plasma
taus_nbi = 6.27e8 .* 2 .* tep .^ (3/2) ./ (nep./ 1e6) ./ 17;
profli.sn0_nbi_th  = nD .* taus_nbi .* pnbi;
profli.sn0_nbi_th  = profli.sn0_nbi_th .* ((ndd_nbi_th ./ max(1,trapz(x,profli.sn0_nbi_th .* vpr,2))) * ve);

% source de neutron beam-beam
profli.sn0_nbi_nbi  =  (taus_nbi .* pnbi ) .^ 2;
profli.sn0_nbi_nbi = profli.sn0_nbi_nbi .* ((ndd_nbi_nbi ./ max(1,trapz(x,profli.sn0_nbi_nbi .* vpr,2))) * ve);



  