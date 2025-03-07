% input npar0
npar0 = 2;
zs    = post.zerod
geo = post.z0dinput.geo;
cons = post.z0dinput.cons;
op0d =  post.z0dinput.option;
width = 32 .* 0.011; % pour TS largeur toroidal du grill LH
% masse effective
switch post.z0dinput.option.gaz
case 1
 	agaz    = 1 .* ones(size(zs.temps));
	zgaz    = 1 .* ones(size(zs.temps));
	nmaj    = zs.n1m;
case 2
 	agaz    = 2 .* ones(size(zs.temps));;
	zgaz    = 1 .* ones(size(zs.temps));
	nmaj    = zs.n1m;
case 3
 	agaz   = (3 .* post.z0dinput.cons.iso + 2) ./ (post.z0dinput.cons.iso + 1);
	zgaz   = 1 .* ones(size(zs.temps));
	nmaj   = zs.n1m;
case 5
 	agaz   = (3 .* real(post.z0dinput.cons.iso) + 2) ./ (real(post.z0dinput.cons.iso) + 1);
	zgaz   = (2 .* real(post.z0dinput.cons.iso) + 1) ./ (real(post.z0dinput.cons.iso) + 1);
	nmaj   = zs.n1m;
    warning('nHe3onD & nTonD not yet implemented !');
case 11
 	agaz   = (11 .* post.z0dinput.cons.iso + 1) ./ (post.z0dinput.cons.iso + 1);
	zgaz   = (5 .* post.z0dinput.cons.iso + 1) ./ (post.z0dinput.cons.iso + 1);
	nmaj   = zs.n1m;
otherwise
 	agaz    = 4 .* ones(size(cons.temps));
	zgaz    = 2 .* ones(size(zs.temps));
	nmaj    = zs.nim;
end
% profil de courant et de puissance
swlh = abs(zs.plh ./ max(1,zs.pohm)) > 0.2;
pfweh = max(0,min(1,op0d.fwcd ~=0)) .* zs.picrh_th;
picrh = zs.picrh_th - pfweh;

[void.iboot,void.iohm,void.pohm,void.RR,void.vloop,void.qa,void.q95,void.qmin,void.q0,void.betap,void.piqj, void.wbp,void.dwbpdt, ...
 void.asser,void.tauj,void.li,void.tauip,void.hitb,void.xitb,void.ate,void.aitb,void.hmhd,voidete0,woidfwcorr,xli,profli] = ...
                   zboot0(op0d.vloop,op0d.modeh,op0d.vref,cons.temps,zs.nem,zs.tem,zs.zeff,geo.R,geo.a,geo.K,geo.d,geo.b0, ...
                   zs.ip,zs.icd,zs.li,zs.wth,zs.ane,zs.ate,zs.vp,zs.sp,zs.sext,zs.hitb,1,zs.taue,zs.tauip, ...
                   op0d.limode,zs.peri,zs.inbicd,zs.xnbi,zs.piqnbi,zs.ieccd,zs.xeccd,zs.ifwcd,zs.ilh,zs.xlh,zs.dlh,zs.ifus, ...
                   zs.xfus,zs.jxfus,zs.j0fus,1,zs.tebord,zs.nebord,swlh,zs.pped,pfweh,picrh, ...
                   zs.plh_th,zs.pecrh,zs.pfus_th,zs.pnbi_th,zs.tite,zs.w,zs.nim);

[x,flh,lc,hc,acc,landau] = z0lhacc(3.7e9,npar0,width,agaz,zgaz,zs.temps,xli,profli.nep,profli.tep,profli.qjli,profli.Raxe,profli.rmx,profli.spr,profli.vpr,profli.fdia,profli.bpol);
