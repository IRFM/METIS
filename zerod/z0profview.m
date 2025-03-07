% visulaisation 3D des profils 0D
zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo = post.z0dinput.geo;
t    = zs.temps;
% profil de courant et de puissance
swlh = abs(zs.plh ./ max(1,zs.pohm)) > 0.2;
pfweh = max(0,min(1,op0d.fwcd ~=0)) .* zs.picrh_th;
picrh = zs.picrh_th - pfweh;
if op0d.zeff >2
	zeffin = zs.zeffsc;
else
	zeffin = cons.zeff;
end
if ~isfield(op0d,'Sn_fraction')
    Sn_fraction = 0;
else
    Sn_fraction = op0d.Sn_fraction;
end


profli = post.profil0d;
xli    = profli.xli;
jfus   = profli.jfus;
sfus   = profli.salf;
pfus   = profli.pfus;
t      = profli.temps;

% indicateur Greenwald en haut du piedestal
ne_ped        = interp1(profli.temps, profli.nep(:,end-1),cons.temps,'nearest','extrap');
flag_negr_ped = zs.negr;
flag_negr_ped(ne_ped < zs.negr) = NaN;

h = findobj(0,'type','figure','tag','z0profview');
if isempty(h)
       h=figure('tag','z0profview');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')


flag = 0;
k = 3;
subplot(k,k-1,1)
zplotprof(gca,t,xli,profli.tep./1e3,'color','r');
zplotprof(gca,t,xli,profli.tip./1e3,'color','b');
title('T_e (r) & T_i (b) [keV]')
subplot(k,k-1,2)
zplotprof(gca,t,xli,profli.nep./1e19,'color','r');
zplotprof(gca,t,xli,profli.nip./1e19,'color','b');
zplotprof(gca,t,xli,profli.nhep./1e19,'color','g');
zplotprof(gca,t,xli,profli.n1p./1e19,'color','m');
zplotprof(gca,t,xli,profli.nzp./1e19,'color','k');
if Sn_fraction > 0
    zplotprof(gca,t,xli,(1 - Sn_fraction ) .* profli.nwp./1e19,'color','c');
    zplotprof(gca,t,xli,Sn_fraction .* profli.nwp./1e19,'color','c','linestyle','--');
else
    zplotprof(gca,t,xli,profli.nwp./1e19,'color','c');
end
zplotprof(gca,cons.temps,ones(size(cons.temps)) * cat(2,xli(end-1), xli(end-1)),cat(2,flag_negr_ped,flag_negr_ped) ./ 1e19,'color','r','marker','o','linestyle','none');
if Sn_fraction > 0
    title('n_e (r), n_i (b), n_H_e (g), n_H_D_T (m), n_i_m_p (k), n_W (c), n_{Sn} (c--) & Greenwald limit (o) [10^1^9 m^-^3]')
else
    title('n_e (r), n_i (b), n_H_e (g), n_H_D_T (m), n_i_m_p (k), n_W (c) & Greenwald limit (o) [10^1^9 m^-^3]')
end
subplot(k,k-1,3)
zplotprof(gca,t,xli,profli.qjli,'color','b');
s  = pdederive(xli,profli.qjli,0,2,2,1) .* (ones(size(t)) * xli) ./ profli.qjli;
zplotprof(gca,t,xli,s,'color','r');
hold on
hhp = plot([0,1],[1,1],'g',[0,1],[1.5,1.5],'g-.',[0,1],[2,2],'g:');
set(hhp,'linewidth',0.5);
title('q (b) & s (r)')
set(gca,'ylim',[-1,Inf]);
subplot(k,k-1,4)
zplotprof(gca,t,xli,profli.jboot./1e6,'color','b');
zplotprof(gca,t,xli,profli.jlh./1e6,'color','r');
zplotprof(gca,t,xli,profli.jnbicd./1e6,'color','m');
zplotprof(gca,t,xli,profli.jeccd./1e6,'color','c');
zplotprof(gca,t,xli,profli.jfwcd./1e6,'color','g');
zplotprof(gca,t,xli,profli.jfus./1e6,'color','k');
zplotprof(gca,t,xli,profli.jrun./1e6,'color','b','linestyle','none','marker','.');
title('J_N_I (b:boot, r:LH, m:NBICD, c:ECCD, g:FWCD, k:fus & dot-b:jrun) [MA m^-^2]')
subplot(k,k-1,5)
%zplotprof(gca,t,xli,profli.jli./1e6,'color','b');
zplotprof(gca,t,xli,profli.jeff./1e6,'color','b');
zplotprof(gca,t,xli,(profli.jboot+profli.jlh+profli.jnbicd+profli.jeccd+profli.jfwcd+profli.jfus+profli.jrun)./1e6,'color','r');
title('<J> (b) & J_N_I (r) [MA m^-^2]')
xlabel('x (normalized radius)');
subplot(k,k-1,6)
padd = profli.pohm + profli.pnbi + profli.plh + profli.pecrh  + profli.picrh + profli.pfweh + profli.pfus;
zplotprof(gca,t,xli,padd./1e6,'color','r');
zplotprof(gca,t,xli,profli.pfus./1e6,'color','b');
zplotprof(gca,t,xli,profli.pohm./1e6,'color','k');
zplotprof(gca,t,xli,profli.prad./1e6,'color','g');
zplotprof(gca,t,xli,profli.pbrem./1e6,'color','m');
zplotprof(gca,t,xli,profli.pcyclo./1e6,'color','c');
title('Padd (r), Pfus (b), Pohm (k), Pline (g), Pbrem (m) & Pcyclo (c) [MW m^-^3]')
xlabel('x (normalized radius)');




