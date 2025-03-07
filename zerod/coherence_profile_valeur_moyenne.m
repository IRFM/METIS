% coherence valeur moyenne profile
nim  = trapz(post.profil0d.xli,post.profil0d.nip .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);
nhem = trapz(post.profil0d.xli,post.profil0d.nhep .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);
nimpm = trapz(post.profil0d.xli,post.profil0d.nzp .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);
nwm = trapz(post.profil0d.xli,post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);
n1m = trapz(post.profil0d.xli,post.profil0d.n1p .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);


fullscreen = get(0,'ScreenSize');
h = findobj(0,'type','figure','tag','cavep');
if isempty(h)
	h=figure('tag','cavep');
else
	figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen);


subplot(2,1,1)
semilogy(post.profil0d.temps,nim,'b',post.zerod.temps,post.zerod.nim,'b-.',post.profil0d.temps,nhem,'r',post.zerod.temps,post.zerod.nhem,'r-.', ...
         post.profil0d.temps,nimpm,'c',post.zerod.temps,post.zerod.nimpm,'c-.',post.profil0d.temps,nwm,'m',post.zerod.temps,post.zerod.nwm,'m-.', ...
         post.profil0d.temps,n1m,'k',post.zerod.temps,post.zerod.n1m,'k-.');

xlabel('time (s)');
ylabel('m^{-3}');
if post.z0dinput.option.Sn_fraction > 0
    legend('<n_i> from profile','<n_i> 0D METIS','<n_{He}> from profile','<n_{He}> 0D METIS','<n_{imp}> from profile','<n_{imp}> 0D METIS', ...
        '<n_{W + Sn}> from profile','<n_{W + Sn}> 0D METIS','<n_{HDT}> from profile','<n_{HDT}> 0D METIS');
else
    legend('<n_i> from profile','<n_i> 0D METIS','<n_{He}> from profile','<n_{He}> 0D METIS','<n_{imp}> from profile','<n_{imp}> 0D METIS', ...
        '<n_{W}> from profile','<n_{W}> 0D METIS','<n_{HDT}> from profile','<n_{HDT}> 0D METIS');
end


nim0d = interp1(post.zerod.temps,post.zerod.nim,post.profil0d.temps,'pchip',NaN);
nem0d = interp1(post.zerod.temps,post.zerod.nem,post.profil0d.temps,'pchip',NaN);
nhem0d = interp1(post.zerod.temps,post.zerod.nhem,post.profil0d.temps,'pchip',NaN);
nimpm0d = interp1(post.zerod.temps,post.zerod.nimpm,post.profil0d.temps,'pchip',NaN);
nwm0d = interp1(post.zerod.temps,post.zerod.nwm,post.profil0d.temps,'pchip',NaN);
n1m0d = interp1(post.zerod.temps,post.zerod.n1m,post.profil0d.temps,'pchip',NaN);
subplot(2,1,2)

plot(post.profil0d.temps,(nim - nim0d) ./ nem0d,'b',post.profil0d.temps,(nhem - nhem0d) ./ nem0d,'r', ...
     post.profil0d.temps,(nimpm - nimpm0d) ./ nem0d,'c',post.profil0d.temps,(nwm - nwm0d) ./ nem0d,'m', ...
     post.profil0d.temps,(n1m - n1m0d) ./ nem0d,'k');


