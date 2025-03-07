% script de test de l'electroneutarlite
zs  = post.zerod;
cons   = post.z0dinput.cons;
times = cons.temps;
temps = post.profil0d.temps;
% modification des donnees
noms = fieldnames(cons);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(cons,nomc);
	valn  = interp1(times,val,temps,'linear');
	indbad      = find(any(~isfinite(valn),2));
	if ~isempty(indbad)
		valn(indbad,:) = ones(length(indbad),1) * val(end,:);
	end
	cons = setfield(cons,nomc,valn);
end
cons.temps = temps;
noms = fieldnames(zs);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(zs,nomc);
	if length(val) == length(times)
		val  = interp1(times,val,temps,'nearest');
		zs = setfield(zs,nomc,val);
	end
end

if isfield(post.z0dinput.option,'Sn_fraction') && (post.z0dinput.option.Sn_fraction > 0)
    zw1 = (1 - post.z0dinput.option.Sn_fraction) .* z0wavez(post.profil0d.tep) + post.z0dinput.option.Sn_fraction .*  z0snavez(post.profil0d.tep);
    zw2 = (1 - post.z0dinput.option.Sn_fraction) .* z0wavez(post.profil0d.tep) .^ 2 + post.z0dinput.option.Sn_fraction .*  z0snavez(post.profil0d.tep) .^2;
else
    zw1 = z0wavez(post.profil0d.tep);
    zw2 = z0wavez(post.profil0d.tep) .^ 2; 
end

profli = post.profil0d;
x      = profli.xli;
option = post.z0dinput.option;

zu2w = trapz(x,zw2 .*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x, post.profil0d.nzp .* post.profil0d.vpr,2);
zu1w = trapz(x,zw1 .*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x, post.profil0d.nzp .* post.profil0d.vpr,2);
zu1 = (option.zimp + option.rimp .* option.zmax);
%zu1 = zu1 + zu1w;
zu2 = (option.zimp .^ 2 + option.rimp .* option.zmax .^ 2);
%zu2 = zu2 + zu2w;

del = zs.nem - zs.n1m - 2 .* zs.nhem - zu1 .* zs.nimpm -  ...
      trapz(x,zw1 .*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x,post.profil0d.vpr,2);

zeffvol = zs.zeff;
dzeff = zs.nem  .* zeffvol - zs.n1m - 4 .* zs.nhem - zu2 .* zs.nimpm -  ...
        trapz(x,zw2.*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x,post.profil0d.vpr,2);
dni   = zs.nim - zs.n1m - zs.nhem - zs.nimpm - zs.nwm;


h = findobj(0,'type','figure','tag','electro1');
if isempty(h)
       h=figure('tag','electro1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(2,2,1);
plot(cons.temps,del./zs.nem);
ylabel('relative error on n_{e,averaged})');
subplot(2,2,2);
plot(cons.temps,dzeff./zs.nem./zs.zeff);
ylabel('relative error on Zeff_{averaged})');
subplot(2,2,3);
plot(cons.temps,dni./zs.nem);
ylabel('relative error on n_{ions,averaged})');
%nhemp = trapz(profli.xli,profli.nhep .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2);
subplot(2,2,4);
plot(cons.temps,(zs.ne0 - profli.nep(:,1))./zs.ne0);
ylabel('relative error on n_{e,0})');
xlabel('time (s)');
edition2

h = findobj(0,'type','figure','tag','electro11');
if isempty(h)
       h=figure('tag','electro11');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(2,1,1)
plot(temps,zs.nem,'r',temps,zs.n1m + 2 .* zs.nhem + zu1 .* zs.nimpm +  ...
      trapz(x,zw1 .*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x,post.profil0d.vpr,2),'b-.');
ylabel('m^-3');
legend('n_{e,avaraged}','sum(Z_i * n_{i,averaged})');

subplot(2,1,2)
plot(temps,zs.zeff,'r',temps,(zs.n1m + 4 .* zs.nhem + zu2 .* zs.nimpm +  ...
        trapz(x,zw2 .*  post.profil0d.nwp .* post.profil0d.vpr,2) ./ trapz(x,post.profil0d.vpr,2)) ./ zs.nem,'b-.');
legend('Z_{eff}','sum(Z_i^2 * n_{i,averaged}) / n_{e,avaraged}');
edition2

dne = (profli.nep - profli.n1p - 2.* profli.nhep - zu1 .* profli.nzp - zw1 .*  post.profil0d.nwp ) ./ profli.nep;
dzeff = (profli.nep .* profli.zeff - profli.n1p - 4.* profli.nhep - zu2 .* profli.nzp - zw2 .*  post.profil0d.nwp) ./ profli.nep;
h = findobj(0,'type','figure','tag','electro2');
if isempty(h)
       h=figure('tag','electro2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(2,2,1)
plot(profli.xli,dne);
ylabel('relative error on n_{e,profile})');
subplot(2,2,2)
plot(profli.xli,dzeff);
ylabel('relative error on Zeff_{e,profile})');
xlabel('r/a');
edition2



h = findobj(0,'type','figure','tag','electro3');
if isempty(h)
       h=figure('tag','electro3');
else
       figure(h);
end   
clf
subplot(2,3,1)
plot(cons.temps,(zs.nem-trapz(profli.xli,profli.nep .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2))./ zs.nem,'o');
ylabel('relative error on n_{e,profile})');
subplot(2,3,2)
plot(cons.temps,(zs.n1m-trapz(profli.xli,profli.n1p .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2)) ./ zs.nem,'o');
ylabel('relative error on n_{HDT,profile})');
subplot(2,3,3)
plot(cons.temps,(zs.nhem-trapz(profli.xli,profli.nhep .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2)) ./ zs.nem,'o');
ylabel('relative error on n_{He,profile})');
subplot(2,3,4)
plot(cons.temps,(zs.nimpm-trapz(profli.xli,profli.nzp .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2)) ./ zs.nem,'o');
ylabel('relative error on n_{imp,profile})');
subplot(2,3,5)
plot(cons.temps,(zs.nwm-trapz(profli.xli,profli.nwp .* profli.vpr,2)./trapz(profli.xli,profli.vpr,2)) ./ zs.nem,'o');
ylabel('relative error on n_{W,profile})');
edition2

h = findobj(0,'type','figure','tag','electro31');
if isempty(h)
       h=figure('tag','electro31');
else
       figure(h);
end   
clf

ne_rec   = profli.n1p + 2 .* profli.nhep + zu1 .* profli.nzp +  zw1 .*  post.profil0d.nwp;
zeff_rec = profli.n1p + 4 .* profli.nhep + zu2 .* profli.nzp +  zw2  .*  post.profil0d.nwp;
zeff_rec = zeff_rec ./ ne_rec;

subplot(2,2,1)
zplotprof(gca,temps,x,profli.nep,'color','b');
zplotprof(gca,temps,x,ne_rec,'color','r');
xlabel('r/a')
ylabel('m^{-3}');
legend('n_{e,profile}','sum(Z_i * n_{i,profile})');

subplot(2,2,2)
zplotprof(gca,temps,x,profli.zeff,'color','b');
zplotprof(gca,temps,x,zeff_rec,'color','r');
xlabel('r/a')
ylabel('Z_{eff}');
legend('Z_{eff,profile}','sum(Z_i ^ 2 * n_{i,profile}) / n_{e,profile}');
%set(gca,'ylim',[1,Inf])
edition2
