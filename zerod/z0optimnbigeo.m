% function for the optimisation of NBI geomery used in sycomore2metis
% shine_throught_limite = maximum power lost in shine throught (W).
% for testing:
% x0 = cat(2,post.z0dinput.option.rtang,post.z0dinput.option.rtang2,post.z0dinput.option.zext,post.z0dinput.option.zext);
% [objective,rtang,rtang2,zext,zext2,falign,fni,soq,shine] = z0optimnbigeo(x0,post,1e6)
function [objective,rtang,rtang2,zext,zext2,ialign,fini,soq,shine,qnot,qmin] = z0optimnbigeo(x0,post,shine_through_limit,mode_ni)


% decodage
rtang  = x0(1);
rtang2 = x0(2);
zext   = x0(3);
zext2  = x0(4);
% for steady state shot 
if length(x0) > 4
  fpnbi      = x0(5);
  if length(x0) > 5
    fpnbi      = fpnbi + sqrt(-1) .* x0(6);  
  end
else
  fpnbi      = 1;
end
if imag(fpnbi) == 0
  fpnbi_mul = fpnbi .* (1 + sqrt(-1));
else
  fpnbi_mul = fpnbi;
end
if nargin < 4
  mode_ni = 2;
end
%disp(x0)
% extraction
zs = post.zerod;
profli = post.profil0d;
option = post.z0dinput.option;
geo    = post.z0dinput.geo;
cons   = post.z0dinput.cons;
ve     = ones(size(profli.xli));
% resampling if useful
temps = post.profil0d.temps;
vt = ones(size(temps));
if length(zs.temps) ~= length(temps)
    noms = fieldnames(zs);
    t    = zs.temps;
    for k=1:length(noms)
	var = zs.(noms{k});
	if size(var,1) == length(t)
	  zs.(noms{k}) = interp1(t,var,temps,'nearest','extrap');  
	end
    end
    zs.temps = temps;
    t    = cons.temps;
    noms = fieldnames(cons);
    for k=1:length(noms)
	var = cons.(noms{k});
	if size(var,1) == length(t)
	  cons.(noms{k}) = interp1(t,var,temps,'nearest','extrap');  
	end
    end
    cons.temps = temps;
    noms = fieldnames(geo);
    for k=1:length(noms)
	var = geo.(noms{k});
	if size(var,1) == length(t)
	  geo.(noms{k}) = interp1(t,var,temps,'nearest','extrap');  
	end
    end
    
end

% change in parameters
option.rtang  = rtang;
option.rtang2 = rtang2;
option.zext   = zext;
option.zext2  = zext2;

% useful data
zu1 = (option.zimp + option.rimp .* option.zmax);
zu2 = (option.zimp .^ 2 + option.rimp .* option.zmax .^ 2);

if (option.W_effect == 1) && isfield(profli,'nep') && isfield(profli,'nzp')
    [profli.nwp,zs.nwm,zu1w,zu2w] = z0acctungsten(option,geo,cons,zs,profli);
    zu1 = zu1 + zu1w; 
    zu2 = zu2 + zu2w;
else
    profli.nwp = zeros(length(cons.temps),21);
    zs.nwm     = zeros(length(cons.temps),1);
end

% factor for bootstrap current
pnbitot  = (real(cons.pnbi) + imag(cons.pnbi));
fpnbi    = trapz(temps,real(fpnbi_mul) .* real(cons.pnbi) + imag(fpnbi_mul) .* imag(cons.pnbi),1) ./ max(1,trapz(temps,pnbitot,1));
pnbi_new = real(fpnbi_mul) .* real(cons.pnbi) + sqrt(-1) .* imag(fpnbi_mul) .* imag(cons.pnbi);

% recomputation of current sources
[ilh,ifwcd,ieccd,inbicd,ecrit_nbi,taus_nbi,etalh0,etalh1,zs.etalh, ...
             xnbi,piqnbi,frnbi,mu0_nbi,xeccd,xlhout,dlhout,profli] =  ...
             zicd0(temps,zs.plh,zs.picrh,zs.pecrh,pnbi_new,zs.te0 ./ (1 + zs.ate),zs.nem,zs.zeff, ...
                   geo.R,geo.a,geo.K,zs.RR,zs.vloop,zs.ip,zs.ate,zs.qa,zs.betaptot + zs.li/2, ...
                   option,zs.hmhd,zs.ane,zs.nebord,zs.tebord,cons.ftnbi,zs.meff,zs.taue,geo.b0,zs.nbar,zs.qmin,zs.d0,zs.modeh, ...
		   zs.nDm,zs.nTm,zs.n1m,zs.nhem,cons.xece,zs.efficiency,zs.nimpm,zu1,zu2,zs.inbicd(1),profli,option.Sn_fraction);

% normalised
jnbicd   = real(profli.jnbishape);
jnbicd   = jnbicd .* ((real(inbicd)./  max(1,trapz(post.profil0d.xli,post.profil0d.spr .* jnbicd,2))) * ve);
if any(imag(inbicd))
  jnbicd2   = imag(profli.jnbishape);
  jnbicd2   = jnbicd2 .* ((imag(inbicd)./  max(1,trapz(post.profil0d.xli,post.profil0d.spr .* jnbicd2,2))) * ve);
  jnbicd   = jnbicd  + jnbicd2;
end       

% recompute modified current and safety factor profile
jeff    = post.profil0d.jeff - post.profil0d.jnbicd + jnbicd + (sqrt(fpnbi) - 1) .*  post.profil0d.jboot; % linearised approximation
jni     = post.profil0d.jni  - post.profil0d.jnbicd + jnbicd + (sqrt(fpnbi) -1 ) .*  post.profil0d.jboot;
ini     = max(0,zs.ini - real(zs.inbicd) - imag(zs.inbicd) + real(inbicd) + imag(inbicd) + (sqrt(fpnbi) -1) .*  zs.iboot);
jeff_blind = post.profil0d.jeff;
jeff_blind(jeff_blind < 1) = 1;
iota    = max(0.01,1./ post.profil0d.qjli .* (1 + tanh((jeff - post.profil0d.jeff) ./jeff_blind)));
q       = 1 ./ iota;
q(:,end)= post.profil0d.qjli(:,end);
s       = pdederive(post.profil0d.xli,q,2,2,2,1) ./ q .* (vt  * post.profil0d.xli);
soq     = trapz(post.profil0d.xli,abs(s)./q.* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2);

%figure(21);clf;plot(soq);drawnow
% current alignment 
ialign = 1 - trapz(profli.xli,profli.nep ./ max(1,profli.tep) .* abs(profli.jeff - jni),2) ./ ...
      			max(eps,trapz(profli.xli,profli.nep ./ max(1,profli.tep) .* abs(profli.jeff),2));
ialign = min(1,max(0,ialign));

% non inductive fraction 
fini = ini ./ max(1,zs.ipar);

% shine throught
shine = (1 - real(frnbi)) .* real(cons.pnbi) + (1 - imag(frnbi)) .* imag(cons.pnbi);
% stiffness of the penality (gives 1 for shine through exceeding the mimit by 1%).
kshine = 10;
% computation of objective to be minimized
% main objective is fini with not to bad alignement. Confinement is a secondary goal to prevent unwanted solutions
% the shine through limit is mandatory
%objective = poid_f_ini .* abs(fini - 1) - min(1,ialign) - (1/3) .* soq + exp(kshine .* xshine) .* (xshine > 0);

% ponderation
if (length(x0) > 4) && (sum(pnbitot < (0.9 .* max(pnbitot))) > 2)
  % for steady state only flattot is taken into account
  pnbitot(pnbitot < (0.9 .* max(pnbitot))) = 0;
end
%objective = trapz(temps,objective .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
fini   = trapz(temps,fini .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
ialign = trapz(temps,ialign .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
shine  = trapz(temps,shine .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
soq    = trapz(temps,soq .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
qnot   = trapz(temps,0.5 .* (q(:,1) + q(:,2)) .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
%qmin   = trapz(temps,0.5 .* min(q,[],2) .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
qmin   = trapz(temps,min(q,[],2) .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
xshine = (shine - shine_through_limit) ./ shine_through_limit;
%objective = - min(fini,1) - min(1,ialign - max(0,fini -1)) - min(1,0.5 .* (fini + ialign)) .* (min(1,soq) + min(2.1,qmin) .* (length(x0) > 4)) +  ...
switch mode_ni
case 1
  pen       = 3 .* max(0,fini -(2/3));
  objective = pen - max(0,1 - pen) .* (min(1,ialign) +  3 .* min(1,soq) + min(1.2,qmin)) +  ...
              exp(kshine .* xshine) .* (xshine > 0) + (qnot - 1.2) .* (qnot > 1.2);
otherwise
  objective = - min(fini,1) - min(1,ialign - max(0,fini - 1)) - min(1,0.5 .* (fini + ialign)) .* (min(1,soq) + min(1.2 + (length(x0) > 4),qmin)) +  ...
                exp(kshine .* xshine) .* (xshine > 0) + (qnot - 3) .* (qnot > 3);
end

% penality for  inefficient power increase
%obdisp = objective; 
%fpnbi_mul
%objective = objective + 0.5 .* (real(fpnbi_mul) + imag(fpnbi_mul)) .* abs(2 - ialign - min(1,fini));
fini_mem   = zs.ini ./ max(1,zs.ipar);
fini_mem   = trapz(temps,fini_mem .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));%fini_mem - fini
ialign_mem = trapz(temps,zs.ialign .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));%ialign_mem  - ialign
%disp([objective,max(0,real(fpnbi_mul) + imag(fpnbi_mul) - 2) .* (max(0,ialign_mem  - ialign)  +  max(0,fini_mem - fini)), ...
%      objective + max(0,real(fpnbi_mul) + imag(fpnbi_mul) - 2) .* (max(0,ialign_mem  - ialign)  +  max(0,fini_mem - fini))]);
objective  = objective + max(0,real(fpnbi_mul) + imag(fpnbi_mul) - 2) .* (max(0,ialign_mem  - ialign)  +  max(0,fini_mem - fini));

% penality for current drive in pedestal
inbi_int  = cumtrapz(post.profil0d.xli,post.profil0d.spr .* jnbicd,2);
iboot_int = cumtrapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jboot,2);
fiped     = (inbi_int(:,end) - inbi_int(:,end - 2)) ./ max(1,iboot_int(:,end) - iboot_int(:,end - 2));
%figure(21);clf;plot(fiped);drawnow
fiped     = trapz(temps,fiped .* pnbitot,1) ./ max(1,trapz(temps,pnbitot,1));
objective = objective + 3 .* max(0,fiped);
%  disp(x0););
%  disp([obdisp,objective,fini,ialign,xshine,soq]);
%  keyboard
%  if ~isfinite(objective)
%      keyboard
%  end
  