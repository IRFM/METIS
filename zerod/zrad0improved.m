% function to compute improved radiative power near the edge of the plasma
% [pradhep,pradimpp,pradmaxp,pradwp,pradhd] =  ...
% zrad0improved(post.profil0d,post.z0dinput.option.zimp,post.z0dinput.option.zmax,post.z0dinput.option.rimp,post.z0dinput.option.noncoronal, ...
% interp1(post.zerod.temps,post.zerod.taue,post.profil0d.temps,'linear','extrap'))
function [pradhep,pradimpp,pradmaxp,pradwp,pradhdt,pradb] = zrad0improved(profil0d,zimp,zmax,rimp,noncoronal,taue,Sn_fraction,te_max_interp,nboron)
	
% adaptation pour le vrai calcul du rayonnement non coronal
taup      = imag(taue);
taue      = real(taue);
taup(taup == 0) = taue(taup == 0);
% calcul du temps local de residence
% calcul du temps local de residence
if isfield(profil0d,'nep') && isfield(profil0d,'tep')
    pen = profil0d.nep .* profil0d.tep;
    pen = pen ./ (max(eps,max(pen,[],2)) * ones(size(profil0d.xli)));
    tau_resid = (taue * ones(size(profil0d.xli))) .* pen + (min(taue,taup) * ones(size(profil0d.xli))) .* (1 - pen);
else
  tau_resid = taue * ones(size(profil0d.xli));
end

%cooling rate	
[ua,uz,post]=z0coefpost;

% indice des impuretes
indhe  =  find(ua==4 & uz == 2);
dz     = abs(uz - zimp);
indimp =  min(find(dz == min(dz)));
dz     = abs(uz - zmax);
indmax =  min(find(dz == min(dz)));
dz     = abs(uz - 5);
indb =  min(find(dz == min(dz)));


% profils normalises
tep    = profil0d.tep ./ 1e3;
nep    = profil0d.nep ./ 1e6;
nhep   = profil0d.nhep ./ 1e6;
nzp    = profil0d.nzp ./ 1e6;
nwp    = profil0d.nwp ./ 1e6;

% for boron
ne         = trapz(profil0d.xli,profil0d.nep,2);
nboronp    = profil0d.nep .* ((nboron ./ max(1e13,ne)) * ones(1,size(profil0d.nep,2))) / 1e6;

% nouvelle coordonnee radiale 
nbp  = 301;
beta = 3;
xl   = linspace(1,0,nbp);
xnew = 1 - (1 - exp(beta .* xl)) ./ (1 - exp(beta));
xnew = cat(2, - xnew(2), xnew,2.* xnew(end) - xnew(end-1));   

% interpolation
tep  = max(1e-3,interp1(profil0d.xli',tep',xnew,'pchip','extrap')');
nep  = max(0,interp1(profil0d.xli',nep',xnew,'pchip','extrap')');
nhep = max(0,interp1(profil0d.xli',nhep',xnew,'pchip','extrap')');
nzp  = max(0,interp1(profil0d.xli',nzp',xnew,'pchip','extrap')');
nwp  = max(0,interp1(profil0d.xli',nwp',xnew,'pchip','extrap')');
n1p  = max(0,interp1(profil0d.xli',profil0d.n1p',xnew,'pchip','extrap')');
nboronp   = max(0,interp1(profil0d.xli',nboronp',xnew,'pchip','extrap')');
vpr  = max(0,interp1(profil0d.xli',profil0d.vpr',xnew,'pchip','extrap')');
tau_resid  = max(0,interp1(profil0d.xli',tau_resid',xnew,'pchip','extrap')');


lzhe   = post(indhe).lz .* 0.1;
tzhe   = post(indhe).te;
lzhe   = cat(2,1e-38,lzhe,lzhe(end) + 5.355e3 .* 4  .*(sqrt(te_max_interp)-sqrt(tzhe(end))) .* 1e-28);
tzhe   = cat(2,1e-4,tzhe,te_max_interp);

lzimp    = post(indimp).lz .* 0.1;
tzimp    = post(indimp).te;
lzimp    = cat(2,1e-38,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max_interp)-sqrt(tzimp(end))) .* 1e-28);
tzimp    = cat(2,1e-4,tzimp,te_max_interp);

lzmax    = post(indmax).lz .* 0.1;
tzmax    = post(indmax).te;
lzmax    = cat(2,1e-38,lzmax,lzmax(end) + 5.355e3 .* uz(indmax)  .*(sqrt(te_max_interp)-sqrt(tzmax(end))) .* 1e-28);
tzmax    = cat(2,1e-4,tzmax,te_max_interp);

lzb   = post(indb).lz .* 0.1;
tzb   = post(indb).te;
lzb   = cat(2,1e-38,lzb,lzb(end) + 5.355e3 .* 4  .*(sqrt(te_max_interp)-sqrt(tzb(end))) .* 1e-28);
tzb   = cat(2,1e-4,tzb,te_max_interp);


if noncoronal == -1
	rzhe  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,'He') .* 1e12;
	if any(~isfinite(rzhe(:))) 
	    rzhe_p    = reshape(10 .^ pchip(log10(tzhe),log10(lzhe),log10(tep(:))),size(tep));	
	    rzhe(~isfinite(rzhe)) = rzhe_p(~isfinite(rzhe));
	end
	rzimp  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,zimp) .* 1e12;
	if any(~isfinite(rzimp(:))) 
	    rzimp_p   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp),log10(tep(:))),size(tep));
	    rzimp(~isfinite(rzimp)) = rzimp_p(~isfinite(rzimp));
	end
	rzmax  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,zmax) .* 1e12;
	if any(~isfinite(rzmax(:))) 
	    rzmax_p   = reshape(10 .^ pchip(log10(tzmax),log10(lzmax),log10(tep(:))),size(tep));		
	    rzmax(~isfinite(rzmax)) = rzmax_p(~isfinite(rzmax));
    end
    
    if any(nbp(:) > sqrt(eps))
        fprintf('b');  %noncoronal radiative data are missing for boron
    end
    rzb    = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tep(:))),size(tep));
    
elseif noncoronal == 1
	rzhe   = z0noncronal(tzhe,lzhe,xnew,tep,nep,taue);
	rzimp  = z0noncronal(tzimp,lzimp,xnew,tep,nep,taue);
	rzmax  = z0noncronal(tzmax,lzmax,xnew,tep,nep,taue);
    rzb    = z0noncronal(tzb,lzb,xnew,tep,nep,taue);
else
	rzhe    = reshape(10 .^ pchip(log10(tzhe),log10(lzhe),log10(tep(:))),size(tep));
	rzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp),log10(tep(:))),size(tep));
	rzmax   = reshape(10 .^ pchip(log10(tzmax),log10(lzmax),log10(tep(:))),size(tep));
    rzb     = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tep(:))),size(tep));
end

% dans ce cas on calcul que le rayonnement dans le plasma de coeur
pradhep   = nep .* nhep .* rzhe;
pradimpp  = nep .* nzp  .* rzimp;
pradmaxp  = nep .* nzp  .* rzmax .*  rimp;
pradb = nep .* nboronp .* rzb;

% cas du tungstene traite part
indw =  find(uz == 74,1);
lzw    = post(indw).lz .* 0.1;
tzw    = post(indw).te;
lzw    = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(te_max_interp)-sqrt(tzw(end))) .* 1e-28);
tzw    = cat(2,1e-4,tzw,te_max_interp);
if noncoronal == 1
	rzw  = z0noncronal(tzw,lzw,xnew,tep,nep,taue);
else
	rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
end
pradwp   = nep .* nwp .* rzw;
if Sn_fraction > 0
    % Sn case
    indsn =  find(uz == 50,1);
    lzsn    = post(indsn).lz .* 0.1;
    tzsn    = post(indsn).te;
    lzsn    = cat(2,1e-38,lzsn,lzsn(end) + 5.355e3 .* 50  .*(sqrt(te_max_interp)-sqrt(tzsn(end))) .* 1e-28);
    tzsn    = cat(2,1e-4,tzsn,te_max_interp);
    if noncoronal == 1
        rzsn  = z0noncronal(tzsn,lzsn,profli.xli,tep,nep,taue);
    else
        rzsn  = reshape(10 .^ pchip(log10(tzsn),log10(lzsn),log10(tep(:))),size(tep));
    end
    pradwp =  (1 - Sn_fraction) .* pradwp + Sn_fraction .*  nep .* nwp .* rzsn;
    
end



% rayonnement de HDT
% ajout de la recombinaison
% recombinaison (Halpha + ...)
trec    = [0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
srec    = [7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 1e-6;% m^3/s;
rrec    = reshape(exp(pchip(log(trec),log(srec),log(tep(:)))),size(tep));
%pradhdt_old = 5.355e3 .* (1e6 .* nep./ 1e20) .* (n1p ./ 1e20) .* sqrt(tep) + ...
%	      1e6 .*nep .* n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3

% amelioration du calcul du rayonnement de l'hydorgene
indhdt  = min(find(uz == 1));
lzhdt   = post(indhdt).lz .* 0.1;
tzhdt   = post(indhdt).te;
lzhdt   = cat(2,1e-38,lzhdt,lzhdt(end) + 5.355e3 .* 4  .*(sqrt(te_max_interp)-sqrt(tzhdt(end))) .* 1e-28);
tzhdt   = cat(2,1e-4,tzhdt,te_max_interp);
if noncoronal == 1
	rzhdt   = z0noncronal(tzhdt,lzhdt,xnew,tep,nep,taue);
else
	rzhdt    = reshape(10 .^ pchip(log10(tzhdt),log10(lzhdt),log10(tep(:))),size(tep));
end
pradhdt = nep .* (n1p ./ 1e6) .* rzhdt + (nep .* 1e6) .* n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3

%figure(21);hold on;plot(xnew,mean(pradhdt_old,1),'c',xnew,mean(pradhdt,1),'m');drawnow
%keyboard
% post traitement concervant l'energie lors du reechantionnage
pradhep   = resample_energy(xnew,pradhep,vpr,profil0d.xli);
pradimpp  = resample_energy(xnew,pradimpp,vpr,profil0d.xli);
pradmaxp  = resample_energy(xnew,pradmaxp,vpr,profil0d.xli);
pradwp    = resample_energy(xnew,pradwp,vpr,profil0d.xli);
pradhdt   = resample_energy(xnew,pradhdt,vpr,profil0d.xli);
pradb     = resample_energy(xnew,pradb,vpr,profil0d.xli);


% rechantillonage conservant l'energie
function [pout,err] = resample_energy(xnew,pnew,vpr,xout)

p_inte      = cumtrapz(xnew,pnew .* vpr,2);
v_inte      = cumtrapz(xnew,vpr,2);
p_inte_out  = interp1(xnew',p_inte',xout','nearest')';
v_inte_out  = interp1(xnew',v_inte',xout','nearest')';
pout_s      = interp1(xnew',pnew',xout','nearest')';
%pout        = (cat(2,pout_s(:,1),diff(p_inte_out,1,2) ./ diff(v_inte_out,1,2)) +  ...
%              cat(2,diff(p_inte_out,1,2) ./ diff(v_inte_out,1,2),pout_s(:,end))) ./ 2;
pdiff        = diff(p_inte_out,1,2) ./ diff(v_inte_out,1,2);
pdiff(~isfinite(pdiff)) = 0;
xdiff        = (xout(1:end-1)+ xout(2:end)) ./ 2;
xdiff        = cat(2,xnew(1:2),xdiff,xnew(end-1:end));
%pdiff        = cat(2,pnew(:,1:2),pdiff,pnew(:,end-1:end));
%pdiff        = cat(2,pnew(:,1),(pnew(:,2) + pdiff(:,1)) ./ 2 ,pdiff,(pnew(:,end - 1) + pdiff(:,end)) ./ 2,pnew(:,end));
pdiff        = cat(2,pnew(:,1),pdiff(:,1) ,pdiff,pdiff(:,end),pnew(:,end));
pout         = interp1(xdiff',pdiff',xout','pchip')';

vpr_out     = interp1(xnew',vpr',xout','nearest')';
err = abs(trapz(xout,vpr_out .* pout,2) - (p_inte_out(:,end) -  p_inte(:,2))) ./ (p_inte_out(:,end) - p_inte(:,2));


% no plot
return
% don't work
% renormalisation finale 

if any(~isfinite(pout(:))) || any(pout(:)< 0)
  keyboard
end


figure(21);clf
subplot(3,1,1)
plot(xnew,pnew,'r',xout,pout,'.b')
subplot(3,1,2)
t=1:length(p_inte_out(:,end));
plot(t,trapz(xout,vpr_out .* pout,2),'b',t,p_inte_out(:,end) -  p_inte(:,2),'r')
subplot(3,1,3)
plot(xnew,p_inte - p_inte(:,2) * ones(1,size(p_inte,2)),'r',xout,cumtrapz(xout,pout .* vpr_out,2),'.b')
drawnow

%keyboard
