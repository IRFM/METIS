% script pour tracer les contributions du rayonnement
improved = 0;
if post.z0dinput.option.matthews < 0
  improved = 1;
end
if ~isfield(post.z0dinput.option,'Sn_fraction')
    Sn_fraction = 0;
else
    Sn_fraction = post.z0dinput.option.Sn_fraction;
end
if ~isfield(post.z0dinput.option,'te_max')
    te_max = 105;
else
    te_max = max(100,post.z0dinput.option.te_max /1e3 +5);
end

taue = interp1(post.zerod.temps,post.zerod.taue,post.profil0d.temps,'linear','extrap');
taup = interp1(post.zerod.temps,post.zerod.taup,post.profil0d.temps,'linear','extrap');

zimp = post.z0dinput.option.zimp;
zmax = post.z0dinput.option.zmax;
rimp = post.z0dinput.option.rimp;
[ua,uz,postcoef]=z0coefpost;

indhe  =  find(ua==4 & uz == 2);
dz     = abs(uz - zimp);
indimp =  min(find(dz == min(dz)));
dz     = abs(uz - zmax);
indmax =  min(find(dz == min(dz)));
dz     = abs(uz - 5);
indb =  min(find(dz == min(dz)));

tep    = post.profil0d.tep ./ 1e3;
nep    = post.profil0d.nep ./ 1e6;
nhep   = post.profil0d.nhep ./ 1e6;
n1p    = post.profil0d.n1p ./ 1e6;
nzp    = post.profil0d.nzp ./ 1e6;
nwp    = post.profil0d.nwp ./ 1e6;

switch post.z0dinput.option.gaz
    case 11
        fnboronm = interp1(post.zerod.temps,post.zerod.nTm./ max(1e13,post.zerod.nem),post.profil0d.temps,'linear','extrap');
        nboronm = interp1(post.zerod.temps,post.zerod.nTm,post.profil0d.temps,'linear','extrap');
    otherwise
        fnboronm = zeros(size(tep,1),1);
        nboronm = zeros(size(tep,1),1);
end
% for boron
nbp    = post.profil0d.nep .* (fnboronm * ones(1,size(post.profil0d.nep,2))) / 1e6;


lzhe   = postcoef(indhe).lz .* 0.1;
tzhe   = postcoef(indhe).te;
lzhe   = cat(2,1e-38,lzhe,lzhe(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzhe(end))) .* 1e-28);
tzhe   = cat(2,1e-4,tzhe,te_max);

lzimp    = postcoef(indimp).lz .* 0.1;
tzimp    = postcoef(indimp).te;
lzimp    = cat(2,1e-38,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-28);
tzimp    = cat(2,1e-4,tzimp,te_max);

lzmax    = postcoef(indmax).lz .* 0.1;
tzmax    = postcoef(indmax).te;
lzmax    = cat(2,1e-38,lzmax,lzmax(end) + 5.355e3 .* uz(indmax)  .*(sqrt(te_max)-sqrt(tzmax(end))) .* 1e-28);
tzmax    = cat(2,1e-4,tzmax,te_max);

lzb   = postcoef(indb).lz .* 0.1;
tzb   = postcoef(indb).te;
lzb   = cat(2,1e-38,lzb,lzb(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzb(end))) .* 1e-28);
tzb   = cat(2,1e-4,tzb,te_max);


%  rzhe    = reshape(10 .^ pchip(log10(tzhe),log10(lzhe),log10(tep(:))),size(tep));
%  rzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp),log10(tep(:))),size(tep));
%  rzmax   = reshape(10 .^ pchip(log10(tzmax),log10(lzmax),log10(tep(:))),size(tep));

% calcul du temps local de residence
pen = post.profil0d.nep .* post.profil0d.tep;
pen = pen ./ (max(eps,max(pen,[],2)) * ones(size(post.profil0d.xli)));
tau_resid = (taue * ones(size(post.profil0d.xli))) .* pen + (min(taue,taup) * ones(size(post.profil0d.xli))) .* (1 - pen);

if post.z0dinput.option.noncoronal == -1
	rzhe  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,'He') .* 1e12;
	if any(~isfinite(rzhe(:)))
	    rzhe_p    = reshape(10 .^ pchip(log10(tzhe),log10(lzhe),log10(tep(:))),size(tep));	
	    %figure(21);subplot(2,2,1);loglog(tep(:),rzhe(:),'or',tep(:),rzhe_p(:),'b.');
	    rzhe(~isfinite(rzhe)) = rzhe_p(~isfinite(rzhe));
	end
	rzimp  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,zimp) .* 1e12;
	if any(~isfinite(rzimp(:))) 
	    rzimp_p   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp),log10(tep(:))),size(tep));
	    %figure(21);subplot(2,2,2);loglog(tep(:),rzimp(:),'or',tep(:),rzimp_p(:),'b.');
	    rzimp(~isfinite(rzimp)) = rzimp_p(~isfinite(rzimp));
	end
	rzmax  = zlightznoncoronal(tep .* 1e3,nep .* 1e6,tau_resid,zmax) .* 1e12;
	if any(~isfinite(rzmax(:))) 
	    rzmax_p   = reshape(10 .^ pchip(log10(tzmax),log10(lzmax),log10(tep(:))),size(tep));		
	    %figure(21);subplot(2,2,3);loglog(tep(:),rzmax(:),'or',tep(:),rzmax_p(:),'b.');
	    rzmax(~isfinite(rzmax)) = rzmax_p(~isfinite(rzmax));
	end
	%figure(21);subplot(2,2,4);loglog(tep(:),nep(:).* tau_resid(:) .* 1e6,'or',tep(:),1e16,'g',tep(:),1e19,'g');drawnow
    if any(nbp(:) > sqrt(eps))
        wraning('noncoronal radiative data are missing for boron');
    end
    rzb    = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tep(:))),size(tep));

elseif post.z0dinput.option.noncoronal == 1
	rzhe   = z0noncronal(tzhe,lzhe,post.profil0d.xli,tep,nep,taue);
	rzimp  = z0noncronal(tzimp,lzimp,post.profil0d.xli,tep,nep,taue);
	rzmax  = z0noncronal(tzmax,lzmax,post.profil0d.xli,tep,nep,taue);
    rzmax  = z0noncronal(tzmax,lzmax,post.profil0d.xli,tep,nep,taue);
    rzb     = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tep(:))),size(tep));

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
pradb = nep .* nbp .* rzb;

% cas du tungsten traite a part
indw =  find(uz == 74,1);
lzw    = postcoef(indw).lz .* 0.1;
tzw    = postcoef(indw).te;
lzw    = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(te_max)-sqrt(tzw(end))) .* 1e-28);
tzw    = cat(2,1e-4,tzw,te_max);
%rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
if post.z0dinput.option.noncoronal == 1
	rzw  = z0noncronal(tzw,lzw,post.profil0d.xli,tep,nep,taue);
else
	rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
end
pradwp   = nep .* nwp .* rzw;
if Sn_fraction > 0
    pradwp = (1 - Sn_fraction) .* pradwp;
    % Sn
    indsn =  find(uz == 50,1);
    lzsn    = postcoef(indsn).lz .* 0.1;
    tzsn    = postcoef(indsn).te;
    lzsn    = cat(2,1e-38,lzsn,lzsn(end) + 5.355e3 .* 74  .*(sqrt(te_max)-sqrt(tzsn(end))) .* 1e-28);
    tzsn    = cat(2,1e-4,tzsn,te_max);
    %rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
    if post.z0dinput.option.noncoronal == 1
        rzsn  = z0noncronal(tzsn,lzsn,post.profil0d.xli,tep,nep,taue);
    else
        rzsn   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
    end
    pradsnp   = Sn_fraction .* nep .* nwp .* rzsn;
else
    pradsnp = zeros(size(pradwp));
end

% rayonnement de HDT
% ajout de la recombinaison
% recombinaison (Halpha + ...)
trec    = [0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
srec    = [7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 1e-6;% m^3/s;
rrec    = reshape(exp(pchip(log(trec),log(srec),log(post.profil0d.tep(:)))),size(post.profil0d.tep));
%pradhdt = 5.355e3 .* (post.profil0d.nep./ 1e20) .* (post.profil0d.n1p ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3) + ...
%		post.profil0d.nep .* post.profil0d.n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3

% amelioration du calcul du rayonnement de l'hydrogene
indhdt  = min(find(uz == 1));
lzhdt   = postcoef(indhdt).lz .* 0.1;
tzhdt   = postcoef(indhdt).te;
lzhdt   = cat(2,1e-38,lzhdt,lzhdt(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzhdt(end))) .* 1e-28);
tzhdt   = cat(2,1e-4,tzhdt,te_max);
if post.z0dinput.option.noncoronal == 1
	rzhdt   = z0noncronal(tzhdt,lzhdt,post.profil0d.xli,tep,nep,taue);
else
	rzhdt    = reshape(10 .^ pchip(log10(tzhdt),log10(lzhdt),log10(tep(:))),size(tep));
end
pradhdt = nep .* n1p .* rzhdt + post.profil0d.nep .* post.profil0d.n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3

switch improved
case 1
    [~,~,~,pradsnp,~] = zrad0improved(post.profil0d,zimp,zmax,rimp,post.z0dinput.option.noncoronal,taue,1,te_max,nboronm);    
    %
    [pradhep,pradimpp,pradmaxp,pradwp,pradhdt,pradb] = zrad0improved(post.profil0d,zimp,zmax,rimp,post.z0dinput.option.noncoronal,taue,0,te_max,nboronm);
    pradwp   = (1 - Sn_fraction) .* pradwp;
    pradsnp   = Sn_fraction .* pradsnp;
end

% resample 0d data
clear zs
noms = fieldnames(post.zerod);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(post.zerod,nomc);
	if length(val) == length(post.zerod.temps)
		val  = interp1(post.zerod.temps,val,post.profil0d.temps,'nearest');
		zs.(nomc) = val;
	end
end

if post.z0dinput.option.gaunt == 1
    [gg1,Cbrem] = z0gaunt_brem(post.profil0d.tep ./ 1e3,1);
    gg2         = z0gaunt_brem(post.profil0d.tep ./ 1e3,2);
    gg5         = z0gaunt_brem(post.profil0d.tep ./ 1e3,5);
    gg_zimp     = z0gaunt_brem(post.profil0d.tep ./ 1e3,zimp);
    gg_zmax     = z0gaunt_brem(post.profil0d.tep ./ 1e3,zmax);
    gg_w        = z0gaunt_brem(post.profil0d.tep ./ 1e3,z0wavez(post.profil0d.tep));
    Cbrem       = Cbrem .* 1e40;

    dpbrem1  = Cbrem .* gg1 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.n1p ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremhe = Cbrem .* gg2 .* 4 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nhep ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremmax  = rimp .* Cbrem .* gg_zmax .* zmax .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nzp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremimp  = Cbrem .* gg_zimp .* zimp .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nzp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremb    = Cbrem .* gg5 .* 25 .* (post.profil0d.nep ./ 1e20) .*  ...
                 (post.profil0d.nep .* (fnboronm * ones(1,size(post.profil0d.nep,2))) ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    
    
    % cas du tungstene
    dpbremw  = Cbrem .* gg_w .* z0wavez(post.profil0d.tep) .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nwp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    if Sn_fraction > 0
        dpbremw = (1 - Sn_fraction) .* dpbremw;
        dpbremsn  = Sn_fraction .* Cbrem .* gg_w .* z0snavez(post.profil0d.tep) .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nwp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    else
        dpbremsn  = zeros(size(dpbremw));
    end
else
    dpbrem1  = 5.355e3 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.n1p ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremhe = 5.355e3 .* 4 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nhep ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremmax  = rimp .* 5.355e3 .* zmax .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nzp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremimp  = 5.355e3 .* zimp .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nzp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    dpbremb  = 5.355e3 .* 5 .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nep .* (fnboronm * ones(1,size(post.profil0d.nep,2))) ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);

    % cas du tungsten
    dpbremw  = 5.355e3 .* z0wavez(post.profil0d.tep) .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nwp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    if Sn_fraction > 0
        dpbremw = (1 - Sn_fraction) .* dpbremw;
        dpbremsn  = Sn_fraction .* 5.355e3 .* z0snavez(post.profil0d.tep) .^ 2 .* (post.profil0d.nep ./ 1e20) .* (post.profil0d.nwp ./ 1e20) .* sqrt(post.profil0d.tep ./ 1e3);
    else
        dpbremsn  = zeros(size(dpbremw));
    end
end
% correction relativiste de Pbrem
switch post.z0dinput.option.cor_rel_brem
    case 'improved'
        xrel = brem_rel_cor(post.profil0d.tep,post.profil0d.zeff);
        
    otherwise
        % P.E. Stott, Plasma. Physics and Controlled Fusion 47 (2005) 1305-1338 ref [15]
        tec2   = 511e3;% me c^2
        xrel   = (1 + 2 .* post.profil0d.tep ./ tec2) .* (1 + (2 ./ post.profil0d.zeff) .* (1 - (1 + post.profil0d.tep ./ tec2) .^ (-1)));
end
dpbrem1   = real(dpbrem1 .* xrel);
dpbremhe  = real(dpbremhe .* xrel);
dpbremmax = real(dpbremmax .* xrel);
dpbremimp = real(dpbremimp .* xrel);
dpbremw   = real(dpbremw .* xrel);
dpbremsn  = real(dpbremsn .* xrel);
dpbremb   = real(dpbremb .* xrel);



% for line brem alone
pbrem_tot = dpbrem1 + dpbremhe + dpbremmax + dpbremimp + dpbremw + dpbremsn + dpbremb;
warning off
coef_norm = post.profil0d.pbrem ./ pbrem_tot;
warning on
coef_norm(~isfinite(coef_norm)) = 1;
coef_norm(coef_norm < 0) = 1;
dpbrem1   = dpbrem1   .* coef_norm;
dpbremhe  = dpbremhe  .* coef_norm;
dpbremmax = dpbremmax .* coef_norm;
dpbremimp = dpbremimp .* coef_norm;
dpbremw   = dpbremw   .* coef_norm;
dpbremsn  = dpbremsn  .* coef_norm;
dpbremb   = dpbremb  .* coef_norm;

for k = 1:11
    % final cut to be compatible with zrad0
    pradhep   = max(1.01 .* dpbremhe,pradhep);
    pradimpp  = max(1.01 .* dpbremimp,pradimpp);
    pradmaxp  = max(1.01 .* dpbremmax,pradmaxp);
    pradwp    = max(1.01 .* dpbremw,pradwp);
    pradsnp   = max(1.01 .* dpbremsn,pradsnp);
    pradhdt   = max(1.01 .* dpbrem1,pradhdt);
    pradb   = max(1.01 .* dpbremb,pradb);

    % for line radiation + brem
    prad_tot = pradhep + pradimpp + pradmaxp + pradwp + pradsnp + pradhdt + pradb;
    warning off
    coef_norm = (post.profil0d.prad + post.profil0d.pbrem) ./ prad_tot;
    warning on
    coef_norm(~isfinite(coef_norm)) = 1;
    coef_norm(coef_norm < 0) = 1;
    %%%figure(21);hold on ;plot(post.profil0d.temps,coef_norm);drawnow
    pradhep   = coef_norm .* pradhep;
    pradimpp  = coef_norm .* pradimpp;
    pradmaxp  = coef_norm .* pradmaxp;
    pradwp    = coef_norm .* pradwp;
    pradsnp   = coef_norm .* pradsnp;
    pradhdt   = coef_norm .* pradhdt;
    pradb     = coef_norm .* pradb;
end

if post.z0dinput.option.matthews >= 1
  % warning
  warndlg('Radiative power in this simulation is computed using scaling law and not cooling rate.', 'Contribution from each ion species maybe inaccurate !');
end



fullscreen = get(0,'ScreenSize');
h = findobj(0,'type','figure','tag','z0plotrad');
if isempty(h)
       h=figure('tag','z0plotrad');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen)


zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradhep ./ 1e6,'color','g');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremhe ./ 1e6,'color','g','linestyle','-.');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradimpp ./ 1e6,'color','b');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremimp ./ 1e6,'color','b','linestyle','-.');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradmaxp ./ 1e6,'color','m');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremmax ./ 1e6,'color','m','linestyle','-.');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradwp ./ 1e6,'color','r');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremw ./ 1e6,'color','r','linestyle','-.');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradsnp ./ 1e6,'color','r','linestyle',':');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremsn ./ 1e6,'color','r','linestyle','--');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradhdt ./ 1e6,'color','c');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbrem1 ./ 1e6,'color','c','linestyle','-.');
switch post.z0dinput.option.gaz
    case 11
        zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradb ./ 1e6,'color','c','linestyle',':');
        zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremb ./ 1e6,'color','c','linestyle','-.');
end
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.pcyclo ./ 1e6,'color','k','linestyle',':');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,(pradhep + pradimpp + pradmaxp + pradwp + pradsnp + pradhdt + pradb + post.profil0d.pcyclo)./ 1e6,'color','k');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,(dpbremhe + dpbremimp + dpbremmax + dpbremw + dpbremsn + dpbrem1 + dpbremb)./ 1e6,'color','k','linestyle','-.');
%zplotprof(gca,post.profil0d.temps,post.profil0d.xli,(post.profil0d.prad + post.profil0d.pbrem + post.profil0d.pcyclo) ./ 1e6,'color','y','linestyle','-.');
switch post.z0dinput.option.gaz
    case 11
        legend('He (line + brem)','He (brem)',sprintf('Z = %d (line + brem)',zimp),sprintf('Z = %d (brem + brem)',zimp), ...
            sprintf('Z = %d (line + brem)',zmax),sprintf('Z = %d (brem)',zmax),'W (line + brem)','W (brem)','Sn (line + brem)','Sn (brem)', ...
            'HDT (line + brem)','HDT (brem)','Boron (line + brem)','Boron (brem)','Cyclo','Total (line + brem + cyclo)','Total (brem)');
    otherwise
       legend('He (line + brem)','He (brem)',sprintf('Z = %d (line + brem)',zimp),sprintf('Z = %d (brem + brem)',zimp), ...
            sprintf('Z = %d (line + brem)',zmax),sprintf('Z = %d (brem)',zmax),'W (line + brem)','W (brem)','Sn (line + brem)','Sn (brem)', ...
            'HDT (line + brem)','HDT (brem)','Cyclo','Total (line + brem + cyclo)','Total (brem)'); 
end
set(gca,'yscale','log');
xlabel('x')
ylabel('MW/m^3')
title(sprintf('METIS : %s@%d core plasma radiation', ...
          post.z0dinput.machine,post.z0dinput.shot));

z0loglin(gca);
edition2

h = findobj(0,'type','figure','tag','z0plotrad0d');
if isempty(h)
       h=figure('tag','z0plotrad0d');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen)

pradhep_tot = trapz(post.profil0d.xli,pradhep .* post.profil0d.vpr,2);
pradbremhe_tot = trapz(post.profil0d.xli,dpbremhe .* post.profil0d.vpr,2);
pradimpp_tot = trapz(post.profil0d.xli,pradimpp .* post.profil0d.vpr,2);
pradbremimp_tot = trapz(post.profil0d.xli,dpbremimp .* post.profil0d.vpr,2);
pradmaxp_tot = trapz(post.profil0d.xli,pradmaxp .* post.profil0d.vpr,2);
pradbremmax_tot = trapz(post.profil0d.xli,dpbremmax .* post.profil0d.vpr,2);
pradwp_tot = trapz(post.profil0d.xli,pradwp .* post.profil0d.vpr,2);
pradsnp_tot = trapz(post.profil0d.xli,pradsnp .* post.profil0d.vpr,2);
pradbremw_tot = trapz(post.profil0d.xli,dpbremw .* post.profil0d.vpr,2);
pradbremsn_tot = trapz(post.profil0d.xli,dpbremsn .* post.profil0d.vpr,2);
pradhdtp_tot = trapz(post.profil0d.xli,pradhdt .* post.profil0d.vpr,2);
pradbrem1_tot = trapz(post.profil0d.xli,dpbrem1 .* post.profil0d.vpr,2);
pradbp_tot = trapz(post.profil0d.xli,pradb .* post.profil0d.vpr,2);
pradbremb_tot = trapz(post.profil0d.xli,dpbremb .* post.profil0d.vpr,2);
pcyclo_tot = trapz(post.profil0d.xli,post.profil0d.pcyclo .* post.profil0d.vpr,2);
switch post.z0dinput.option.gaz
    case 11
        plot(post.profil0d.temps,pradhep_tot ./ 1e6,'g',post.profil0d.temps,pradbremhe_tot ./ 1e6,'g-.', ...
            post.profil0d.temps,pradimpp_tot ./ 1e6,'b',post.profil0d.temps,pradbremimp_tot ./ 1e6,'b-.', ...
            post.profil0d.temps,pradmaxp_tot ./ 1e6,'m',post.profil0d.temps,pradbremmax_tot ./ 1e6,'m-.', ...
            post.profil0d.temps,pradwp_tot ./ 1e6,'r',post.profil0d.temps,pradbremw_tot ./ 1e6,'r-.', ...
            post.profil0d.temps,pradsnp_tot ./ 1e6,'r--',post.profil0d.temps,pradbremsn_tot ./ 1e6,'r:', ...
            post.profil0d.temps,pradhdtp_tot ./ 1e6,'c',post.profil0d.temps,pradbrem1_tot ./ 1e6,'c-.', ...
            post.profil0d.temps,pradbp_tot ./ 1e6,'c--',post.profil0d.temps,pradbremb_tot ./ 1e6,'c:', ...
            post.profil0d.temps,pcyclo_tot ./ 1e6,'k:', ...
            post.profil0d.temps,(pradhep_tot + pradimpp_tot + pradmaxp_tot + pradwp_tot + pradhdtp_tot + pcyclo_tot) ./ 1e6,'k', ...
            post.profil0d.temps,(pradbremhe_tot + pradbremimp_tot + pradbremmax_tot + pradbremw_tot + pradbrem1_tot) ./ 1e6,'k-.');
    otherwise
        plot(post.profil0d.temps,pradhep_tot ./ 1e6,'g',post.profil0d.temps,pradbremhe_tot ./ 1e6,'g-.', ...
            post.profil0d.temps,pradimpp_tot ./ 1e6,'b',post.profil0d.temps,pradbremimp_tot ./ 1e6,'b-.', ...
            post.profil0d.temps,pradmaxp_tot ./ 1e6,'m',post.profil0d.temps,pradbremmax_tot ./ 1e6,'m-.', ...
            post.profil0d.temps,pradwp_tot ./ 1e6,'r',post.profil0d.temps,pradbremw_tot ./ 1e6,'r-.', ...
            post.profil0d.temps,pradsnp_tot ./ 1e6,'r--',post.profil0d.temps,pradbremsn_tot ./ 1e6,'r:', ...
            post.profil0d.temps,pradhdtp_tot ./ 1e6,'c',post.profil0d.temps,pradbrem1_tot ./ 1e6,'c-.', ...
            post.profil0d.temps,pcyclo_tot ./ 1e6,'k:', ...
            post.profil0d.temps,(pradhep_tot + pradimpp_tot + pradmaxp_tot + pradwp_tot + pradhdtp_tot + pcyclo_tot) ./ 1e6,'k', ...
            post.profil0d.temps,(pradbremhe_tot + pradbremimp_tot + pradbremmax_tot + pradbremw_tot + pradbrem1_tot) ./ 1e6,'k-.');
end
     %, ...
     %zs.temps,(zs.prad + zs.pbrem + zs. pcyclo) / 1e6,'.y');
ylabel('MW');
xlabel('time (s)');
switch post.z0dinput.option.gaz
    case 11
         legend('He (line + brem)','He (brem)',sprintf('Z = %d (line + brem)',zimp),sprintf('Z = %d (brem)',zimp), ...
            sprintf('Z = %d (line + brem)',zmax),sprintf('Z = %d (brem)',zmax),'W (line + brem)','W (brem)','Sn (line + brem)','Sn (brem)', ...
            'HDT (line + brem)','HDT (brem)','Boron (line + brem)','Boron (brem)','Cyclo','Total (line + brem + cyclo)','Total (brem)');
   otherwise
        legend('He (line + brem)','He (brem)',sprintf('Z = %d (line + brem)',zimp),sprintf('Z = %d (brem)',zimp), ...
            sprintf('Z = %d (line + brem)',zmax),sprintf('Z = %d (brem)',zmax),'W (line + brem)','W (brem)','Sn (line + brem)','Sn (brem)', ...
            'HDT (line + brem)','HDT (brem)','Cyclo','Total (line + brem + cyclo)','Total (brem)');
end
title(sprintf('METIS : %s@%d core plasma radiation', ...
          post.z0dinput.machine,post.z0dinput.shot));


edition2

% if 0
%     fullscreen = get(0,'ScreenSize');
%     h = findobj(0,'type','figure','tag','z0plotrad_line');
%     if isempty(h)
% 	  h=figure('tag','z0plotrad_line');
%     else
% 	  figure(h);
%     end   
%     clf
%     set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
% 	    'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen)
% 
% 
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradhep - dpbremhe)  ./ 1e6,'color','g');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradimpp - dpbremimp) ./ 1e6,'color','b');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradmaxp - dpbremmax)./ 1e6,'color','m');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradwp - dpbremw) ./ 1e6,'color','r');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradhdt -dpbrem1) ./ 1e6,'color','c');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,max(0,pradhep - dpbremhe)  ./ 1e6 + max(0,pradimpp - dpbremimp) ./ 1e6 + max(0,pradmaxp - dpbremmax) ./ 1e6  + ...
%     max(0,pradwp - dpbremw) ./ 1e6 + max(0,pradhdt -dpbrem1) ./ 1e6 ,'color','k');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.prad./ 1e6,'color','y','linestyle','-.');
%     legend('He (line)',sprintf('Z = %d (line)',zimp), ...
%     sprintf('Z = %d (line)',zmax),'W (line)', ...
%     'HDT (line + rec)','Total (line)');
%     set(gca,'yscale','log');
%     xlabel('x')
%     ylabel('MW/m^3')
%     z0loglin(gca);
%     edition2
% 
%     fullscreen = get(0,'ScreenSize');
%     h = findobj(0,'type','figure','tag','z0plotrad_brem');
%     if isempty(h)
% 	  h=figure('tag','z0plotrad_brem');
%     else
% 	  figure(h);
%     end   
%     clf
%     set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
% 	    'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen)
% 
% 
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremhe  ./ 1e6,'color','g');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremimp ./ 1e6,'color','b');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremmax./ 1e6,'color','m');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbremw ./ 1e6,'color','r');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dpbrem1 ./ 1e6,'color','c');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,(dpbremhe + dpbremimp + dpbremmax + dpbremw + dpbrem1)./ 1e6,'color','k');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.pbrem./ 1e6,'color','y','linestyle','-.');
%     legend('He (brem)',sprintf('Z = %d (brem)',zimp), ...
%     sprintf('Z = %d (brem)',zmax),'W (brem)', ...
%     'HDT (brem)','Total (brem)');
%     set(gca,'yscale','log');
%     xlabel('x')
%     ylabel('MW/m^3')
%     z0loglin(gca);
%     edition2
% 
% 
%     % verification prad_w
%     indw =  find(uz == 74,1);
%     lzw    = postcoef(indw).lz .* 0.1;
%     tzw    = postcoef(indw).te;
%     lzw    = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(te_max)-sqrt(tzw(end))) .* 1e-28);
%     tzw    = cat(2,1e-4,tzw,te_max);
%     %rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
%     if post.z0dinput.option.noncoronal == 1
%         rzw  = z0noncronal(tzw,lzw,post.profil0d.xli,tep,nep,taue);
%     else
%         rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
%     end
%     prad_w = nep .* nwp .* rzw;
%     figure(23);
%     hold on
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,prad_w./ 1e6,'color','r');
%     zplotprof(gca,post.profil0d.temps,post.profil0d.xli,pradwp./ 1e6,'color','b');
% 
% end
