% calcul le rayonnement 0d
% ane = piquage de ne
% ate = piquage de te
% rw  = coefficient de reflexion des parois
% Vp = volume du plasma
% S  = surface externe du plasma
% nboron is the volume averaged boron density for fuel pB11. It should be
% set to zero for other mixture. It is a time dependent vector [nbt,1].
function [prad,pbrem,pcyclo,pradsol,profli] = zrad0(ne,te,zeff,zmax,zimp,rimp,zgaz,R,a,K,Bt,ane,ate,rw,Vp,S,modeh,frad, ...
                                                    te0,nbar,norme,telim,nesol_ave,taue,dsol,ya0_w,z_prad_option,sol_rad, ...
                                                    gaunt,noncoronal,temps,profli,Sn_fraction,te_max,cor_rel_brem,nboron)

% adaptation pour le vrai calcul du rayonnement non coronal
taup      = imag(taue);
taue      = real(taue);
taup(taup == 0) = taue(taup == 0);
% calcul du temps local de residence
if isfield(profli,'nep') && isfield(profli,'tep')
    pen = profli.nep .* profli.tep;
    pen = pen ./ (max(eps,max(pen,[],2)) * ones(size(profli.xli)));
    tau_resid = (taue * ones(size(profli.xli))) .* pen + (min(taue,taup) * ones(size(profli.xli))) .* (1 - pen);
else
  tau_resid = taue;
end
%figure(21);plot(profli.xli,tau_resid);drawnow

% limit for extrapolation (in keV)
te_max = max(105,te_max / 1e3 + 5);


if isfield(profli,'n1p')
	% calcul de la charge a mettre dans le scaling ; formule de Stangeby 
	% ref : The Plasma Boundary of Magnetic Fusion Devices,Stangeby, p 627, IOP 2000
	z_stangeby = min(114,(profli.n1p  ./ profli.nep - profli.zeff) ./ min(-0.001,profli.n1p  ./ profli.nep - 1));
        z_stangeby = trapz(profli.xli,profli.nep .* z_stangeby .* profli.vpr,2) ./ trapz(profli.xli,profli.nep  .* profli.vpr,2);
%  figure(21);
%  plot(temps,z_stangeby);
%  drawnow
else
	z_stangeby = zmax;
end

switch  z_prad_option
case 'zimp'
	z_prad = zimp;
case 'Stangeby'
	z_prad = z_stangeby;
otherwise
	z_prad = zmax;
end
improved = 0;
if norme < 0
  norme = 0;
  improved = 1;
  % G-F Matthews et al , Journal of Nuclear Materials, 1997 revue 1999
  % G. F. Matthews Nuclear Fusion vol 39  1999, p 27
  pradmat = real(frad .*  1e6 ./ 4.5  .* (zeff - zgaz) .* (nbar./1e20) .^ 1.89 .* S .^ 0.94 ./ z_prad .^ 0.12);
  pradsol = zeros(size(ne));
elseif norme == 2
  norme = 1;
  % Highly radiating type-III ELMy H-mode with low plasma core pollution
  % J. Rapp et al
  % Journal of Nuclear Materials 390???391 (2009) 238???241
  pradsol = zeros(size(ne));
  pradmat = real(frad .*  1e6 ./ 40  .* (zeff - zgaz) .* (nbar./1e20) .^ 1.5 .* S .^ 0.94 ./ z_prad .^ 0.12 .* a .* R ./ max(1e-6,taue));
else
  % G-F Matthews et al , Journal of Nuclear Materials, 1997 revue 1999
  % G. F. Matthews Nuclear Fusion vol 39  1999, p 27
  pradmat = real(frad .*  1e6 ./ 4.5  .* (zeff - zgaz) .* (nbar./1e20) .^ 1.89 .* S .^ 0.94 ./ z_prad .^ 0.12);
  pradsol = zeros(size(ne));
end

if isfield(profli,'nzp')
    [ua,uz,post]=z0coefpost;
    
    indhe  =  find(ua==4 & uz == 2);
    dz     = abs(uz - zimp);
    indimp =  min(find(dz == min(dz)));
    dz     = abs(uz - zmax);
    indmax =  min(find(dz == min(dz)));
    dz     = abs(uz - 5);
    indb =  min(find(dz == min(dz)));
      
    tep    = profli.tep ./ 1e3;
    nep    = profli.nep ./ 1e6;
    nhep   = profli.nhep ./ 1e6;
    n1p    = profli.n1p ./ 1e6;
    nzp    = profli.nzp ./ 1e6;
    nwp    = profli.nwp ./ 1e6;
    
    % for boron
    nbp    = profli.nep .* ((nboron ./ max(1e13,ne)) * ones(1,size(profli.nep,2))) / 1e6;
    
    lzhe   = post(indhe).lz .* 0.1;
    tzhe   = post(indhe).te;
    lzhe   = cat(2,1e-38,lzhe,lzhe(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzhe(end))) .* 1e-28);
    tzhe   = cat(2,1e-4,tzhe,te_max);
    
    lzimp    = post(indimp).lz .* 0.1;
    tzimp    = post(indimp).te;
    lzimp    = cat(2,1e-38,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-28);
    tzimp    = cat(2,1e-4,tzimp,te_max);
    
    lzmax    = post(indmax).lz .* 0.1;
    tzmax    = post(indmax).te;
    lzmax    = cat(2,1e-38,lzmax,lzmax(end) + 5.355e3 .* uz(indmax)  .*(sqrt(te_max)-sqrt(tzmax(end))) .* 1e-28);
    tzmax    = cat(2,1e-4,tzmax,te_max);
    
    lzb   = post(indb).lz .* 0.1;
    tzb   = post(indb).te;
    lzb   = cat(2,1e-38,lzb,lzb(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzb(end))) .* 1e-28);
    tzb   = cat(2,1e-4,tzb,te_max);
 
    if noncoronal == -1
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
            fprintf('b');  %noncoronal radiative data are missing for boron
        end
        rzb    = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tep(:))),size(tep));
        
    elseif noncoronal == 1
        rzhe   = z0noncronal(tzhe,lzhe,profli.xli,tep,nep,taue);
        rzimp  = z0noncronal(tzimp,lzimp,profli.xli,tep,nep,taue);
        rzmax  = z0noncronal(tzmax,lzmax,profli.xli,tep,nep,taue);
        rzb   = z0noncronal(tzb,lzb,profli.xli,tep,nep,taue);
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
    
    % cas du tungstene traite part
    indw =  find(uz == 74,1);
    lzw    = post(indw).lz .* 0.1;
    tzw    = post(indw).te;
    lzw    = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(te_max)-sqrt(tzw(end))) .* 1e-28);
    tzw    = cat(2,1e-4,tzw,te_max);
    if noncoronal == 1
        rzw  = z0noncronal(tzw,lzw,profli.xli,tep,nep,taue);
    else
        rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tep(:))),size(tep));
    end
    pradwp   = nep .* nwp .* rzw;
    if Sn_fraction > 0
        % Sn case
        indsn =  find(uz == 50,1);
        lzsn    = post(indsn).lz .* 0.1;
        tzsn    = post(indsn).te;
        lzsn    = cat(2,1e-38,lzsn,lzsn(end) + 5.355e3 .* 50  .*(sqrt(te_max)-sqrt(tzsn(end))) .* 1e-28);
        tzsn    = cat(2,1e-4,tzsn,te_max);
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
    rrec    = reshape(exp(pchip(log(trec),log(srec),log(profli.tep(:)))),size(tep));
    %pradhdt_old = 5.355e3 .* (profli.nep./ 1e20) .* ( profli.n1p ./ 1e20) .* sqrt(profli.tep ./ 1e3) + ...
    %          profli.nep .* profli.n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3
    
    % amelioration du calcul du rayonnement de l'hydorgene
    indhdt  = min(find(uz == 1));
    lzhdt   = post(indhdt).lz .* 0.1;
    tzhdt   = post(indhdt).te;
    lzhdt   = cat(2,1e-38,lzhdt,lzhdt(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzhdt(end))) .* 1e-28);
    tzhdt   = cat(2,1e-4,tzhdt,te_max);
    if noncoronal == 1
        rzhdt   = z0noncronal(tzhdt,lzhdt,profli.xli,tep,nep,taue);
    else
        rzhdt    = reshape(10 .^ pchip(log10(tzhdt),log10(lzhdt),log10(tep(:))),size(tep));
    end
    pradhdt = nep .* n1p .* rzhdt + profli.nep .* profli.n1p .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3
    
    %figure(21);clf;plot(profli.xli,mean(pradhdt_old,1),'b',profli.xli,mean(pradhdt,1),'r');drawnow
    
    switch improved
        case 1
            [pradhep,pradimpp,pradmaxp,pradwp,pradhdt,pradb] = zrad0improved(profli,zimp,zmax,rimp,noncoronal,taue + sqrt(-1) .* taup,Sn_fraction,te_max,nboron);
    end
    
    profli.fprad  = real(frad .* (pradhep + pradimpp + pradmaxp + pradhdt + pradwp + pradb));
    prad          = trapz(profli.xli,profli.fprad .* profli.vpr,2);
    
    if norme == 1
        fact = min(100,max(0.01,max(1,pradmat) ./ max(1,prad)));
        prad = max(1,fact .* prad);
        profli.fprad  = profli.fprad .* (fact * ones(size(profli.xli)));
        pradhep       = pradhep .* (fact * ones(size(profli.xli)));
        pradimpp      = pradimpp .* (fact * ones(size(profli.xli)));
        pradmaxp      = pradmaxp .* (fact * ones(size(profli.xli)));
        pradhdt       = pradhdt .* (fact * ones(size(profli.xli)));
        pradwp        = pradwp .* (fact * ones(size(profli.xli)));
        pradb         = pradb .* (fact * ones(size(profli.xli)));
    else
        % le rayonnement dans le divertor est donne par la difference
        switch sol_rad
            case 'coupled'
                pradsol = max(0,pradmat - prad);
        end
    end
else
    % G-F Matthews et al , Journal of Nuclear Materials, 1997 revue 1999
    % G. F. Matthews Nuclear Fusion vol 39  1999, p
    prad    = frad .* pradmat
    switch sol_rad
        case 'coupled'
            pradsol = max( 2./ 3 .* pradmat, (1-frad) .* pradmat);
    end
end
% securite 1er temps
prad(1) = 0;

% references :
%
%  *  Eletron cyclotron radiative transfer in fusion plasma, F. Albajar et all, Nucl. Fus.,
%     vol 42, 2002, pp 670-678.
%
%  *  Improved calculation of synchrotron radiation losses in relistic tokamak plamas, F. Albajar et al,
%     Nucl. Fus.,vol 41, 2001, pp 665-678.
%
%
%  * LATF mdoel : 
%    F. Albajar et al, NF 49, 2009 p 115017
%
% facteur de forme des profils
%
if imag(rw) ~= 0
   % le calcul n'est pas requis
   pcyclo = NaN .* te0;
elseif (rw < 0) && isfield(profli,'tep') &&  isfield(profli,'nep')
   % LATF model : 
   ctrub = 8.2e-6.* sqrt(1-abs(rw));
   mu    = 9.10938188e-31 .* 2.99792458e8 .^ 2 ./ (profli.tep .* 1.602176462e-19);
   m2t   = (1 + 1.93030 ./ mu) ./ (1- 0.58167 ./ mu);
   dpdv  = ctrub ./ (2 .* pi) .^ 2  ./  sqrt(profli.rmx(:,end) * ones(size(profli.xli))) .* sqrt(profli.nep ./ 1e20) .*  ...
           (profli.tep ./ 1e3 .* sqrt(profli.bpol .^ 2 +  (profli.fdia  ./ profli.Raxe) .^ 2)) .^ (5/2) .* m2t;
   Am1   = profli.epsi(:,end);
   ga    = sqrt(1 + 9 .* Am1 ./ sqrt(te ./ 1e3));

   pcyclo = ga .* trapz(profli.xli,dpdv .* profli.vpr,2) .* 1e6;   

else
  if isfield(profli,'tep')
	  [alphaf,betaf,teout] = z0expote(profli.xli(1:end-1),profli.tep(:,1:end-1));
	  te0                  = profli.tep(:,1) ./ 1e3;
	  ne0                  = profli.nep(:,1) ./ 1e20;
	  %figure(61);clf
	  %plot(profli.xli,profli.tep,'b',profli.xli,teout,'or');
	  %pause
  else
	  betaf  = 2; % exposant de r/a dans Te(r/a)
	  alphaf = ate;
	  te0     = te0 ./ 1e3;
	  ne0   = ne .* (1 + ane) ./ 1e20;

  end
  if isfield(profli,'nep')
	% pour utiliser le piquage vrai et non pas le piquage demande
	gamma    = profli.nep(:,1) ./ trapz(profli.xli,profli.vpr .* profli.nep,2) .* trapz(profli.xli,profli.vpr,2) - 1;
	%figure(21);
	%clf
	%plot(temps,ane,'.r',temps,gamma);
	%drawnow
  else
  	gamma  = ane;
  end
  ka = (gamma + 3.87 .* alphaf + 1.46 ) .^ (-0.79) .* ...
	  (1.98 + alphaf) .^ 1.36 .* betaf .^ 2.14 .* ...
	  (betaf .^ 1.53 + 1.87 .* alphaf - 0.16) .^ (-1.33);

  % facteur de rapport d'aspect
  rap = max(1.5,min(15,R ./ a));
  ga  = 0.93 .* (1 + 0.85 .* exp( - 0.82 .* rap));


  % pa0
  pa0   = 6.04e3 .* a .* ne0 ./ Bt;

  % puissance totale par rayon
  pcyclo  = real(1e6 .* 3.84e-8 .* (1-abs(rw)) .^ 0.62 .* R .* a .^ 1.38 .* K .^ 0.79 .* ...
	    Bt .^ 2.62 .* ne0 .^ 0.38 .* te0 .* (16 + te0) .^ 2.61 .* ...
	    (1 + 0.12 .* te0 .* ((1-max(0,rw)) ./ pa0) .^ 0.41) .^ (-1.51) .* ...
	    ka .* ga);

end

% Pbrem  (rapport dimensionnement)
if (gaunt == 1) && isfield(profli,'tep')
    [gg1,Cbrem] = z0gaunt_brem(profli.tep ./ 1e3,1);
    gg2         = z0gaunt_brem(profli.tep ./ 1e3,2);
    gg5         = z0gaunt_brem(profli.tep ./ 1e3,5);
    gg_zimp     = z0gaunt_brem(profli.tep ./ 1e3,zimp);
    gg_zmax     = z0gaunt_brem(profli.tep ./ 1e3,zmax);
    gg_w        = z0gaunt_brem(profli.tep ./ 1e3,z0wavez(profli.tep));
    gg_sn        = z0gaunt_brem(profli.tep ./ 1e3,z0snavez(profli.tep));
    Cbrem       = Cbrem .* 1e40;
    
    dpbrem1  = Cbrem .* gg1 .* (profli.nep ./ 1e20) .* (profli.n1p ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremhe = Cbrem .* gg2 .* 4 .* (profli.nep ./ 1e20) .* (profli.nhep ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremmax  = rimp .* Cbrem .* gg_zmax .* zmax .^ 2 .* (profli.nep ./ 1e20) .* (profli.nzp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremimp  = Cbrem .* gg_zimp .* zimp .^ 2 .* (profli.nep ./ 1e20) .* (profli.nzp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremb    = Cbrem .* gg5 .* 25 .* (profli.nep ./ 1e20) .*  ...
                 (profli.nep .* ((nboron ./ max(1e13,ne)) * ones(1,size(profli.nep,2))) ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    
    % cas du tungstene
    dpbremw  = Cbrem .* gg_w .* z0wavez(profli.tep) .^ 2 .* (profli.nep ./ 1e20) .* (profli.nwp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    if Sn_fraction > 0
        dpbremw = (1- Sn_fraction) .* dpbremw + ...
                  Sn_fraction .* Cbrem .* gg_sn .* z0snavez(profli.tep) .^ 2 .* (profli.nep ./ 1e20) .* (profli.nwp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    end
    
    if isfield(profli,'fprad')
        % on spepare les contribution de raies et le brem
        % fpradmem = profli.fprad;
        profli.fprad  = real(frad .* (max(0,pradhep - dpbremhe) + max(0,pradmaxp - dpbremmax) + ...
            max(0,pradimpp - dpbremimp) + max(0,pradhdt - dpbrem1) + max(0,pradwp - dpbremw) + max(0,pradb - dpbremb)));
        
        prad   = trapz(profli.xli,profli.fprad .* profli.vpr,2);
    end
    dpbrem   = dpbrem1 + dpbremhe + dpbremmax + dpbremimp + dpbremw + pradb;
    %pbremnc  = trapz(profli.xli,profli.vpr .*dpbrem,2);
    % correction relativiste de Pbrem
    switch cor_rel_brem
        case 'improved'
            xrel = brem_rel_cor(profli.tep,profli.zeff);
            
        otherwise
            % P.E. Stott, Plasma. Physics and Controlled Fusion 47 (2005) 1305-1338 ref [15]
            tec2   = 511e3;% me c^2
            xrel   = (1 + 2 .* profli.tep ./ tec2) .* (1 + (2 ./ profli.zeff) .* (1 - (1 + profli.tep ./ tec2) .^ (-1)));
    end
    dpbrem = real(dpbrem .* xrel);
    profli.pbrem = dpbrem ;
    % integration
    pbrem   = trapz(profli.xli,profli.vpr .*dpbrem,2);
    if ~isfield(profli,'fprad')
        % prad contient deja pbrem
        % attention la separation brem/rad est problematique
        prad     = max(1.01 .* pbrem,prad);
        prad     = prad - pbrem;
    end
    
elseif isfield(profli,'tep')
    dpbrem1  = 5.355e3 .* (profli.nep ./ 1e20) .* (profli.n1p ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremhe = 5.355e3 .* 4 .* (profli.nep ./ 1e20) .* (profli.nhep ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremmax  = rimp .* 5.355e3 .* zmax .^ 2 .* (profli.nep ./ 1e20) .* (profli.nzp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremimp  = 5.355e3 .* zimp .^ 2 .* (profli.nep ./ 1e20) .* (profli.nzp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    dpbremb  = 5.355e3 .* 5 .^ 2 .* (profli.nep ./ 1e20) .* (profli.nep .* ((nboron ./ max(1e13,ne)) * ones(1,size(profli.nep,2))) ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    
    % cas du tungstene
    dpbremw  = 5.355e3 .* z0wavez(profli.tep) .^ 2 .* (profli.nep ./ 1e20) .* (profli.nwp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    if Sn_fraction > 0
        dpbremw = (1- Sn_fraction) .* dpbremw + ...
            Sn_fraction .* 5.355e3 .* z0snavez(profli.tep) .^ 2 .* (profli.nep ./ 1e20) .* (profli.nwp ./ 1e20) .* sqrt(profli.tep ./ 1e3);
    end
    
    if isfield(profli,'fprad')
        % on spepare les contribution de raies et le brem
        %fpradmem = profli.fprad;
        profli.fprad  = real(frad .* (max(0,pradhep - dpbremhe) + max(0,pradmaxp - dpbremmax) + ...
            max(0,pradimpp - dpbremimp) + max(0,pradhdt - dpbrem1) + max(0,pradwp - dpbremw) + max(0,pradb - dpbremb)));
        
        prad   = trapz(profli.xli,profli.fprad .* profli.vpr,2);
    end
    dpbrem   = dpbrem1 + dpbremhe + dpbremmax + dpbremimp + dpbremw + dpbremb;
    %pbremnc  = trapz(profli.xli,profli.vpr .*dpbrem,2);
    % correction relativiste de Pbrem
    switch cor_rel_brem
        case 'improved'
            xrel = brem_rel_cor(profli.tep,profli.zeff);
            
        otherwise
            % P.E. Stott, Plasma. Physics and Controlled Fusion 47 (2005) 1305-1338 ref [15]
            tec2   = 511e3;% me c^2
            xrel   = (1 + 2 .* profli.tep ./ tec2) .* (1 + (2 ./ profli.zeff) .* (1 - (1 + profli.tep ./ tec2) .^ (-1)));
    end
    dpbrem = real(dpbrem .* xrel);
    profli.pbrem = dpbrem ;
    % integration
    pbrem   = trapz(profli.xli,profli.vpr .*dpbrem,2);
    if ~isfield(profli,'fprad')
        % prad contient deja pbrem
        % attention la separation brem/rad est problematique
        prad     = max(1.01 .* pbrem,prad);
        prad     = prad - pbrem;
    end
else
    gbr      = (1 + ane) .^ 2 .* sqrt(1 + ate) ./ (1 + 2 .* ane + 0.5 .* ate);
    pbrem    = real(1e6 .* 1.68e-2 .* zeff .* (ne./1e20) .^ 2 .* sqrt(te./1e4) .* Vp .* gbr);
    %pbremnc  =  pbrem;
    % prad contient deja pbrem
    % attention la separation brem/rad est problematique
    prad     = max(1.01 .* pbrem,prad);
    prad     = real(prad - pbrem);
end

% calcul du rayonnement de la sol
% continuite sur la sepratrice
%  lambda   = mean(a./100);
%  rsol    = linspace(0,lambda .* 5,21);
%  vr      = ones(size(rsol));
%  fsol    = exp(-rsol ./ lambda);
%  vprsol  = (profli.vpr_tor(:,end) * vr)  .* ( 1 + (1 ./ profli.rmx(:,end)) * rsol);
%  psol    =  profli.fprad(:,end) * fsol;
%  pradsol = max(trapz(rsol,psol .* vprsol,2),max(0,pradmat - prad - pbrem));
% reference  : Post  et al, PoP vol 2 june 1995, p 2328-
if isfield(profli,'nzp')
    [ua,uz,post]=z0coefpost;
    
    indhe  =  find(ua==4 & uz == 2);
    dz     = abs(uz - zimp);
    indimp =  min(find(dz == min(dz)));
    dz     = abs(uz - zmax);
    indmax =  min(find(dz == min(dz)));
    dz     = abs(uz - 5);
    indb =  min(find(dz == min(dz)));

    %lambda   = mean(a./50);
    lambda   = mean(2 .* dsol);
    rsol    = linspace(0,lambda .* 5,21);
    vr      = ones(size(rsol));
    fsol    = exp(-rsol ./ lambda);
    
    vprsol  = (profli.vpr_tor(:,end) * vr)  .* ( 1 + (1 ./ profli.rmx(:,end)) * rsol);
    
    ltoln = 1.4;
    % temperature dans la sol
    teform = fsol .^ (1./ ltoln);
    
    % densite moyenne le long de la ligne de champ
    fact_ave = nesol_ave ./ profli.nep(:,end);
    switch sol_rad
        case 'coupled'
            tesol    = max(telim * ones(size(teform)),profli.tep(:,end)  * teform) ./ 1e3;
        otherwise
            tesol    = max(0,profli.tep(:,end)  * teform) ./ 1e3;
    end
    % eviter log 0
    tesol = max(eps,tesol);
    
    nesol    = (fact_ave .* profli.nep(:,end)) * fsol ./ 1e6;
    n1sol    = (fact_ave .* profli.n1p(:,end))  * fsol ./ 1e6;
    nhesol   = (fact_ave .* profli.nhep(:,end)) * fsol ./ 1e6;
    nzsol    = (fact_ave .* profli.nzp(:,end)) * fsol ./ 1e6;
    % for boron
    nbsol    = fact_ave .* profli.nep(:,end) .* (nboron ./ max(1e13,ne))  * fsol / 1e6;
    
    
    lzhe   = post(indhe).lz .* 0.1;
    tzhe   = post(indhe).te;
    lzhe   = cat(2,1e-38,lzhe,lzhe(end) + 5.355e3 .* 4  .* (sqrt(te_max)-sqrt(tzhe(end))) .* 1e-28);
    tzhe   = cat(2,1e-4,tzhe,te_max);
    
    lzimp    = post(indimp).lz .* 0.1;
    tzimp    = post(indimp).te;
    lzimp    = cat(2,1e-38,lzimp,lzimp(end) + 5.355e3 .* uz(indimp)  .*(sqrt(te_max)-sqrt(tzimp(end))) .* 1e-28);
    tzimp    = cat(2,1e-4,tzimp,te_max);
    
    lzmax    = post(indmax).lz .* 0.1;
    tzmax    = post(indmax).te;
    lzmax    = cat(2,1e-38,lzmax,lzmax(end) + 5.355e3 .* uz(indmax)  .*(sqrt(te_max)-sqrt(tzmax(end))) .* 1e-28);
    tzmax    = cat(2,1e-4,tzmax,te_max);
    
    lzb   = post(indb).lz .* 0.1;
    tzb   = post(indb).te;
    lzb   = cat(2,1e-38,lzb,lzb(end) + 5.355e3 .* 4  .*(sqrt(te_max)-sqrt(tzb(end))) .* 1e-28);
    tzb   = cat(2,1e-4,tzb,te_max);
 

    if noncoronal == 1
        rzhe     = z0noncronal(tzhe,lzhe,rsol,tesol,nesol,taue,dsol,a);
        rzimp    = z0noncronal(tzimp,lzimp,rsol,tesol,nesol,taue,dsol,a);
        rzmax    = z0noncronal(tzmax,lzmax,rsol,tesol,nesol,taue,dsol,a);
        if any(nbp(:) > sqrt(eps))
            fprintf('b');  %noncoronal radiative data are missing for boron
        end
        rzb    = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tesol(:))),size(tesol));
    else
        rzhe    = reshape(10 .^ pchip(log10(tzhe),log10(lzhe),log10(tesol(:))),size(tesol));
        rzimp   = reshape(10 .^ pchip(log10(tzimp),log10(lzimp),log10(tesol(:))),size(tesol));
        rzmax   = reshape(10 .^ pchip(log10(tzmax),log10(lzmax),log10(tesol(:))),size(tesol));
        rzb    = reshape(10 .^ pchip(log10(tzb),log10(lzb),log10(tesol(:))),size(tesol));
    end
    % dans ce cas on calcul que le rayonnement dans le plasma de coeur
    pradhesol   = nesol .* nhesol .* rzhe;
    pradimpsol  = nesol .* nzsol  .* rzimp;
    pradmaxsol  = nesol .* nzsol  .* rzmax .* rimp;
    pradbsol    = nesol .* nbsol  .* rzb;
    
    
    % cas du tungstene traite a part
    % nwsol  = (ya0_w * vr) .* nesol;
    % le calcul du rayonnement dans le divertor est fait a part maintenant
    nwsol    = (fact_ave .* profli.nwp(:,end)) * fsol ./ 1e6;
    indw   = find(uz == 74,1);
    lzw    = post(indw).lz .* 0.1;
    tzw    = post(indw).te;
    lzw    = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(te_max)-sqrt(tzw(end))) .* 1e-28);
    tzw    = cat(2,1e-4,tzw,te_max);
    if noncoronal == 1
        rzw     = z0noncronal(tzw,lzw,rsol,tesol,nesol,taue,dsol,a);
    else
        rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(tesol(:))),size(tesol));
    end
    pradwsol  = nesol .* nwsol  .* rzw;
    %
    if Sn_fraction > 0
        % Sn case
        indsn =  find(uz == 50,1);
        lzsn    = post(indsn).lz .* 0.1;
        tzsn    = post(indsn).te;
        lzsn    = cat(2,1e-38,lzsn,lzsn(end) + 5.355e3 .* 50  .*(sqrt(te_max)-sqrt(tzsn(end))) .* 1e-28);
        tzsn    = cat(2,1e-4,tzsn,te_max);
        if noncoronal == 1
            rzsn  = z0noncronal(tzsn,lzsn,rsol,tesol,nesol,taue,dsol,a);
        else
            rzsn  = reshape(10 .^ pchip(log10(tzsn),log10(lzsn),log10(tesol(:))),size(tesol));
        end
        pradwsol  = (1 - Sn_fraction) .* pradwsol + Sn_fraction .* nesol .* nwsol  .* rzsn;
    end
    
    
    
    % rayonnement de HDT
    % ajout de la recombinaison
    % recombinaison (Halpha + ...)
    trec    = [0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
    srec    = [7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 1e-6;% m^3/s;
    rrec    = reshape(exp(pchip(log(trec),log(srec),log(tesol(:).*1e3))),size(tesol));
    if noncoronal == 1
        pradhdtsol = 5.355e3 .* (nesol./ 1e20 .* 1e6) .* ( n1sol ./ 1e20 .* 1e6) .* sqrt(tesol) + ...
            nesol .* 1e6 .* n1sol .* 1e6 .* rrec .* (13.6 + (3/2) .* (1 + ((profli.tip(:,end) ./ profli.tep(:,end)) * ones(size(teform))) .* tesol)) .* 1.6022e-19; % passage en W/m^3
    else
        pradhdtsol = 5.355e3 .* (nesol./ 1e20 .* 1e6) .* ( n1sol ./ 1e20 .* 1e6) .* sqrt(tesol) + ...
            nesol .* 1e6 .* n1sol .* 1e6 .* rrec .* 13.6 .* 1.6022e-19; % passage en W/m^3
    end
    
    
    switch sol_rad
        case 'coupled'
            fpradsol      = real(frad .* (pradhesol + pradimpsol + pradmaxsol + pradhdtsol + pradwsol + pradbsol));
            pradsol       = pradsol + dsol .* trapz(rsol ./ lambda,fpradsol .* vprsol,2);
        otherwise
            fpradsol      = real(pradhesol + pradimpsol + pradmaxsol + pradhdtsol + pradwsol + pradbsol);
            pradsol       = dsol .* trapz(rsol ./ lambda,fpradsol .* vprsol,2);
    end
    
end


% source externe pour PRAD
if isappdata(0,'PLINE_EXP') & isfield(profli,'qjli');
	pline_exp  = getappdata(0,'PLINE_EXP');
	pline      = max(1,interp1_ex(pline_exp.temps,pline_exp.prad,temps,'nearest','extrap'));
	pline      = pchip(pline_exp.x,pline,profli.xli);

	% recopie
	profli.fprad  = real(pline);
end

% renormalisation finale pour amelioree la precision
prad_in       = trapz(profli.xli,profli.fprad .* profli.vpr,2);
profli.fprad  = profli.fprad .* ((prad ./max(eps,prad_in)) * ones(1,length(profli.xli)));
