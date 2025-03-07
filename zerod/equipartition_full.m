function [qeib0,qei] = equipartition_full(te,ti,ne,nH,nD,nT,nHe,nimp1,nimp2,nW,nSn,Zimp1,Zimp2,nbp,nHe3)


% list of atomic number
% search for impurities
[A,Z] = chargemasse;
ind1 = find(Z == Zimp1,1);
if isempty(ind1)
    Aimp1 = 7/3 * Zimp1;
else
    Aimp1 = A(ind1);
end
ind2 = find(Z == Zimp2,1);
if isempty(ind2)
    Aimp2 = 7/3 * Zimp2;
else
    Aimp2 = A(ind2);
end
ZW  = z0wavez(te);
ZSn = z0snavez(te);


% Coulomb logarithm
warning off
lnei          =  15.2 - 0.5 .* log(ne ./ 1e20) + log(te ./ 1e3);
warning on
ind = find(~isfinite(lnei) | (lnei <10));
if ~isempty(ind)
    lnei(ind) = 10 .* ones(1,length(ind));
end

qeib0_sp = rel_ex(te,ti,1) .* spitzer_1(te,ti,ne,nH,1,1,lnei) + rel_ex(te,ti,2) .* spitzer_1(te,ti,ne,nD,1,2,lnei) + ...
           rel_ex(te,ti,3) .* spitzer_1(te,ti,ne,nT,1,3,lnei) + rel_ex(te,ti,4) .* spitzer_1(te,ti,ne,nHe,2,4,lnei) + ...
           rel_ex(te,ti,Aimp1) .* spitzer_1(te,ti,ne,nimp1,Zimp1,Aimp1,lnei) + rel_ex(te,ti,Aimp2) .* spitzer_1(te,ti,ne,nimp2,Zimp2,Aimp2,lnei) + ...
           rel_ex(te,ti,183.84) .* spitzer_1(te,ti,ne,nW,ZW,183.84,lnei) + rel_ex(te,ti,118.71) .* spitzer_1(te,ti,ne,nSn,ZSn,118.71,lnei) + ...
           rel_ex(te,ti,11) .* spitzer_1(te,ti,ne,nbp,5,11,lnei) + rel_ex(te,ti,3.02) .* spitzer_1(te,ti,ne,nHe3,2,3.02,lnei);
    
% relativistic correction on collision time 
qeib0 = rel_cor(te) .* qeib0_sp;         

% correction for large Ti/Te
%from:
% Modification of classical Spitzer ionelectron energy transfer rate
% for large ratios of ion to electron temperatures, T.H. Rider and P.J. Catto,
% Phys. Plasmas 2, 18731885 (1995), https://doi.org/10.1063/1.871274
% formula 56
% la somme doit etre entre 1 et 50 (Rider Phd) à verifier !
factor_tiote = term(te,ti,ne,nH,1,1) + term(te,ti,ne,nD,1,2) +  ...
               term(te,ti,ne,nT,1,3) + term(te,ti,ne,nHe,2,4) +  ....
               term(te,ti,ne,nimp1,Zimp1,Aimp1) + term(te,ti,ne,nimp2,Zimp2,Aimp2) +  ....
               term(te,ti,ne,nW,ZW,183.84) + term(te,ti,ne,nSn,ZSn,118.71) + ...
               term(te,ti,ne,nbp,5,11) + term(te,ti,ne,nbp,2,3.02);
           
% should be between less than 50
factor_tiote = min(50,factor_tiote);
qeib0 = exp(-(3 .* sqrt(pi) ./ 4 .* factor_tiote)  .^ (2/3)) .* qeib0;
if nargout > 1
    qei = phys.e .* (te - ti) .* qeib0;
end


% for testingr
return
% physical constantes
phys = cphys;


taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
    ((phys.e .* te) .^ (3/2) ./ ne ./ lnei);
factor_w_sn = nW .* ZW  .^ 2 ./ 183.84  +  nSn .* ZSn .^ 2 ./ 118.71 ;
qeib0_old    = 3 .* phys.me ./ phys.mp ./ taues .* (nH + nD ./ 2 + nT ./ 3 + nHe  +  ...
        (Zimp1 ./ (7/3)) .* nimp1 + (Zimp2./ (7/3)) .* nimp2 + factor_w_sn);
qei_old     = qeib0_old .* phys.e .* (te - ti);

figure(21);
time=(1:size(te,1))';
xx = linspace(0,1,size(te,2));
subplot(2,2,1)
zplotprof(gca,time,xx,qeib0_sp,'color','r');
hold on
zplotprof(gca,time,xx,qeib0_old,'color','b');
subplot(2,2,2);
qei = phys.e .* (te - ti) .* qeib0;
zplotprof(gca,time,xx,qei,'color','r');
hold on
zplotprof(gca,time,xx,qei_old,'color','b');
subplot(2,2,3)
zplotprof(gca,time,xx,exp(-(3 .* sqrt(pi) ./ 4 .* factor_tiote)  .^ (2/3)));
subplot(2,2,4)
zplotprof(gca,time,xx,rel_cor(te),'color','b');
hold on 
zplotprof(gca,time,xx,rel_ex(te,ti,1),'color','r');

drawnow




function rc = rel_cor(te)
% relativistic corection
% valid at all energy
% from Thesis:  Susan Stepney, Relativistic thermal plasma, Institute of
% Astronomy, University of Cambridge, September 1983
% formula 2.19 
% this is the Klein-Nishina corrections decrease scattering cross section
% term

% physical constantes
phys = cphys;

tstar = phys.e .* te ./ phys.me ./ phys.c .^ 2; 
%
%epsi_low = 3/2 + 15/8 .* tstar - 15/8 .* tstar .^ 2;
%epsi_high= 3 - 1./ tstar + 1 ./ 2 ./ tstar .^ 2;
epsi  = (besselk(1, 1./ tstar) ./ besselk(2, 1 ./ tstar) + 3 .* tstar - 1) ./ tstar;

% normalisation to be used with qei
% il y a sans doute une inversion
% le temps d'equipartition augmente avec te relativist dans le papier 
% donc la puissance echangée diminue avec ce facteur.
rc = (3./2) ./ epsi;
rc(~isfinite(rc)) = 1;


function ec = rel_ex(te,ti,Ai)
% relativistic correction to energy exchange compare to Spitzer formula 
% valid at all energy
% from Thesis:  Susan Stepney, Relativistic thermal plasma, Institute of
% Astronomy, University of Cambridge, September 1983
% formulas 2.23 and 2.34
% this is equivalent to the Cordey correction in previous cited work but
% without approximation
% in not take into account the Klein-Nishina corrections decrease scattering cross section

% physical constantes
phys = cphys;

te = phys.e .* te ./ phys.me ./ phys.c .^ 2;
ti = phys.e .* ti ./ (Ai .* phys.ua) ./ phys.c .^ 2;

ex_sp  = (te + ti) .^ -(3/2);
ex_rel = exp(-1 ./ te) ./ besselk(2,1 ./ te) .* sqrt(te ./ (te + ti) .^ 3) .* ...
         (2 .* (te + ti) .^ 2 + 2 .* (te + ti) +1);
     
ec = sqrt(2 .* pi) ./ 2 .* ex_rel ./ ex_sp; 
ec(~isfinite(ec)) = 1;
% verifier l'effet final avec l'approximation de Cordey = OK

function spt = spitzer_1(te,ti,ne,ni,Zi,Ai,lnei)
%from:
% Modification of classical Spitzer ionelectron energy transfer rate
% for large ratios of ion to electron temperatures, T.H. Rider and P.J. Catto,
% Phys. Plasmas 2, 18731885 (1995), https://doi.org/10.1063/1.871274

% physical constantes
phys = cphys;

te = phys.e .* te;
ti = phys.e .* ti;
% some units problem or missing numerical factor of value phys.ee / pi
%spt = 4 .* sqrt(2 .* pi.* Ai .* phys.ua .* phys.me) .* Zi .^ 2 .* phys.e .^ 4 .*  ...
%      ni .* ne .* lnei ./ (Ai .* phys.ua .* te + phys.me .* ti) .^ (3/2); 
spt = 3 .* phys.me ./ (Ai .* phys.ua) .*  ...
      (sqrt(2) .* phys.e .^ 4 ./ (12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 .* sqrt(phys.me))) .*  ...
      Zi .^ 2  .* ne .* ni .* lnei./ ...
      (te + phys.me ./ (Ai .* phys.ua) .* ti) .^ (3/2);


% term in large Ti/Te correction formula  
function out = term(te,ti,ne,ni,Zi,Ai)
% correction for large Ti/Te
%from:
% Modification of classical Spitzer ionelectron energy transfer rate
% for large ratios of ion to electron temperatures, T.H. Rider and P.J. Catto,
% Phys. Plasmas 2, 18731885 (1995), https://doi.org/10.1063/1.871274
% formula 56

% physical constantes
phys = cphys;

out = Zi .^ 2 .*  (ni ./ ne) .* (phys.me ./ phys.ua ./ Ai) .* (ti ./te);






