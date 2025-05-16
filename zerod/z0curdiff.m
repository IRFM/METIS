% pour test :
% [psiv,dpsidtv,qv,jresv,eparv,Fv,bpolv,liv,factv] = ...
%   z0curdiff(temps,x,eta,jli - jres,ptot,spr,kx,Raxe, ...
%   fdia,a,d,ipv,vloop,lif);
function [psi,dpsidt,q,jmoy,jres,epar,F,bpol,li, ...
          rmx,vpr,grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr,phi,dphidx, ...
          difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dpdpsi] = ...
	   z0curdiff(t,x,eta,jni,ptot,Sp,kx,...
	  Raxe,b0,a,dx,ip,vloop,liin,asser,Vp,Fin,qdds,peri,Rsepa,Zsepa, ...
	  amorti,drmdt,qeff,psi_old,phi_old,difcurconv,lao_change,evolution,edge_flux, ...
	  mode_expo_inte,cronos_regul,ddsmode,w1,epsq,flux_ip,L_ext,q0_dds_trig,betap1crit)

%save('loc1');	  
% utilisation du mexfile si disponible	 
if isappdata(0,'MEXSOLVER_IN_METIS')
	mexsolver = getappdata(0,'MEXSOLVER_IN_METIS');
else	
	mexsolver =  [];
end
if isempty(mexsolver)
	repmex=which(strcat('mexpde1dsolver.',mexext));
	if ~isempty(repmex)
		mexsolver = 1;
	else
		mexsolver = 0;	
	end
	setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
end  

% just for backward compatibility
if iscomplex(ptot(:))
    hollow = 1;
else
    hollow = 0;
end
ptot = real(ptot);

% decodage
freebie    = imag(evolution);
evolution  = real(evolution);


% pour test 
%mode_expo_inte =1;
%mode_expo_inte = 0;

% on travail a ip donne pour le 0D
% definition
ve = ones(size(x));
vt = ones(size(t));

% decodage qdds
s1crit = imag(qdds);
qdds   = real(qdds);

% r0
r0  =  Raxe(:,end);
Raxe = max(r0 * ve,min((r0 + a ./ 4) * ve,Raxe));
if isempty(Fin)
	Fin =  (r0 .* b0) * ve;
end
% element de geometrie
expo  = max(difcurconv);
rho   = sqrt(phi_old ./ pi ./ (b0 * ve));

if any(max(rho,[],2) > max(Raxe(:,1), a .* sqrt(1 + kx(:,end) .^ 2)))
    indbad = find(max(rho,[],2) > max(Raxe(:,1), a .* sqrt(1 + kx(:,end) .^ 2)));
    fprintf('B');
    rho(indbad,:) = (a(indbad) * x) .* sqrt(1 + kx(indbad,:) .^ 2);
end
%save('loc2');	
if isappdata(0,'EQUILIBRIUM_EXP') 
    equi_ext    = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
    rmx         = equi_ext.rmx;
    rm          = equi_ext.rmx(:,end);
    rho         = rmx;
    phi         = equi_ext.phi_tor;
    dphidx      = equi_ext.dphidx_tor;
    vpr         = equi_ext.vpr_tor;
    spr         = equi_ext.spr_tor;
    grho2r2     = equi_ext.grho2r2;
    r2i         = equi_ext.r2i;
    ri          = equi_ext.ri;
    C2          = equi_ext.C2;
    C3          = equi_ext.C3;
    grho        = equi_ext.grho;
    grho2       = equi_ext.grho2;
    phiplasma   = equi_ext.phiplasma;
    badjac      = zeros(size(rmx));
    indbad      = [];
else
    % le probleme d'imprecision numerique arrive apres cette ligne
    [phi,dphidx,vpr,grho2r2,r2i,ri,C2,C3,grho,grho2,phiplasma,badjac] = gg0d(x,Raxe,a,kx,dx,Vp,Sp,Fin,Rsepa,Zsepa,rho,expo);
    %  save(tempname,'x','Raxe','a','kx','dx','Vp','Sp','Fin','Rsepa','Zsepa','rho','expo', ...
    %       'phi','dphidx','vpr','grho2r2','r2i','ri','C2','C3','grho','grho2','phiplasma');
    % save('loc3');
    % le probleme d'imprecision numerique arrive avant cette ligne
    % securite
    if any(C2(:,end) <= 0) | any(any(vpr(:,2:end) <= 0))| any(any(phi(:,2:end) <= 0))| any(any(badjac(:,2:end)))
        indbad = find((C2(:,end) <= 0)| any(vpr(:,2:end) <= 0,2) | any(phi(:,2:end) <= 0,2)| any(badjac(:,2:end),2));
        %fprintf('number of time slices corrected for geometry = %d\n', length(indbad))
        if ~isempty(Rsepa) & ~isempty(Zsepa)
            indsepabad = indbad;
        else
            indsepabad = [];
        end
        [phi(indbad,:),dphidx(indbad,:),vpr(indbad,:),grho2r2(indbad,:),r2i(indbad,:),ri(indbad,:),C2(indbad,:),C3(indbad,:),grho(indbad,:),grho2(indbad,:),phiplasma(indbad)] = ...
            gg0d(x,Raxe(indbad,:),a(indbad),kx(indbad,:),dx(indbad,:),Vp(indbad),Sp(indbad),Fin(indbad,:), ...
            Rsepa(indsepabad,:),Zsepa(indsepabad,:),rho(indbad,:),0);
        %fprintf('MORPHING : data inconsistency @ %g s \n',t(indbad));
        fprintf('m');
    end
    if any(C2(:,end) <= 0) | any(any(vpr(:,2:end) <= 0))| any(any(phi(:,2:end) <= 0))
        fprintf('M');
        indbad = find((C2(:,end) <= 0)| any(vpr(:,2:end) <= 0,2) | any(phi(:,2:end) <= 0,2));
        indok  = find((C2(:,end) > 0) & all(vpr(:,2:end)> 0,2) & all(phi(:,2:end) > 0,2));
        if ~isempty(indok)
            phi(indbad,:)     = ones(length(indbad),1) * mean(phi(indok,:),1);
            dphidx(indbad,:)  = ones(length(indbad),1) * mean(dphidx(indok,:),1);
            vpr(indbad,:)     = ones(length(indbad),1) * mean(vpr(indok,:),1);
            grho2r2(indbad,:) = ones(length(indbad),1) * mean(grho2r2(indok,:),1);
            r2i(indbad,:)     = ones(length(indbad),1) * mean(r2i(indok,:),1);
            ri(indbad,:)      = ones(length(indbad),1) * mean(ri(indok,:),1);
            C2(indbad,:)      = ones(length(indbad),1) * mean(C2(indok,:),1);
            C3(indbad,:)      = ones(length(indbad),1) * mean(C3(indok,:),1);
            grho(indbad,:)    = ones(length(indbad),1) * mean(grho(indok,:),1);
            grho2(indbad,:)   = ones(length(indbad),1) * mean(grho2(indok,:),1);
        end
    end
    % c'est la coordonnee de flux toroidal rho de cronos
    rmx   = sqrt(phi ./ pi ./ (b0 * ve));
    rm    = rmx(:,end);
    % element de surface pour le calcul du courant total
    spr = vpr .* ri ./ 2 ./ pi;
end
% constante
mu0  = 4 .* pi .* 1e-7;
%save('loc4');	

% la valeur de C2 au bord est numeriquement difficile a calculee avec peu de points
% renormalisation en utilisant qeff
psid1edge   = -(2*pi) .* mu0 .* rm .* ip ./ C2(:,end);
qedge       =  max(1,Fin(:,end) .* C3(:,end) ./ 4 ./ pi .^ 2 ./ abs(psid1edge) .* rm);	
% trop complique -pas stable numeriquement
%gq          = mean(dx(:,end)) .^ 3; % fixe avec des runs helena (c'est un gain de boucle via li)
%qfit        = max(qedge,(1 + gq) .* max(1,qeff) - gq .* qedge);
qfit        = max(qedge,(qeff + qedge) ./ 2);
%C2(:,end)   = max(0.3,min(3,max(1,qfit) ./ qedge)) .* C2(:,end);
if ~isappdata(0,'EQUILIBRIUM_EXP') 
    C2(:,end)   = max(0.3,min(3,max(1,qfit) ./ qedge)) .* C2(:,end);
end

%save('loc5');	

if cronos_regul == 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modif DM 16/09/16: regularisation du facteur de securite au centre
    % kreg >= 2, kreg = 2 -> regularisation faible, kreg = 3 -> meilleure regularisation
    % Modif DM 18/11/16: N'utiliser kreg = 3 que si jni est suffisamment regulier pres de l'axe magnetique.
    kjni=find(mean(abs(jni(:,1:5)),2) > 0.75*max(mean(abs(jni(:,1:5)),2)));
    if mean((max(jni(kjni,1:5),[],2)-min(jni(kjni,1:5),[],2))./(1+abs(mean(jni(kjni,1:5),2)))) > 0.5
        kreg=2;
    else
        kreg=3;
    end
    krm1=kreg-1;
    krp1=kreg+1;
    krp2=kreg+2;
    Areg=eye(2)/[x(1,krp1:krp2);x(1,krp1:krp2).^2];  
end

if isempty(psi_old)
	% calcul de psi initial (cf. maple)
	lipj = liin(1);
	lipj(~isfinite(lipj)) = 1;
	ip(ip <= 0) =1;
	piqj = max(0.1,min(10,(exp(lipj)-1.65)./0.89));
	jini = (1 - x.^2) .^ piqj;
	iini = rm(1) .* trapz(x,spr(1,:) .* jini,2);
	jini = jini .* ip(1) ./ iini;
	jini = max(jni(1,:),jini);
	iini = rm(1) .* trapz(x,spr(1,:) .* jini,2);
	jini = jini .* ip(1) ./ iini;
	interi = -  rm(1) .* cumtrapz(x,mu0 .* ri(1,:) .* vpr(1,:) .*jini,2);
	inter = interi ./ max(eps,C2(1,:));
	inter(1) = 0;
	psi_init =  rm(1) .* cumtrapz(x,inter,2);
	psi_init =  psi_init - psi_init(end);
	% normalisation precise
	% dpsidx s'annule au centre
	psid1    = pdederive(x,psi_init,0,2,2,1);
	% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
	psid2    = pdederive(x,psi_init,1,0,2,2);
	% c2cd1 ne s'annule pas  au centre (c2c est nul au centre)
	% il ne doit pas changer de signe au bord
	% pour que ip soit precisement defini, sa derive seconde s'annulle au bord
	c2cd1    = pdederive(x,C2(1,:),2,0,2,1);
	% c2cd1 doit toujours etre > 0
	c2cd1(c2cd1 < 0)   = 0;
	% calcul de jmoy global
	jmoy1      =  - psid1 .* c2cd1       ./   (mu0 .* max(eps,vpr(1,:)) .*ri(1,:) .* (rm(1) .*ve) .^ 2);
	jmoy1(1)   = 0;
	jmoy2     =  - grho2r2(1,:)  .* psid2   ./   (mu0  .*ri(1,:) .* (rm(1) .* ve).^ 2) ;
	jmoy      =  jmoy1 + jmoy2;
	% calcul de j0 est maintenant ok, il faut regularsier j(2)
    if cronos_regul ~= 5
        jmoy(1:2)  = pchip(x(3:11),jmoy(3:11),x(1:2));
    end
	ipinit   = (rm(1) .* trapz(x,spr(1,:) .* jmoy,2) - ...
		C2(1,end) ./ (2 .* pi .* mu0 .* rm(1)) .* psid1(end)) ./ 2;
	psi_init = psi_init .* ip(1) ./ ipinit;
else
	psi_init = psi_old(1,:);
end
%save('loc6');	
if evolution == 1
        psi_old(length(t) - 1,:) = psi_old(length(t) - 2,:);
        psi_old(length(t),:) = psi_old(length(t) - 2,:);
end


% correction du au chagement de coordonnees Lao -> flux toroidal
% uniquement pour le soveur de diffusion 
% les donnees calculees avec gg0d reste en coordonnee de Lao
if lao_change == 1
	eta_tor      = exp(tsplinet(rmx,log(eta),rmx(:,end) * x));
	jni_tor      = tsplinet(rmx,jni,rmx(:,end) * x);
	Fin_tor      = tsplinet(rmx,Fin,rmx(:,end) * x);
	grho2r2_tor  = tsplinet(rmx,grho2r2,rmx(:,end) * x);
	r2i_tor      = tsplinet(rmx,r2i,rmx(:,end) * x);
	ri_tor       = tsplinet(rmx,ri,rmx(:,end) * x);
	C2_tor       = tsplinet(rmx,C2,rmx(:,end) * x);
	C3_tor       = tsplinet(rmx,C3,rmx(:,end) * x);
	vpr_tor      = tsplinet(rmx,vpr,rmx(:,end) * x);
	spr_tor      = tsplinet(rmx,spr,rmx(:,end) * x);	
	phi_tor      = tsplinet(rmx,phi,rmx(:,end) * x);	
	ptot_tor     = tsplinet(rmx,ptot,rmx(:,end) * x);	
	x_tor        = tsplinet(rmx,vt * x ,rmx(:,end) * x);	

    if any(any(C2_tor < 0)) || any(any(vpr_tor < 0))|| any(any(r2i_tor < 0))|| any(any(grho2r2_tor < 0)) || any(any(Fin_tor < 0))
        fprintf('W');
        indbad = find(any(C2_tor(:,2:end) <= 0,2)| any(vpr_tor(:,2:end) <= 0,2) |  ...
                      any(phi_tor(:,2:end) <= 0,2) | any(grho2r2_tor(:,2:end) <= 0,2) | ...
                      any(Fin_tor(:,2:end) <= 0,2));
        ptot_tor(indbad,:)    =  ptot(indbad,:);
        Fin_tor(indbad,:)     =  Fin(indbad,:);
        eta_tor(indbad,:)     =  eta(indbad,:);
        jni_tor(indbad,:)     =  jni(indbad,:);
        phi_tor(indbad,:)     =  phi(indbad,:);
        vpr_tor(indbad,:)     =  vpr(indbad,:);
        spr_tor(indbad,:)     =  spr(indbad,:);
        grho2r2_tor(indbad,:) =  grho2r2(indbad,:);
        r2i_tor(indbad,:)     =  r2i(indbad,:) ;
        ri_tor(indbad,:)      =  ri(indbad,:);
        C2_tor(indbad,:)      =  C2(indbad,:);
        C3_tor(indbad,:)      =  C3(indbad,:);
        vtx                   =  vt*x; 
        x_tor(indbad,:)       =  vtx(indbad,:);
    end
    
	if evolution == 1
		psi      = tsplinet(rmx,psi_old,rmx(:,end) * x);
	end
else
	eta_tor      = eta;
	jni_tor      = jni;
	Fin_tor      = Fin;
	grho2r2_tor  = grho2r2;
	r2i_tor      = r2i;
	ri_tor       = ri;
	C2_tor       = C2;
	C3_tor       = C3;
	vpr_tor      = vpr;
	spr_tor      = spr;
	phi_tor      = phi;
	ptot_tor     = ptot;	
	x_tor        = vt * x;	
	
	if evolution == 1
		psi      = psi_old;
	end
end
	
% allignement du flux au bord pour le premier temps
if evolution == 0
	edge_flux  = edge_flux - edge_flux(1);
end

% valeur au bord
dpsidt_edge  = - vloop ./ 2 ./ pi;
tt           = (min(t):min(diff(t)):max(t))';
if length(tt) > 1e6
    tt = linspace(min(t),max(t),1e6+1)';
end
psi_edge     = interp1(tt,cumtrapz(tt,interp1(t,dpsidt_edge,tt,'linear','extrap')),t,'linear','extrap');
% aservissement sur le flux
psi_edge(asser == -1) = edge_flux(asser == -1);
psi_edge(asser == -2) = edge_flux(asser == -2);

% critical curvature (base on critical shear)
term_a         = Fin_tor.* C3_tor ./ 4 ./ pi .^ 2 .* (rmx(:,end) * ve);
term_a_d1      = pdederive(x,term_a,2,2,2,1); 
volume         = (rmx(:,end) * ve) .* cumtrapz(x,vpr_tor,2);
d2psidpsi_crit = term_a_d1 ./ term_a - s1crit .* (rmx(:,end) * ve) .* vpr_tor ./ 2 ./ volume;
d2psidpsi_crit(:,1) = 1 ./ eps;

% creation des matrices pour la resolution de l'equation de diffusion
%coeficient de l'equation du flux
Ajj      =  grho2r2_tor ./ r2i_tor .* eta_tor ./ mu0 ./ (rm * ve) .^2;
warning off
lnc2f    =  log( C2_tor ./ Fin_tor);
lnc2f(:,1)  = 2 .* lnc2f(:,2) -  lnc2f(:,3);
warning on
%Bjj      =  (grho2r2_tor ./ r2i_tor .* eta_tor ./ mu0 ./ (rm * ve) .^ 2).*  ...
%             rpdederive(x,lnc2f,2,2,2,1)  ...
%            + (vt * x) ./ (rm * ve) .* (drmdt * ve);
% ajout du terme en dB0dt
if (evolution == 1) && (freebie == 1)
    bpcor    =  (vt * x) ./ 2 .* ((z0dxdt_freebie(b0,t) ./ b0) * ve);
else
    bpcor    =  (vt * x) ./ 2 .* ((z0dxdt(b0,t) ./ b0) * ve);
end
Bjj      =  (grho2r2_tor ./ r2i_tor .* eta_tor ./ mu0 ./ (rm * ve) .^ 2).*  ...
             pdederive(x,lnc2f,2,2,2,1)  ...
            + (vt * x) ./ (rm * ve) .* (drmdt * ve) + bpcor;
Cjj      = 0 .* ve;
Djj      = (b0 * ve) .* eta_tor ./ Fin_tor ./ r2i_tor .* jni_tor;


%  figure(51);
%  clf
%  plot(t,Bjj,'b',t,bpcor,'r');
%  drawnow


% dimension
nbrho = length(x);
dx    = mean(diff(x));
% boucle sur le temps
A           = zeros(nbrho,1,1);
B           = zeros(nbrho,1,1);
C           = zeros(nbrho,1,1);
D           = zeros(nbrho,1);
AP          = zeros(nbrho,1,1);
BP          = zeros(nbrho,1,1);
CP          = zeros(nbrho,1,1);
DP          = zeros(nbrho,1);
% conditions au limites
V0          = zeros(2,1);
V1          = zeros(2,1);
T0          = zeros(2,1,3);
T1i          = zeros(2,1,3);
T1v          = zeros(2,1,3);
T1h          = zeros(2,1,3);

T0(1,1,2)   = 1;
T0(2,1,2)   = 1;

T1i(2,1,2)   = 1;
T1i(1,1,2)   = 1;

T1v(2,1,1)   = 1;
T1v(1,1,1)   = 1;

%T1h(1,1,1)   = 1;
%T1h(2,1,2)   = 1;

if evolution == 0
	psi         = 0 .* (vt * ve);
	psi(1,:)    = psi_init;
end

dpsidt      = zeros(size(psi));
maskdds      = 0 .* vt;
indice_inv   =  0 .* vt;

% flag convergence
convf       = 0 .* vt;
convf(1)    = 1;
warnex      = 0;

% boucle sur le temps
if evolution == 1
	inddeb = length(t) -1;
	convf(1:(length(t) -2)) = 1;
else
	inddeb = 2 ;
end

%save('loc7');	

for k = inddeb:length(t)
	A(:,1,1)    = zreshape(Ajj(k-1,:),nbrho,1,1);
	B(:,1,1)    = zreshape(Bjj(k-1,:),nbrho,1,1);
	D(:,1)      = zreshape(Djj(k-1,:),nbrho,1,1);
	AP(:,1,1)   = zreshape(Ajj(k,:),nbrho,1,1);
	BP(:,1,1)   = zreshape(Bjj(k,:),nbrho,1,1);
	DP(:,1)     = zreshape(Djj(k,:),nbrho,1,1);
        F           = psi(k-1,:)';
	dtin        = t(k) - t(k - 1);
	delta_psi   = psi_edge(k) - psi_edge(k - 1);
	        if asser(k) == -2
		FP          = F - F(end) + psi_edge(k);
   		V1(1,1)     = flux_ip(k-1);
   		V1(2,1)     = flux_ip(k);

		T1h(1,1,2)   =  L_ext(k-1) ./ (2.* pi .* mu0 .* rm(k-1)) .* C2_tor(k-1,end);
		T1h(2,1,2)   =  L_ext(k)  ./  (2.* pi .* mu0 .* rm(k))   .* C2_tor(k,end);
		T1h(1,1,1)   = 1;
		T1h(2,1,1)   = 1;
		T1           = T1h;

	elseif asser(k) == -1
               % allignement du flux au bord en cas de commutation de l'asservissement (t_switch finit)
		if evolution == 0
		      if asser(k-1) ~= -1
			    psi_edge(k:end) = psi_edge(k:end) - edge_flux(k-1) + psi(k-1,end);
		      end
		end

		FP          = F - F(end) + psi_edge(k);
   		V1(1,1)     = F(end);
   		V1(2,1)     = psi_edge(k);
		T1          = T1v;	
	elseif asser(k) == 1
		FP          = F + delta_psi;
   		V1(1,1)     = F(end);
   		V1(2,1)     = FP(end);
		T1          = T1v;
	else
		FP          = F;
   		V1(1,1)     = -(2*pi) .* mu0 .* rm(k-1) .* ip(k-1) ./ C2_tor(k-1,end);
   		V1(2,1)     = -(2*pi) .* mu0 .* rm(k) .* ip(k) ./ C2_tor(k,end);
		T1          = T1i;
		
	end
	
	if evolution == 1	
                if mode_expo_inte == 1
			fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,-1,dx,dtin);
			if any(~isfinite(fp))
			    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
			end
		elseif  isdeployed && (mexsolver == 1)
                        try
			    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
                        catch
                            disp(lasterr)
			    warning('mexpde1dsolver is not working; switch to mfile solver');
			    setappdata(0,'MEXSOLVER_IN_METIS',0);
                            mexsolver   = 0;
			    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
                        end
                elseif mexsolver == 1
			fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
		else
			fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
		end
		% control de la condition a limite
		psid1k    = pdederive(x,fp',0,2,2,1);
	else
		
                if mode_expo_inte == 1
			fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,-1,dx,dtin);
			if any(~isfinite(fp))
			      fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0.5,dx,dtin);
			end
		elseif  isdeployed && (mexsolver == 1)
                        try
			    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0.5,dx,dtin);
                        catch
                            disp(lasterr)
			    warning('mexpde1dsolver is not working; switch to mfile solver');
			    setappdata(0,'MEXSOLVER_IN_METIS',0);
                            mexsolver   = 0;
			    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0.5,dx,dtin);
                       end
		elseif mexsolver == 1
			fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0.5,dx,dtin);
		else
			fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0.5,dx,dtin);
		end
		% control de la condition a limite
		psid1k    = pdederive(x,fp',0,2,2,1);
		
	
		dF        = abs(fp - F);
		if (abs(dF(end)) > abs(2 .* delta_psi)) & (asser(k) == 1)
		        if  isdeployed && (mexsolver == 1)
				try
				    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
				catch
				    disp(lasterr)
				    warning('mexpde1dsolver is not working; switch to mfile solver');
				    setappdata(0,'MEXSOLVER_IN_METIS',0);
                                    mexsolver   = 0;
				    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
				end
			elseif mexsolver == 1
				fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
			else
				fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
			end
			% control de la condition a limite
			psid1k    = pdederive(x,fp',0,2,2,1);
		elseif ((abs(psid1k(end) + (2*pi) .* mu0 .* rm(k) .* ip(k) ./ C2_tor(k,end)) ./ max(eps,abs(psid1k(end)))) > 1e-2)  & (asser(k) == 0)
			if  isdeployed && (mexsolver == 1)
				try
				    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);	
				catch
				    disp(lasterr)
				    warning('mexpde1dsolver is not working; switch to mfile solver');
				    setappdata(0,'MEXSOLVER_IN_METIS',0);
                                    mexsolver   = 0;
				    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);	
				end
			elseif mexsolver == 1
				fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);	
			else
				fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);	
			end
			% control de la condition a limite
			psid1k    = pdederive(x,fp',0,2,2,1);
		elseif any(psid1k>0)
			if  isdeployed && (mexsolver == 1)
				try
				    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
				catch
				    disp(lasterr)
				    warning('mexpde1dsolver is not working; switch to mfile solver');
				    setappdata(0,'MEXSOLVER_IN_METIS',0);
                                    mexsolver   = 0;
				    fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
				end
			elseif mexsolver == 1
				fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
			else
				fp          = pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,dx,dtin);
			end
			if warnex == 0
				fprintf('!');
				warnex = 1;
			end
			% control de la condition a limite
			psid1k    = pdederive(x,fp',0,2,2,1);
				
		end
	end	
	
	dF        = abs(fp - F);
	% control de la condition a limite
	%psid1k    = pdederive(x,fp',0,2,2,1);
	convf(k)     =1;
	
	if (abs(dF(end)) > abs(2 .* delta_psi)) & (asser(k) == 1)
		fp          = fp - fp(end) + F(end) + delta_psi;
		convf(k)    = 0;
	elseif any(psid1k>0) & (psid1k(end) > 0)
		psid1k      = pdederive(x,F',0,2,2,1);
		psid1k(end)  =  - (2*pi) .* mu0 .* rm(k) .* ip(k) ./ max(C2_tor(k,:));
		fp          = cumtrapz(x,psid1k)';
		fp          = fp - fp(end) + F(end) + delta_psi;
		convf(k)     =0;
	end	
	psi(k,:)    = fp';

    if cronos_regul == 5
        psid1k    = pdederive(x,fp',0,2,2,1);
        % Modif DM 16/09/16: regularisation du facteur de securite au centre.
        % En coordonnees toroidales normalisees, phi = phimax * x^2 = pi * rm^2 * B0 * x^2
        % => q(x) = -(dphi/dx)/(2*pi*dpsi/dx) = - rm^2 * B0 * x/(dpsi/dx).
        % q(0) fini => on impose dpsi/dx = a1*x + a2*x^2 au voisinage de x = 0
        % => q(0) = -B0 * rm^2/a1 pour x = 0.
        % On choisit a1 et a2 tels que dpsi/dx soit conserv?? pour x=x(kreg+1) et x=x(kreg+2).
        % Ceci conserve Bpol, donc le courant toroidal ?? l'int??rieur de x = x(kreg+1).
        % On int??gre le dpsi/dx regularise pour obtenir le nouveau psi en x=x(1:kreg-1).
        % Les nouveaux dpsi/dx et psi satisfont dpsi/dx = pdederive(x,psi,0,2,2,1).
        
        % Calcul des coefficients du d??veloppement limite
        Breg=psid1k(krp1:krp2);
        Xreg=Breg*Areg;
        % Calcul de dpsi/dx
        psid1reg=Xreg*[x(1,1:kreg);x(1,1:kreg).^2];
        
        % Correction de psi(x) au centre: Integration par inversion de pdederive(x,psi_init,0,2,2,1)
        psi(k,1:krm1) = psi(k,3:krp1)-psid1reg(2:kreg).*(x(1,3:krp1)-x(1,1:krm1));
    end
	
	% correction de psi pour les dds
	%psid1k    = pdederive(x,fp',0,2,2,1);
	if (qdds > 0) & (convf(k) == 1) 
        psid1dds  =  - Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve)./ qdds;
        psid1dds(1) = 0; % symetrie
        indcor    = max(find(psid1k < psid1dds));
        if (indcor > 1) & (indcor < 18)
            psicor = cumtrapz(x,max(psid1k,psid1dds),2);
            psi(k,1:indcor) = psicor(1:indcor) - psicor(indcor+1) + psi(k,indcor+1);
            maskdds(k) = 1;
            indice_inv(k) = indcor;
        end
    elseif (qdds < 0 ) & (convf(k) == 1)
        psid1trig  =  - Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve)./ abs(qdds);
        psid1trig(1) = 0; % symetrie
        indtrig    = max(find(psid1k < psid1trig));
        % case with additional trigger on critical magnetic shear
        if ~isempty(indtrig) && ((s1crit > 0) || (betap1crit > 0)) && (indtrig > 1) && (indtrig < 18)
            indtrig_mem_ = indtrig;
            d2psidpsi = pdederive(x,psi(k,:),2,2,2,2) ./ psid1k;
            
            % betap1
            bpol1k  = sqrt(grho2r2_tor(k,:)) .* abs(psid1k) ./ (rm(k)*ve);
            emp1k   = trapz(x(1:indtrig),bpol1k(1:indtrig) .^ 2 .* vpr_tor(k,1:indtrig),2) ./ 2 ./ mu0;
            if emp1k(1) > 0
                betap1k  = trapz(x(1:indtrig),(ptot(k,1:indtrig) - ptot(k,indtrig)*ones(1,indtrig)) .* vpr_tor(k,1:indtrig),2) ./ emp1k;
            else
                betap1k = 0;
            end
            
            if ((d2psidpsi(indtrig) > d2psidpsi_crit(k,indtrig)) || (s1crit == 0))  && ((betap1k < betap1crit) || (betap1crit == 0))
                if q0_dds_trig > 0
                    psid1trig0  =  - Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve)./ abs(q0_dds_trig);
                    % central value is 0
                    if psid1k(2) > psid1trig0(2)
                        % no ST
                        indtrig = [];
                        %fprintf('n');
                        %else
                        %fprintf('u');
                    end
                else
                    % no ST
                    indtrig = [];
                end
            end
        end
        % accrochage en frequence
        if isappdata(0,'TE_EXP')
            TE_EXP =getappdata(0,'TE_EXP');
            if isfield(TE_EXP,'dte_trig')
                te0 = TE_EXP.te(:,1);
                dte = diff(te0) ./ (te0(1:end-1) + te0(2:end)) .* 2;
                ind_dte = find(TE_EXP.temps >= t(k),1);
                if ~isempty(ind_dte)
                    if isfield(TE_EXP,'ton_trig')
                        ton_trig = TE_EXP.ton_trig;
                    else
                        ton_trig = t(1);
                    end
                    if isfield(TE_EXP,'toff_trig')
                        toff_trig = TE_EXP.toff_trig;
                    else
                        toff_trig = t(end);
                    end
                    if (dte(ind_dte) < -abs(TE_EXP.dte_trig)) && (t(k) >= ton_trig) && (t(k) <= toff_trig)
                        fprintf('c');
                    else
                        indtrig = [];
                    end
                end
            end
        end
        
        if ~isempty(indtrig) && (indtrig > 1) && (indtrig < 18)
            if ddsmode > 0
                q         = - Fin_tor(k,2:end) .* C3_tor(k,2:end) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve(2:end)) ./ psid1k(2:end);
                %% q         = interpos(x(2:end),q,x,-1,[1 1],[1e32 1e32],[length(x) ones(size(q(2:end)))]);
                %%figure(21);plot(x(2:end),q,x,interpos(x(2:end),q,x,-1,[1 1],[1e32 1e32],[length(x) ones(size(q(2:end)))]),x,max(0,pchip(x(2:end),q,x)));pause(0.1);
                q         = max(0,pchip(x(2:end),q,x));
                %[qnew,indcor] = z0dds(x,psi(k,:),q,ddsmode,w1,epsq);
                indcor = 0;
                try
                    [qnew,indcor] = z0dds(x,psi(k,:),q,ddsmode,w1,epsq);
                catch
                    fprintf('s');
                end
                if ~isempty(indcor) && (indcor > 0)
                    indcor = max(1,min(17,indcor));
                    inte   = -  Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve) ./ qnew;
                    psicor = cumtrapz(x,inte,2);
                    % raccordement a rmix :
                    psicor(1:indcor)       = psicor(1:indcor) - psicor(indcor) + psi(k,indcor);
                    psicor((indcor+1):end) = psi(k,(indcor+1):end);
                    maskdds(k) = 1;
                    indice_inv(k) = indcor;
                    psi(k,:) = psicor;
                end
            else
                psid1dds  =  - Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve)./ max(1,ceil(abs(qdds)));
                psid1dds(1) = 0; % symetrie
                indcor    = max(1,min(17,max(indtrig,max(find(psid1k <= psid1dds)) + 1)));
                if ~isempty(indcor) && (indcor> 0)
                    psicor = cumtrapz(x,max(psid1k,psid1dds),2);
                    psi(k,1:indcor) = psicor(1:indcor) - psicor(indcor+1) + psi(k,indcor+1);
                    maskdds(k) = 1;
                    indice_inv(k) = indcor;
                end
            end
        end
    end

	% correction instabilite trou de courant
	psid1hole  =  - Fin_tor(k,:) .* C3_tor(k,:) ./ 4 ./ pi .^ 2 .* (rm(k,:) * ve)./ (pi .* (qeff(k) * ve));
	indcor    = max(find(psid1k > psid1hole));
	if (indcor > 1) & (indcor < 18) & (convf(k) == 1)
		psicor = cumtrapz(x,min(psid1k,psid1hole),2);
		psi(k,1:indcor) = psicor(1:indcor) - psicor(indcor+1) + psi(k,indcor+1);		
	end
 
end
%save('loc8');	

% le mode asser est 0 ou 1
asser = abs(asser);

% indicateur de non convergence 
if any(convf == 0)
	idx = round(length(find(convf == 0)) ./ length(t) .* 10);
	if idx > 3
		fprintf('X');
	else
		fprintf('x');	
	end
	
	difcurconv  = convf .* difcurconv + (~convf) .* max(0,expo - 4);
end	

% dpsidx s'annule au centre
psid1    = pdederive(x,psi,0,2,2,1);
psid1cons         = -(2*pi) .* mu0 .* rm .* ip ./ C2_tor(:,end);
psid1(asser == 0,end) = psid1cons(asser ==0);
if any(asser == 2)
  psid1flux = (flux_ip - psi(:,end)) .* (2.* pi .* mu0 .* rm) ./ (L_ext .* C2_tor(:,end));
  psid1(asser == 2,end) = psid1flux(asser == 2);
end
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
psid2    = pdederive(x,psi,1,0,2,2);
%
% natural regulation thank to a remark from D. Moreau
%
if cronos_regul == 4
    % near magnetic axis
    q_not = b0 ./ (psi(:,1) - psi(:,2)) .* rm .^ 2 .* x(2) .^ 2 ./ 2;
    psid2(:,1) = - b0 ./ q_not .* rm .^ 2;
end
% c2cd1 ne s'annule pas  au centre (c2c est nul au centre)
% il ne doit pas changer de signe au bord
% pour que ip soit precisement defini, sa derive seconde s'annulle au bord
c2cd1    = pdederive(x,C2_tor,2,0,2,1);
% c2cd1 doit toujours etre > 0
c2cd1(c2cd1 < 0)   = 0;
% calcul de jmoy global
jmoy1      =  - psid1 .* c2cd1       ./   (mu0 .* max(eps,vpr_tor) .*ri_tor .* (rm *ve) .^ 2);
jmoy1(:,1) = 0;
jmoy2     =  - grho2r2_tor  .* psid2   ./   (mu0 .* ri_tor .* (rm * ve).^ 2) ;
jmoy      =  jmoy1 + jmoy2;
% calcul de j0 est maintenant ok
if cronos_regul == 5
    % Modif DM 16/09/16: calcul de jmoy(:,1)
    % jmoy(:,1:2)  = pchip(x(3:11),jmoy(:,3:11),x(1:2));
    % figure(9);clf;plot(x,jmoy(end,:),'b');
    jmoy(:,1)  = interp1(x(2:3)',jmoy(:,2:3)',x(1),'linear','extrap')';
    % figure(9);hold on;plot(x,jmoy(end,:),'r');
elseif cronos_regul == 4
  % nothing
elseif cronos_regul == 0
  jmoy(:,1:2)  = pchip(x(3:11),jmoy(:,3:11),x(1:2));
else
  jmoy(:,2)  = pchip(x([1,3:11]),jmoy(:,[1,3:11]),x(2));
end
%jmoy(:,1:2)  = pchip(x(3:11),jmoy(:,3:11),x(1:2));
%jmoy(:,2)  = pchip(x([1,3:11]),jmoy(:,[1,3:11]),x(2));

% jmoy pour FREEBIE
c2cd1_    = pdederive(x,C2_tor,2,2,2,1);
psid2_    = pdederive(x,psi,1,2,2,2);
jmoy1       =  - psid1 .* c2cd1_       ./   (mu0 .* max(eps,vpr_tor) .*ri_tor .* (rm *ve) .^ 2);
if cronos_regul == 5
    % Modif DM 16/09/16: calcul de jmoy1(:,1)
    % jmoy1(:,1)  = pchip(x(2:end),jmoy1(:,2:end),x(1));
    % figure(10);clf;plot(x,jmoy1(end,:),'b');
    jmoy1(:,1)  = interp1(x(2:3)',jmoy1(:,2:3)',x(1),'linear','extrap')';
    % figure(10);hold on;plot(x,jmoy1(end,:),'r');
else
    jmoy1(:,1)  = pchip(x(2:end),jmoy1(:,2:end),x(1));
end
jmoy2       =  - grho2r2_tor  .* psid2_   ./   (mu0 .* ri_tor .* (rm * ve).^ 2) ;
jgs         =  jmoy1 + jmoy2;

% ip en sortie 
ipout   = rm .* trapz(x,spr_tor .* jmoy,2);
%ipout((asser == 0) | ((asser == 1) & (convf == 0))) = ip((asser == 0) | ((asser == 1) & (convf == 0)));
ipout((asser == 0) | ((asser >= 1) & (convf == 0))) = ip((asser == 0) | ((asser >= 1) & (convf == 0)));
if any(asser == 2)
  ipout(asser == 2) = - (flux_ip(asser == 2) - psi(asser == 2,end))  ./ L_ext(asser == 2);
  %disp([ipout,rm .* trapz(x,spr_tor .* jmoy,2),ip,- (flux_ip(asser == 2) - psi(asser == 2,end))  ./ L_ext(asser == 2)]);
end

% calcul de F
warning off
dpdpsi   = pdederive(x,ptot_tor,0,1,2,1) ./ min(-eps,psid1);
if hollow ~= 1
    dpdpsi(psid1 > -eps) = 0;
end
% estimation simple de la valeur sur l'axe : (dP/dpsi ~= 0 sur l'axe meme si dP/drho = 0)
% peu d'incidence sur l'integrale
if cronos_regul == 4
    % fit of P on order 2 polynomial in center to have deriavtive on magnetic axis
    % dpdx =  2 * dpdx_1 * x and p = p0 + dpdx_1 * x^2
    dpdx_1 = (ptot_tor(:,2) - ptot_tor(:,1)) ./ x(2) .^ 2;
    dpdpsi(:,1) = - 2 .* dpdx_1 ./ b0 .* q_not ./ rm .^ 2;
else
    dpdpsi(:,1) = dpdpsi(:,2);
end
df2dpsi  = 2 .* mu0 .* (jmoy .* ri_tor - dpdpsi)./ r2i_tor;
df2dx    = df2dpsi .* psid1;
F2       = (r0 .^ 2 .* b0 .^ 2) *ve + cumtrapz(x(end:-1:1),df2dx(:,end:-1:1),2);
F2       = min((r0 .^ 2 .* b0 .^ 2) *ve .* 2, max((r0 .^ 2 .* b0 .^ 2) *ve ./ 2 ,F2));
F        = sqrt(F2(:,end:-1:1));
warning on


% derivee temporelle
dpsidt =z0dxdt(psi,t);
dpsidt(1,:) = dpsidt(2,:) - dpsidt(2,end) + dpsidt_edge(1);
% calcul de epar
epar = - F ./ (b0 * ve).* r2i_tor .* dpsidt;

% calcul de q
if cronos_regul == 5
    % Modif DM 16/09/16: regularisation du facteur de securite au centre
    warning off
    Breg=psid1(:,krp1:krp2);
    Xreg=Breg*Areg;
    dphidpsi0=2*pi*((b0.*rm.*rm))./(Xreg*[1;x(1)]);
    q = abs(pdederive(psi,phi_tor,4,2,2,1,dphidpsi0)) ./ 2 ./ pi;
    warning on
elseif (cronos_regul == 3) || (cronos_regul == 4)
    warning off
    q = abs(pdederive(psi,phi_tor,1,2,2,1)) ./ 2 ./ pi;
    warning on
    if cronos_regul == 4
        q(:,1) = q_not;
    end
else
    warning off
    q = abs(pdederive(psi,phi_tor,0,2,2,1)) ./ 2 ./ pi;
    warning on
    % regularisation au centre 
    q(:,1:2)  = pchip(x(3:11),q(:,3:11),x(1:2));

    if cronos_regul >= 1 && cronos_regul < 3
      % debut test calcul q0 comme dans cronos
      inter = mu0 .* kx(:,1) .* Raxe(:,1) .* jmoy(:,1) ./ 2 ./ (1+ kx(:,1) .^ 2);
      warning off
      q0    = max(0.1,F(:,1) ./Raxe(:,1) ./ inter ./ 2);
      warning on
      q0(inter == 0) = q(inter == 0,1);
      if cronos_regul == 2
	    q0  = 0.5 .* (q0 + q(:,1));
      end
      q(:,1:2)  = pchip(cat(2,0,x(3:11)),cat(2,q0,q(:,3:11)),x(1:2));
      % fin test calcul q0 comme dans cronos
    end
end

% calcul de la valeur au bord
dx_tordx    = pdederive(x,x_tor,2,2,2,1);
% calcul precis pour le bord
q(:,end)  = abs(dphidx(:,end) ./ psid1(:,end) .*  dx_tordx(:,end) ./ 2 ./ pi);

% estimation precise de dptotdspi sur l'axe et de ffprim
if cronos_regul ~= 4
    d2pdx2 = pdederive(x,ptot,1,0,2,2);
    dpdpsi(:,1)   = - (d2pdx2(:,1) ./ rm) .* q(:,1) ./ 4 ./ pi ./ phi_tor(:,end); 
    df2dpsi(:,1)  = 2 .* mu0 .* (jmoy(:,1) .* ri_tor(:,1) - dpdpsi(:,1))./ r2i_tor(:,1);
end

% calcul de jeff
Fd1    = pdederive(x,F,0,2,2,1);
jeff       = - F .^ 2 ./ (mu0 .* max(eps,vpr_tor)) ./  (rm * ve) .^ 2 ./ (b0 * ve) .* ...
               (c2cd1 ./ F .* psid1  +  C2_tor ./ F .* psid2 - C2_tor .* Fd1  ./ F .^ 2 .* psid1);
if cronos_regul == 5    
    % Modif DM 16/09/16: calcul de jeff(:,1)
    jeff(:,1)  = interp1(x(2:3)',jeff(:,2:3)',x(1),'linear','extrap')';
elseif cronos_regul == 4
    vpr_tor_o_rho = vpr_tor(:,2) ./ rm ./ x(2);
    psid1_o_rho   = -b0 ./ q_not .* rm;
    C2_tor_o_rho  = vpr_tor_o_rho .* grho2r2_tor(:,1);
    jeff(:,1) = - F(:,1) .^ 2 ./ (mu0 .* vpr_tor_o_rho) ./  rm .^ 2 ./ b0 .* ...
                  (c2cd1(:,1) ./ F(:,1) .* psid1_o_rho  +  C2_tor_o_rho ./ F(:,1) .* psid2(:,1));
else
    jeff(:,1:2)  = pchip(x(3:11),jeff(:,3:11),x(1:2));
end
jres       = jeff - jni_tor;


% calcul de epar pour le premier temps
% contrainte : courant ohmique est dpsi/dt au bord
epar_1 = eta_tor(1,:) .* jres(1,:); 
epar_1(end) = epar(1,end);
%figure(21);clf;plot(x,epar(1,:),'r',x,epar_1,'b');drawnow
if all(isfinite(epar_1)) && (evolution == 0)
    epar(1,:) = epar_1;
    dpsidt(1,:) = - epar_1 ./ F(1,:) .* (b0(1) * ve) ./ r2i_tor(1,:);
end

% calcul de ej
terme  = c2cd1 .* psid1 + C2_tor .* psid2;
ejp    = - epar .* (b0 * ve) ./ ((rm * ve) .^2 .* mu0 .* F .*max(eps,vpr_tor) .* r2i_tor ) .* terme;
%drmdt  = pdederive(t,rm,2,2,1,1);
% ce terme n'est pas exprime correctement 
%ref : V. D. Pustovitov PPCF 54 (2012) 124036
%ejc    = pdederive(x,ptot,0,2,2,1) .* ((drmdt ./ rm) * ve);
drmxdt = z0dxdt(rmx,t);
if evolution == 1
     mat = cat(2,t,ones(size(t)));
     mat_inv = pinv(mat);
     for kml = 1:size(rmx,2)
        %pp = polyder(polyfit(t,rmx(:,kml),1))
        pp = (mat_inv(1,:) * rmx(:,kml))';
        drmxdt(:,kml) = pp;
    end
    %figure(138);clf;plot(x,drmxdt,'.',x,z0dxdt(rmx,t));drawnow
end

%figure(121);clf;plot(t,drmdt,'r',t,drmxdt(:,end),'b',t,bpcor(:,end),'k');drawnow

if lao_change == 1
	drmxdt       = tsplinet(rmx,drmxdt,rmx(:,end) * x);
end
ejc    = pdederive(x,ptot,0,2,2,1) .* (drmxdt ./ (rm * ve));
ejp(:,1)  = 2 .* ejp(:,2) - ejp(:,3);
ej        = ejp - ejc;

% calcul de li
bpol  = sqrt(grho2r2_tor) .* abs(psid1) ./ (rm*ve);
li     = 2 .* rm .* trapz(x,bpol .^ 2  .* vpr_tor,2) ./ ((mu0 .*ip ) .^ 2  .* r0);

% poyting 
% ref :Ohmic ???ux consumption during initial operation of the NSTX spherical torus
% J.E. Menard et al, NF 2001 41 1197
% ici on prend en compte que la composante est utile pour la consommation de flux poloidal.
poynting = rm .* trapz(x,vpr_tor .* (z0dxdt(bpol .^ 2 ./ 2 ./ mu0,t) + ej),2);


%save('loc9');	
% changement de coordonnees inverse flux toroidal vers Lao
if lao_change == 1
	
	psi         = tsplinet(rmx(:,end) * x,psi,rmx);
	dpsidt      = tsplinet(rmx(:,end) * x,dpsidt,rmx);
	qnew        = tsplinet(rmx(:,end) * x,q,rmx);
	q           = qnew .* (qnew > 0) + q .* (qnew <= 0);
	jmoy        = tsplinet(rmx(:,end) * x,jmoy,rmx);
	jgs         = tsplinet(rmx(:,end) * x,jgs,rmx);
	jres        = tsplinet(rmx(:,end) * x,jres,rmx);
	jeff        = tsplinet(rmx(:,end) * x,jeff,rmx);
	epar        = tsplinet(rmx(:,end) * x,epar,rmx);
	F           = tsplinet(rmx(:,end) * x,F,rmx);
	bpol        = tsplinet(rmx(:,end) * x,bpol,rmx);
	ej          = tsplinet(rmx(:,end) * x,ej,rmx);
	df2dpsi	    = tsplinet(rmx(:,end) * x,df2dpsi,rmx);
	dpdpsi      = tsplinet(rmx(:,end) * x,dpdpsi,rmx);	
    

end

%save('loc10');	
%keyboard

% --------------------------------------------------------------
function s=zreshape(e,nx,ny,nz)

% continuite au centre pour les NaN
if ~isfinite(e(1))
	e(1) = (61/46) .* e(2) - (9/23) .* e(3) + (3/46) .* e(4);
end

% suppression des NaN
ind = find(~ isfinite(e));
if ~isempty(ind)
	e(ind) = zeros(1,length(ind));
end

% mise en forme du vecteur
s = reshape(e,nx,ny,nz);






