% ZEQUI_HELENA calcul de l'equilibre (marche pour toutes les machines)
%-------------------------------------------------------------------------
% fichier zequi_helena.m ->  zequi_helena
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule l'equilibre a un temps donne. Elle remplie la structure "datak.equi".
% La normalisation des flux utilise la valeur de psi sur la derniere surface magnetique, sauf si ip est donne
% en entree.
%
% syntaxe  :
%
%     [equi,mhd_cd,straightline] = zequi_helena(cons,geo,equi,prof,phys,ip,x,change);
%
% entrees :
%
%     cons    =  structure des consignes de calcul de la fonction (pas de parametre pour cette fonction)
%     geo     =  datak.geo
%     equi    =  datak.equi
%     prof    =  datak.prof
%     phys    =  param.phys
%     ip      =  courant plasma pour la normalisation (A)
%     x       =  param.gene.x
%     change  =  flag indiquant un changement de geometrie
%     exact   =  flag pour le reclacul rapide et exact de helena en utilisant P' et FF'
%
% sorties :
%
%     equi          = structure equi de datak (cf. zinit.m)
%     mhd_cd        = modulation par la mhd de l'efficacite de generation de courant
%     straightline  = structure de donnees pour les coordonnees tel que les lignes de champs magnetiques sont des geodesiques
%
%  a l'initialisation
%     equi.valeur    = structure des parametres
%     equi.nbrho   = nombre de point en rho de la grille RZ
%     equi.nbtheta = nombre de point en theta de la grille RZ
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 31/05/2005.
%
%
% liste des modifications :
%
%  * 02/03/2001 -> Z des surface compatible avec la separatrice
%  * 07/06/2001 -> mode 'ecrit' dans helena si pb au 1er appel
%  * 15/06/2001 -> remplacement de interp1 par tsplinet
%  * 10/07/2001 -> ajout d'une correction sur les courants creux
%  * 10/07/2001 -> ajout du calcul de rhog
%  * 12/07/2001 -> ajout du calcul de raxe
%  * 23/07/2001 -> correction de la valeur au centre de grho2
%  * 20/09/2001 -> ajout de equi.shear
%  * 28/09/2001 -> modification de rhog (mise en concordance avce sa definition)
%  * 09/10/2001 -> suppression des zregular pas vraiment justifies
%  * 08/11/2001 -> correction de rho dans le cas ou helena plante
%  * 28/05/2002 -> probleme lorsque le profil de courant est negatif au centre
%                  et que helena n'a pas converge (min(ind) -1 = 0)
%  * 24/07/2002 -> ajout de la securite sur la dsmf
%  * 27/08/2002 -> modification de la securite sur la dsmf
%  * 04/06/2003 -> version pour utiliser amix si fail ~=0
%  * 04/09/2003 -> correction d'un bug si geo.R et geo.Z contiennent des NaNs alors que geo.mode <> 2
%  * 17/09/2003 -> ajout de la possibilite d'appeler le module deltap
%  * 13/10/2003 -> ajout de psiRZ
%  * 03/11/2003 -> ajout des criteres mhd analytiques
%  * 07/01/2004 -> modifcation de la formule corrective du gradient de pression au centre pour les DDS
%  * 29/03/2004 -> correction passage parametres en mode fail ~0 : li a la mauvaise place.
%  * 31/03/2004 -> ajout du calcul de c2c.
%  * 03/09/2004 -> retour a l'ancienne valeur de c2c = vpr * grho2r2.
%  * 15/11/2004 -> modification du calcul de la valeur de c2c au bord
%  * 22/11/2004 -> suppression de la modification du 15/11
%  * 25/11/2004 -> utilisation du calcul inverse pour c2c
%  * 26/11/2004 -> ajout des parametres XR1 et SIG1
%  * 06/12/2004 -> changement de la valeur par defaut de xr1
%  * 06/12/2004 -> changement de la valeur par defaut de xr1
%  * 20/01/2005 -> changement du clacul de c2c
%  * 09/02/2005 -> ajout des sortie pour le reclacul rapide et precis de helena
%  * 01/03/2005 -> ajout des sortie pour les coordonnees de champs droits
%  * 01/03/2005 -> option d'appel avec reclacul exact
%  * 18/05/2005 -> modification de la securite sur la derive de ptot au centre
%  * 25/05/2005 -> suppression du try/catch : probleme en matlab 7
%  * 26/05/2005 -> correction du troisieme appel a helena en mode fail = 1
%  * 31/05/2005 -> securite pour les temps a convergence difficile : nouvelle valeur de ptot au centre si ifail =1
%  * 06/04/2006 -> modification dpdpsi + ajout nouvelles entrees
%--------------------------------------------------------------
%
function [equi,mhd_cd,memoire,straightline] = ...
    zequi_helena(memoire,cons,pfcur,geo,equi,prof,phys,ipin,x,change,exact)

% mode initialisation
% pour economise la memoire
nbmode  = 1026;
% fonction auto declarante
if nargin <=1

    valeur.xr1        = 1.1;
    type.xr1          = 'real';
    borne.xr1         = [0,1.1];
    defaut.xr1        = 1.1;
    info.xr1          = 'Gaussian center to perform an HELENA input finer grid; si >1 :no effect';

    valeur.sig1        = 0.1;
    type.sig1          = 'real';
    borne.sig1         = [1e-2,0.3];
    defaut.sig1        = 0.1;
    info.sig1          = 'Gaussian widthto perform an HELENA input finer grid; si >1 :no effect';

    valeur.dponoff        = 0;
    type.dponoff          = 'integer';
    borne.dponoff         = {0,1};
    defaut.dponoff        = 0;
    info.dponoff          = 'delta_prime module call if = 1';

    valeur.mmax           = 3;
    type.mmax             = 'integer';
    borne.mmax            = [1,7];
    defaut.mmax           = 3;
    info.mmax             = 'maximum poloidal number for the deltap module ';

    valeur.nmax           = 2;
    type.nmax             = 'integer';
    borne.nmax            = [1,7];
    defaut.nmax           = 3;
    info.nmax             = 'maximum toroidal number for the deltap module';

    valeur.ng            = 50;
    type.ng              = 'integer';
    borne.ng             = [10,200];
    defaut.ng            = 50;
    info.ng              = 'number of finite element used [50]';

    valeur.np            = 1;
    type.np              = 'integer';
    borne.np             = {0,1};
    defaut.np            = 1;
    info.np              = 'finite element evaluation [1]';

    valeur.sigs          = 0.1;
    type.sigs            = 'integer';
    borne.sigs           = [0.01,1];
    defaut.sigs          = 0.1;
    info.sigs            = 'concentration of finite element around the  resonance [0.1]';

    valeur.mhd_bal        = 0.5;           % facteur de modulation de l'efficacite de generation de courant par le ballooning
    type.mhd_bal          = 'float';
    borne.mhd_bal         = [0,1];
    defaut.mhd_bal        = 0.5;
    info.mhd_bal          = 'modulation factor of the current drive efficiency for the ballooning mode';

    valeur.errcur        = 1e-4;           % facteur de modulation de l'efficacite de generation de courant par le ballooning
    type.errcur          = 'float';
    borne.errcur         = [1e-3,1e-5];
    defaut.errcur        = 1e-4;
    info.errcur          = 'error allowed on the current profile';

    valeur.errit        = 1e-7;           % facteur de modulation de l'efficacite de generation de courant par le ballooning
    type.errit          = 'float';
    borne.errit         = [1e-6,1e-8];
    defaut.errit        = 1e-6;
    info.errit          = 'internal parameter for the iteration';

    valeur.modep        = 'Ptot';           % facteur de modulation de l'efficacite de generation de courant par le ballooning
    type.modep          = 'list';
    borne.modep         = {'Ptot','Pth','Pth norm'};
    defaut.modep        = 'Ptot';
    info.modep          = 'Pression used in HELNA (Pth norm = Pth * Wtot / Wth)';

    interface.ts    = '';
    interface.jet   = '';

    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
    sortie.interface=interface;

    sortie.description = 'Equilibrium computation';   % description (une ligne) de la fonction
    sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
    sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
    sortie.resume ='';                           % nom de la fonction qui ecrit le resume du parametrage


    sortie.nbrho   = 101;
    sortie.nbtheta = 65;
    sortie.nbmode  = nbmode;

    equi = sortie;
    return

end


% compatibilite
if ~isfield(cons,'xr1')
    xrsig = [99,99];
elseif cons.xr1 > 1
    xrsig = [99,99];
else
    xrsig = [cons.xr1,cons.sig1];
end
if ~isfield(cons,'errcur')
    errcur = 10-3;
else
    errcur = cons.errcur;
end
if ~isfield(cons,'errit')
    errit = 10-7;
else
    errit = cons.errit;
end

% selon le mode d'appel
if nargin < 11
    exact = 0;
elseif isempty(exact)
    exact = 0;
elseif exact == 1
    % test de compatibilite avec le jeu de donnees
    if ~isfield(equi,'df2RZ') | ~isfield(equi,'dprRZ') |  ~isfield(equi,'frmode') |  ~isfield(equi,'bnorme')
        warning('ZEQUI_HELENA : the exact mode can not be use with this old data set');
        exact = 0;
    elseif any(~isfinite(double(equi.df2RZ))) | any(~isfinite(double(equi.dprRZ)))| any(~isfinite(double(equi.frmode)))| any(~isfinite(equi.bnorme))
        warning('ZEQUI_HELENA : the exact mode can not be use with this old data set');
        exact = 0;
    elseif all(double(equi.df2RZ) == 0) | all(double(equi.dprRZ) == 0)  | all(equi.bnorme ==0)
        warning('ZEQUI_HELENA : the exact mode can not be use with this old data set');
        exact = 0;
    end
end

% flag de retry si echec premier appel
retry =0;
langue = getappdata(0,'langue_cronos');
% mise en forme des parametres
% psi est une corordonnee et doit etre corrigee des fluctuations numeriques
% correction du centre pourquoi ?
%psipin = zregular(x,prof.psi);
psipin = prof.psi;
% correction monotonie et signe
absdpsidx  = abs(pdederive(x,psipin,0,2,2,1));
ind        = find(absdpsidx == 0);
if ~isempty(ind)
    mindpsidx      = min(absdpsidx(absdpsidx > 0));
    absdpsidx(ind) = mindpsidx;
end
psipin = cumtrapz(x,absdpsidx);

%  % schema numerique conservateur
%  psipin   = - (prof.psi - prof.psi(1));
%  absdpsi  = abs(diff(psipin));
%  ind        = find(absdpsi == 0);
%  if ~isempty(ind)
%     mindpsi      = min(absdpsi(absdpsi > 0));
%     absdpsi(ind) = mindpsi;
%  else
%     ind = 1;
%  end
%  psipin = cat(2,0,cumsum(absdpsi));

% changement d'echelle (normalisation)
psipin   =  2 .* pi .* (psipin -psipin(1)) ;
if ~isfield(cons,'modep')
    cons.modep = 'Ptot';
end
% changement d'unite de ptot
if strcmp(cons.modep,'Ptot')
    ptot     = prof.ptot ./ 1e3 ./ phys.e;  % Pa -> keV.*m^-3
end
if strcmp(cons.modep,'Pth')
    ptot     = (prof.pe+prof.pion) ./ 1e3 ./ phys.e;  % Pa -> keV.*m^-3
end
if strcmp(cons.modep,'Pth norm')
    if all(isfinite(equi.vpr))
        we      = (3/2) .* zintvol(prof.pe,x,equi.vpr,equi.rhomax);
        wi      = (3/2) .* zintvol(prof.pion,x,equi.vpr,equi.rhomax);
        wd     = (3/2) .* zintvol(prof.ptot,x,equi.vpr,equi.rhomax);
        wth      = we+wi;
        ptot     = (prof.pe+prof.pion) ./ 1e3 ./ phys.e * wd / wth;  % Pa -> keV.*m^-3
    else
        we       = trapz(x,prof.pe .* x,2);
        wi       = trapz(x,prof.pion .* x,2);
        wd       = trapz(x,prof.pion .* x,2);
        wth      = we+wi;
        ptot     = (prof.pe+prof.pion) ./ 1e3 ./ phys.e * wd / wth;  % Pa -> keV.*m^-3;
    end
end


% condition de convergence de helena
% 1 - monotone
dpdpsi   = pdederive(psipin,ptot,0,2,2,1);
if any(dpdpsi > 0)
    ind            = find(dpdpsi > 0);
    indm           = find(dpdpsi <0);
    indin          = min(indm);
    dpdpsiok       = dpdpsi(indm);
    dpdpsimin      = dpdpsiok(min(find(dpdpsi(indm) == max(dpdpsi(indm)))));
    dpdpsi(ind)    = dpdpsimin .* ones(1,length(ind));
    if indin >1
        dpdpsi(1:indin) = dpdpsimin .* ones(1,indin);
    end
    % petit ajout pour ameliorer la convergence
    ptotnew  = cumtrapz(psipin,dpdpsi);
    ptot     = ptotnew - ptotnew(end) + ptot(end) .* (ptot(end)>0);
end
% 2 - borne
dpdpsi    = pdederive(psipin,ptot,0,2,2,1);
mu        = psipin ./ psipin(end);
indmu       = min(find(mu >0.01));
dpdpsimin = mean(-abs(dpdpsi(ind:end)));
%  probleme de signe ?
%  JU -> oui
ind       = find((dpdpsi > dpdpsimin) & (mu < mu(indmu)));
if ~isempty(ind)
    % correction effet d'arrondi d'indice pour DDS
    if length(ind)  >3
        ind(end) = [];
    end
    dpdpsi(ind) = dpdpsimin .* ones(1,length(ind));
    ptotnew  = cumtrapz(psipin,dpdpsi);
    ptot     = ptotnew - ptotnew(end) + ptot(end) .* (ptot(end)>0);
end

% calcul de ptotcor
dpdpsicor = dpdpsi;
dpdpsi_c  = - abs((max(ptot) - min(ptot)) ./ (max(psipin) - min(psipin)));
dpdpsicor(1:5) = 0.1 .* dpdpsi_c ;
ptotcor  = cumtrapz(psipin,dpdpsicor);
ptotcor     = ptotcor - ptotcor(end) + ptot(end) .* (ptot(end)>0);


% calcul de ptot corrig? au centre (avec interpos)
%sigma    = [length(psipin) ones(size(psipin(2:end)))];
%ptot     = interpos(psipin,ptot_ori,psipin,-1,[1 1],[dpdpsimin 1e32],sigma);

% jmoy pour l'equilibre
jmoy = prof.jmoy;

% verification de l'existence de z0
if isfinite(geo.z0)
    z0 = geo.z0;
else
    z0 = 0;
end

% securite sur la dsmf
if geo.mode == 2       % correction 04/09/2003 :  on evite que ca plante si geo.R et geo.Z ne sont pas remplis alors que geo.mode <> 2
    
    % re ordered LCFS points
    R_ = geo.R;
    Z_ = geo.Z;
    r_ = [];
    z_ = [];
    % initial point
    ind_ = find(R_ == max(R_));
    r_ = R_(ind_(1));
    z_ = Z_(ind_(1));
    R_(ind_(1)) = [];
    Z_(ind_(1)) = [];
    while ~isempty(R_)
      dd_ = sqrt((R_-r_(end)) .^ 2 + (Z_-z_(end)) .^ 2);
      ind_ = find(dd_ == min(dd_));
      if dd_(ind_(1) > 0)
        r_ = cat(2,r_,R_(ind_(1)));
        z_ = cat(2,z_,Z_(ind_(1)));
      end
      R_(ind_) = [];
      Z_(ind_) = [];
    end
    if z_(2) < z_(1)
          r_ = r_(end:-1:1);
          z_ = z_(end:-1:1);
    end
    geo.R = r_;
    geo.Z = z_;
 
    
    geo.R   = double(geo.R(1:(end-1)));
    %geo.Z   = double(geo.Z(1:(end-1))) - z0;
    % modification du 06/12/2004
    geo.Z   = double(geo.Z(1:(end-1)));


    % suprression des points double
    indbad = find((diff(geo.R) == 0) & (diff(geo.Z) == 0));
    while ~isempty(indbad)
        geo.R(indbad) = [];
        geo.Z(indbad) = [];
        indbad = find((diff(geo.R) == 0) & (diff(geo.Z) == 0));
    end
    ind1 = min(find(geo.R == min(geo.R)));
    ind2 = min(find(geo.R == max(geo.R)));

    r0_  = 0.5 .* (geo.R(ind1) + geo.R(ind2));
    z0_  = geo.Z(ind2);
    % in case ind2>1, otherwise theta is not monotonic
    if ind2 ~= 1
        geo.R = [geo.R(ind2:end) geo.R(ind2-1:-1:1)];
        geo.Z = [geo.Z(ind2:end) geo.Z(ind2-1:-1:1)];
    end
    %aa   = unwrap(angle((geo.R - r0) + sqrt(-1) .* (geo.Z -z0)));
    aa    = atan2(geo.Z -z0_,geo.R - r0_);
    aa    = unwrap(aa .* (aa >=0) + (aa + 2*pi) .* (aa< 0));
    if sum(diff(aa)  < 0) > sum(diff(aa)  > 0)
        geo.R = geo.R(end:-1:1);
        geo.Z = geo.Z(end:-1:1);
        aa = aa(end:-1:1);
    end
    if any(diff(aa) < 0)
        indbad = find(diff(aa)<=0) + 1;
        while(~isempty(indbad)) && (length(geo.R) > 2)
            geo.R(indbad) = [];
            geo.Z(indbad) = [];
            aa(indbad)    = [];
            indbad = find(diff(aa)<=0) + 1;
        end
    end
    if length(geo.R) < 3
        geo.mode = 1;
    end
    %    % correction de la phase
    %    ind0   = min(find(aa >= 0));
    %   ind    = [ind0:(length(geo.R)-1),1:ind0];
    %   geo.R  = geo.R(ind);
    %   geo.Z  = geo.Z(ind);
    %
    %end
    [aa,ind_aa] = sort(aa);
    geo.R  = geo.R(ind_aa);
    % modification du 06/12/2004
    geo.Z  = geo.Z(ind_aa) - z0_;
    geo.R(end+1) = geo.R(1);
    geo.Z(end+1) = geo.Z(1);
end


% initialisation de l'equilibre et normalisation
% indique au code d'equilibre qu'il faut utilise psi pour la normalisation

ip = -1;
if ~isempty(ipin)
    if isfinite(ipin)
        ip =abs(ipin);
    end
end
%fprintf('ip = %g\n',ip);
if isappdata(0,'equilibre_helena_init')
    mem.init    = getappdata(0,'equilibre_helena_init');
    mem.b       = getappdata(0,'equilibre_helena_b');
    mem.df2     = getappdata(0,'equilibre_helena_df2');
    mem.kkbig   = getappdata(0,'equilibre_helena_kkbig');
    mem.ptotold = getappdata(0,'equilibre_helena_ptotold');
    mem.corj    = getappdata(0,'equilibre_helena_corj');
else
    mem.init = 0;
end
if isempty(mem.init)
    mem.init = 0;
end

% si change corj -> 0
if change ~=0
    mem.corj =0;
elseif ~isfield(mem,'corj')
    mem.corj =0;
end
% si probleme lors de l'appel precedent
if jmoy(1)./max(jmoy)*100 < 5
    mem.corj = 1;
end
if mem.corj == 1
    % correction du courant au bord
    if any(jmoy <0)
        ind    = find(jmoy <0);
        %      min(ind) - 1 doit etre > 1
        if x(max(min(ind)-1,1)) > 0.95
            ind = max(find(x <=0.95)):length(x);
        end
        if min(ind) > 1
            jmoymem   = jmoy;
            jcor      = jmoy(min(ind)-1);
            xcor      = x(min(ind)-1);
            jmoy(ind) = jcor .* (1-x(ind).^2)./(1-xcor^2);
            % renormalisation de jmoy pour conserver ip
            rapjmoy   = trapz(x,jmoymem .* x) ./ trapz(x,jmoy .* x);
            jmoy      = jmoy .* rapjmoy;
        else
            disp('Jmoy can not be modified ...')
        end
    end
    % correction des profils creux
    % modification du 10/07/2001
    % correction du point x central
    creux =0.5;
    if (jmoy(1).*mean(jmoy)) <= (creux .* mean(jmoy) .^ 2)
        jmoymem   = jmoy;
        indx     = max(find(x <= 0.05));
        indm     = min(find(abs(jmoy) == max(abs(jmoy))));
        indj     = max(find (((jmoy.*mean(jmoy)) <= (creux .* mean(jmoy) .^ 2)) & (x <= x(indm))))+1;
        if isempty(indj)
            indj=1;
        end
        if indj > length(x)
            inj= length(x);
        end
        ind      = max(indx,indj);
        j0       = creux .* mean(jmoy);
        alpha    = (jmoy(ind) - j0) ./ x(ind) .^ 2;
        jmoy(1:ind) = j0 + alpha .* x(1:ind) .^ 2;
        % renormalisation de jmoy pour conserver ip
        rapjmoy   = trapz(x,jmoymem .* x) ./ trapz(x,jmoy .* x);
        jmoy      = jmoy .* rapjmoy;
    end
    % fin modification du 10/07/2001
end


% selon le mode d'appel
if (ip < 0) & (equi.fail == 0) & (mem.init == 1) & ~isempty(mem.ptotold) & ...
        ~isempty(mem.b) & ~isempty(mem.df2) & ~isempty(mem.kkbig) & (change == 0) & (exact == 0)

    % l'equilibre precedent est valide
    init(1)    = mem.init;
    b       = mem.b;
    df2     = mem.df2;
    kkbig   = mem.kkbig;
    ptotold = mem.ptotold;

    % recalcul le b
    dpnewdpsi     = pdederive(equi.psi,ptot,2,2,2,1);
    dpolddpsi     = pdederive(equi.psi,ptotold,2,2,2,1);
    bnew          = b .*  abs(mean(dpnewdpsi(1:3)) ./ mean(dpolddpsi(1:3)));
    disp([b,bnew]);
    if ~isfinite(bnew) | (bnew == 0)
        init(1) = 0;
    elseif abs(b-bnew) > abs(0.1 .* b)
        init(1) = 0;
    elseif abs(bnew) > 10
        init(1) = 0;
    else
        disp('fast mode of HELENA')
        b=bnew;
    end

else
    % l'equilibre precedent n'est pas valide
    init(1)    = 0 ;
    b       = 0;
    df2     = 0 .* ones(size(equi.df2RZ));
    kkbig   = [];
    ptotold = [];
end

% appel de la fonction de calcul de l'equilibre
disp('calling HELENA (magnetic equilibrium) ->')
%disp('from zequi_helena ->')
%
% probleme du 29 octobre 2001
%
%
% geo.Z(end)-z0 proche de zero ce qui dans le calcul de theta effectue par helena
% donne 0 a la place de 2*pi
% correction effectue dans helmex77
% V. Basiuk
%
if exact == 0
    init(2) = errcur;
    init(3) = errit;
    [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
        r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail, ...
        nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
        helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig,0,xrsig);
				% (19 args on RHS)

else

    % l'equilibre precedent n'est pas valide
    init(1)    = 0 ;
    b       = 0;
    df2     = 0 .* ones(size(equi.df2RZ));
    kkbig   = [];
    ptotold = [];
    init(2) = errcur;
    init(3) = errit;

    [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
        r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail, ...
        nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
        helmex77(-2 .* pi .* (equi.psi - equi.psi(1)),equi.ptot ./ 1e3 ./ phys.e,equi.jmoy,geo.r0,geo.b0,equi.ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z, ...
        init,equi.bnorme,double(equi.df2RZ),kkbig,0,xrsig,double(equi.dprRZ),double(equi.frmode));


end
%disp('helena returned')

%save last_helena
if ifail ~=0

    % petite correction pour choc iter delicat ( a cause de la forme de jmoy(psi))
    jmoy(1) = jmoy(indmu);

    %save equi_helena_ifail1
    if exact == 0


        % deuxieme appel systematique avec amix = 0.9 et nbiter = 1000
        % l'equilibre precedent n'est pas valide
        init(1)  = 0 ;
        b       = 0;
        df2     = 0 .* ones(size(equi.df2RZ));
        kkbig = [];
        init(2) = errcur;
        init(3) = errit;

        [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
            r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail, ...
            nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
            helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig,2,xrsig);
    else

        % l'equilibre precedent n'est pas valide
        init(1)    = 0 ;
        b       = 0;
        df2     = 0 .* ones(size(equi.df2RZ));
        kkbig   = [];
        ptotold = [];
        init(2) = errcur;
        init(3) = errit;

        [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
            r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail, ...
            nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
            helmex77(-2 .* pi .* (equi.psi - equi.psi(1)),equi.ptot ./ 1e3 ./ phys.e,equi.jmoy,geo.r0,geo.b0,equi.ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z, ...
            init,equi.bnorme,double(equi.df2RZ),kkbig,2,xrsig,double(equi.dprRZ),double(equi.frmode));

    end

    if ifail ~=0
        %save equi_helena_ifail2
    end


end

% si appel mode rapide ou jmoy <0 au bord
if ((init(1) == 1) | (jmoy(end) <0)) & (ifail ~=0) & (exact == 0)
    if init(1) == 1
        disp('HELENA did not converged in fast mode')
        disp('passing in standard mode');
    else
        disp('HELENA did not converged  ')
        disp('standard mode with current profile modification ');
    end
    retry =1;
    % correction du courant au bord
    if any(jmoy <0)
        % flag pour l'appel suivant
        mem.corj =1;
        % correction proprement dite
        ind    = find(jmoy <0);
        %  correction si min(ind) = 1
        if x(max(min(ind)-1,1)) > 0.95
            ind = max(find(x <=0.95)):length(x);
        end
        if min(ind) > 1
            jmoymem   = jmoy;
            jcor      = jmoy(min(ind)-1);
            xcor      = x(min(ind)-1);
            jmoy(ind) = jcor .* (1-x(ind).^2)./(1-xcor^2);
            % renormalisation de jmoy pour conserver ip
            rapjmoy   = trapz(x,jmoymem .* x) ./ trapz(x,jmoy .* x);
            jmoy      = jmoy .* rapjmoy;
        else
            if strcmp(langue,'anglais')
                disp('Jmoy can not be modified ...')
            else
                disp('Jmoy ne peut pas etre corrige ...')
            end
        end
    end
    % correction des profils creux
    % modification du 10/07/2001
    % correction du point x central
    creux =0.5;
    if (jmoy(1).*mean(jmoy)) <= (creux .* mean(jmoy) .^ 2)
        % flag pour l'appel suivant
        mem.corj =1;
        % correction proprement dite
        jmoymem  = jmoy;
        indx     = max(find(x <= 0.05));
        indm     = min(find(abs(jmoy) == max(abs(jmoy))));
        indj     = max(find (((jmoy.*mean(jmoy)) <= (creux .* mean(jmoy) .^ 2)) & (x <= x(indm))))+1;
        if isempty(indj)
            indj  = 1;
        end
        if indj > length(x)
            indj  = length(x);
        end
        ind      = max(indx,indj);
        j0       = creux .* mean(jmoy);
        alpha    = (jmoy(ind) - j0) ./ x(ind) .^ 2;
        jmoy(1:ind) = j0 + alpha .* x(1:ind) .^ 2;
        % renormalisation de jmoy pour conserver ip
        rapjmoy   = trapz(x,jmoymem .* x) ./ trapz(x,jmoy .* x);
        jmoy      = jmoy .* rapjmoy;
    end
    % fin modification du 10/07/2001


    % utilisation de ptotcor
    ptot = ptotcor;

    % l'equilibre precedent n'est pas valide
    init(1)  = 0 ;
    b       = 0;
    df2     = 0 .* ones(size(equi.df2RZ));
    kkbig = [];
    ptotold = [];
    init(2) = errcur;
    init(3) = errit;

    [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
        r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail, ...
        nchi,cpsurf,radius,gem11,gem12,gem33,rspot,zspot,brspot,bzspot]= ...
        helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig,2,xrsig);
    % [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
    % r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,dpr,kkbig,li,dimerc,drmerc,balcrit,frmode,ifail]= ...
    % helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig,2,xrsig);

    if ifail ~=0
        %save equi_helena_ifail3
    end
end

% Appel depuis Matlab :
% [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp, ...
%    rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,pout, ...
%    r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,kkbig,li,ifail]= ...
%    helmex77(psipin,ptot,jin,R0,B0,Ip,lsc,a,elong,triangh,triangb,Rext,Zext,init,b,df2,kkbig);
%
% Inputs :
% psipin : flux poloidal (V*s, doit etre nul au centre) [Xpts,1]
% ptot : ne Te + ni Ti (keV * m-3) [Xpts,1]
% jin : densite de courant (A/m2) [Xpts,1]
% R0 : grand rayon (m) [1,1]
% B0 : champ magnetique (T) [1,1]
% Ip : courant plasma (A) [1,1], Ip < 0 -> utilisation de psip pour le calcul du courant
% ls    % : type d'input pour la derniere surface magnetique fermee [1,1]
%        0 = plasma symetrique, donnee de a, elong, triang
%        1 = plasma asym???rique, donnee de a, elong, triangh, triangb
%        2 = plasma asym???rique, donnee de la LS    % en (R,Z)
% a  : petit rayon (m) [1,1]
% elong : ellipticite de la derniere surface magnetique fermee (1 pour TS)
% triangh : triangularite haute de la derniere surface magnetique fermee (0 pour TS)
% triangb : triangularite basse de la derniere surface magnetique fermee (0 pour TS)
% Rext  : grand rayon de la derniere surface de flux [250,1]
% Zext  : altitude de la derniere surface de flux [250,1]
% init(1) : si =1, alors on part de l'equilibre initial specifie par df2,b,kkbig [1,1]
%        si =0, alors on part d'un equilibre initial quelconque
% b : parametre b de l'equilibre initial (-) [1,1]
% df2 : vecteur df2 de l'equilibre initial (-) [101,1]
% kkbig : matrice kkbig deja calculee (-) [KKLDA,4*MAXNODE]
%
%
% Outputs :
% <> designe la moyenne de la variable sur une surface de flux poloidal
% C2     : Vp*<gradient(rho)^2/R^2 > [Xptsout,1]
% C3     : Vp*< 1/R^2 > [Xptsout,1]
% rho2m  : <|gradient(rho)|^2> [Xptsout,1]
% rho    : coordonee de flux toroidal (racine du flux toroidal) (m) [Xptsout,1]
% [R,Z]  : coordonnee geometrique associe a chaque petit rayon [100,Xptsout]
% Vp     : V' (dV/drho a rho(i)-m**2) [Xptsout,1]
% drhoor : <|grad(rho)|/R> [Xptsout,1]
% qpsi   : facteur de securite [Xptsout,1]
% psip   : flux poloidal (V*s) [Xptsout,1]
% ftra   : fraction de piegees [Xptsout,1]
% Fdia   : fonction diamagnetique [Xptsout,1]
% psiT   : flux toroidal (V*s) [Xptsout,1]
% B2     : <B^2> (T^2) [Xptsout,1]
% invB2  : moyenne de l'inverse de B^2 sur une surface de flux (T^-2) [Xptsout,1]
% Sp     : dS/drho (<tau>) [Xptsout,1]
% rav    : <R> (m) [Xptsout,1]
% oor    : <1/R> (m^-1) [Xptsout,1]
% drhoav : <|grad(rho)|> [Xptsout,1]
% xshift : decentrement (m) [Xptsout,1]
% xrad   : petit rayon geometrique (m) [Xptsout,1]
% xell   : ellipticite (-) [Xptsout,1]
% xtriapos : triangularite superieure (Z>0) (-) [Xptsout,1]
% xtrianeg : triangularite inferieure (Z<0) (-) [Xptsout,1]
% jout   : densite de courant recalculee d'apres l'equilibre trouve (A/m2) [Xptsout,1]
% ipout  : Ip = Ipout*max(xrad)*B0/mu0
% pout   : pression totale recalculee d'apres l'equilibre trouve (keV * m-3) [Xptsout,1]
% r2 : <R^2> (m^2) [Xptsout,1]
% oor2 : <1/R2> (m^-2) [Xptsout,1]
% r2tau2 : <R^2/tau^2> (m^-6) [Xptsout,1]
% drho2ob2 : < |grad(rho)|^2/B^2 > (m^2.T^-2) [Xptsout,1]
% r3otau3 : <R^3/tau^3> (m^-3) [Xptsout,1]
% r3otau : <R^3/tau> (m) [Xptsout,1]
% BPR     : Coordonnee R du champ magnetique sur un maillage [R Z]
% BPZ     : Coordonnee Z du champ magnetique sur un maillage [R Z]
% b       : parametre b de l'equilibre final (-) [1,1]
% df2     : vecteur df2 de l'equilibre initial (-) [101,1]
% kkbig   : matrice kkbig deja calculee (-) [KKLDA,4*MAXNODE]
% li      : inductance interne
% DIMERC : Ideal Mercier criterion
% DRMERC : resistive Mercier criterion
% BALCRIT : ballooning stability criterion (stable if >1, marginally stable if =1, unstable if <1, except if s<0 : balcrit=0 (always stable with negative magnetic shear))
% ifail  : code d'erreur (0 si pas d'erreur) [1,1]
if any(~isfinite(psipout))
    ifail = 118218;
end

equi.fail           = ifail;
equi.conv           = 1;
equi.bnorme         = b;
equi.psi0           = psipout(1) ./ 2 ./ pi;

% calcul de rhomax
if rho(1) ~= 0
      if rho(1) > eps
        disp('ZEQUI_HELENA: Problem in HELENA rho output. rho(1) above eps');
      end
      rho(1)=0;
end
equi.rhomax         = rho(end);
rhoin               = equi.rhomax .* x';

% calcul des grandeurs derivees
if ifail == 0
    Vpb         =  Vp;
    Vpb(1)      =  eps;
    grho2r2     =  C2 ./ Vpb;
    grho2r2(1)  = 0;
    grho2r2     =  zcentre(grho2r2);   % probleme
    r2i         =  C3 ./ Vpb;
    r2i         =  zcentre(r2i);
    oor(1)=0;
    oor         =  zcentre(oor);
    drhoor(1)=0;
    drhoor      =  zcentre(drhoor);
    drhoav(1)=0;
    drhoav      =  zcentre(drhoav);
    rav         =  zcentre(rav);
    r2          =  zcentre(r2);
    r2tau2      =  zcentre(r2tau2);
    drho2ob2    =  zcentre(drho2ob2);
    r3otau3     =  zcentre(r3otau3);
    r3otau      =  zcentre(r3otau);

    % reechantillonage des sorties
    % debut modification du 15/06/2001
    % code original :
    % equi.q              = interp1(rho,qpsi,rhoin,'spline')';
    % equi.psi            = - interp1(rho,psipout,rhoin,'spline')' ./ 2 ./ pi + prof.psi(1) ;
    % equi.phi            = interp1(rho,psiT,rhoin,'spline')' ;
    % equi.ftrap          = interp1(rho,ftra,rhoin,'spline')';
    % equi.F              = interp1(rho,Fdia,rhoin,'spline')';
    % equi.b2i            = interp1(rho,invB2,rhoin,'spline')';
    % equi.b2             = interp1(rho,B2,rhoin,'spline')';
    % equi.vpr            = interp1(rho,Vp,rhoin,'spline')';
    % equi.spr            = interp1(rho,Sp,rhoin,'spline')';
    % equi.d              = interp1(rho,xshift,rhoin,'spline')';
    % equi.ind            = zeros(size(x));
    % equi.trh            = interp1(rho,xtriapos,rhoin,'spline')';
    % equi.trl            = interp1(rho,xtrianeg,rhoin,'spline')';
    % equi.e              = interp1(rho,xell,rhoin,'spline')';
    % equi.a              = interp1(rho,xrad,rhoin,'spline')';
    % equi.grho2r2        = interp1(rho,grho2r2,rhoin,'spline')';
    % equi.grhor          = interp1(rho,drhoor,rhoin,'spline')';
    % equi.ri             = interp1(rho,oor,rhoin,'spline')';
    % equi.grho2          = interp1(rho,rho2m,rhoin,'spline')';
    % equi.grho           = interp1(rho,drhoav,rhoin,'spline')';
    % equi.rmoy           = interp1(rho,rav,rhoin,'spline')';
    % equi.jmoy           = interp1(rho,jout,rhoin,'spline')';
    % equi.ptot           = interp1(rho,pout .* 1e3 .* phys.e,rhoin,'spline')';
    % equi.r2             = interp1(rho,r2,rhoin,'spline')';
    % equi.r2i            = interp1(rho,r2i,rhoin,'spline')';
    % equi.r2tau2         = interp1(rho,r2tau2,rhoin,'spline')';
    % equi.grho2b2        = interp1(rho,drho2ob2,rhoin,'spline')';
    % equi.r3tau3         = interp1(rho,r3otau3,rhoin,'spline')';
    % equi.r3tau          = interp1(rho,r3otau,rhoin,'spline')';

    % debut du code modifie
    if all(rho == 0)
        rho = linspace(0,1,length(rho))';
    end
    equi.q              = helena_interp(rho',qpsi',rhoin');
    equi.shear          = x ./ equi.q .* rpdederive(x,equi.q,0,2,2,1);
    equi.psi            = - helena_interp(rho',psipout',rhoin') ./ 2 ./ pi + prof.psi(1) ;
    %equi.psi            = - interp1(rho',psipout',rhoin','linear') ./ 2 ./ pi + prof.psi(1) ;
    equi.phi            = helena_interp(rho',psiT',rhoin') ;
    equi.ftrap          = helena_interp(rho',ftra',rhoin');
    equi.ftrap(1)       = 0;
    equi.F              = helena_interp(rho',Fdia',rhoin');
    equi.b2i            = helena_interp(rho',invB2',rhoin');
    equi.b2             = helena_interp(rho',B2',rhoin');
    equi.vpr            = helena_interp(rho',Vp',rhoin');
    if equi.vpr(1) ~= 0
        disp('ZEQUI_HELENA: Problem in helena_interp');
    end
    equi.spr            = helena_interp(rho',Sp',rhoin');
    equi.d              = helena_interp(rho',xshift',rhoin');
    % compatibilite ascendante
    if isfield(equi,'ind')
        equi.ind            = zeros(size(x));
    else
        equi.indh           = NaN .* x;
        equi.indl           = NaN .* x;
    end
    equi.trh            = helena_interp(rho',xtriapos',rhoin');
    equi.trl            = helena_interp(rho',xtrianeg',rhoin');
    equi.e              = helena_interp(rho',xell',rhoin');
    equi.a              = helena_interp(rho',xrad',rhoin');
    equi.grho2r2        = helena_interp(rho',grho2r2',rhoin');
    equi.grhor          = helena_interp(rho',drhoor',rhoin');
    equi.ri             = helena_interp(rho',oor',rhoin');
    equi.grho2          = helena_interp(rho',rho2m',rhoin');
    % correction au centre de grho2 (23/07/2001)
    equi.grho2          = zregular(x,equi.grho2);
    equi.c2c            = helena_interp(rho',C2',rhoin');
    equi.grho           = helena_interp(rho',drhoav',rhoin');
    equi.rmoy           = helena_interp(rho',rav',rhoin');
    equi.jmoy           = helena_interp(rho',jout',rhoin');
    equi.ptot           = helena_interp(rho',pout' .* 1e3 .* phys.e,rhoin');
    equi.r2             = helena_interp(rho',r2',rhoin');
    equi.r2i            = helena_interp(rho',r2i',rhoin');
    equi.r2tau2         = helena_interp(rho',r2tau2',rhoin');
    equi.grho2b2        = helena_interp(rho',drho2ob2',rhoin');
    equi.r3tau3         = helena_interp(rho',r3otau3',rhoin');
    equi.r3tau          = helena_interp(rho',r3otau',rhoin');
    % fin du code modifie

    % ip et li
    equi.ip = ipout .* max(xrad) .* geo.b0 ./ phys.mu0;
    equi.li = li;

    % calcul de c2c
    %equi.c2c            = - equi.rhomax .* cumtrapz(x,equi.jmoy .* phys.mu0 .* equi.ri .* equi.vpr) ./  (rpdederive(x,equi.psi,2,2,2,1) ./ equi.rhomax);
    equi.c2c       = equi.vpr .* equi.grho2r2;
    equi.c2c(1)        = eps;
    %equi.c2c       = zcentre(equi.c2c);
    % indice de  l'avant dernier point de la grille helena sur la
    %indfin         = min(max(find(rhoin <= rho(end-1))),length(x) - 1) - 1;
    % calcul de la derniere valeur de c2c avec la formule de ip pour coherence
    % datak.gene.ip        = - 1 ./ ( 2 * pi * phys.mu0 .* datak.equi.rhomax) .* datak.equi.c2c(end) .*  datak.prof.psid1(end);
    %psid1          = rpdederive(x,equi.psi,0,2,2,1);
    %equi.c2c(end) = - equi.ip .* ( 2 * pi * phys.mu0 .* equi.rhomax) ./ psid1(end);
    %equi.c2c       = spline(x([1:indfin,end]),equi.c2c([1:indfin,end]),x);


    % nouveau calccul inverse de c2c
    %psid1         = rpdederive(x,equi.psi,0,2,2,1) ./ equi.rhomax;
    %psid1(1)     = eps;
    %inte         = - equi.rhomax .* phys.mu0 .* cumtrapz(x,equi.vpr .* equi.ri .* equi.jmoy,2);
    %equi.c2c     = inte ./ psid1;
    %equi.c2c(1) = eps;

    %figure(19);clf;plot(x,equi.c2c,x,c2c);drawnow;disp('in zequi_helena');
    % controle
    %warning off
    %psid1     = rpdederive(x,equi.psi,0,2,2,1);
    %terme1                = equi.c2c ./ equi.rhomax .* psid1;
    %jmoyc                  = - 1 ./ (phys.mu0 .* equi.vpr .* equi.rhomax) .* rpdederive(x,terme1,2,2,2,1) ./ equi.ri;
    %jmoyc                  = zcentre(jmoyc);
    %warning on
    %keyboard


    % debut du code de control
    % equic.q              = interp1(rho,qpsi,rhoin,'spline')';
    % equic.psi            = - interp1(rho,psipout,rhoin,'spline')' ./ 2 ./ pi + prof.psi(1) ;
    % equic.phi            = interp1(rho,psiT,rhoin,'spline')' ;
    % equic.ftrap          = interp1(rho,ftra,rhoin,'spline')';
    % equic.F              = interp1(rho,Fdia,rhoin,'spline')';
    % equic.b2i            = interp1(rho,invB2,rhoin,'spline')';
    % equic.b2             = interp1(rho,B2,rhoin,'spline')';
    % equic.vpr            = interp1(rho,Vp,rhoin,'spline')';
    % equic.spr            = interp1(rho,Sp,rhoin,'spline')';
    % equic.d              = interp1(rho,xshift,rhoin,'spline')';
    % equic.ind            = zeros(size(x));
    % equic.trh            = interp1(rho,xtriapos,rhoin,'spline')';
    % equic.trl            = interp1(rho,xtrianeg,rhoin,'spline')';
    % equic.e              = interp1(rho,xell,rhoin,'spline')';
    % equic.a              = interp1(rho,xrad,rhoin,'spline')';
    % equic.grho2r2        = interp1(rho,grho2r2,rhoin,'spline')';
    % equic.grhor          = interp1(rho,drhoor,rhoin,'spline')';
    % equic.ri             = interp1(rho,oor,rhoin,'spline')';
    % equic.grho2          = interp1(rho,rho2m,rhoin,'spline')';
    % equic.grho           = interp1(rho,drhoav,rhoin,'spline')';
    % equic.rmoy           = interp1(rho,rav,rhoin,'spline')';
    % equic.jmoy           = interp1(rho,jout,rhoin,'spline')';
    % equic.ptot           = interp1(rho,pout .* 1e3 .* phys.e,rhoin,'spline')';
    % equic.r2             = interp1(rho,r2,rhoin,'spline')';
    % equic.r2i            = interp1(rho,r2i,rhoin,'spline')';
    % equic.r2tau2         = interp1(rho,r2tau2,rhoin,'spline')';
    % equic.grho2b2        = interp1(rho,drho2ob2,rhoin,'spline')';
    % equic.r3tau3         = interp1(rho,r3otau3,rhoin,'spline')';
    % equic.r3tau          = interp1(rho,r3otau,rhoin,'spline')';
    % save control_equi equi equic;
    %
    % % test des differences
    % names  = fieldnames(equic);
    % for k = 1:length(names)
    %   new = getfield(equi,names{k});
    %   old = getfield(equic,names{k});
    %   err = sqrt(sum((new-old).^2))./sqrt(sum((new+old).^2)).*4;
    %   fprintf('delta(%s) = %g\n',names{k},err);
    % end
    %
    %
    % fin du code de control

    % fin de la modification du 15/06/2001

    % le R et le Z
    equi.R     = single(shiftdim(R,-1));
    % z est compatible avec la separatrice
    % correction du 06/12/2004
    if geo.mode == 2
        equi.Z     = single(shiftdim(Z,-1) + z0_);
    else
        equi.Z     = single(shiftdim(Z,-1) + z0);
    end
    equi.rhoRZ = single(shiftdim(rho,-1));
    equi.psiRZ = single(shiftdim(-psipout./ 2 ./ pi + prof.psi(1),-1));
    equi.df2RZ = single(shiftdim(df2,-1));
    equi.dprRZ = single(shiftdim(dpr,-1));
    equi.frmode= single(shiftdim(frmode(1:nbmode),-1));

    % La carte de champ
    BPHI       = (Fdia * ones(1,size(R,2))) ./ R;
    equi.BR    = single(shiftdim(BPR,-1));
    equi.BZ    = single(shiftdim(BPZ,-1));
    equi.BPHI  = single(shiftdim(BPHI,-1));


    % calcul de rhog
    % ajout du 10/07/2001
    equi.rhog = equi.a ./ equi.a(end);
    % ajout du 12/07/2001
    raxe      = (min(R') + max(R')) ./ 2;
    equi.raxe = helena_interp(rho',raxe,rhoin');
    equi.zaxe = NaN .* equi.raxe;

    % donnees MHD
    equi.mhd.ballooning     = helena_interp(rho',(1- balcrit') .* (balcrit' ~= 0),rhoin');
    equi.mhd.mercier  = helena_interp(rho',- drmerc',rhoin');
    equi.mhd.ideal    = helena_interp(rho',0.25 - dimerc',rhoin');
    % si mhd_cd n'est pas calculer par mishka
    if any(~isfinite(prof.mhd_cd))
        v                       = equi.mhd.ballooning .* (equi.mhd.ballooning > 0);
        v                       = 1 - v ./ max(eps,max(v));
        mhd_cd                  = 1 - cons.mhd_bal .* v;
    else
        mhd_cd                  = prof.mhd_cd;
    end
    % memorisation des donnees pour le calcul rapide
else
    disp('HELENA failure, last equilibrium is kept')
    mhd_cd                  = prof.mhd_cd;
end
if (ifail == 0) & ~isempty(b) & ~isempty(df2) & ~isempty(kkbig) &(retry == 0)

    % l'equilibre precedent est valide
    mem.init    = 1;
    mem.b       = b;
    mem.df2     = df2;
    mem.kkbig   = kkbig;
    mem.ptotold = ptot;
else
    % l'equilibre precedent n'est pas valide
    mem.init  = 0 ;
    mem.b       = 0;
    mem.df2     = 0 .* ones(size(equi.df2RZ));
    mem.kkbig = [];
    mem.ptotold = [];
end

setappdata(0,'equilibre_helena_init',mem.init);
setappdata(0,'equilibre_helena_b',mem.b);
setappdata(0,'equilibre_helena_df2',mem.df2);
setappdata(0,'equilibre_helena_kkbig',mem.kkbig);
setappdata(0,'equilibre_helena_ptotold',mem.ptotold);
setappdata(0,'equilibre_helena_corj',mem.corj);


% geodesique (nchi,cpsurf,radius,gem11,gem12,gem33)
straightline.nchi   = nchi;
straightline.cpsurf = cpsurf;
straightline.radius = radius;
straightline.psi    = psipout;
straightline.rho    = rho;
straightline.q      = qpsi;
straightline.Fdia   = Fdia;

straightline.gem11  = reshape(gem11,nchi,length(psipout))';
straightline.gem11(:,end+1) = straightline.gem11(:,1);

straightline.gem12  = reshape(gem12,nchi,length(psipout))';
straightline.gem12(:,end+1) = straightline.gem12(:,1);

straightline.gem33  = reshape(gem33,nchi,length(psipout))';
straightline.gem33(:,end+1) = straightline.gem33(:,1);

straightline.rspot = reshape(rspot,nchi,length(psipout))';
straightline.rspot(:,end+1) = straightline.rspot(:,1);


if geo.mode == 2
    straightline.zspot = reshape(zspot,nchi,length(psipout))' + z0_;
else
    straightline.zspot = reshape(zspot,nchi,length(psipout))' + z0;
end
straightline.zspot(:,end+1) = straightline.zspot(:,1);

straightline.brspot = reshape(brspot,nchi,length(psipout))';
straightline.brspot(:,end+1) = straightline.brspot(:,1);

straightline.bzspot = reshape(bzspot,nchi,length(psipout))';
straightline.bzspot(:,end+1) = straightline.bzspot(:,1);


% si le calcul de deltaprime est demande
if (cons.dponoff  == 1)
    small  = 1e-10;
    rw     = 1;
    m_deltap = NaN .* ones(size(x));
    deltap = NaN .* ones(size(x));
    % recherche des nombre m et n tel que q = m/n;
    for m = 1:cons.mmax
        for n=1:cons.nmax
            ql= m ./n;
            d = abs(equi.q - ql);
            indmn = min(find(d == min(d)));
            if (indmn > 1) & (ql >= min(equi.q)) & (ql <= max(equi.q))
                if ~isfinite(m_deltap(indmn))
                    % on calcul deltap pour le m le plus petit
                    param1           = [n m cons.ng cons.np 1 cons.sigs 0.5 0];
                    tabtr            = dpfem_1(x,equi.q,equi.jmoy,param1,small,rw);
                    deltap(indmn)   = tabtr(4,1);
                    m_deltap(indmn) = m;
                end
            end
        end
    end

    % ecriture dans la structure equi.mhd
    equi.mhd.deltap   = deltap;
    equi.mhd.m_deltap = m_deltap;

end


% integralle volumique
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax
function s=zintvol(e,x,vpr,rhomax)

s = rhomax.*trapz(x,vpr .* e,2);



% interpolation dans helena
function yout = helena_interp(xin,yin,xout)

yout = pchip(xin,yin,xout);

