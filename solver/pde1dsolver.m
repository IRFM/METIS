% PDE1DSOLVER resoud un systeme de pde 1d couplees
%--------------------------------------------------------------
% fichier pde1dsolver.m ->  pde1dsolver
%
%
% fonction Matlab 5 : 
%
% Cette fonction resoud un systeme d'equation aux derivees partielles
% a une dimension par la methode des trapezes (methode utilisant des 
% differnces finies). Les equations doivent etre de la forme :
%
% diff(F,t) = A *  diff(F,x,x) + B * diff(F,x) + C * F + D
%
% La resolution par la methode des trapezes est une resolution au 
% sens des moindre carres utilisant la methode explicite et implicite.
% Les conditions aux limites doivent etre donnees au deux extremites
% de l'intervalle x -> [ x0 x1]. Initialement, la matrice F doit etre 
% connue.
%
% Il y a un probleme lorsque les condition aux limites sont nulles. Le calcul
% des indices des elements a mettre a 0 est fausse (changement de type de 
% conditions au limites). la correction consiste a ajouter eps aux valeurs nulles.
% Cela peut fausser les valeurs si le solveur est utilise avec des donnees petites.
%
%
%
% syntaxe :
% 
% [fp,{dfpdx,d2fpdx2,FS,F,ALPHA,ALPHAP,IDENTITE}] = ...
% pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt)
%
%
% Les indices des matrices :
%
% Le premier indice (k) des matrices de coefficients est l'indice 
% d'espace :  
%                 x=x0+(k-1)*dx  , k->[1..K]
% 
% Le deuxieme indice (n) des matrices de coefficients est l'indice 
% qui donne le coefficient de l'equation m qui relie l'equation m
% a l'equation n (n -> [1 M])
% 
% Le troisieme indice (m) des matrices de coefficients est l'indice 
% qui indique le rang (ou numero) de l'equation (m -> [1..M])
% 
% Les entrees :
% 
%  A  = matrice des coefficients lies a la derivee seconde de F 
%       au temps t: A=A(k,n,m)
%
%  B  = matrice des coefficients lies a la derivee premiere de F 
%       au temps t: B=B(k,n,m)
%
%  C  = matrice des coefficients lies a  F au temps t: 
%       C=C(k,n,m)
%
%  D  = matrice des coefficients independant de  F 
%       au temps t: D=D(k,m)
%
% 
%  AP  = matrice des coefficients lies a la derivee seconde de F 
%        au temps t+dt: AP=AP(k,n,m)
%
%  BP  = matrice des coefficients lies a la derivee premiere de F 
%        au temps t+dt. BP=BP(k,n,m)
%
%  CP  = matrice des coefficients lies a  F au temps t+ dt:
%        C=C(k,n,m)
%
%  DP  = matrice des coefficients independant de  F 
%        au temps t+dt: DP=DP(k,m)
%
%  F   = matrice des valeurs au temps t    : F=F(k,m)
%
%  FP  = matrice des valeurs au temps t+dt : FP=FP(k,m)
%        (utiliser pour les equations en mode interpretatif,
%         les valeurs inconnues doivent etre remplace par des 0)
%
% Les conditions aux limites sont de la forme :
% 
%   a * diff(F,x,x)(xk) + b * diff(F,x)(xk) + c * F(xk) = d
%   avec xk = x0 ou x1
%   
%  V0  = matrice des valeurs des conditions aux limites en x=x0:
%        V0=V0(l,m) avec l=1 pour le temps t et l=2 pour le temps
%        t+dt.
%
%  T0  = vecteur indiquant le type de condition au limite 
%        en x=x0 : T0=T0(l,m,j) avec l=1 pour le temps t et 
%        l=2 pour le temps t+dt. L'indice j est relatif la 
%        derivee de F en x0 a l'ordre j-1 :
%               j = 1  -> coefficient portant sur F        
%               j = 2  -> coefficient portant sur diff(F,x)        
%               j = 3  -> coefficient portant sur diff(F,x,x)        
%        
%
%  V1  = matrice des valeurs des conditions aux limites en x=x1:
%        V1=V1(l,m) avec l=1 pour le temps t et l=2 pour le temps
%        t+dt.
%
%  T1  = vecteur indiquant le type de condition wu limite 
%        en x=x1 : T1=T1(l,m,j) avec l=1 pour le temps t et 
%        l=2 pour le temps t+dt. L'indice j est relatif la 
%        derivee de F en x0 a l'odre j-1 :
%               j = 1  -> coefficient portant sur F        
%               j = 2  -> coefficient portant sur diff(F,x)        
%               j = 3  -> coefficient portant sur diff(F,x,x)
%
% mode = vecteur indiquant, pour chaque equation, si elle est en 
%        mode preditif ou interpretatif : mode =mode(m). 
%        Les valeurs possible sont :
%                    0  -> predictif
%                    1  -> interpretatif
%
%  f   = scalaire indiquant le mode de fonctionnement du solver
%        f -> [0.0..1.0] & -1 . Valeurs particulieres :
%                0   -> implicite pure (pb de recherche d'equilibre)
%                0.5 -> Crank-Nicolson 
%                0.3 -> preditif stable
%                1   -> explicite pure
%                -1  -> exponential integrator order 2
%
%        remarque : le solveur est stable pour f -> [0.0..0.5]
%
%  dx   =  pas d'espace (par defaut 1)
%  
%  dt   =  pas de temps (par defaut 1)  
% 
% sortie :
%
%  fp        =  matrices de la nouvelle estimation des valeurs au temps
%               t+dt : fp=fp(k,m)    
%
%  dfpdx     =  matrices de la nouvelle estimation de la derivee 
%               au temps t+dt : dfpdx=dfpdx(k,m)    
%
%  d2fpdx2   =  matrices de la nouvelle estimation de la derive 
%               2ieme au temps t+dt : d2fpdx2=d2fpdx2(k,m)    
%
%  les autres sorties ne servent que lors des tests de la fonction
%
% exemple : lancer testpde1dsolver -> test complet + exemple 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.9, du 29/05/2002.
% 
% ordre de compilation : mcc -x pde1dsolver
%                        ( chemin vers le compilateur -> addpath /usr/local/matlab5/toolbox/compiler)
% 
% liste des modifications : 
% 
% * 06/08/99 -> calcul des derivees en sortie uniquement si necessaire
% * 23/08/99 -> mise a zero des coefficients des equations interpretatives 
% * 23/08/99 -> retrait des Nan et Inf dans les donnees d'entrees -> securite
% * 14/06/2001 -> utilisation des matrices creuses pour tous les elements de l'equation
% * 14/06/2001 -> compilation avec mcc (V2.0)
% * 15/06/2001 -> suppression de la compilatrion (gain trop petit ~10%)
% * 29/05/2002 -> ajout de eps pour les conditions aux limites nulles
%--------------------------------------------------------------
%
function [fp,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= ...
         pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt)

% mode expo_inte
mode_expo_inte = 0;

% appel du mexfile si disponible
if (nargout == 1) && (f  >= 0)
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
			disp('Using mexfile version of pde1dsolver');
		else
			mexsolver = 0;	
		end
		setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
	end  
	%
	if mexsolver == 1
		if (f == 1) | all(mode == 1)
			warning('the explicit Euler algorithm must not be used with mexpde1dsolver; switch to expo_euler solver');
            		mode_expo_inte = 1;
            		f              = 0;
                        setappdata(0,'MEXSOLVER_IN_METIS',0);
       		elseif isdeployed
                        try
			    fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt);
			    return

                        catch
			    warning('mexpde1dsolver is not working; switch to expo_euler solver');
			    mode_expo_inte = 1;
			    f              = 0;
			    setappdata(0,'MEXSOLVER_IN_METIS',0);
                        end
                else
            		fp          = mexpde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt);
            		return
        	end
	end
elseif f < 0
    mode_expo_inte = 1;
    f              = 0;
end

% pas de verifications des entrees - gain de temps

% correction des conditions aux limites si  == 0 (ajout de eps)
% debut correction
V1(V1 == 0) = eps;
V0(V0 == 0) = eps;
% fin correction 

% les dimensions
K=size(F,1);
M=size(F,2);


% mise en forme de FP (nul pour le mode predictif)
% modification du 07/08/2000 -> securite
n=find(~mode);
if ~isempty(n)
	FP(:,n) = zeros(K,length(n));
end
% fin modification du 07/08/2000

% modification du 23/08/2000
% gain de temps
n=find(mode);
if ~isempty(n)
	A(:,:,n)   = zeros(K,M,length(n));
	B(:,:,n)   = zeros(K,M,length(n));
	C(:,:,n)   = zeros(K,M,length(n));
	D(:,n)     = zeros(K,length(n));
	AP(:,:,n)  = zeros(K,M,length(n));
	BP(:,:,n)  = zeros(K,M,length(n));
	CP(:,:,n)  = zeros(K,M,length(n));
	DP(:,n)    = zeros(K,length(n));
end

% retrait des NaN et des inf dans les donnees d'entrees -> securite
ind =find(~isfinite(F));
if ~isempty(ind)
	F(ind) = zeros(1,length(ind));
	zverbose('=>>>> NaN ou Inf dans les donnees des PDE !!!\n');
end
ind =find(~isfinite(FP));
if ~isempty(ind)
	FP(ind) = zeros(1,length(ind));
	zverbose('=>>>> NaN ou Inf dans les donnees des PDE !!!\n');
end
% fin modification du 23/08/2000 


% constante dimentionnelle
dtx2 = dt ./ dx .^ 2;
dt2x = dt ./ 2 ./ dx;

% 1 - combinaison des coeficients
a = dtx2 .* A - dt2x .* B;
b = - 2 .* dtx2 .* A  + dt .* C;
c = dtx2 .* A + dt2x .* B;
d = dt .* D;

ap = dtx2 .* AP - dt2x .* BP;
bp = - 2 .* dtx2 .* AP  + dt .* CP;
cp = dtx2 .* AP + dt2x .* BP;
dp = dt .* DP;

s = f .* d + (1-f) .* dp;

% 2 - condition aux limites en x=x0
% a - conditionnement des donnees
u = (T0(:,:,3) ./ (dx.^2) + T0(:,:,2) ./ (2.*dx));
v = (T0(:,:,1) - 2 .* T0(:,:,3) ./ (dx.^2));
w = (T0(:,:,3) ./ (dx.^2) - T0(:,:,2) ./ (2*dx));

warning off
U = - u ./ w;
V = - v ./ w;
W = V0 ./ w;
X = - u ./ v;
Y = V0 ./ v;
Z = V0 ./ u;
warning on

comp = ones(1,M);

% b - cas  w~=0
i1 = find( w(1,:) ~= 0 );
i2 = find( w(2,:) ~=0) ;
if ~ (isempty(i1) & isempty(i2))
	b(1,i1,:) = b(1,i1,:) + a(1,i1,:) .* permute(V(1,i1)' * comp,[3 1 2]);
	bp(1,i2,:) = bp(1,i2,:) + ap(1,i2,:) .* permute(V(2,i2)' * comp,[3 1 2]);

	c(1,i1,:) = c(1,i1,:) + a(1,i1,:) .* permute(U(1,i1)' * comp,[3 1 2]);
	cp(1,i2,:) = cp(1,i2,:) + ap(1,i2,:) .* permute(U(2,i2)' * comp,[3 1 2]);

	s(1,:) = s(1,:) + squeeze( f .* sum( a(1,i1,:) .* permute(W(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( ap(1,i2,:) .* permute(W(2,i2)' * comp,[3 1 2]),2))';
	                     
	a(1,i1,:)   = zeros(1,length(i1),M);
	ap(1,i2,:)  = zeros(1,length(i2),M);
end

% c - cas w ==0 et v ~=0
i1 = find(( w(1,:) == 0 ) & (v(1,:) ~= 0));
i2 = find(( w(2,:) == 0 ) & (v(2,:) ~= 0));
if ~ (isempty(i1) & isempty(i2))
	c(1,i1,:) = c(1,i1,:) + b(1,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]);
	cp(1,i2,:) = cp(1,i2,:) + bp(1,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]);

	b(2,i1,:) = b(2,i1,:) + a(2,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]);
	bp(2,i2,:) = bp(2,i2,:) + ap(2,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]);

	s(1,:) = s(1,:) + squeeze( f .* sum( b(1,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( bp(1,i2,:) .* permute(Y(2,i2)' * comp,[3 1 2]),2))';

	s(2,:) = s(2,:) + squeeze( f .* sum( a(2,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( ap(2,i2,:) .* permute(Y(2,i2)' * comp,[3 1 2]),2))';
	
	F(1,i1)  = X(1,i1) .* F(2,i1)  + Y(1,i1);
	FP(1,i2) = X(2,i2) .* FP(2,i2) + Y(2,i2);
                     
	a(2,i1,:)   = zeros(1,length(i1),M);
	ap(2,i2,:)  = zeros(1,length(i2),M);
	b(1,i1,:)   = zeros(1,length(i1),M);
	bp(1,i2,:)  = zeros(1,length(i2),M);
end

% d - cas w ==0, v ==0 et u~=0
i1 = find(( w(1,:) == 0 ) & (v(1,:) == 0) & (u(1,:) ~= 0));
i2 = find(( w(2,:) ==0) & (v(2,:) == 0) & (u(2,:) ~= 0));
if ~ (isempty(i1) & isempty(i2))
	s(1,:) = s(1,:) + squeeze( f .* sum( c(1,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( cp(1,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';

	s(2,:) = s(2,:) + squeeze( f .* sum( b(2,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( bp(2,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';
                     
	s(3,:) = s(3,:) + squeeze( f .* sum( a(3,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( ap(3,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';

	F(2,i1)  = Z(1,i1);                     
	FP(2,i2) = Z(2,i2);                     
	                     
	a(3,i1,:)   = zeros(1,length(i1),M);
	ap(3,i2,:)  = zeros(1,length(i2),M);
	b(2,i1,:)   = zeros(1,length(i1),M);
	bp(2,i2,:)  = zeros(1,length(i2),M);
	c(1,i1,:)   = zeros(1,length(i1),M);
	cp(1,i2,:)  = zeros(1,length(i2),M);
end

% 3 - condition aux limites en x=x1
% a - conditionnement des donnees
u = (T1(:,:,3) ./ (dx.^2) + T1(:,:,2) ./ (2.*dx));
v = (T1(:,:,1) - 2 .* T1(:,:,3) ./ (dx.^2));
w = (T1(:,:,3) ./ (dx.^2) - T1(:,:,2) ./ (2*dx));

warning off
U = - w ./ u;
V = - v ./ u;
W = V1 ./ u;
X = - w ./ v;
Y = V1 ./ v;
Z = V1 ./ w;
warning on 

comp = ones(1,M);

% b - cas  u~=0
i1 = find( u(1,:) ~= 0 );
i2 = find( u(2,:) ~=0) ;
if ~ (isempty(i1) & isempty(i2))
	a(K,i1,:) = a(K,i1,:) + c(K,i1,:) .* permute(U(1,i1)' * comp,[3 1 2]);
	ap(K,i2,:) = ap(K,i2,:) + cp(K,i2,:) .* permute(U(2,i2)' * comp,[3 1 2]);

	b(K,i1,:) = b(K,i1,:) + c(K,i1,:) .* permute(V(1,i1)' * comp,[3 1 2]);
	bp(K,i2,:) = bp(K,i2,:) + cp(K,i2,:) .* permute(V(2,i2)' * comp,[3 1 2]);

	s(K,:) = s(K,:) + squeeze( f .* sum( c(K,i1,:) .* permute(W(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( cp(K,i2,:) .* permute(W(2,i2)' * comp,[3 1 2]),2))';

	c(K,i1,:)   = zeros(1,length(i1),M);
	cp(K,i2,:)  = zeros(1,length(i2),M);
end                     
% c - cas u ==0 et v ~=0
i1 = find(( u(1,:) == 0 ) & (v(1,:) ~= 0));
i2 = find(( u(2,:) ==0) & (v(2,:) ~= 0));
if ~ (isempty(i1) & isempty(i2))
	a(K,i1,:) = a(K,i1,:) + b(K,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]);
	ap(K,i2,:) = ap(K,i2,:) + bp(K,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]);
	
	% bug corrige le 07/08/2000
	b(K-1,i1,:) = b(K-1,i1,:) + c(K-1,i1,:) .* permute(X(1,i1)' * comp,[3 1 2]);
	bp(K-1,i2,:) = bp(K-1,i2,:) + c(K-1,i2,:) .* permute(X(2,i2)' * comp,[3 1 2]);

	s(K,:) = s(K,:) + squeeze( f .* sum( b(K,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( bp(K,i2,:) .* permute(Y(2,i2)' * comp,[3 1 2]),2))';

	s(K-1,:) = s(K-1,:) + squeeze( f .* sum( c(K-1,i1,:) .* permute(Y(1,i1)' * comp,[3 1 2]),2) + ...
	                         (1 - f) .* sum( cp(K-1,i2,:) .* permute(Y(2,i2)' * comp,[3 1 2]),2))';
                     
	F(K,i1)  = X(1,i1) .* F(K-1,i1)  + Y(1,i1);
	FP(K,i2) = X(2,i2) .* FP(K-1,i2) + Y(2,i2);
	                         	                         
	c(K-1,i1,:)   = zeros(1,length(i1),M);
	cp(K-1,i2,:)  = zeros(1,length(i2),M);
	b(K,i1,:)   = zeros(1,length(i1),M);
	bp(K,i2,:)  = zeros(1,length(i2),M);
end
% d - cas u ==0, v ==0 et w ~=0
i1 = find(( u(1,:) == 0 ) & (v(1,:) == 0) & (w(1,:) ~= 0));
i2 = find(( u(2,:) ==0) & (v(2,:) == 0) & (w(2,:) ~= 0));
if ~ (isempty(i1) & isempty(i2))
	s(K,:) = s(K,:) + squeeze( f .* sum( a(K,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                     (1 - f) .* sum( ap(K,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';

	s(K-1,:) = s(K-1,:) + squeeze( f .* sum( b(K-1,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                         (1 - f) .* sum( bp(K-1,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';
                     
	s(K-2,:) = s(K-2,:) + squeeze( f .* sum( c(K-2,i1,:) .* permute(Z(1,i1)' * comp,[3 1 2]),2) + ...
	                         (1 - f) .* sum( cp(K-2,i2,:) .* permute(Z(2,i2)' * comp,[3 1 2]),2))';

	F(K-1,i1)  = Z(1,i1);                     
	FP(K-1,i2) = Z(2,i2);   
	
	a(K,i1,:)     = zeros(1,length(i1),M);
	ap(K,i2,:)    = zeros(1,length(i2),M);
	b(K-1,i1,:)   = zeros(1,length(i1),M);
	bp(K-1,i2,:)  = zeros(1,length(i2),M);
	c(K-2,i1,:)   = zeros(1,length(i1),M);
	cp(K-2,i2,:)  = zeros(1,length(i2),M);
end
% 4 - equation predictive et interpretative
n = find(mode);            % c'est du matlab 5 !
if ~isempty(n)

   % decalage d'indice pour ce terme -1
   Fa=F(1:(K-1),n);  
   FPa=FP(1:(K-1),n);  
   % on cree un matrice de la forme [F,F,F...F] (M fois) 
   % sur la troisieme dimension
   Fa=reshape(Fa(:)*ones(1,M),K-1,length(n),M);     
   FPa=reshape(FPa(:)*ones(1,M),K-1,length(n),M);     

   % decalage d'indice pour ce terme 0
   Fb=F(:,n);  
   FPb=FP(:,n);  
   % on cree un matrice de la forme [F,F,F...F] (M fois) 
   % sur la troisieme dimension
   Fb=reshape(Fb(:)*ones(1,M),K,length(n),M);     
   FPb=reshape(FPb(:)*ones(1,M),K,length(n),M);     
   
   % decalage d'indice pour ce terme +1
   Fc=F(2:K,n);  
   FPc=FP(2:K,n);  
   % on cree un matrice de la forme [F,F,F...F] (M fois) 
   % sur la troisieme dimension
   Fc=reshape(Fc(:)*ones(1,M),K-1,length(n),M);     
   FPc=reshape(FPc(:)*ones(1,M),K-1,length(n),M);     
   
   % modification du terme source    
   s     = s + squeeze( sum( ...
           f  .*  (  cat(1 , zeros(1,length(n),M), a(2:K,n,:).*Fa) + ...
                     b(:,n,:) .* Fb + ...
                     cat(1 , c(1:(K-1),n,:) .*Fc , zeros(1,length(n),M))) + ...
       (1 - f) .* (  cat(1 , zeros(1,length(n),M), ap(2:K,n,:).*FPa) + ...
                     bp(:,n,:) .* FPb + ...
                     cat(1 , cp(1:(K-1),n,:) .*FPc , zeros(1,length(n),M))) , 2));
                                    
    % annulation des coefficients deja utilises
    comp=zeros(K,length(n),M);
    a(:,n,:)=comp;
    ap(:,n,:)=comp;
    b(:,n,:)=comp;
    bp(:,n,:)=comp;
    c(:,n,:)=comp;
    cp(:,n,:)=comp;
    
end

% 5 - creation des matrices creuses et des vecteurs pour
%     la resolution du systeme lineaire (au sens des moindres carres)

% a - Les vecteurs de donnees au temps t
FS = F(:);
S  = s(:) + FS;

% vecteur d'indice
ia  = 2:K;
ib  = 1:K;
ic  = 1:(K - 1);

% initialisation des buffer a vide
alpha  = [];
alphap = [];
p      = [];
q      = [];

% boucle de remplissage des bufffers
for m = 1:M
	for n = 1:M
		
		% pour a
		alpha    = cat(1,alpha,a(ia,n,m));
		alphap   = cat(1,alphap,ap(ia,n,m));
		p        = cat(1,p,ia' + K .* (m - 1)); 
		q        = cat(1,q,ia' - 1 + K .* (n - 1)); 
		
		% pour b
		alpha    = cat(1,alpha,b(ib,n,m));
		alphap   = cat(1,alphap,bp(ib,n,m));
		p        = cat(1,p,ib' + K .* (m - 1)); 
		q        = cat(1,q,ib' + K .* (n - 1)); 
		
		% pour c
		alpha    = cat(1,alpha,c(ic,n,m));
		alphap   = cat(1,alphap,cp(ic,n,m));
		p        = cat(1,p,ic' + K .* (m - 1)); 
		q        = cat(1,q,ic' + 1 + K .* (n - 1)); 
		
	end
end
	
% duplication des indice avnat retrait des elements non nul	
pp = p;
qp = q;

% g - retrait des elements nuls
iz = find(alpha  ~= 0);
alpha = alpha(iz);
p  =p(iz);
q = q(iz);

izp = find(alphap ~= 0);
alphap = alphap(izp);
pp = pp(izp);
qp = qp(izp);


% 6 - resolution du systeme 
%fpv =  (IDENTITE - ALPHAP) \ ( ALPHA * FS  + S);
% modification du 14/06/2001 
if mode_expo_inte == 0
	% h - creation des matrices creuses 
	km = K .* M ;
	ALPHA     = sparse(p,q, f .* alpha, km , km);
	ALPHAP    = sparse(pp,qp, (1 - f) .*  alphap, km, km);
	IDENTITE  = sparse(1:km, 1:km,1, km, km);
    	fpv =  full((IDENTITE - ALPHAP) \ ( ALPHA * sparse(FS)  + sparse(S)));    
 
else
	% exponential integrator
	% h - creation des matrices creuses 
	km = K .* M ;
	ALPHAP    = sparse(pp,qp, (1 - f) .*  alphap, km, km);
	[fpv,conv] = expo_inte(ALPHAP, S(:) - FS(:), FS(:), M);
	
	%fpvn =  full((IDENTITE - ALPHAP) \ ( ALPHA * sparse(FS)  + sparse(S)));
	%indp = 1:21;
	%figure(21);hold on ;plot(indp,fpvn(indp) -  FS(indp),'g');drawnow
	
        % cas de non convergence, implicite pure
	if conv == 0
		% h - creation des matrices creuses 
		ALPHA     = sparse(p,q, f .* alpha, km , km);
		IDENTITE  = sparse(1:km, 1:km,1, km, km);
		fpv =  full((IDENTITE - ALPHAP) \ ( ALPHA * sparse(FS)  + sparse(S)));
	end
end



%isa(( ALPHA * sparse(FS)  + sparse(S)),'sparse')
%isa((IDENTITE - ALPHAP),'sparse')

%fprintf('Rang = %d, rcond = %g\n',sprank(IDENTITE - ALPHAP),1 ./ condest(IDENTITE - ALPHAP));

% 7 - mise en forme de la sortie
fp = reshape(fpv,[K,M]);

% 8  - on recupere les donnees interpretees (modifier par les conditions aux limites)
n = find(FP);            
if ~isempty(n)	
  fp(n)=FP(n);
end
 
% 9 - calcul des deriveees
if nargout >1 

	% a - prolongation de la derivee 2 
	%     par continiute de la derive 3.
	fp0   = 4 .* fp(1,:)  - 6 .* fp(2,:)   + 4 .* fp(3,:)    - fp(4,:);
	fpkp1 = 4 .* fp(K,:)  - 6 .* fp(K-1,:) + 4 .* fp(K-2,:)  - fp(K-3,:);

	% b - matrice pour le calcul des derivees
	fpp = cat(1,fp0,fp,fpkp1);

	% c - calcule des derivees
	im = 1:K;
	i0 = 2:(K+1);
	ip = 3:(K +2);

	dfpdx   = (fpp(ip,:) - fpp(im,:)) ./ (2 .* dx);
	d2fpdx2 = (fpp(ip,:) - 2 .* fpp(i0,:) + fpp(im,:)) ./ (dx .^ 2);  
	
end


function [y_np1,conv] = expo_inte(L,G,y_n,nbeq)

% L = operateur lineaire (apporximation du jacobien), 
% G =  valeur de la partie non lineaire (source)
% y_n = valeur precedente

%L0 = diag(zeros(1,size(L,2)));
L1 = speye(size(L));
L0 = 0 .* L1;
%issparse(cat(1,cat(2,L, L1),cat(2,L0,L0)))
LF = zspexpm(cat(1,cat(2,L, L1),cat(2,L0,L0)));
%LFa = expm(cat(1,cat(2,L, L1),cat(2,L0,L0)));
%max(abs(LF(:) -LFa(:)))
phi0 = LF(1:size(L,1),1:size(L,2));
phi1  = LF(1:size(L,1),(1 + size(L,2)):end);
%  LF = expm(cat(1,cat(2,L, L1,L0),cat(2,L0,L0,L1),cat(2,L0,L0,L0)));
%  phi0 = LF(1:size(L,1),1:size(L,2));
%  phi1 = LF(1:size(L,1),(1 + size(L,2)):(2 .* size(L,2)));
%  phi2 = LF(1:size(L,1),(1 + 2 .* size(L,2)):end);
y_np1 = phi0 * y_n + phi1 * G;

% methode expokit (reference deja testee mais la precision est moindre)
%[y_np1_n, err] = phiv( 1, L, G, y_n);

% cette methode introduit un offset : la cause est inconnue !!!
%  % vecteur propres et valeurs propres
%  if nbeq == 1
%      [V,D]      = eig(full(L)); % attention le preconditionneur ajoute du bruit numerique si les donnees ont des echelles differentes!
%  else
%      [V,D]      = eig(full(L),'nobalance');
%  end
%  DD         = diag(D);
%  Vm1        = inv(V);
%  
%  % reduction du bruit
%  bruit = Vm1*V - diag(diag(Vm1*V));
%  bruit = pi .* max(abs(bruit(:)));
%  rap = abs(DD) ./ max(realmin,max(abs(DD)));
%  DD(rap <= bruit) = 0;
%  
%  % diagonalisation du probleme
%  U          = Vm1 * G;
%  z_n        = Vm1 * y_n(:);
%  
%  % facteur integrant
%  texp    = exp(DD);
%  % amelioration de la precision pour le petit pas de temps
%  texpm1  = expm1(DD);
%  
%  % EDT + collocation in t_n + dt
%  
%  % prise en compte exacte des valeurs propres nulles
%  DDc = DD;
%  ind0 = find(DD== 0);
%  DDc(ind0) = 1;
%  texpm1(ind0) = 0;
%  texp(ind0) = 1;
%  
%  % integration du probleme diagonal
%  % EDT + collocation in t_n + dt/2
%  z_np1 = z_n .* texp + U ./ DDc .* texpm1;
%  
%  % retour dans l'espace direct
%  y_np1_n = V * z_np1;
%  y_np1 = V * z_np1;

% securite si valeur  imaginaires
if any(imag(y_np1(:)))
   if any((abs(imag(y_np1)) ./ max(eps,abs(real(y_np1)))) > sqrt(eps))
      warning('EXPO_EULER: imag ');
   end
   y_np1 = real(y_np1);
   conv = 0;
else
   conv = 1;
end

%indp = 1:21;
%figure(21);clf;plot(indp,y_np1(indp) - y_n(indp),'or',indp,y_np1_n(indp) - y_n(indp),'+b',indp,y_np1_n(indp) -y_np1(indp),'m');
