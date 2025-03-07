function [T,c]=test2eq_pde1dsolver(auto)

% sript de test #1 du solver : test des dimensions
disp('Test de la fonction pde1dsolver')
disp('test 2 equations')
disp('Probleme non lineaire de la diffusion d''un produit concentre dans un liquide')
disp('en presence d''un gradient de temperature et d''une source de chauffage.')
disp('La conductance thermique du liquide est modifiee par la concentration du produit.')

dbstop if error
K = 21; % nombre de points d'espaces
M = 2;  % nombre d'equations
L=200;  % nombre de temps

if nargin ==0
	alpha=input('Le coefficient de diffusion du concentre dans le liquide ? '); 
	xie=input('Le coefficient de diffusion de la chaleur dans le liquide ? ');
	c0=input('la concentration au centre du tube ? ');
	t0=input('la temperature en x=0 ? ');
	t1=input('la temperature en x=1 ? ');
	q=input('le flux de chaleur au centre du tube ? ');
else
	alpha=0.1; 
	xie=0.01;
	c0=pi;
	t0=-1;
	t1=1;
	q=1;	
end
t=zeros(L+1,1);                 % vecteur des temps 
x=linspace(0,1,K);              % vecteur d'expace 0..1
dx=mean(diff(x));               % pas d'espace
T=zeros(L+1,K);                 % matrice des temperatures
T(1,:)=linspace(t0,t1,K);       % initialisation du profile initial (lineaire)
c=zeros(L+1,K);                 % matrice des concentrations
c(1,fix(K/2)+1)=c0;             % le concentre est mis au centre du tube

% i,itialisation des matrices a zero (la plupart des elements sont nuls pour ce probleme)
A=zeros(K,M,M);
AP=A;
B=A;
BP=A;
C=A;
CP=A;
D=zeros(K,M);
DP=D;

% conditions au limites (elle ne dependent pas du temps)
% initialisation a zeros (la plupart des valeurs sont nulles)
V0=zeros(2,M);
T0=zeros(2,M,3);
V1=V0;
T1=T0;
% condition en x=0
T0(:,1,2)=ones(2,1);         % condition sur la derivee premiere pour l'equation sur c 
T0(:,2,1)=ones(2,1);         % condition sur la valeur pour T
V0(:,2) = t0 .* ones(2,1);   % valeur de T au bord 

% condition en x=0
T1(:,1,2)=ones(2,1);         % condition sur la derivee premiere pour l'equation sur c 
T1(:,2,1)=ones(2,1);         % condition sur la valeur pour T
V1(:,2) = t1 .* ones(2,1);   % valeur de T au bord 

% calcul de dt
dt = (dx.^2 /max(alpha,xie))/6;

% equation 1 = concentration
%equation 2 = temperture
for l=1:L
	t(l+1)=t(l)+dt;
	
	% matrice A au temps l
	A(:,1,1)=alpha.*ones(1,K);      % derivee 2 de c equation sur c
	A(:,1,2)=alpha.*c(l,:);         % derivee 2 de c eqution sur T
	A(:,2,2)=xie.*ones(1,K);        % derivee 2 de T eqution sur T
	% matrice A au temps l+1
	AP(:,1,1)=alpha.*ones(1,K);
	AP(:,1,2)=alpha.*c(l+1,:);
	AP(:,2,2)=xie.*ones(1,K);
	
	% matrice B au temps l
	B(:,1,2)=alpha.*pdederive(x,c(l,:),2,2,2);  % derivee 1 de c equation sur T 
	% matrice B au temps l+1
	BP(:,1,2)=alpha.*pdederive(x,c(l+1,:),2,2,2);  
	
	% terme source pour la chaleur
	D(fix(K/2)+1,2)=q;
	
	% creation des matrice F
	F=cat(2,c(l,:)',T(l,:)');
	% toutes les eqautions sont predictives
	FP0=zeros(size(F));
	mode=zeros(2,1);
	
	% memoire pour la boucle de convergence
	FMEM=cat(2,c(l+1,:)',T(l+1,:)');
	
	
	% calcul explicite pour la premiere evaluation des coefficient non lineaire couple
   [FPE,FS,S,ALPHA,ALPHAP,IDENTITE]= ...
   pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP0,V0,T0,V1,T1,mode,1,dx,dt); 
   FP=FPE;
   % spy des matrice
   %figure(1)
   %spy(ALPHA)
   %figure(2)
   %spy(ALPHAP)
   %drawnow
	% boucle de convergence  
	n=5;
	while (sum(sum((FMEM-FP).^2))>1e-3)&(n>0)
		FMEM=FP;
		n=n-1;
		
      % recalcul des coeffcicients
 	   AP(:,1,2)=alpha.*FP(:,1)';
 	   BP(:,1,2)=alpha.*pdederive(x,FP(:,1)',1,1,2);  
 	   [FP,FS,S,ALPHA,ALPHAP,IDENTITE]= ...
      pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP0,V0,T0,V1,T1,mode,1/2,dx,dt); 
		
      % spy des matrice
      %figure(1)
      %spy(ALPHA)
      %figure(2)
      %spy(ALPHAP)
      %drawnow
   end     
   % copie des valeurs dans c et T
   c(l+1,:)=FP(:,1)';
   T(l+1,:)=FP(:,2)';
     
   % plot des profils au temps l, l+1  premiere et deuxieme approximation   
   if nargin ==0
	   figure(1);
	   hold off 
	   plot(x,c(l,:),'c:',x,c(l+1,:),'c-',x,T(l,:),'m:',x,T(l+1,:),'m-');
	   title(num2str(t(l+1)));
	   drawnow
   else
   	fprintf('.');
   	if rem(l,80)==0
   		fprintf('\n');
   	end
   end
end
