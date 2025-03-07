function T=test3_pde1dsolver(auto)

% sript de test #1 du solver : test des dimensions
disp('Test de la fonction pde1dsolver')
disp('test une equation')
disp('Valeurs aux bornes de la derivee 2 donnees')

dbstop if error
K = 21; % nombre de points d'espaces
M = 1;  % nombre d'equations
L=200;  % nombre de temps

if nargin ==0
	xie=input('Le coefficient de diffusion de la chaleur dans le liquide ? ');
	t0=input('la derivee 2 de la temperature en x=0 ? ');
	t1=input('la derivee 2 de la temperature en x=1 ? ');
	q=input('le flux de chaleur au centre du tube ? ');
	f=input('f ? ');
else
	xie=0.01;
	t0=-1;
	t1=1;
	q=pi;
	f=0.487;
end


t=zeros(L+1,1);                 % vecteur des temps 
x=linspace(0,1,K);              % vecteur d'expace 0..1
dx=mean(diff(x));               % pas d'espace
T=zeros(L+1,K);                 % matrice des temperatures
TD1 = zeros(L+1,K);
TD2 = zeros(L+1,K);

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
T0(:,1,3)=ones(2,1);         % condition sur la derivee 2 de   T
V0(:,1) = t0 .* ones(2,1);   % valeur de T au bord 

% condition en x=0
T1(:,1,3)=ones(2,1);         % condition sur la derivee 2 de T
V1(:,1) = t1 .* ones(2,1);   % valeur de T au bord 

% calcul de dt
dt = (dx.^2 /xie)/6;

% equation 1 = concentration
%equation 2 = temperture
for l=1:L
	t(l+1)=t(l)+dt;
	
	% matrice A au temps l
	A(:,1,1)=xie.*ones(1,K);        % derivee 2 de T eqution sur T
	% matrice A au temps l+1
	AP(:,1,1)=xie.*ones(1,K);
	
	% terme source pour la chaleur
	D(fix(K/2)+1,1)=q;
	DP(fix(K/2)+1,1)=q;
	
	% creation des matrice F
	F=T(l,:)';
	% toutes les eqautions sont predictives
	FP=zeros(size(F));
	mode=0;
	
  	[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= ...
      pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt); 
		
   % copie des valeurs dans c et T
   T(l+1,:)=FP(:,1)';
   TD1(l+1,:)=dfpdx(:,1)';
   TD2(l+1,:)=d2fpdx2(:,1)';
  
   % plot des profils au temps l, l+1  premiere et deuxieme approximation   
   if nargin ==0
	   figure(1);
	   hold off 
	   plot(x,T(l,:),'m:',x,T(l+1,:),'m-');
	   hold on
   	td=TD1(l,:);
   	tdp=TD1(l+1,:);
	   plot(x,td,'c:',x,tdp,'c-');
   	td=TD2(l,:);
   	tdp=TD2(l+1,:);
	   plot(x,td,'g:',x,tdp,'g-');
	   title(num2str(t(l+1)));
	   pause(1)
   else
   	fprintf('.');
   	if rem(l,80)==0
   		fprintf('\n');
   	end
   end
end
