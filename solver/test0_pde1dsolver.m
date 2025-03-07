function [A1,AP1,A2,AP2,A3,AP3,A4,AP4,A5,AP5,A6,AP6,A7,AP7]=test0_pde1dsolver

A1=[];AP1=[];A2=[];AP2=[];A3=[];AP3=[];A4=[];AP4=[];A5=[];AP5=[];A6=[];AP6=[];A7=[];AP7=[];

% script de test #1 du solver : test des operateurs
disp('Test de la fonction pde1dsolver')
disp('Premier test : test des operateurs')
dbstop if error
disp(' test pour une equation seule')

K = 5; % nombre de points d'espaces
M = 1;  % nombre d'equations

% initialisation des matrices a zero (la plupart des elements sont nuls pour ce probleme)
A=zeros(K,M,M);
AP=A;
B=A;
BP=A;
C=A;
CP=A;
D=zeros(K,M);
DP=D;
F=D;
FP=D;
% conditions au limites (elle ne dependent pas du temps)
% initialisation a zeros (la plupart des valeurs sont nulles)
V0=zeros(2,M);
T0=zeros(2,M,3);
V1=V0;
T1=T0;
	
disp('operateur derivee 2')
A=ones(K,M,M);
AP=A;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A1=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP1=ALPHAP;
disp('operateur derivee 1')
A=0.*A;
AP=A;
B=2.*ones(K,M,M);
BP=B;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A2=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP2=ALPHAP;
disp('operateur * ')
B=0.*B;
BP=B;
C=ones(K,M,M);
CP=C;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A3=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP3=ALPHAP;

K = 5; % nombre de points d'espaces
M = 2;  % nombre d'equations

% initialisation des matrices a zero (la plupart des elements sont nuls pour ce probleme)
A=zeros(K,M,M);
AP=A;
B=A;
BP=A;
C=A;
CP=A;
D=zeros(K,M);
DP=D;
F=D;
FP=D;
% conditions au limites (elle ne dependent pas du temps)
% initialisation a zeros (la plupart des valeurs sont nulles)
V0=zeros(2,M);
T0=zeros(2,M,3);
V1=V0;
T1=T0;
	
disp('operateur derivee 2 croisees')
A(:,1,2)=ones(K,1,1);
A(:,2,1)=ones(K,1,1);
AP=A;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A5=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP5=ALPHAP;
disp('operateur derivee 1 croisee ')
A=0.*A;
AP=A;
B(:,1,2)=2.*ones(K,1,1);
B(:,2,1)=2.*ones(K,1,1);
BP=B;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A6=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP6=ALPHAP;
disp('operateur * croise')
B=0.*B;
BP=B;
C(:,1,2)=ones(K,1,1);
C(:,2,1)=ones(K,1,1);
CP=C;
% calcul explicite 
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,1,1,1); 
disp(full(ALPHA))
A7=ALPHA;
% calcul implicite
[FP,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,0,0,1,1); 
disp(full(ALPHAP))
AP7=ALPHAP;
