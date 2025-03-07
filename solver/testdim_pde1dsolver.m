function testdim_pde1dsolver(auto)
% sript de test #1 du solver : test des dimensions
disp('Test de la fonction pde1dsolver')
disp('test des dimensions')
dbstop if error
K = 5+fix(20*rand(1));
M = 1+fix(5*rand(1));
A=randn(K,M,M).*(rand(K,M,M)>0.5);
AP=randn(K,M,M).*(rand(K,M,M)>0.5);
B=randn(K,M,M).*(rand(K,M,M)>0.5);
BP=randn(K,M,M).*(rand(K,M,M)>0.5);
C=randn(K,M,M).*(rand(K,M,M)>0.5);
CP=randn(K,M,M).*(rand(K,M,M)>0.5);
D=randn(K,M).*(rand(K,M)>0.5);
DP=randn(K,M).*(rand(K,M)>0.5);
F=randn(K,M).*(rand(K,M)>0.5);
FP=randn(K,M).*(rand(K,M)>0.5);
V0=randn(2,M).*(rand(2,M)>0.5);
V1=randn(2,M).*(rand(2,M)>0.5);
T0=randn(2,M,3).*(rand(2,M,3)>0.5);
T1=randn(2,M,3).*(rand(2,M,3)>0.5);
mode=(rand(M,1)>0.75);
f=rand(1);
dx=rand(1);
dt=rand(1);

[fp,dfpdx,d2fpdx2,FS,S,ALPHA,ALPHAP,IDENTITE]= ...
pde1dsolver(A,B,C,D,AP,BP,CP,DP,F,FP,V0,T0,V1,T1,mode,f,dx,dt);

fprintf('K = %d, M = %d, f = %f\n',K,M,f);
fprintf('dx = %f, dt = %f \n',dx,dt);
disp('mode :')
disp(mode)
disp(' ');
whos
disp(' ');
if nargin ==0
	
	figure(1)
	spy(ALPHA)
	title('Alpha')

	figure(2)
	spy(ALPHAP)
	title('Alphap')

	figure(3)
	spy(IDENTITE)
	title('identite')

	figure(4)
	plot(FS,'or')
	hold on
	plot(S,'+g')
	
	pause(1);
end
