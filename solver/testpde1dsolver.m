% script de test du solveur
pwdmem=pwd;
cdt=which('pde1dsolver');
cdt=cdt(1:(max(findstr(cdt,'/'))-1));
cd(cdt);
!rm test_pde1dsolver.txt
diary test_pde1dsolver.txt

disp(time)
disp(pwd)

[s,t]=unix('ls -alg');
disp(t);

disp('Script de test automatique du solveur de PDE 1D')
disp(' ')
disp('-------------------------------------------------------------------')
disp('Premier test : test informatique sur les dimensions des matrices');
for k=1:10
	fprintf('Passage # %d\n',k);
	testdim_pde1dsolver;
end
disp('Le test des dimensions a ete effectue avec succes');
disp(' ')
disp('-------------------------------------------------------------------')
disp('Test des operateurs :');
[A1t,AP1t,A2t,AP2t,A3t,AP3t,A4t,AP4t,A5t,AP5t,A6t,AP6t,A7t,AP7t]=test0_pde1dsolver;
load reference1
if all(all(A1t==A1))
	disp('L''operateur derivee 2 explicite ok');
else
	disp('Probleme sur l''operateur derivee 2 explicite');
end
if all(all(AP1t==AP1))
	disp('L''operateur derivee 2 implicite ok');
else
	disp('Probleme sur l''operateur derivee 2 implicite');
end
if all(all(A2t==A2))
	disp('L''operateur derivee 1 explicite ok');
else
	disp('Probleme sur l''operateur derivee 1 explicite');
end
if all(all(AP2t==AP2))
	disp('L''operateur derivee 1 implicite ok');
else
	disp('Probleme sur l''operateur derivee 1 implicite');
end
if all(all(A3t==A3))
	disp('L''operateur * explicite ok');
else
	disp('Probleme sur l''operateur  * explicite');
end
if all(all(AP3t==AP3))
	disp('L''operateur * implicite ok');
else
	disp('Probleme sur l''operateur  * implicite');
end

if all(all(A5t==A5))
	disp('L''operateur derivee 2 explicite croise pour 2 equations ok');
else
	disp('Probleme sur l''operateur derivee 2 explicite croise pour 2 equations ');
end
if all(all(AP5t==AP5))
	disp('L''operateur derivee 2 implicite  croise pour 2 equations  ok');
else
	disp('Probleme sur l''operateur derivee 2 implicite croise pour 2 equations ');
end
if all(all(A6t==A6))
	disp('L''operateur derivee 1 explicite  croise pour 2 equations  ok');
else
	disp('Probleme sur l''operateur derivee 1 explicite croise pour 2 equations ');
end
if all(all(AP6t==AP6))
	disp('L''operateur derivee 1 implicite  croise pour 2 equations  ok');
else
	disp('Probleme sur l''operateur derivee 1 implicite croise pour 2 equations ');
end
if all(all(A7t==A7))
	disp('L''operateur * explicite croise pour 2 equations  ok ');
else
	disp('Probleme sur l''operateur  * explicite croise pour 2 equations ');
end
if all(all(AP7t==AP7))
	disp('L''operateur * implicite croise pour 2 equations  ok ');
else
	disp('Probleme sur l''operateur  * implicite croise pour 2 equations ');
end
disp(' ')
disp('-------------------------------------------------------------------')
disp('Test des conditions aux limites -> valeur donnees')
load reference2
Tt=test1_pde1dsolver(1);
fprintf('\n');
if all(all(Tt==T))
	disp('Le test des conditions aux limites, valeurs donnees, est ok');
else
	disp('Probleme lors du test des conditions aux limites, valeurs donnees');
end
disp(' ')
disp('-------------------------------------------------------------------')
disp('Test des conditions aux limites -> derivees donnees')
load reference3
Tt=test2_pde1dsolver(1);
fprintf('\n');
if all(all(Tt==T))
	disp('Le test des conditions aux limites, derivees donnees, est ok');
else
	disp('Probleme lors du test des conditions aux limites, derivees donnees');
end
disp(' ')
disp('-------------------------------------------------------------------')
disp('Test des conditions aux limites -> derivees 2 donnees')
load reference4
Tt=test3_pde1dsolver(1);
fprintf('\n');
if all(all(Tt==T))
	disp('Le test des conditions aux limites, derivees 2 donnees, est ok');
else
	disp('Probleme lors du test des conditions aux limites, derivees  2 donnees');
end
disp(' ')
disp('-------------------------------------------------------------------')
disp('Test de 2 equations non lineaires couplees')
load reference5
[Tt,ct]=test2eq_pde1dsolver(1);
fprintf('\n');
if all(all(Tt==T)) & all(all(ct==c))
	disp('Le test des equations non lineaires couplees est ok');
else
	disp('Probleme lors du test des equations non lineaires couplees');
end

disp(' ')
disp('-------------------------------------------------------------------')
disp('Test intpretatif/predictif')
[Tm,cm]=test2eq_pde1dsolver(c);
fprintf('\n');
if all(all(Tm==Tt)) & all(all(cm==ct))
	disp('Le test intpretatif/predictif  est ok');
else
	disp('Probleme lors du test intpretatif/predictif');
end
disp(' ')
disp('-------------------------------------------------------------------')
disp('Fin')


diary off
cd(pwdmem);