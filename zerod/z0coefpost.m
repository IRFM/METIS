% COEFPOST coefficient de rayonnement de POST
%-------------------------------------------------------------------------
% fichier coefpost.m ->  coefpost
%
%
% fonction Matlab 5 :
%
% Cette fonction retourne le coefficients "radiatives cooling rates"
% encore appeles les coefficients de Post pour les principales
% impuretees du plasma. Ces donnees sont valables si l'equilibre
% coronal est realise. 
% Les valeurs donn�s sont des points caracteristiques des courbes
% il faut les utiliser en log/log a travers des splines cubiques.
% la temperature electronique est en KeV.
% les densit�en cm^-3 !!!!
% les Lz sont en egr !!!!
% Z est la charge des ions completement ionise
% A est la nombre atomique des ions
%
% La densite volumique de puissance rayonnee est donnee par :
%
% Pz = Ne*Nz*Lz
%
% 1 erg = 1e-7 J 
%
% ref : Steady_state radiative cooling rates for low-density high-temperature plasmas,
%       D. E. Post et al, Atomic data and nuclear tables, 20, 397-439, 1977.
%
% syntaxe  :
%
%     [a,z,post]=coefpost;
%
% sorties :
%
%     a,z  =  tableaux des nombres atomique et des charges des elements
%     post =  tableaux des valeurs carateristiques des coefficients de Post
%             (au extremite de l'intervalle et les extrema)
%
%      a a(k) et z(k) correspond post(k).
%      post(k).te  = temperature electronique en keV
%      post(k).lz  = coefficient de post (erg * s^-1 * cm^-3)
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 25/10/2000.
%
%
% liste des modifications :
%
%  * ajout de Li
%  
%
%--------------------------------------------------------------
%
function [a,z,post]=z0coefpost(old)

persistent tabmat
% nouvelle donnees ADAS
if nargin > 0
  tabmat = [];
elseif isempty(tabmat)
  try 
    load('Lz_zave.mat')   
  end
end
%! il faut recalculer les points a 100 kev a partir du brem seul
k=1;
if ~isempty(tabmat)
    nb=32;
else
    nb=26;
end
sg.te=[];
sg.lz=[];
a=NaN.*ones(1,nb);
z=NaN.*ones(1,nb);
post(nb)=sg;
% pour He3 et He4
post(k).te = [1e-3, 3e-3,    7e-3,    2.5e-1,  1e2];
post(k).lz = [1e-24, 3.5e-23, 7.5e-21, 2.2e-23, 2.3e-22];
z(k)=2;
a(k)=3;
if isfield(tabmat,'He')
    post(k).te = tabmat.He.data(:,1)';
    post(k).lz = tabmat.He.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3, 3e-3,    7e-3,    2.5e-1,  1e2];
post(k).lz = [1e-24, 3.5e-23, 7.5e-21, 2.2e-23, 2.3e-22];
z(k)=2;
a(k)=4;
if isfield(tabmat,'He')
    post(k).te = tabmat.He.data(:,1)';
    post(k).lz = tabmat.He.data(:,2)' * 1e13;
end
k=k+1;

% pour Li6 et Li7
post(k).te = [1e-3,  5e-3,    7e-3,    2e-2,    5e-1,   1e2];
post(k).lz = [1e-24, 2.7e-23, 4.8e-23, 1.8e-21, 7e-23,  5e-22];
z(k)=3;
a(k)=6;
if isfield(tabmat,'Li')
    post(k).te = tabmat.Li.data(:,1)';
    post(k).lz = tabmat.Li.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3, 3e-3,    7e-3,    2.5e-1,  1e2];
post(k).lz = [1e-24, 3.5e-23, 7.5e-21, 2.2e-23, 2.3e-22];
z(k)=3;
a(k)=7;
if isfield(tabmat,'Li')
    post(k).te = tabmat.Li.data(:,1)';
    post(k).lz = tabmat.Li.data(:,2)' * 1e13;
end
k=k+1;

% pour Be8 et Be9impBe.Z= 4.;
%post(k).te = [1e-3,  2.5e-3, 5e-3, 1.3e-2, 4e-2,  1e0,     1e2];
%post(k).lz = [1e-24, 4e-19, 6e-22, 7e-23,  3e-21, 1.7e-22, 9.5e-22];
post(k).te =[1.000e+00,1.500e+00,2.000e+00,3.000e+00,5.000e+00,7.000e+00,1.000e+01,1.500e+01,2.000e+01,3.000e+01,5.000e+01,7.000e+01,1.000e+02,1.500e+02,2.000e+02,3.000e+02,5.000e+02,7.000e+02,1.000e+03,2.000e+03,5.000e+03,1.000e+04,2.000e+04]./1e3;
post(k).lz =[6.593e-27,2.167e-26,1.257e-26,1.008e-27,8.386e-29,2.625e-29,1.119e-29,6.329e-29,2.764e-28,4.857e-28,3.013e-28,1.453e-28,7.450e-29,4.054e-29,2.836e-29,1.873e-29,1.304e-29,1.154e-29,1.106e-29,1.243e-29,1.789e-29,2.484e-29,3.493e-29]./1e-7;
z(k)=4;
a(k)=8;
if isfield(tabmat,'Be')
    post(k).te = tabmat.Be.data(:,1)';
    post(k).lz = tabmat.Be.data(:,2)' * 1e13;
end
k=k+1;
%post(k).te = [1e-3,  2.5e-3, 5e-3, 1.3e-2, 4e-2,  1e0,     1e2];
%post(k).lz = [1e-24, 4e-19, 6e-22, 7e-23,  3e-21, 1.7e-22, 9.5e-22];
post(k).te =[1.000e+00,1.500e+00,2.000e+00,3.000e+00,5.000e+00,7.000e+00,1.000e+01,1.500e+01,2.000e+01,3.000e+01,5.000e+01,7.000e+01,1.000e+02,1.500e+02,2.000e+02,3.000e+02,5.000e+02,7.000e+02,1.000e+03,2.000e+03,5.000e+03,1.000e+04,2.000e+04]./1e3;
post(k).lz =[6.593e-27,2.167e-26,1.257e-26,1.008e-27,8.386e-29,2.625e-29,1.119e-29,6.329e-29,2.764e-28,4.857e-28,3.013e-28,1.453e-28,7.450e-29,4.054e-29,2.836e-29,1.873e-29,1.304e-29,1.154e-29,1.106e-29,1.243e-29,1.789e-29,2.484e-29,3.493e-29]./1e-7;
z(k)=4;
a(k)=9;
if isfield(tabmat,'Be')
    post(k).te = tabmat.Be.data(:,1)';
    post(k).lz = tabmat.Be.data(:,2)' * 1e13;
end
k=k+1;

% pour B10 et B11
post(k).te = [1e-3,  2e-3,  2e-2,  6.8e-2, 2e0,     1e2];
post(k).lz = [1e-24, 2e-19, 5e-22, 5e-21,  3.2e-22, 1.3e-21];
z(k)=5;
a(k)=10;
if isfield(tabmat,'B')
    post(k).te = tabmat.B.data(:,1)';
    post(k).lz = tabmat.B.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3,  2e-3,  2e-2,  6.8e-2, 2e0,     1e2];
post(k).lz = [1e-24, 2e-19, 5e-22, 5e-21,  3.2e-22, 1.3e-21];
z(k)=5;
a(k)=11;
if isfield(tabmat,'B')
    post(k).te = tabmat.B.data(:,1)';
    post(k).lz = tabmat.B.data(:,2)' * 1e13;
end
k=k+1;

% pour C12 C13 C14
post(k).te = [1e-3,  2.5e-3,  3e-3,  7e-3,  4e-2,    1e-1,    2e0,     1e2];
post(k).lz = [1e-24, 3.5e-20, 1e-19, 8e-19, 1.1e-21, 8.3e-21, 5.3e-22, 2.2e-21];
z(k)=6;
a(k)=12;
if isfield(tabmat,'C')
    post(k).te = tabmat.C.data(:,1)';
    post(k).lz = tabmat.C.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3,  2.5e-3,  3e-3,  7e-3,  4e-2,    1e-1,    2e0,     1e2];
post(k).lz = [1e-24, 3.5e-20, 1e-19, 8e-19, 1.1e-21, 8.3e-21, 5.3e-22, 2.2e-21];
z(k)=6;
a(k)=13;
if isfield(tabmat,'C')
    post(k).te = tabmat.C.data(:,1)';
    post(k).lz = tabmat.C.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3,  2.5e-3,  3e-3,  7e-3,  4e-2,    1e-1,    2e0,     1e2];
post(k).lz = [1e-24, 3.5e-20, 1e-19, 8e-19, 1.1e-21, 8.3e-21, 5.3e-22, 2.2e-21];
z(k)=6;
a(k)=14;
if isfield(tabmat,'C')
    post(k).te = tabmat.C.data(:,1)';
    post(k).lz = tabmat.C.data(:,2)' * 1e13;
end
k=k+1;

% pour N14 
post(k).te = [1e-3,  2e-3,    4e-3,    1e-2,  6e-2,    2e-1   2.3e0,   1e2];
post(k).lz = [1e-24, 1.2e-20, 7.5e-20, 9e-19, 2.2e-21, 1e-20, 8.5e-22, 3e-21];
z(k)=7;
a(k)=14;
if isfield(tabmat,'N')
    post(k).te = tabmat.N.data(:,1)';
    post(k).lz = tabmat.N.data(:,2)' * 1e13;
end
k=k+1;

% pour O16 O17 o18
post(k).te = [1e-3,  2.5e-3,  5e-3 ,   2e-2,  1e-1,    2e-1,    3e0,     1e2];
post(k).lz = [1e-24, 2.8e-21, 8.5e-20, 9e-19, 3.5e-21, 1.3e-20, 1.3e-21, 3.9e-21];
z(k)=8;
a(k)=16;
if isfield(tabmat,'O')
    post(k).te = tabmat.O.data(:,1)';
    post(k).lz = tabmat.O.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3,  2.5e-3,  5e-3 ,   2e-2,  1e-1,    2e-1,    3e0,     1e2];
post(k).lz = [1e-24, 2.8e-21, 8.5e-20, 9e-19, 3.5e-21, 1.3e-20, 1.3e-21, 3.9e-21];
z(k)=8;
a(k)=17;
if isfield(tabmat,'O')
    post(k).te = tabmat.O.data(:,1)';
    post(k).lz = tabmat.O.data(:,2)' * 1e13;
end
k=k+1;
post(k).te = [1e-3,  2.5e-3,  5e-3 ,   2e-2,  1e-1,    2e-1,    3e0,     1e2];
post(k).lz = [1e-24, 2.8e-21, 8.5e-20, 9e-19, 3.5e-21, 1.3e-20, 1.3e-21, 3.9e-21];
z(k)=8;
a(k)=18;
if isfield(tabmat,'O')
    post(k).te = tabmat.O.data(:,1)';
    post(k).lz = tabmat.O.data(:,2)' * 1e13;
end
k=k+1;

% pour Ne
post(k).te = [1e-3,  4e-3,    7e-3,  4e-2,  2e-1,    4e-1,    5e0,     1e2];
post(k).lz = [1e-24, 1.8e-21, 6e-20, 8e-19, 6.5e-21, 1.8e-20, 2.7e-21, 6.1e-21];
z(k)=10;
a(k)=20;
if isfield(tabmat,'Ne')
    post(k).te = tabmat.Ne.data(:,1)';
    post(k).lz = tabmat.Ne.data(:,2)' * 1e13;
end
k=k+1;

% pour Si
post(k).te = [1e-3,  2e-2,  1.3e-1,  5e-1,    9e-1,  1e1,   1e2];
post(k).lz = [1e-24, 4e-21, 6.5e-19, 1.5e-20, 3e-20, 7e-21, 1.3e-20];
if isfield(tabmat,'Si')
    post(k).te = tabmat.Si.data(:,1)';
    post(k).lz = tabmat.Si.data(:,2)' * 1e13;
end
z(k)=14;
a(k)=28;
k=k+1;

% pour Ar
%post(k).te = [1e-3,  3e-2,    5.5e-2, 2.5e-1, 1e0,     2e0,     1.5e1    1e2];
%post(k).lz = [1e-24, 1.1e-18, 9e-20,  5e-19,  2.3e-20, 3.8e-20, 1.3e-20, 2.3e-20];
post(k).te = [2.164e+00,3.430e+00,5.437e+00,8.617e+00,1.366e+01,2.164e+01,3.430e+01,5.437e+01,8.617e+01,1.366e+02,2.164e+02,3.430e+02,5.437e+02,8.617e+02,1.366e+03,2.164e+03,3.430e+03,5.437e+03,8.617e+03,1.366e+04,2.164e+04,3.430e+04,5.437e+04,8.617e+04]./1e3;
post(k).lz = [3.020e-29,5.623e-28,4.786e-28,1.738e-26,4.571e-26,8.511e-26,7.762e-26,1.318e-26,1.995e-26,3.715e-26,3.631e-26,2.455e-26,1.660e-26,5.623e-27,4.074e-27,4.074e-27,3.090e-27,2.089e-27,1.585e-27,1.445e-27,1.479e-27,1.585e-27,1.778e-27,2.089e-27]./1e-7;
z(k)=18;
a(k)=40;
if isfield(tabmat,'Ar')
    post(k).te = tabmat.Ar.data(:,1)';
    post(k).lz = tabmat.Ar.data(:,2)' * 1e13;
end
k=k+1;


% pour Ti
post(k).te = [1e-3,  2e-2,    2e-1,  5e-1,  2e0,     4e0,     2.8e1,   1e2];
post(k).lz = [1e-24, 2.5e-18, 3e-19, 5e-19, 3.2e-20, 4.5e-20, 2.8e-20, 3.5e-20];
z(k)=22;
a(k)=48;
if isfield(tabmat,'Ti')
    post(k).te = tabmat.Ti.data(:,1)';
    post(k).lz = tabmat.Ti.data(:,2)' * 1e13;
end
k=k+1;

% pour V
post(k).te = [1e-3,  2e-2     2.5e-1   5.2e-1, 2e0,     4.8e0, 3e1,   2e2];
post(k).lz = [1e-24, 2.8e-18, 3.3e-19, 5e-19,  3.8e-20, 5e-20, 3e-20, 4.5e-20];
z(k)=23;
a(k)=51;
if isfield(tabmat,'V')
    post(k).te = tabmat.V.data(:,1)';
    post(k).lz = tabmat.V.data(:,2)' * 1e13;
end
k=k+1;

% pour Cr
post(k).te = [1e-3,  2e-2,    4e-2,    3.5e-1,  6e-1,  2.8e0, 6e0,   3e1,     1e2];
post(k).lz = [1e-24, 2.1e-18, 3.2e-18, 3.8e-19, 5e-19, 4e-20, 5e-20, 3.3e-20, 4.3e-20];
z(k)=24;
a(k)=52;
if isfield(tabmat,'Cr')
    post(k).te = tabmat.Cr.data(:,1)';
    post(k).lz = tabmat.Cr.data(:,2)' * 1e13;
end
k=k+1;

% pour Fe
post(k).te = [1e-3,  2e-2,  7e-2,    2.2e-1, 5e-1,  8e-1,    3e0,     8e0,   3.5e1,   1e2];
post(k).lz = [1e-24, 1e-18, 3.5e-18, 1e-18,  4e-19, 5.1e-19, 4.5e-20, 6e-20, 4.3e-20, 5.1e-20];
z(k)=26;
a(k)=56;
if isfield(tabmat,'Fe')
    post(k).te = tabmat.Fe.data(:,1)';
    post(k).lz = tabmat.Fe.data(:,2)' * 1e13;
end
k=k+1;

% pour Ni
post(k).te = [1e-3,  3e-2,  1e-1,  8e-1,    1e0,   4e0,   1e1,   3e1,   1e2];
post(k).lz = [1e-24, 1e-18, 3e-18, 4.3e-19, 5e-19, 5e-20, 6e-20, 5e-20, 6e-20];
z(k)=28;
a(k)=59;
if isfield(tabmat,'Ni')
    post(k).te = tabmat.Ni.data(:,1)';
    post(k).lz = tabmat.Ni.data(:,2)' * 1e13;
end
k=k+1;

% pour Cu
post(k).te = [1e-3,  3e-2,  1.3e-1, 8.5e-1,  1.1e0, 5e0,   1e1,     3.5e1,   1e2];
post(k).lz = [1e-24, 7e-19, 3e-18,  4.5e-19, 5e-19, 5e-20, 6.7e-20, 5.5e-20, 7e-20];
z(k)=29;
a(k)=64;
if isfield(tabmat,'Cu')
    post(k).te = tabmat.Cu.data(:,1)';
    post(k).lz = tabmat.Cu.data(:,2)' * 1e13;
end
k=k+1;

% pour Cs
post(k).te = [1e-3,  1e-1,    9e-1,    2e0,   7e0,   3e1,   1e2];
post(k).lz = [1e-24, 7.5e-18, 1.8e-18, 2e-18, 5e-19, 2e-19, 2.6e-19];
z(k)=55;
a(k)=133;
if isfield(tabmat,'Cs')
    post(k).te = tabmat.Cs.data(:,1)';
    post(k).lz = tabmat.Cs.data(:,2)' * 1e13;
end
k=k+1;

% pour W
post(k).te = [1e-3,  1e-1,  2e-1,    1e0,   3e0,     7e0,     1.9e1,   3e1,     6e1,     1e2];
post(k).lz = [1e-24, 3e-18, 1.6e-18, 6e-18, 1.9e-18, 1.9e-18, 5.3e-19, 5.3e-19, 4.3e-19, 4.7e-19];
z(k)=74;
a(k)=184;


% nouvelle valeur pour W
% ref : T. Putterich et al, NF 50 (2010) p 025012
post(k).te = [
1
30
40
50
60
70
100
150
200
300
400
500
600
800
1000
1500
2000
2300
2700
3000
3500
4000
5000
6000
7000
10000
12000
15000
20000
25000
30000
40000
100000]' ./ 1e3;

post(k).lz =  [
1e-37
2.32E-31
2.56E-31
2.55E-31
2.08E-31
1.73E-31
1.67E-31
1.57E-31
1.65E-31
1.99E-31
2.42E-31
2.91E-31
3.23E-31
3.79E-31
4.21E-31
4.59E-31
4.38E-31
3.94E-31
2.88E-31
2.28E-31
1.96E-31
1.88E-31
1.87E-31
1.79E-31
1.67E-31
1.33E-31
1.15E-31
9.47E-32
7.30E-32
6.11E-32
5.47E-32
4.95E-32 
4.95E-32
]' .* 1e6 ./ 1e-7;

if isfield(tabmat,'W')
    post(k).te = tabmat.W.data(:,1)';
    post(k).lz = tabmat.W.data(:,2)' * 1e13;
end

k=k+1;

% Kr if available
if isfield(tabmat,'Kr')
    post(k).te = tabmat.Kr.data(:,1)';
    post(k).lz = tabmat.Kr.data(:,2)' * 1e13;
    z(k)=36;
    a(k)=83.8;
    k=k+1;
    
end
% Xenon if available
if isfield(tabmat,'Xe')
    post(k).te = tabmat.Xe.data(:,1)';
    post(k).lz = tabmat.Xe.data(:,2)' * 1e13;
    z(k)=54;
    a(k)=131.3;
    k=k+1;
end

% Chlorine if available
if isfield(tabmat,'Cl')
    post(k).te = tabmat.Cl.data(:,1)';
    post(k).lz = tabmat.Cl.data(:,2)' * 1e13;
    z(k)=17;
    a(k)=34.453;
    k=k+1;
end

% HDT if available
if isfield(tabmat,'H')
    post(k).te = tabmat.H.data(:,1)';
    post(k).lz = tabmat.H.data(:,2)' * 1e13;
    z(k)=1;
    a(k)=1;
    k=k+1;
    post(k).te = tabmat.H.data(:,1)';
    post(k).lz = tabmat.H.data(:,2)' * 1e13;
    z(k)=1;
    a(k)=2;
    k=k+1;
    post(k).te = tabmat.H.data(:,1)';
    post(k).lz = tabmat.H.data(:,2)' * 1e13;
    z(k)=1;
    a(k)=3;
    k=k+1;
end

% complete with additional tabmat data
%zadded = [];
if ~isempty(tabmat)
   noms = fieldnames(tabmat);
   for k=1:length(noms)
       if ~any(tabmat.(noms{k}).Z == z)
           %zadded(end+1) = tabmat.(noms{k}).Z;
           z(end+1) = tabmat.(noms{k}).Z;
           a(end+1) = tabmat.(noms{k}).A;
           post(end+1).te = tabmat.(noms{k}).data(:,1)';
           post(end).lz = tabmat.(noms{k}).data(:,2)' * 1e13;          
       end          
   end
   %disp(zadded)
end

for k = 1:length(post)
    indok      = find((post(k).te > 1e-4) & (post(k).te <= 1e2));
    post(k).te = post(k).te(indok);
    post(k).lz = post(k).lz(indok);
end    
    


    