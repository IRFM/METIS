% function [Z,A,s] = interpret_symb(symb);
%
% Traduction des symboles de Mendeleev en numéro atomique et nombre de masse
%
% R. Guirlet, 31/10/2003
%
%
function [Z,A,s] = interpret_symb(symb);
%
if iscellstr(symb)
	symb = char(symb);
end
if length(symb)>1 & symb(2)==' '
	symb = symb(1);
end
%
if strcmp(symb,'H') | strcmp(symb,'h')
	Z = 1;
	A = 1;
	s = 'H';
elseif strcmp(symb,'D') | strcmp(symb,'d')
	Z = 1;
	A = 2;
	s = 'D';
elseif strcmp(symb,'T') | strcmp(symb,'t')
	Z = 1;
	A = 3;
	s = 'T';
elseif symb=='HE' | symb=='he' | symb=='He'
	Z = 2;
	A = 4;
	s = 'He';
elseif symb=='LI' | symb=='li' | symb=='Li'
	Z = 3;
	A = 7;
	s = 'Li';
elseif symb=='BE' | symb=='be' | symb=='Be'
	Z = 4;
	A = 9;
	s = 'Be';
elseif symb=='B' | symb=='b' | symb == 'BO' | symb == 'bo' | symb == 'Bo'
	Z = 5;
	A = 11;
	s = 'B';
elseif symb=='C' | symb=='c' | symb=='CA' | symb=='ca' | symb=='Ca'
	Z = 6;
	A = 12;
	s = 'C';
elseif symb=='N' | symb=='n' | symb == 'AZ' | symb == 'az' | symb == 'Az'
	Z = 7;
	A = 14;
	s = 'N';
elseif symb=='O' | symb=='o' | symb=='OX' | symb=='ox' | symb=='Ox'
	Z = 8;
	A = 16;
	s = 'O';
elseif symb=='F' | symb=='f' | symb=='FL' | symb=='fl' | symb=='Fl'
	Z = 9;
	A = 19;
	s = 'F';
elseif symb=='NE' | symb=='ne' | symb=='Ne'
	Z = 10;
	A = 20;
	s = 'Ne';
elseif symb=='NA' | symb=='na' | symb=='Na'
	Z = 11;
	A = 23;
	s = 'Na';
elseif symb=='MG' | symb=='mg' | symb=='Mg'
	Z = 12;
	A = 24;
	s = 'Mg';
elseif symb=='AL' | symb=='al' | symb=='Al'
	Z = 13;
	A = 27;
	s = 'Al';
elseif symb=='SI' | symb=='si' | symb=='Si'
	Z = 14;
	A = 28;
	s = 'Si';
elseif symb=='P' | symb=='p'
	Z = 15;
	A = 31;
	s = 'P';
elseif symb=='S' | symb=='s'
	Z = 16;
	A = 32;
	s = 'S';
elseif symb=='CL' | symb=='cl' | symb=='Cl'
	Z = 17;
	A = 35;
	s = 'Cl';
elseif symb=='AR' | symb=='ar' | symb=='Ar'
	Z = 18;
	A = 40;
	s = 'Ar';
elseif symb=='K' | symb=='k'
	Z = 19;
	A = 39;
	s = 'K';
elseif symb=='CA' | symb=='ca' | symb=='Ca'
	Z = 20;
	A = 40;
	s = 'Ca';
elseif symb=='SC' | symb=='sc' | symb=='Sc'
	Z = 21;
	A = 445;
	s = 'Sc';
elseif symb=='TI' | symb=='ti' | symb=='Ti'
	Z = 22;
	A = 48;
	s = 'Ti';
elseif symb=='V' | symb=='v'
	Z = 23;
	A = 51;
	s = 'V';
elseif symb=='CR' | symb=='cr' | symb=='Cr'
	Z = 24;
	A = 52;
	s = 'Cr';
elseif symb=='MN' | symb=='mn' | symb=='Mn'
	Z = 25;
	A = 55;
	s = 'Mn';
elseif symb=='FE' | symb=='fe' | symb=='Fe'
	Z = 26;
	A = 56;
	s = 'Fe';
elseif symb=='CO' | symb=='co' | symb=='Co'
	Z = 27;
	A = 59;
	s = 'Co';
elseif symb=='NI' | symb=='ni' | symb=='Ni'
	Z = 28;
	A = 59;
	s = 'Ni';
elseif symb=='CU' | symb=='cu' | symb=='Cu'
	Z = 29;
	A = 64;
	s = 'Cu';
elseif symb=='ZN' | symb=='zn' | symb=='Zn'
	Z = 30;
	A = 65;
	s = 'Zn';
elseif symb=='GA' | symb=='ga' | symb=='Ga'
	Z = 31;
	A = 70;
	s = 'Ga';
elseif symb=='GE' | symb=='ge' | symb=='Ge'
	Z = 32;
	A = 73;
	s = 'Ge';
elseif symb=='AS' | symb=='as' | symb=='As'
	Z = 33;
	A = 75;
	s = 'As';
elseif symb=='SE' | symb=='se' | symb=='Se'
	Z = 34;
	A = 79;
	s = 'Se';
elseif symb=='BR' | symb=='br' | symb=='Br'
	Z = 35;
	A = 80;
	s = 'Br';
elseif symb=='KR' | symb=='kr' | symb=='Kr'
	Z = 36;
	A = 84;
	s = 'Kr';
elseif symb=='ZR' | symb=='zr' | symb=='Zr'
	Z = 40;
	A = 91;
	s = 'Zr';
elseif symb=='MO' | symb=='mo' | symb=='Mo'
	Z = 42;
	A = 96;
	s = 'Mo';
elseif symb=='PD' | symb=='pd' | symb=='Pd'
	Z = 46;
	A = 106;
	s = 'Pd';
elseif symb=='AG' | symb=='ag' | symb=='Ag'
	Z = 47;
	A = 108;
	s = 'Ag';
elseif symb=='XE' | symb=='xe' | symb=='Xe'
	Z = 54;
	A = 131;
	s = 'Xe';
elseif symb=='W' | symb=='w'
	Z = 74;
	A = 184;
	s = 'W';
else
	disp(['interpret_symb.m - the symbol ' symb ' cannot be interpreted by this routine'])
	a = input('Please enter the atomic number and mass number [Z A] : ');
	Z = min(a); A = max(a);
end
