% function s = interpret_symb(Z);
%
% Traduction des numéros atomiques en symboles de Mendeleev 
%
% R. Guirlet, 5/12/2006
%
%
function s = interpret_symb(Z);
%
if iscellstr(Z) | ischar(Z)
	disp('interpret_Z - Z should be an integer - Stop')
	return
end
%
if Z == 1;
	A = 1;
	s = 'H';
elseif Z == 2;
	A = 4;
	s = 'He';
elseif Z == 3;
	A = 6;
	s = 'Li';
elseif Z == 4;
	A = 9
	s = 'Be';
elseif Z == 5;
	A = 11;
	s = 'B';
elseif Z == 6;
	A = 12;
	s = 'C';
elseif Z == 7;
	A = 14;
	s = 'N';
elseif Z == 8;
	A = 16;
	s = 'O';
elseif Z == 9;
	A = 19;
	s = 'F';
elseif Z == 10;
	A = 20;
	s = 'Ne';
elseif Z == 11;
	A = 23;
	s = 'Na';
elseif Z == 12;
	A = 24;
	s = 'Mg';
elseif Z == 13;
	A = 27;
	s = 'Al';
elseif Z == 14;
	A = 28;
	s = 'Si';
elseif Z == 15;
	A = 31;
	s = 'P';
elseif Z == 16;
	A = 32;
	s = 'S';
elseif Z == 17;
	A = 35;
	s = 'Cl';
elseif Z == 18;
	A = 40;
	s = 'Ar';
elseif Z == 19;
	A = 39;
	s = 'K';
elseif Z == 20;
	A = 40;
	s = 'Ca';
elseif Z == 21;
	A = 445;
	s = 'Sc';
elseif Z == 22;
	A = 48;
	s = 'Ti';
elseif Z == 23;
	A = 51;
	s = 'V';
elseif Z == 24;
	A = 52;
	s = 'Cr';
elseif Z == 25;
	A = 55;
	s = 'Mn';
elseif Z == 26;
	A = 56;
	s = 'Fe';
elseif Z == 27;
	A = 59;
	s = 'Co';
elseif Z == 28;
	A = 59;
	s = 'Ni';
elseif Z == 29;
	A = 64;
	s = 'Cu';
elseif Z == 30;
	A = 65;
	s = 'Zn';
elseif Z == 31;
	A = 70;
	s = 'Ga';
elseif Z == 32;
	A = 73;
	s = 'Ge';
elseif Z == 33;
	A = 75;
	s = 'As';
elseif Z == 34;
	A = 79;
	s = 'Se';
elseif Z == 35;
	A = 80;
	s = 'Br';
elseif Z == 36;
	A = 84;
	s = 'Kr';
elseif Z == 40;
	A = 91;
	s = 'Zr';
elseif Z == 42;
	A = 96;
	s = 'Mo';
elseif Z == 46;
	A = 106;
	s = 'Pd';
elseif Z == 47;
	A = 108;
	s = 'Ag';
elseif Z == 54;
	A = 131;
	s = 'Xe';
elseif Z == 74;
	A = 184;
	s = 'W';
else
	disp(['interpret_symb.m - the symbol ' symb ' cannot be interpreted by this routine'])
	a = input('Please enter the atomic number and mass number [Z A] : ');
	Z == min(a); A = max(a);
end
