function [A,Z,name] = chargemasse(element)

liste = {
{'Ag',47,107.8682, 1.21},
{'Al',13,26.98154, 1.09},
{'Au',79,196.96654},
{'B',  5,10.81},
{'Be', 4,9.01},
{'C', 6,12.011},
{'Cd',48,112.41},
{'Co',27,58.93},
{'Cr',24,51.996},
{'Cu',29,63.546},
{'Fe',26,55.845},
{'Ge',32,72.61},
{'Hf',72,178.49},
{'Ir',77,192.217},
{'Mn',25,54.94},
{'Mo',42,95.94},
{'Nb',41,92.90638},
{'Ni',28,58.6934},
{'Os',76,190.2},
{'Pb',82,207.2},
{'Pd',46,106.42},
{'Pt',78,195.08},
{'Re',75,186.21},
{'Rh',45,102.9055},
{'Ru',44,101.07},
{'Si',14,28.0855},
{'Sn',50,118.69},
{'Th',90,232.04},
{'Ti',22,47.9},
{'Ta',73,180.9479},
{'U', 92,238.03},
{'V', 23,50.9415},
{'W', 74,183.84},
{'Zn',30,65.38},
{'Zr',40,91.22},
{'Ar', 18,39.95},
{'D',   1,2},
{'T',   1,3},
{'Ga', 31,69.75},
{'H',   1,1},
{'He',  2,4},
{'He4',  2,4},
{'He3', 2,3.02},
{'Kr', 36,83.8},
{'Li',  3,6.94},
{'Na', 11,22.99},
{'Ne', 10,20.18},
{'O',   8,16.00},
{'Xe', 54,131.3},
{'N',7,14},
{'Cl',17,34.453},
};

if nargin == 0
    for k=1:length(liste)
	name{k} = liste{k}{1};
	Z(k)    = liste{k}{2};
	A(k)  = liste{k}{3};
    end
    return
end

ind_impur =  NaN;
for k = 1:length(liste)
    if ischar(element)
	if strmatch(element,liste{k}{1},'exact')
	      ind_impur = k;
	      break
	end
    else
       if liste{k}{2} == element
	    ind_impur = k;
	    break
       end
    end
end

name = liste{ind_impur}{1};
Z = liste{ind_impur}{2};
A = liste{ind_impur}{3};


