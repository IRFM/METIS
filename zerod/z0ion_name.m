function [rep,A_true] = z0ion_name(Z,A)

% impurities
liste = {
{'Ag',47,107.8682 },
{'Al',13,26.98154 },
{'Au',79,196.96654},
{'B',  5,10.81   },
{'Be', 4,9.01    },
{'C', 6,12.011   },
{'Cd',48,112.41},
{'Co',27,58.93},
{'Cr',24,51.996},
{'Cu',29,63.546},
{'Fe',26,55.845},
{'Ge',32,72.61},
{'Hf',72,178.49},
{'Ir',77,192.217},
{'Mn',25,54.94 },
{'Mo',42,95.94 },
{'Nb',41,92.90638},
{'Ni',28,58.6934},
{'Os',76,190.2  },
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
{'U', 92,238.03 },
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
{'He3', 2,3},
{'Kr', 36,83.8},
{'Li',  3,6.94},
{'Na', 11,22.99},
{'Ne', 10,20.18},
{'O',   8,16.00},
{'Xe', 54,131.3},
{'N',7,14},
};
for k = 1:length(liste)
    name_liste{k} = liste{k}{1};
    num_liste(k)  = liste{k}{2};
    a_liste(k)  = liste{k}{3};
end

dd = abs(Z - num_liste);
ind_name   = find(dd == min(dd));
if ~isempty(ind_name)
  name_liste = name_liste(ind_name);
  num_liste  = num_liste(ind_name);
  a_liste    = a_liste(ind_name);
  if (length(ind_name) > 1) && (nargin > 1)
        dd = abs(a_liste - A);
        ind_ok = find(dd == min(dd));
	rep    = name_liste{ind_ok};
	A_true = a_liste(ind_ok);
  else
	rep = name_liste{1};
	A_true = a_liste(1);
  end
elseif nargin == 2
  rep = sprintf('Ion_{%d}^{%d}',A,Z);
  A_true = A;
else
  rep = sprintf('Ion^{%d}',Z);
  A_true = Z .* (7/3);
end