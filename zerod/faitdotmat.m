% cree un fichier matlab avec les donnees
error('do not use');
list = dir('*Lzdata.dat')
tabmat = [];
for k=1:length(list)
    nom = list(k).name;
    if nom(2) == '_'
	element = nom(1);
    else
	element = nom(1:2);
    end	
    tabmat.(element) = importdata(nom,' ',2)
    switch element
    case 'Ar'
      tabmat.(element).Z = 18;
      tabmat.(element).A = 39.95;
    case 'C'
      tabmat.(element).Z = 6;
      tabmat.(element).A = 12.011;
    case 'H'
      tabmat.(element).Z = 1;
      tabmat.(element).A = 1;
      tabmat.('D') = tabmat.(element)
      tabmat.('D').A = 2;
      tabmat.('T') = tabmat.(element)
      tabmat.('T').A = 3;
    case 'Kr'
      tabmat.(element).Z = 36;
      tabmat.(element).A = 83.8;
    case 'Ne'
      tabmat.(element).Z = 10;
      tabmat.(element).A = 20.18;
    case 'O'
      tabmat.(element).Z = 8;
      tabmat.(element).A = 16;
    case 'W'
      tabmat.(element).Z = 74;
      tabmat.(element).A = 183.84;
    case 'Be'
      tabmat.(element).Z = 4;
      tabmat.(element).A = 9.01;
    case 'Fe'
      tabmat.(element).Z = 26;
      tabmat.(element).A = 55.845;
    case 'He'
      tabmat.(element).A = 4;
      tabmat.(element).Z = 3;
      tabmat.('He3') = tabmat.(element)
      tabmat.('He3').A = 3;
      tabmat.('He4') = tabmat.(element)     
    case 'N'
      tabmat.(element).Z = 7;
      tabmat.(element).A = 14;
    case 'Ni'
      tabmat.(element).Z = 28;
      tabmat.(element).A = 58.6934;
    case 'Si'
      tabmat.(element).Z = 14;
      tabmat.(element).A = 28.0855;
    case 'Xe'
      tabmat.(element).Z = 54;
      tabmat.(element).A = 131.3;
  
    otherwise
      error('missing element');
    end
    
end
    
    
 save Lz_zave tabmat   