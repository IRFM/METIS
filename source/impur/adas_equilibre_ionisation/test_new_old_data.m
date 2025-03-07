function test_new_old_data

tabmat = [];
load('Lz_zave.mat','tabmat');
tabmat_new = tabmat;
tabmat = [];
load('/Applications/software/metis/zerod/Lz_zave.mat','tabmat');


% list of elements
list = fieldnames(tabmat_new);


% boucle sur les éléments
for k=1:length(list)
     [A,Z,name] = chargemasse(list{k});
     try
         figure('WindowStyle','docked');
         subplot(2,2,1)
         semilogx(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,3),tabmat_new.(list{k}).data(:,1)*1e3,tabmat_new.(list{k}).data(:,3),'.');
         title(list{k});
         ylabel('<Z>')
         %xlabel('Te (eV)');
         legend('ADAS','Putterrich');
         subplot(2,2,2)
         semilogx(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,4),tabmat_new.(list{k}).data(:,1)*1e3,tabmat_new.(list{k}).data(:,4),'.');
         title(list{k});
         ylabel('<Z^2>')
         xlabel('Te (eV)');
         %legend('ADAS','Putterrich');
         subplot(2,2,3)
         loglog(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,2),tabmat_new.(list{k}).data(:,1)*1e3,tabmat_new.(list{k}).data(:,2),'.');
         title(list{k});
         ylabel('Lz (W*m^3')
         xlabel('Te (eV)');
         %legend('ADAS','Putterrich');
         drawnow
     catch
         disp(name)
     end
end
