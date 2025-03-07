function test_Putterich_data

% chemin
addpath(fullfile(fileparts(which(mfilename)),'IAEA_NIST'));
%
tabmat = [];
load Lz_zave_HDT0.mat

% list of elements
list = fieldnames(tabmat);

% get data
data = lz_zmean_Putterich_2019;

% boucle sur les éléments
for k=1:length(list)
     [A,Z,name] = chargemasse(list{k});
     try
         te = data.te;
         Zave = data.(name).Zave;
         Lz = data.(name).Lz/1e6;
         figure;
         subplot(2,2,1)
         semilogx(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,3),te,Zave);
         title(list{k});
         ylabel('<Z>')
         %xlabel('Te (eV)');
         legend('ADAS','Putterrich');
         subplot(2,2,2)
         semilogx(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,4),te,Zave.^2);
         title(list{k});
         ylabel('<Z^2>')
         xlabel('Te (eV)');
         %legend('ADAS','Putterrich');
         subplot(2,2,3)
         loglog(tabmat.(list{k}).data(:,1)*1e3,tabmat.(list{k}).data(:,2),te,Lz);
         title(list{k});
         ylabel('Lz (W*m^3')
         xlabel('Te (eV)');
         %legend('ADAS','Putterrich');
         drawnow
     catch
         disp(name)
     end
end
