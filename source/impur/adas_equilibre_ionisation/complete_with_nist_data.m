function complete_with_nist_data(matfname,density)

% input
if (nargin < 2) || isempty(density)
    density = 1e19;
end

% chemin
addpath(fullfile(fileparts(which(mfilename)),'IAEA_NIST'));
%
tabmat = [];
load(fullfile(fileparts(which(mfilename)),matfname),'tabmat');

% all possible elements
[A,Z,names] = chargemasse;
% list of elements
list = fieldnames(tabmat);
% missing elements
missing = setdiff(names(:),list(:))


% boucle sur les éléments
for k=1:length(missing)
     [A,Z,name] = chargemasse(missing{k});
     try
         [te,Zave,Lz] = get_iaea_nist_data(Z,density);
         [te_rep,Z_rep,rep] = get_iaea_nist_repartition(Z,density);
     catch
         fprintf('%s(%g,%g) not available\n',name,Z,A);
         disp(lasterr);
         te     = [];
         Zave   = [];
         Lz     = [];
     end
     if ~isempty(te)
         figure;
         subplot(2,2,1)
         semilogx(te,Zave);
         title(missing{k});
         ylabel('<Z>')
         %xlabel('Te (eV)');
         legend('NIST');
         subplot(2,2,2)
         semilogx(te,Zave.^2);
         title(missing{k});
         ylabel('<Z^2>')
         xlabel('Te (eV)');
         %legend('NIST');
         subplot(2,2,3)
         loglog(te,Lz);
         title(missing{k});
         ylabel('Lz (W*m^3)')
         xlabel('Te (eV)');
         %legend('NIST');
         subplot(2,2,4)
         semilogx(te_rep(:),rep');
         title(missing{k});
         ylabel('Abondance')
         xlabel('Te (eV)');
         leg = cell(1,length(Z_rep));
         for ml=1:length(Z_rep)
             leg{ml}= sprintf('Z = %d',Z_rep(ml));
         end
         legend(leg);
         %legend('NIST');
         drawnow
         %
         tabmat.(name).A = A;
         tabmat.(name).Z = Z;
         tabmat.(name).data = cat(2,te(:) ./ 1e3,Lz(:),Zave(:),Zave(:));
         tabmat.(name).Zl = Z_rep(:)';
         tabmat.(name).rep    = rep';
         tabmat.(name).year   = {[]};
         tabmat.(name).source = {'NIST/IAEA'};
         tabmat.(name).quality = NaN * te_rep;
         tabmat.(name).Te     = te_rep(:) ./ 1e3;
         
     end
end

% new name
newname = fullfile(fileparts(which(mfilename)),sprintf('%s_completed_NIST.mat',strrep(matfname,'.mat','')));
save(newname,'tabmat');

