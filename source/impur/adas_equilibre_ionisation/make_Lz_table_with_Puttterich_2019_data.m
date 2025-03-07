% script pour fabriquer la table de Lz pour METIS
source = 'T. Pütterich et al 2019 Nucl. Fusion 59 056013';

% input
density = 1e19;

% chemin
addpath(fullfile(fileparts(which(mfilename)),'IAEA_NIST'));


% reference for repartition
ref = load('Lz_zave_HDT0_completed_Putterich_completed_NIST.mat');

%
tabmat = [];
% all possible elements
[A,Z,names] = chargemasse;
% list of elements
% missing elements
missing = names(:);


% get data
data = lz_zmean_Putterich_2019;

% initialisation de la table
for k=1:length(A)
    [A,Z,name] = chargemasse(missing{k});
    
    try
        te = data.te;
        Zave = data.(name).Zave;
        Lz = data.(name).Lz/1e6;
    catch
        fprintf('%s(%g,%g) not available\n',name,Z,A);
        disp(lasterr);
        te     = [];
        Zave   = [];
        Lz     = [];
    end
    try
        Z_rep   = ref.tabmat.(name).Zl;
        rep     = ref.tabmat.(name).rep;
        quality = ref.tabmat.(name).quality;
        te_rep  = ref.tabmat.(name).Te;
    catch
        te_rep = [];
        Z_rep  = [];
        rep    = [];
    end
     
    if ~isempty(Lz)
         tabmat.(name).A = A;
         tabmat.(name).Z = Z;
         tabmat.(name).data = cat(2,te(:) ./ 1e3,Lz(:),Zave(:),Zave(:));
         tabmat.(name).Zl = Z_rep(:);
         tabmat.(name).rep    = rep';
         tabmat.(name).year   = {[]};
         tabmat.(name).source = {'T. Pütterich et al 2019 Nucl. Fusion 59 056013'};
         tabmat.(name).quality = quality;
         tabmat.(name).Te     = te_rep;
         
        % graphe
        figure('color','w','position',[557         165        1101         717]);
        subplot(2,2,1)
        loglog(tabmat.(name).data(:,1),tabmat.(name).data(:,2),'b');
        sl = 'Putterich';
        title(sprintf('%s @ %s/%d',name,sl,tabmat.(name).year{1}));
        xlabel('Te (keV)');
        ylabel('Lz (W m^3)')
        set(gca,'xlim',[1e-3,50]);
        legend('METIS');
        subplot(2,2,2)
        loglog(tabmat.(name).data(:,1),tabmat.(name).data(:,3),'b');
        xlabel('Te (keV)');
        ylabel('<Z>')
        set(gca,'xlim',[1e-3,50]);
        legend('METIS');
        subplot(2,2,3)
        semilogx(tabmat.(name).Te,tabmat.(name).rep);
        xlabel('Te (keV)');
        ylabel('Abundances');
        set(gca,'xlim',[1e-3,50]);
        eval(sprintf('legend(%s''NaN'')',sprintf('''%g'',',tabmat.(name).Zl)))
        subplot(2,2,4)
        loglog(tabmat.(name).data(:,1),tabmat.(name).data(:,4),'b');
        xlabel('Te (keV)');
        ylabel('<Z^2>')
        set(gca,'xlim',[1e-3,50]);
        hgsave(sprintf('Lz_Putterich_%s',name));
        print(gcf, '-dpng', sprintf('Lz_Putterich_%s',name));
        
    else
        fprintf('No data for %s\n',name);
    end
end
% for compatibility
tabmat.He4 = tabmat.He;
tabmat.He3 = tabmat.He;
save Lz_zave_Putterich_2019.mat tabmat;

