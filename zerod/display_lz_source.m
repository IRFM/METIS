function display_lz_source

% load data
tabmat = [];
load(fullfile(fileparts(which(mfilename)),'Lz_zave.mat'),'tabmat');
%
disp('---------------------------------------------------------------')
disp('Since the publication of the reference paper for METIS (J.F. Artaud et al 2018 Nucl. Fusion 58 105001),')
disp('atomic data in METIS has been updated. Additionnaly to the METIS reference, depending on the plasma composition, following references has to be add in any publication citing METIS:')
disp(' * for ADAS and open ADAS database:  Summers H.P. and O’Mullane M.G. 2011 Atomic data and modelling for fusion: the ADAS project AIP Conf. Proc. 1344 17987')
disp(' * for T. Pütterich data: T. Pütterich et al 2019 Nucl. Fusion 59 056013');
disp(' * for NIST data in IAEA database: NIST data has been dowload from IAEA website: https://www-amdis.iaea.org/FLYCHK/')
disp('     * NIST References:')
disp('          1) FLYCHK: an extension to the K-shell spectroscopy kinetics model FLY, H. -K. Chung, W. L. Morgan and R. W. Lee, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 81, November 2003, Pages 107-115');
disp('          2) FLYCHK: Generalized population kinetics and spectral model for rapid spectroscopic analysis for all elements, H.-K. Chung, M.H. Chen, W.L. Morgan, Y. Ralchenko and R.W. Lee, High Energy Density Physics, Volume 1, Issue 1, December 2005, Pages 3-12');
disp(' ')
disp('Data source for cooling rate coefficient currently used in METIS for each element is:')
noms = sort(fieldnames(tabmat));
for k= 1:length(noms)
    switch tabmat.(noms{k}).source{1}
        case 'adas'
            if ~isempty(tabmat.(noms{k}).year) &&  ~isempty(tabmat.(noms{k}).year{1})
                fprintf('Source for cooling rate coefficient of %3.3s is %s (%d)\n',noms{k},'ADAS database',tabmat.(noms{k}).year{1});
            else
                fprintf('Source for cooling rate coefficient of %3.3s is %s\n',noms{k},'ADAS database');
            end
        case 'web'
            if ~isempty(tabmat.(noms{k}).year) &&  ~isempty(tabmat.(noms{k}).year{1})
                fprintf('Source for cooling rate coefficient of %3.3s is %s (%d)\n',noms{k},'Open ADAS database',tabmat.(noms{k}).year{1});
            else
                fprintf('Source for cooling rate coefficient of %3.3s is %s\n',noms{k},'Open ADAS database');
            end
        case 'NIST/IAEA'
            fprintf('Source for cooling rate coefficient of %3.3s is %s\n',noms{k},'NIST data in IAEA database');
        otherwise
            fprintf('Source for cooling rate coefficient of %3.3s is %s\n',noms{k},tabmat.(noms{k}).source{1});
    end
end
disp('---------------------------------------------------------------')
disp(' ')

