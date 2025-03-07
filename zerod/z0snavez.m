%load('Lz_zave.mat') 
% te = logspace(-1,6,1001);
% semilogx(te,z0snavez(te),'r',tabmat.Sn.data(:,1).*1e3,tabmat.Sn.data(:,3),'ok')
function [zave,zave2] = z0snavez(te)

persistent tabmat
persistent te_more
persistent sndata_ext
persistent sndata_ext2



% nouvelle donnees ADAS
if isempty(tabmat)
  try 
    load('Lz_zave.mat','tabmat')   
  end
end

if isfield(tabmat,'Sn')
    if isempty(te_more)
        pp   = polyfit(log(cat(1,tabmat.Sn.data(end-3:end,1),10000)),log(cat(1,tabmat.Sn.data(end-3:end,3),74)),2);
        te_more = cat(1,200,500,1000,10000);
        sndata_ext = exp(polyval(pp,log(te_more)));
        sndata_ext2 = sndata_ext .^ 2;
    end
    zave = interp1(cat(1,0,tabmat.Sn.data(:,1),te_more) .* 1e3,cat(1,0,tabmat.Sn.data(:,3),sndata_ext),te,'pchip','extrap');
    zave = min(50,max(0,zave));
    zave2 = zave .^ 2;
else
    warning('Sn data not available');
    zave  = 50 * ones(size(te));
    zave2 = zave .^ 2;
end
