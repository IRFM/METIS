%load('Lz_zave.mat') 
%semilogx(te,z0wavez(te),'r',te,z0wavez(te,1),'b',tabmat.W.data(:,1).*1e3,tabmat.W.data(:,3),'ok')
function [zave,zave2] = z0wavez(te,old)

persistent tabmat
persistent te_more
persistent wdata_ext
persistent wdata_ext2



% nouvelle donnees ADAS
if nargin > 1
  tabmat = [];
elseif isempty(tabmat)
  try 
    load('Lz_zave.mat')   
  end
end

if isfield(tabmat,'W')
    if isempty(te_more)
        pp   = polyfit(log(cat(1,tabmat.W.data(end-3:end,1),10000)),log(cat(1,tabmat.W.data(end-3:end,3),74)),2);
        te_more = cat(1,200,500,1000,10000);
        wdata_ext = exp(polyval(pp,log(te_more)));
        pp   = polyfit(log(cat(1,tabmat.W.data(end-3:end,1),10000)),log(cat(1,tabmat.W.data(end-3:end,4),74^2)),2);
        wdata_ext2 = exp(polyval(pp,log(te_more)));
    end
    %zave = interp1(cat(1,0,tabmat.W.data(:,1),te_more) .* 1e3,cat(1,0,tabmat.W.data(:,3),exp(polyval(pp,log(te_more)))),te,'pchip','extrap');
    zave = interp1(cat(1,0,tabmat.W.data(:,1),te_more) .* 1e3,cat(1,0,tabmat.W.data(:,3),wdata_ext),te,'pchip','extrap');
    zave = min(74,max(0,zave));
    
    if nargout >1
      %pp   = polyfit(log(cat(1,tabmat.W.data(end-3:end,1),10000)),log(cat(1,tabmat.W.data(end-3:end,4),74^2)),2);
      %te_more = cat(1,200,500,1000,10000);
      %zave2 = interp1(cat(1,0,tabmat.W.data(:,1),te_more) .* 1e3,cat(1,0,tabmat.W.data(:,4),exp(polyval(pp,log(te_more)))),te,'pchip','extrap');
      zave2 = interp1(cat(1,0,tabmat.W.data(:,1),te_more) .* 1e3,cat(1,0,tabmat.W.data(:,4),wdata_ext2),te,'pchip','extrap');
      zave2 = min(74^2,max(0,zave2));
    end
else
    rep = [
    0   0
    30  6.77
    40  7.95
    50  9.13
    60  10.25
    70  11.46
    100 13.47
    150 15.70
    200 17.71
    300 20.18
    400 21.62
    500 22.81
    600 23.75
    800 25.33
    1000 26.67
    1500 30.47
    2000 33.93
    2300 36.08
    2700 40.39
    3000 42.84
    3500 44.43
    4000 45.30
    5000 46.72
    6000 48.20
    7000 50.00
    10000 54.58
    12000 57.02
    15000 59.42
    20000 62.02
    25000 63.64
    30000 64.99
    40000 66.83
    100000  74];


    zave  = interp1(rep(:,1),rep(:,2),te,'pchip','extrap');
    zave  = min(74,max(0,zave));
    zave2 = zave .^ 2;
end
