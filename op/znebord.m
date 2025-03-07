%tout en m^-3
function [nea,nbar_mes,nea_mes]=znebord(nbar)

nbar(nbar < 3e17)     = 3e17;
nbar(~isfinite(nbar)) = 3e17;

nbar_mes = [1e19,3e19,6e19,1e20,3e20];
nea_mes  = [5e17,3.5e18,1e19,2e19,5e19];

pp  = polyfit(log(nbar_mes),log(nea_mes),2);
nea = exp(polyval(pp,log(nbar)));

