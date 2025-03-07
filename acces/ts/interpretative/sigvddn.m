function sv=sigvddn(T)

%function sv=sigvddn(T)
%
% Taux de reaction DD de production de neutron
%                        ---------------------
% sv en m^3/s, T en keV
%
% reference Glasstone-Lovberg 1960

n=find(T==0);
if ~isempty(n)
  T(n)=0.01*ones(size(n));
end
sv=1.17e-20 ./(T.^(2/3)) .*exp(-18.76./T.^(1/3));
