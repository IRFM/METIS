% calcul  de  la section efficace d'arret des neutre rapides  complete
function sv = z0nbistoptot(A,E,pa,pnbi,ne,te,nhe,nlow,nhigh,zlow,zhigh,nw,Sn_fraction)

% temps de slowing down 
tause = max(1e-6,min(100,6.27e8 .* A .* te .^ (3/2) ./ (ne./ 1e6) ./ 17));
% standard
if nargin >= 12
  sv0  = z0suzuki_crx(A,E,ne,te,nhe,nlow,nhigh,nw,zlow,zhigh,Sn_fraction);
  
  %sv0j = z0nbistop(A,E,ne,te,nhe,nlow,nhigh,zlow,zhigh);  
  %figure(21);clf;plot(ne(:).*te(:),sv0(:),'r.',ne(:).*te(:),sv0j(:),'b.');drawnow
else
  sv0 = z0nbistop(A,E,ne,te,nhe,nlow,nhigh,zlow,zhigh);
end
% effet des ions rapides
svbp = z0nbistopfast(A,E,pa,ne,te,pnbi,tause);
% section
sv  =sv0 + svbp;
