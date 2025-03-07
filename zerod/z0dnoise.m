% fonction de debuitage des consignes pour metis
function z0dinput = z0dnoise(z0dinput)
  
pin = z0dinput.cons.ip + z0dinput.cons.plh + z0dinput.cons.pecrh  + z0dinput.cons.picrh + real(z0dinput.cons.pnbi) + imag(z0dinput.cons.pnbi);  
tau  = 4e-3  .* (z0dinput.cons.ip./1e6) .* sqrt(z0dinput.cons.nbar./1e19) .* (pin ./ 1e6) .^ -(2/3) .*  ...
       z0dinput.geo.R .^2 .* sqrt(z0dinput.geo.K);
tau = mean(tau);
width = ceil(5 .* mean(tau) ./ mean(diff(z0dinput.cons.temps)));
if (width ./ 2) == fix(width ./ 2)
	width = width + 1;
end
width = max(width,3);

if isfield(z0dinput,'denoising')
	if strcmp(z0dinput.denoising,'yes')
		width = width .* 2 + 1;
	end
end

%  % filtrage du bruit temporelle introduit par EFIT
%  z0dinput.geo.a      = zautowfilt(z0dinput.geo.a,z0dinput.cons.temps);
%  z0dinput.geo.R      = zautowfilt(z0dinput.geo.R,z0dinput.cons.temps);
%  z0dinput.geo.z0     = zautowfilt(z0dinput.geo.z0,z0dinput.cons.temps);
%  z0dinput.geo.K      = zautowfilt(z0dinput.geo.K,z0dinput.cons.temps);
%  z0dinput.geo.d      = zautowfilt(z0dinput.geo.d,z0dinput.cons.temps);
%  
%  % filtrage separatrice
%  if isfield(z0dinput.exp0d,'Rsepa')
%  	for k =1:size(z0dinput.exp0d.Rsepa,2)
%  		z0dinput.exp0d.Rsepa(:,k) = zautowfilt(z0dinput.exp0d.Rsepa(:,k),z0dinput.cons.temps);
%  		z0dinput.exp0d.Zsepa(:,k) = zautowfilt(z0dinput.exp0d.Zsepa(:,k),z0dinput.cons.temps);
%  	end
%  end
%  % consigne
%  noms = fieldnames(z0dinput.cons);
%  noms(strmatch('temps',noms,'exact')) =[];
%  for k=1:length(noms)
%  	z0dinput.cons.(noms{k}) = zautowfilt(z0dinput.cons.(noms{k}),z0dinput.cons.temps);
%  end

% filtrage du bruit temporelle introduit par EFIT
z0dinput.geo.a      = sgolayfilt(z0dinput.geo.a,1,width);
z0dinput.geo.R      = sgolayfilt(z0dinput.geo.R,1,width);
z0dinput.geo.z0     = sgolayfilt(z0dinput.geo.z0,1,width);
z0dinput.geo.K      = sgolayfilt(z0dinput.geo.K,1,width);
z0dinput.geo.d      = sgolayfilt(z0dinput.geo.d,1,width);

% filtrage separatrice
if isfield(z0dinput.exp0d,'Rsepa')
	for k =1:size(z0dinput.exp0d.Rsepa,2)
		z0dinput.exp0d.Rsepa(:,k) = sgolayfilt(z0dinput.exp0d.Rsepa(:,k),1,width);
		z0dinput.exp0d.Zsepa(:,k) = sgolayfilt(z0dinput.exp0d.Zsepa(:,k),1,width);
	end
end
% consigne
noms = fieldnames(z0dinput.cons);
noms(strmatch('temps',noms,'exact')) =[];
for k=1:length(noms)
	z0dinput.cons.(noms{k}) = sgolayfilt(z0dinput.cons.(noms{k}),1,width);
end
