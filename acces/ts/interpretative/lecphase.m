function [phase,cert] = lecphase(shot,ti,tf,xfcim,tfcim,freq)
xfcim1 = xfcim(:,1);
xfcim2 = xfcim(:,2);
xfcim3 = xfcim(:,3);
tfcim1 = tfcim;
tfcim2 = tfcim;
tfcim3 = tfcim;
phasebrut = [];
phase = [];
cert = [];

if shot >= 29000

  [phas,tphas,void,cert] = tsbase(shot,'gphfci');
  if isempty(phas)
    phase = [0 0 0];
    return
  end
  phas         = abs(interp1(tphas,phas,tfcim,'nearest'));
  ind          = find(xfcim1 > 0.2);
  if ~isempty(ind)
    phase1     = mean(phas(ind,1));
  else 
   phase1      = 0;
  end
  ind          = find(xfcim2 > 0.2);
  if ~isempty(ind)
    phase2     = mean(phas(ind,2));
  else 
   phase2      = 0;
  end
  ind          = find(xfcim3 > 0.2);
  if ~isempty(ind)
    phase3     = mean(phas(ind,3));
  else 
   phase3      = 0;
  end
  phase = [phase1 phase2 phase3];
  return	
elseif shot>=20100,
   [phasebrut,tphase,void,cert]=tsbase(shot,'gphase');
   if ~isempty(tphase)
     tphase = tphase(:,1);
   end
end

if ~isempty(phasebrut)

pathf = fullfile(fileparts(which('lecphase')),'lombard','param.ant');

nomfic   = 'caltirphi';
load(fullfile(pathf,sprintf('liste%s',nomfic)));
numfic   = min(find(shot<choc));
if ~isempty(numfic)
  nomfic = [nomfic num2str(numfic)];
end
%eval(['load /cgc_data/zineb/.upper/lombard/param.ant/' nomfic]);
load(fullfile(pathf,nomfic));

nomfic   = 'calmesphi';
%eval(['load /cgc_data/zineb/.upper/lombard/param.ant/liste' nomfic])
load(fullfile(pathf,sprintf('liste%s',nomfic)));
numfic   = min(find(shot<choc));
if ~isempty(numfic)
  nomfic = [nomfic num2str(numfic)];
end
%eval(['load /cgc_data/zineb/.upper/lombard/param.ant/' nomfic]);
load(fullfile(pathf,nomfic));

nomfic   = 'longsond';
%eval(['load /cgc_data/zineb/.upper/lombard/param.ant/liste' nomfic])
load(fullfile(pathf,sprintf('liste%s',nomfic)));
numfic   = min(find(shot<choc));
if ~isempty(numfic)
  nomfic = [nomfic num2str(numfic)];
end
load(fullfile(pathf,nomfic));
%eval(['load /cgc_data/zineb/.upper/lombard/param.ant/' nomfic]);

lambda=3e8./(freq.*1e6); 

corr(1)=(calmesphi(1)+(longsond(1)-longsond(3))).*360./lambda(1);
corr(2)=(calmesphi(2)+(longsond(2)-longsond(4))).*360./lambda(1);
corr(3)=(calmesphi(3)+(longsond(5)-longsond(7))).*360./lambda(2);
corr(4)=(calmesphi(4)+(longsond(6)-longsond(8))).*360./lambda(2);
corr(5)=(calmesphi(5)+(longsond(9)-longsond(11))).*360./lambda(3);
corr(6)=(calmesphi(6)+(longsond(10)-longsond(12))).*360./lambda(3);

corr(7:8)=calmesphi(7:8).*360./lambda(1);
corr(9:10)=calmesphi(9:10).*360./lambda(2);
corr(11:12)=calmesphi(11:12).*360./lambda(3);

%Correction de l'erreur introduite par les differences de calibration
%des tiroirs de mesure
for index=1:12,
   Phase(:,index)=(phasebrut(:,index)>0).*...
   (caltirphi(index,1).*phasebrut(:,index)+...
   caltirphi(index,2))+(phasebrut(:,index)<0).*...
   (caltirphi(index,3).*phasebrut(:,index)+caltirphi(index,4));
   Phase(:,index)=Phase(:,index)+corr(index);

%on ramene phi entre -45 et +315
   Phase(:,index)=(Phase(:,index)>315).*(Phase(:,index)-360)+...
   (Phase(:,index)<-45).*(Phase(:,index)+360) +...
   ((Phase(:,index)>=-45)&( (Phase(:,index)<=315))).*Phase(:,index);
end

phase = Phase(:,1:6);

phase1 = mean(phase(:,1:2)')';

phase1 = tsample(phase1,tphase,tfcim1);
ind    = xfcim1 > 0.2;
if ~isempty(ind)
  phase1 = mean(phase1(ind));
else
  phase1 = 0;
end

phase2 = mean(phase(:,3:4)')';

phase2 = tsample(phase2,tphase,tfcim2);
ind    = xfcim2 > 0.2;
if ~isempty(ind)
  phase2 = mean(phase2(ind));
else
  phase2 = 0;
end

phase3 = mean(phase(:,5:6)')';

phase3 = tsample(phase3,tphase,tfcim2);
ind    = xfcim3 > 0.2;
if ~isempty(ind)
  phase3 = mean(phase3(ind));
else
  phase3 = 0;
end

phase = [phase1 phase2 phase3];


else

phase = [0 0 0];

end

