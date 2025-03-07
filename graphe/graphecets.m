function graphecets(param,data,post)
% graphecets(param,data,post.ece);
%
% trace les ecart entre les rayons ECE deduit de CRONOS
% et ceux de TPROF et TSBASE
% rho normalise pour TPROF
% grand rayon pour TSBASE
%
% Auteur :
% V. Basiuk
% date de creation
% 15 janvier 2002
% dernieres modifications : 
%
tece         = post.tece;
rhoece       = post.rhoece;
Rece         = post.Rece;
dRece        = post.dRece;
rsuraece     = post.rsuraece;

occur        = 0;
numchoc      = param.from.shot.num;
t            = data.gene.temps;
numstr       = [int2str(numchoc),'.',int2str(occur)];
eval(['[rhoprof,trhoprof,void1,void2]=tsbase(',numstr,',''gprofrhote'');']);
[Rexp,tRexp,void1,void2] = tsbase(fix(numchoc),'gshr');
nrhoece      = interp1(t,rsuraece,trhoprof,'linear');
rhoprof = rhoprof(:,13:end);
val          = abs(rhoprof-nrhoece)./rhoprof;
val2         = abs(rhoprof-nrhoece);
for k=1:size(rhoece,2)

  nval    = val(:,k);
  nval2   = val2(:,k);
  if ~isempty(find(~isnan(nval)))
    mval(k)  = mean(nval(~isnan(nval)));
    mval2(k) = mean(nval2(~isnan(nval2)));
    sval(k)  = std(nval(~isnan(nval)));
    sval2(k) = std(nval2(~isnan(nval2)));
  else
    mval(k)  = 0;
    sval(k)  = 0;
    mval2(k) = 0;
    sval2(k) = 0;
  end
  
end

subplot(2,2,1)

errorbar(mean(rhoprof),mval*100,sval*100,sval*100,'ro')
axis([0 1 0 100])
title(['ecart relatif rho ece CRONOS et TPROF, #',int2str(numchoc)])
xlabel(['r/a TPROF'])
ylabel('%')

subplot(2,2,2)

errorbar(mean(rhoprof),mval2,sval2,sval2,'ro')
axis([0 1 0 0.1])
title(['ecart absolu moyenne, #',int2str(numchoc)])
xlabel(['r/a TPROF'])


subplot(2,2,3)

plotprof(gca,t,1:size(Rece,2),Rece,'color','r');
dRece(isnan(dRece)) = 0;
plotprof(gca,t,1:size(Rece,2),Rece+dRece,'color','r','linestyle','--');
plotprof(gca,tRexp,1:size(Rexp,2),Rexp,'color','b','linestyle','o','linewidth',2);
title('comparison diag eand CRONOS')
legend('Rece','Rece - dRece','Exp.',4)

nRexp    = tsample(Rexp,tRexp,t);
val      = abs(Rece+dRece-nRexp);

for k=1:size(Rece,2)

  nval    = val(:,k);
  if ~isempty(find(~isnan(nval)))
    mval(k)  = mean(nval(~isnan(nval)));
    sval(k)  = std(nval(~isnan(nval)));
  else
    mval(k)  = 0;
    sval(k)  = 0;
  end
  
end

subplot(2,2,4)
nnRexp=nRexp;

nnRexp(isnan(nnRexp))=0;

mRexp=sum(nnRexp,1)./sum(nnRexp~=0,1);
errorbar(mRexp,mval*1000,sval*1000,sval*1000,'ro')
axis([min(nRexp(:)) max(nRexp(:)) 0 15])
title(['ecart absolu moyenne, #',int2str(numchoc)])
xlabel(['R ece TSBASE'])
ylabel('mm')






