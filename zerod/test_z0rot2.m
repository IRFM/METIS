% script de test de z0rot2
% modification des donnees de sortie
profli = post.profil0d;
t = post.zerod.temps;  
tp  = profli.temps;
noms = fieldnames(profli);
for l=1:length(noms)
   nomc = noms{l};
   val  = getfield(profli,nomc);
   if size(val,1) ~= length(t) & size(val,1) > 3
      val  = interp1(tp,val,t,'nereast');
      profli = setfield(profli,nomc,val);
   end
end

[wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,sturb,factrot,profli] = ...
  z0rot2(post.zerod,post.z0dinput.option,post.z0dinput.cons,post.z0dinput.geo,profli);

t = post.zerod.temps;  
        	
[deplw,tb]=tsbase(post.z0dinput.shot,'SDEPLW');
[edeplw,tb]=tsbase(post.z0dinput.shot,'SEDEPLW');
[dtb,tb]=tsbase(post.z0dinput.shot,'SDTBRAG');
nom=tsbase(post.z0dinput.shot,'SCIONBRAG');
% temps
tempsh=[tb'-dtb';tb'+dtb';ones(size(tb'))*NaN];
tempsh=tempsh(:);
% d�lacement de la raie w
deplwh=[deplw';deplw';ones(size(tb'))*NaN];
deplwh=deplwh(:);
% d�lacement de W
deplwv=[deplw'-edeplw';deplw'+edeplw';ones(size(tb'))*NaN];
deplwv=deplwv(:);
% temps
tempsv=[tb';tb';ones(size(tb'))*NaN];
tempsv=tempsv(:);
 
if strcmp(lower(deblank(nom)),'fer')
                    fact =   3.0e8;
elseif  strcmp(lower(deblank(nom)),'chrome')
                   fact =  2.5e8;
else
                   fact =NaN;
end
if fix(post.z0dinput.shot)>28424
	fact = - fact;
end

figure
subplot(2,1,1)
plot(t,snbi,t,sicrh,t,sfus,t,sripth,t,sriplh,t,sripicrh,t,sturb, ...
     t,snbi+sicrh+sfus+sripth+sriplh+sripicrh+sturb)
legend('snbi','sicrh','sfus','sripth','sriplh','sripicrh','sturb','total')
subplot(2,1,2)
plot(t,wrad .* post.z0dinput.geo.R ,'g', ...
     t,profli.vtor(:,1),'b', ...
     [tempsh;tempsv],[deplwh;deplwv].*fact+1e5,'r')
%set(gca,'ylim',[0.9 .* min( deplw.*fact),1.1 .* max(deplw.*fact)]);    
figure
zplotprof(gca,t * ones(1,size(profli.Raxe,2)),profli.Raxe + profli.rmx,profli.vtor./ profli.Raxe,'color','r');
zplotprof(gca,t * ones(1,size(profli.Raxe,2)),profli.Raxe + profli.rmx,profli.omega,'color','g');
