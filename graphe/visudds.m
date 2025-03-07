%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
if ~exist('tempsdds','var')
    tempsdds=5;
end
tempsdds = inputd(' time (s) : ',tempsdds);
if ~exist('sortieFCE','var')
  sortieFCE=[];
end
if ~isfield(sortieFCE,'shot')
  sortieFCE  = FCE(param);
end
if sortieFCE.shot ~= ceil(param.from.shot.num)
  sortieFCE  = FCE(param);
  sortie2 = [];
end
if ~exist('sortie2','var')
  sortie2=[];
end

if ~isfield(sortie2,'val')
  sortie2    = StCont(0,50,8,0.03,sortieFCE,data);
end
  kg = min(find(data.gene.temps>tempsdds));
  R = double(squeeze(data.equi.R(kg,:,:)));
  Z = double(squeeze(data.equi.Z(kg,:,:)));
  x = double(data.equi.rhoRZ(kg,:))./data.equi.rhomax(kg);
  ind  = min(find(sortie2.tdiffece>tempsdds));
  if ~isnan(R)
    R(:,end) = [];
    Z(:,end) = [];
    R(1,:)   = [];
    Z(1,:)   = [];
    nx  = x' * ones(1,size(R,2));
    nx(1,:)  = [];
    valp(1) = griddata(R,Z,nx,sortie2.InvRadNormalOuter(ind),0);
    valp(2) = griddata(R,Z,nx,sortie2.InvRadNormalInner(ind),0);
    valp(3) = griddata(R,Z,nx,sortie2.InvRadInvertedOuter(ind),0);
    valp(4) = griddata(R,Z,nx,sortie2.InvRadInvertedInner(ind),0);
  end
  
h = findobj(0,'type','figure','tag','sawteeth');
if isempty(h)
       h=figure('tag','sawteeth');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
%
%
%
rinv=zeros(size(data.gene.temps));
rinv1 = rinv;
if ~exist('jeux1')
    jeux1=[];
end


for k=param.gene.kmin:param.gene.k
  indinv=max(find(data.prof.q(k,:) <= 1));
  if ~isempty(indinv)
    rinv(k) = param.gene.x(indinv);
  end
  if isfield('jeux1','data')
    if isfield('jeux1.data','prof')
      if ~isempty(jeux1.data.prof)
        indinv1=min(find(jeux1.data.prof.q(k,:) > 1));
        if ~isempty(indinv1)
          rinv1(k) = jeux1.param.gene.x(indinv1);
        end
      end
    end
  end
end

subplot(2,1,1)
ind = sortie2.val(:,1) > 0;
plot(data.gene.temps(ind),mean(sortie2.val(ind,1:4),2),'o',data.gene.temps,rinv)
if isfield(jeux1,'data')
  if isfield('jeux1.data','prof')
    if ~isempty(jeux1.data.prof)
      hold on
      plot(jeux1.data.gene.temps,rinv1,'--')
      hold off
      legend('Exp.','CRO.','Ref.')
      kg1 = min(find(jeux1.data.gene.temps>tempsdds));
    end
  end
else
  legend('Exp.','CRO.')
end
axis([data.gene.temps(param.gene.kmin) data.gene.temps(param.gene.k) 0 0.4])
xlabel('time (s)')
ylabel('x(q=1)')
title('q=1 location in normalized radius from ECE diagnostic and CRONOS run')
subplot(2,1,2)
plot(param.gene.x,data.prof.q(kg,:),[valp(1) valp(1)],[0 4],'r--',...
                                    [valp(2) valp(2)],[0 4],'b--',...
                                    [valp(3) valp(3)],[0 4],'r.-',...
                                    [valp(4) valp(4)],[0 4],'b.-',...
				    [0 0.4],[1 1],'k')
if isfield(jeux1,'data')
  if isfield('jeux1.data','prof')
    if ~isempty(jeux1.data.prof)
      hold on
      plot(param.gene.x,jeux1.data.prof.q(kg1,:),'k--')
      hold off
    end
  end
end


axis([0 0.4 0.5 1.5])
title([' shot ',int2str(param.from.shot.num),' t=',num2str(tempsdds,4),' s'])
ylabel('safety factor')
xlabel('x (normalized radius)')	    
