function rcouche=plotcouche(data,param,temps)
%  [rcouche]=plotcouche(data,param,temps)
%  position des couches fci H, D et He3
%  rcouche(harmonique,minoritaire)
%

h = findobj(0,'type','figure','tag','couche FCI');
if isempty(h)
       h=figure('tag','couche FCI');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
freq  = param.cons.fci.frequence;
ind   = iround(data.gene.temps,temps);
  
a     = data.geo.a(ind);
R0    = data.geo.r0(ind);
R1    = squeeze(double(data.equi.R(ind,:,:)));
Z     = squeeze(double(data.equi.Z(ind,:,:)));
bR    = squeeze(double(data.equi.BR(ind,:,:)));
bZ    = squeeze(double(data.equi.BZ(ind,:,:)));
bphi  = squeeze(double(data.equi.BPHI(ind,:,:)));
btot  = sqrt(bphi.^2 + bR.^2 + bZ.^2);
bmag  = data.geo.b0(ind);
axe   = R1(1,1);
B0    = btot(1,1);
x     = param.gene.x;
Te    = data.prof.te(ind,:);
Ti    = data.prof.ti(ind,:);
Pel   = data.source.fci.el(ind,:);
Pion  = data.source.fci.ion(ind,:);
Pabs  = data.gene.paddfci/1e6;

%
%
%
mp                   = 1.6726e-27;
el                   = 1.6022e-19;
me                   = 9.1095e-31;
%
% champ B
%
[n,m]  = size(R1);
[R,iR] = sort(reshape(R1,1,n*m));
B      = reshape(btot,1,n*m);
B      = B(iR);
fp     = el*B/2/pi/mp/1e6;
fd     = el*B/2/pi/mp/2/1e6;
fhe    = 2*el*B/2/pi/mp/3/1e6;
theta  = linspace(0,2*pi,200);
subplot(2,2,1)
plot(R1(1:10:end,:)',Z(1:10:end,:)','r')


axis('equal');
hold on
xlabel('R (m)');
ylabel('Z (m)');

pfci=abs(data.cons.fci(ind,:));
[pmax,nant] = max(pfci);

f    = freq(nant);
nh   = 1;
nd   = 1;
nhe  = 1;
tex1 = [' '];
tex2 = [' '];
tex3 = [' '];
tex4 = [' '];

rcouche(4,3) = 0;
for k=1:4
  ind    = find(f/k>=fp);  
  if isempty(ind) == 0
     ind       = ind(1);
     if k == 1
       tex1=[' 1H = ',num2str(round(R(ind)*100)/100),' m'];
     end
     if k == 2
       tex2=[' 2H = ',num2str(round(R(ind)*100)/100),' m'];
     end     
     rcouche(nh,1) = R(ind);
     if R(ind) > 2.5
       X     = linspace(R(ind)-0.5,R(ind)+0.8,500);
     else 
       X     = linspace(R(ind)-0.1,R(ind)+0.5,500);  
     end  
     Y     = zeros(size(X));
	  if strcmp(param.from.machine,'TS')
       [val,val0,val1] = coucherip(X,Y,B0,f,k,1,1);
       [val,val2,val3] = coucherip(X,Y,B0,f,k,-1,1);                         
       plot([rcouche(nh,1) rcouche(nh,1)],[-a a],'b',[val0 val2],[0 0],'b',[val0 val0],[-0.1 0.1],'b',[val2 val2],[-0.1 0.1],'b')
     else
       plot([rcouche(nh,1) rcouche(nh,1)],[-a a],'b')
	  end
     text('units','data','position',[rcouche(nh,1) a-0.05],'string',[int2str(k),' H'])
     nh = nh+1;
   end
   ind    = find(f/k>=fd); 
   if isempty(ind) == 0 
      if ind(1) > 1
        ind    = ind(1); 
        rcouche(nd,2) = R(ind);
        if R(ind) > 2.5
          X     = linspace(R(ind)-0.5,R(ind)+0.8,500);
        else 
          X     = linspace(R(ind)-0.1,R(ind)+0.5,500);  
        end
        Y     = zeros(size(X));        
	     if strcmp(param.from.machine,'TS')
          [val,val0,val1] = coucherip(X,Y,B0,f,k,1,2);
          if isempty(val0)
            val0 = R(ind(1));
          end
          [val,val2,val3] = coucherip(X,Y,B0,f,k,-1,2); 			 
          plot([rcouche(nd,2) rcouche(nd,2)],[-a a],'g',[val0 val2],[0 0],'g',[val0 val0],[-0.1 0.1],'g',[val2 val2],[-0.1 0.1],'g')
        else
          plot([rcouche(nd,2) rcouche(nd,2)],[-a a],'g')
		  end
        text('units','data','position',[rcouche(nd,2) a-.25*nd],'string',[int2str(k),' D'])
        if k == 3
          tex3 =['3D = ',num2str(round(R(ind)*100)/100),' m'];
        end
        nd = nd+1;
      end
    end
   ind    = find(f/k>=fhe); 
   if isempty(ind) == 0 
      if ind(1) > 1
        ind    = ind(1); 
        rcouche(nhe,3) = R(ind);
        if R(ind) > 2.5
          X     = linspace(R(ind)-0.5,R(ind)+0.8,500);
        else 
          X     = linspace(R(ind)-0.1,R(ind)+0.5,500);  
        end
        Y     = zeros(size(X));  
	     if strcmp(param.from.machine,'TS')		        
          [val,val0,val1] = coucherip(X,Y,B0,f,k,1,1.5);
          if isempty(val0)
            val0 = R(ind(1));
          end
          [val,val2,val3] = coucherip(X,Y,B0,f,k,-1,2); 
        end
        plot([rcouche(nhe,3) rcouche(nhe,3)],[-a a],'c')
        if k < 3
          text('units','data','position',[rcouche(nhe,3) a-.25*nhe],'string',[int2str(k),' He3'])
        end
        nhe = nhe+1;
      end
    end
  end

text('units','normalized','position',[0.1 0.95],'string',tex1);
text('units','normalized','position',[0.4 0.95],'string',tex2);
text('units','normalized','position',[0.7 0.95],'string',tex3);
text('units','normalized','position',[1 0.95],'string',tex4);

title(['position des couches FCI temps =',int2str(temps),' s'])
plot(axe,0,'r*')

axis([min(R) max(R)+0.1 -a a])
hold off
%
%
%
subplot(2,2,2)
Pfci = abs(data.cons.fci)/1e6;
if size(Pfci,2) == 3
plot(data.gene.temps,Pabs,'ro',data.gene.temps,Pfci(:,1),...
     data.gene.temps,Pfci(:,2),...
     data.gene.temps,Pfci(:,3))
     legend('Pabs',[int2str(freq(1)) ' MHz'],[int2str(freq(2)) ' MHz'],[int2str(freq(3)) ' MHz'])
end	  
if size(Pfci,2) == 4
plot(data.gene.temps,Pabs,'ro',data.gene.temps,Pfci(:,1),...
     data.gene.temps,Pfci(:,2),...
     data.gene.temps,Pfci(:,3),...
     data.gene.temps,Pfci(:,3))
	  legend('Pabs',[int2str(freq(1)) ' MHz'],[int2str(freq(2)) ' MHz'],[int2str(freq(3)) ' MHz'],[int2str(freq(4)) ' MHz'])
end	
ylabel('MW')
xlabel('temps')
title('Puissance / antenne')
  
subplot(2,2,3)
indc=find(rcouche > 0);
nat = ['H' 'H' 'H' 'H';'D' 'D' 'D' 'D';'X' 'X' 'X' 'X']';
for k=1:length(indc)
  if rcouche(indc(k)  > axe)
    Rcor = R1(:,1);
  else
    Rcor = R1(:,34);
  end
  pos(k) = interp1(Rcor,x,rcouche(indc(k)),'nearest');
  ch(k)  = rem(indc(k),4);
  minf(k) = nat(indc(k));
  if strcmp(minf(k), 'X')
    mino(k,:) = ['He3'];
  else
    mino(k,:) = [minf(k) '  '];
  end
  
end
plot(x,Te,x,Ti,[pos;pos],[0*pos;ones(size(pos))*max(Te)])
legend('Te','Ti')
ylabel('eV')
for k=1:length(indc)
text('units','data','string',[int2str(ch(k)) mino(k,:)],'position',[pos(k)*1.01 max(Te)*0.9])
end
title(['B0=',num2str(bmag,3),' T, R0=',num2str(R0,4),' m, a=',num2str(a,4),' m'])

subplot(2,2,4)

plot(x,Pel/1e6,x,Pion/1e6,[pos;pos],[0*pos;ones(size(pos))*max(Pel/1e6)])
for k=1:length(indc)
text('units','data','string',[int2str(ch(k)) mino(k,:)],'position',[pos(k)*1.01 max(Pel/1e6)*0.9])
end
ylabel('MW/m2')
title('d�p�t')

