function [tokamak] = plot_tokamak(name)

%% Jo Lister and Vladimir Dokuka, September 2004
%% recovers the XML tokamak data and plots the tokamak
%% the figure must be set on entry

tokamak = load_tokamak(name)

hold on
for i=1:length(tokamak.pfcoils.coil)
    if(isfield(tokamak.pfcoils.coil{i},'rzdrdz'))
      coil = eval(tokamak.pfcoils.coil{i}.rzdrdz);
      contourr = coil(1)+coil(3)/2*[-1 1 1 -1 -1];
      contourz = coil(2)+coil(4)/2*[1 1 -1 -1 1];
      elseif (isfield(tokamak.pfcoils.coil{i},'contour'))
      contour = eval(tokamak.pfcoils.coil{i}.contour);
      contourr = contour(1,:);
      contourz = contour(2,:);
      else
        disp('pfcoil contour data not found')
        keyboard
    end
    fill(contourr,contourz,'b')
end

rlim=eval(tokamak.limiter.r);
zlim=eval(tokamak.limiter.z);
size = sqrt(mean((rlim-mean(rlim)).^2+(zlim-mean(zlim)).^2));
plot(rlim([1:end 1]),zlim([1:end 1]),'b')

for i = 1:length(tokamak.vessel.turn)
    if(isfield(tokamak.vessel.turn{i},'rzdrdz'))
      coil = eval(tokamak.vessel.turn{i}.rzdrdz);
      contourr = coil(1)+coil(3)/2*[-1 1 1 -1 -1];
      contourz = coil(2)+coil(4)/2*[1 1 -1 -1 1];
      elseif (isfield(tokamak.vessel.turn{i},'contour'))
      contour = eval(tokamak.vessel.turn{i}.contour);
      contourr = contour(1,:);
      contourz = contour(2,:);
      else
        disp('vessel contour data not found')
        keyboard
    end
    plot(contourr,contourz,'k')
end

for i=1:length(tokamak.magdiag.dflux)
 r12 = eval(tokamak.magdiag.dflux{i}.r12);
 z12 = eval(tokamak.magdiag.dflux{i}.z12);
 plot(r12(1),z12(1),'pb')
end

for i=1:length(tokamak.magdiag.bpol)
 r = eval(tokamak.magdiag.bpol{i}.r);
 z = eval(tokamak.magdiag.bpol{i}.z);
 ang = eval(tokamak.magdiag.bpol{i}.polangle);
 quiver(r,z,cos(ang),sin(ang),0.1*size,'r');
end

daspect([1 1 1])

xlabel('R, m','Fontsize',12,'FontWeight','bold')
ylabel('Z, m','FontSize',12,'FontWeight','bold')
set(gca,'Xgrid','on','Ygrid','on','FontSize',12,'FontWeight','bold',...
'LineWidth',2,'box','on')

title(['DINA-CH definition of ' tokamak.general.title])
