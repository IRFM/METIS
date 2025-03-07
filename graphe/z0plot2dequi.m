function z0plot2dequi(post)
%
try
    post = evalin('base','post');
    [temps,R,Z,Rext,Zext]=z0make2dsurface(post);
catch
    return
end
disp('-------------------------------------------')
fprintf('Equilibrium display \n');
disp('time dependent')

fprintf('data size :')
disp(size(R))
disp(' ')

% figure
hf = findobj(0,'type','figure','tag','z0equi2d_standalone');
if isempty(hf)
	hf = figure;
else
	figure(hf);
end
clf
set(hf,'tag','z0equi2d_standalone','DoubleBuffer','on','color',[1 1 1], ...
            'defaultaxesfontsize',16,'defaultaxesfontweight','bold', ...
	    'defaultaxesfontname','times','defaulttextcolor',[0 0 0],'defaultaxescolor',[1 1 1], ...
	    'defaultaxeszcolor',[0 0 0],'defaultaxesxcolor',[0 0 0],'defaultaxesycolor',[0 0 0],'resize','on');
title('Set the rigth figure size and, after,  strike a key');
[xv,yv] = ginput(1);
set(hf,'resize','off');
pause(1)

% limitationdu nombre de surfaces
pas=1;
indr = 1:pas:size(R,2);

for kl =1:length(indr)
	rr= squeeze(R(:,indr(kl),:));
	zz= squeeze(Z(:,indr(kl),:));
	zplotprof(gca,temps,rr,zz,'color',[0 0.7 1]);
	fprintf('.');
end
zplotprof(gca,temps,Rext,Zext,'color','r');
% plot de la paroi
% limiter
rwall = [];
zwall = [];
if isfield(post.z0dinput.option,'first_wall')
    [pfw,ffw] = fileparts(post.z0dinput.option.first_wall);
    filename_wall = fullfile(pfw,sprintf('%s.mat',ffw));
    if exist(filename_wall,'file')
            wall = load(filename_wall);
            if isfield(wall,'R') && isfield(wall,'Z')
                rwall = wall.R(:);
                zwall = wall.Z(:);
            end
    end
end
if isempty(rwall)
    switch upper(post.z0dinput.machine)
    case 'ITER'
        if  exist('iterwall')
            [rwall,zwall] = iterwall;
        else
            rwall = [];
            zwall = [];
        end    
    otherwise
      if ~isempty(strfind(upper(post.z0dinput.machine),'WEST')) 
	      [rwall,zwall] = west_limiter(post.z0dinput.shot);
	      try
		    walldata = Get_Paroi_WEST(post.z0dinput.shot,temps);
		    rwall    = walldata.Rparoi;
		    zwall    = walldata.Zparoi;
		    rwall(:,end+1) = rwall(:,1);
		    zwall(:,end+1) = zwall(:,1);                  
	      end

      elseif ~isempty(strfind(post.z0dinput.machine,'JT-60SA')) || ~isempty(strfind(post.z0dinput.machine,'JT60-SA'))
            wall = load('jt60sawall');
            rwall = wall.wall.data(1:end,1);
            zwall = wall.wall.data(1:end,2);
      else
        rwall = [];
        zwall = [];
      end
    end
end

if isempty(rwall)
    try
	paroi = evalin('base','param.from.paroi');
	if ~isempty(paroi)
	   rawll = paroi.R;
	   zwall = paroi.Z;
	end

    end
end
hold on
if ~isempty(rwall) && ~isempty(zwall)
      if all(size(rwall) > 1) && (size(rwall,1) == length(temps))
	  zplotprof(gca,temps,rwall,zwall,'color','m','linestyle','-');
      else
	  plot(rwall,zwall,'m')
      end
end

if ~isempty(strfind(upper(post.z0dinput.machine),'WEST')) 
    % Plot also LCFS of NICE
    equi = imas_west_get(post.z0dinput.shot,'equilibrium',0,1);
    if isempty(equi)
	disp('NICE+Polarimetry not available, using NICE magnetic only');
	equi = imas_west_get(post.z0dinput.shot,'equilibrium');	
    end
    if ~isempty(equi)
	tequi = equi.time - tsbase(post.z0dinput.shot,'RIGNITRON');
	zplotprof(gca,tequi,squeeze(equi.boundPlasma(:,1,:)),squeeze(equi.boundPlasma(:,2,:)),'color','k','linestyle',':');

    end
end

fprintf('\n');
hold off
xlabel('R (m)')
ylabel('Z (m)')
axis('equal');
title(sprintf('METIS : %s@%d/ 2D equilibrium', ...
          post.z0dinput.machine,post.z0dinput.shot));
