function equiplot(varargin)

if nargin < 2
  
  disp(['----------------------------------'])
  disp(sprintf('To be used with 2 arguments:'))
  disp(sprintf('equiplot(shot_number,run_number))'))
  disp(['----------------------------------'])
  result = [];

else

  %% SHOT AND RUN NUMBER
  shot_number = varargin{1};
  run_number  = varargin{2};

  %% READ THE EQUILIBRIUM
  idx = euitm_open('euitm',shot_number,run_number);
  equi = euitm_get(idx,'equilibrium');
  euitm_close(idx);
  assignin('base','equi',equi);
  
  nt=length(equi);
  if length(varargin) == 3
    nt=varargin{3};
  end
  %% CARTESIAN (R,Z) GRIDS
  r2d   = equi(nt).profiles_2d(1).r;
  z2d   = equi(nt).profiles_2d(1).z;
  rho2d = sqrt(equi(nt).profiles_2d(1).phi./max(equi(nt).profiles_2d(1).phi(:)));
  psi2d = equi(nt).profiles_2d(1).psi;  
  assignin('base','r2d',r2d);
  assignin('base','z2d',z2d);
  assignin('base','rho2d',rho2d);
  assignin('base','psi2d',psi2d);
  
  %% BOOZER STRAIGHTFIELD LINE GRIDS
  if size(equi(nt).coord_sys.position.r,1) > 0
    r2d_pt = equi(nt).coord_sys.position.r;
    z2d_pt = equi(nt).coord_sys.position.z;
    assignin('base','r2d_pt',r2d_pt);
    assignin('base','z2d_pt',z2d_pt);
  end
  
  %% DISPLAY THE FIGURE
  markersize = 20;
  fontsize = 15;
  linewidth = 3;  
  clf
  h=axes;
  set(h,'FontSize',fontsize)
  hold on ; grid on
  axis equal

  %% (R,Z) GRIDS
  if length(strfind(equi(nt).coord_sys.grid_type,'NER=2, NEGP=0'))==0

    contour(r2d,z2d,rho2d,50)
    contour(r2d,z2d,psi2d,30,'k','linewidth',1)
    set(gca,'xlim',[min(r2d(:)),max(r2d(:))])
    set(gca,'ylim',[min(z2d(:)),max(z2d(:))])
    
  %% (PSI,THETA) GRIDS
  else
    
    plot(r2d_pt,z2d_pt,'b')
    set(gca,'xlim',[min(r2d_pt(:)),max(r2d_pt(:))])
    set(gca,'ylim',[min(z2d_pt(:)),max(z2d_pt(:))])

  end
    
end
