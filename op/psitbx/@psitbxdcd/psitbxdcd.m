function dcd = psitbxdcd(varargin)

% PSITBXDCD	PsiTbx Diagnostic-Chord-Description constructor
%   PSITBXDCD(RD,ZD,PHID,PVD,TVD,ND,T) creates a psitbxdcd object with
%   (RD,ZD,PHID) being the position of a reference point along the chord in
%   cylindrical coordinates, PVD the poloidal viewing angle, measured between
%   the midplane (positive above) and the chord, TVD the toroidal viewing angle,
%   measured between the vertical plane of the chord and the position radius
%   (positive when pointing toward positive phi), ND the number of points for
%   grid discretisation (see PSITBXGRID) and T the time vector for moving
%   chords. The sizes of RD,ZD,PHID,PVD,TVD must be the same unless scalar and
%   their last dimension must correspond to LENGTH(T) if present.

% type casting
if nargin == 1
 switch class(varargin{1})
 case 'psitbxdcd', dcd = varargin{1}; return
 case 'char'
  switch upper(varargin{1})
  case 'DEMO'
   dcd = psitbxdcd(1.2,0,pi/8,[linspace(-pi/3,pi/3,20);zeros(1,20)],...
    [zeros(1,20);linspace(-pi/3,pi/3,20)]);
  case 'FAN'
   dcd = psitbxdcd(repmat([1.2,.9],20,1),repmat([0,1],20,1),pi/8,...
         [linspace(-pi/3,pi/3,20)',linspace(-2*pi/3,-pi/3,20)'],0);
  case 'CCD'
  end
  return
 end
end

% default
if nargin < 7, [varargin{nargin+1:7}] = deal([]); end
nd = varargin{6}; if isempty(nd), nd = 81; end
t = varargin{7}(:)';

% check and adapt size
dim = [1,1];
for k = 1:5
 if isempty(varargin{k}), varargin{k} = NaN; end
 if length(varargin{k}) > 1
  if isequal(dim,[1,1])
   dim = size(varargin{k});
  else
   if ~isequal(size(varargin{k}),dim)
    error('Incoherent parameter size')
   end
  end
 end
end
for k = 1:5
 if length(varargin{k}) == 1
  varargin{k} = repmat(varargin{k},dim);
 end
end
if length(t) > 1 & dim(end) ~= length(t)
 error('Incoherent parameter and time size')
end

[dcd.rd,dcd.zd,dcd.phid,dcd.pvd,dcd.tvd] = deal(varargin{1:5});
dcd.nd = nd;
dcd.t = t;

dcd = class(dcd,'psitbxdcd');
