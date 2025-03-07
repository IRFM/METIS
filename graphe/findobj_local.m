% try to accelerate findobj
function hout = findobj_local(hin,varargin)

persistent dmax

if  exist('findobj','builtin') == 0
    hout = findobj(hin,varargin{:});
    return
end

if isempty(dmax)
  dmax = 3;
end
for d=0:5
  hout = builtin('findobj',hin,varargin{:},'-depth',d);
  if ~isempty(hout)
    break
  end
end
if ~isempty(hout) && (d >= dmax)
 dmax = d + 1;
 fprintf('ZDATAPLOT: incrising research tree depth (@%d); you must retry last action if action have aborted\n',dmax);
end