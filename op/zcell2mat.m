function [m] = zcell2mat(c)
%ZCELL2MAT Combine a cell array of matrices into one matrix.
%
%	Synopsis
%
%	  m = zcell2mat(c)
%
%	Description
%
%	  M = ZCELL2MAT(C)
%	    C - Cell array of matrices: {M11 M12... ; M21 M22... ; ...}
%	  returns
%	    M - Single matrix: [M11 M12 ...; M21 M22... ; ...]
%
%	Examples
%
%	  c = {[1 2] [3]; [4 5; 6 7] [8; 9]};
%	  m = zcell2mat(c)
%
%	See also ZMAT2CELL

elements = prod(size(c));

if elements == 1
  m = c{1};
  return
end

if elements == 0
  m = [];
  return
end

[rows,cols] = size(c);

if (cols == 1)
  m = cell(1,rows);
  for i=1:rows
    m{i} = [c{i}]';
  end
  m = [m{:}]';
  
else
  m = cell(1,rows);
  for i=1:rows
    m{i} = [c{i,[1:cols]}]';
  end
  m = [m{1,[1:rows]}]';
end
