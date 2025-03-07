function [c] = zmat2cell(m,rowSizes,colSizes)
%ZMAT2CELL Break matrix up into a cell array of matrices.
%
%	Synopsis
%
%	  zmat2cell(M,R,C);
%
%	Description
%
%	  ZMAT2CELL(M,R,C) takes three arguments,
%	    M - ROWxCOL Matrix.
%	    R - Vector of row sizes (should sum to ROW).
%	    C - Vector of col sizes (should sum to COL).
%	  and returns a cell array of matrices, found using R and C.
%
%	Example
%
%	   M = [1 2 3 4; 5 6 7 8; 9 10 11 12];
%	   C = zmat2cell(M,[1 2],[2 2])
%	    
%	See also ZCELL2MAT


switch nargin

  case 1,
    c = {m};
  
  case 2,
    rows = length(rowSizes);
    c = cell(rows,1);
    rowStart = 0;
    for i=1:rows
      c{i,1} = m(rowStart+[1:rowSizes(i)],:);
      rowStart = rowStart + rowSizes(i);
    end
	
  case 3,
    rows = length(rowSizes);
	cols = length(colSizes);
    c = cell(rows,cols);
    rowStart = 0;
    for i=1:rows
	  colStart = 0;
	  for j=1:cols
        c{i,j} = m(rowStart+[1:rowSizes(i)],colStart+[1:colSizes(j)]);
        colStart = colStart + colSizes(j);
	  end
      rowStart = rowStart + rowSizes(i);
    end

end
