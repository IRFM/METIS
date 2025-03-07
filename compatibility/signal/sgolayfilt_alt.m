function y=sgolayfilt_alt(x,order,width,weight,dim)

persistent B1 B2 B3 B4 B5
if isempty(B1)
    load filter_matrix_B
end
 
if isempty(x)
     y = [];
     return
elseif sum(size(x) > 1) > 1
  if nargin < 5
    dim = 1;
  end
  y = NaN * ones(size(x));
  if dim == 1
      for k =1:size(x,2)
	y(:,k) = sgolayfilt_alt(x(:,k),order,width);
      end
  else
      for k =1:size(x,1)
	y(k,:) = sgolayfilt_alt(x(k,:),order,width);
      end  
  end
  return
end
if order > 5
   error('order > 5 is not implemanted in this version'); 
end

if width < (order +1);
  error('width must be > order + 1');
elseif width > 101
    waring('maximum width is 101 points in this version; width have been change to 101');
end

ss = size(x);
x = x(:);

% selection of the projection matrix B
switch order
    case 1
        B = B1{width};
    case 2
        B = B2{width};       
    case 3
        B = B3{width};
    case 4
        B = B4{width};
    case 5
        B = B5{width};
end
% Preallocate output
y = zeros(size(x));
F = width;
% Compute the transient on
y(1:(F+1)/2-1,:) = flipud(B((F-1)/2+2:end,:))*flipud(x(1:F,:));
% Compute the steady state output
ytemp = filter(B((F-1)./2+1,:),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);
% Compute the transient off
y(end-(F+1)/2+2:end,:) = flipud(B(1:(F-1)/2,:))*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
y = reshape(y,ss); 
