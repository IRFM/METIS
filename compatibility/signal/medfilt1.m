function y=medfilt1_alt(x,width)

 
if isempty(x)
     y = [];
     return
elseif length(size(x)) > 2
    error('this version works only with vertor');
elseif size(x,1) >1 && size(x,2) > 1
    error('this version works only with vertor');    
end


ss = size(x);
x = x(:);



nx = length(x);
blksz = nx;
if rem(width,2)~=1  
    m = width/2;
else
    m = (width-1)/2;
end
X = [zeros(m,1); x; zeros(m,1)];
y = zeros(nx,1);

% Work in chunks to save memory
indr = (0:width-1)';
indc = 1:nx;
for i=1:blksz:nx
    ind = indc(ones(1,width),i:min(i+blksz-1,nx)) + ...
          indr(:,ones(1,min(i+blksz-1,nx)-i+1));
    xx = reshape(X(ind),width,min(i+blksz-1,nx)-i+1);
    y(i:min(i+blksz-1,nx)) = median(xx,1);
end
y = reshape(y,ss); 
