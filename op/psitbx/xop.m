function c = xop(a,ka,op,b,kb)

% XOP   Expand operand operation
%   C=XOP(A,':','-',B,[2,4]) C(I1,I2,I3,I4,...)=A(I1,I2,I3,I4,...)-B(I2,I4)
%   C=XOP(A,[4,2],'-',B,':') C(I1,I2,I3,I4,...)=A(I4,I2)-B(I1,I2,I3,I4,...)

if length(a) <= 1 | length(b) <= 1 | (strcmp(ka,':') & strcmp(kb,':'))
 c = feval(op,a,b);
elseif strcmp(ka,':')
 ka = 1:max([ndims(a) kb]);
 dima = [size(a) ones(1,length(ka)-ndims(a))];
 kper = ka; kper(kb) = [];
 c = feval(op,permute(a,[kb kper]),repmat(b,[ones(1,length(kb)) dima(kper)]));
 c = ipermute(c,[kb kper]);
else
 kb = 1:max([ndims(b) ka]);
 dimb = [size(b) ones(1,length(kb)-ndims(b))];
 kper = kb; kper(ka) = [];
 c = feval(op,repmat(a,[ones(1,length(ka)) dimb(kper)]),permute(b,[ka kper]));
 c = ipermute(c,[ka kper]);
end
