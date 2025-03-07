function alpha = salpha_lookup(s)

if nargin == 0
    s= linspace(-1,10,1001);
end

%alpha critique from Freidberg book
% with two branches ans second stability

%lookup table
al = linspace(0,4,1001);
sp = 0.78 .* al + sqrt(0.72 .* al - 0.36 - 0.11 .* al .^ 2);
sp(iscomplex(sp)) = NaN;
sp = real(sp);
al(~isfinite(sp)) = [];
sp(~isfinite(sp)) = [];


% sm = 0.78 .* al - sqrt(0.72 .* al - 0.36 - 0.11 .* al .^ 2);
% sm(iscomplex(sm)) = NaN;
% sm = real(sm);
% 
%
alpha = min(4,max(sqrt(eps),interp1(sp,al,s,'pchip','extrap')));

if nargout ~= 0
    return
end

% alternative root(W(s,alpha) = 0)
[al1,al2] = zpolyroot(0.5 + 1.39 .* s .^ 2,-2.17 .* s -1,ones(size(s))); 
al1(imag(al1)~= 0) = NaN;
al2(imag(al2)~= 0) = NaN;


alpha_one = (0.5 + 1.39 .* s .^ 2) ./ (1 + 0.83 .* s);

figure
plot(al,sp,'.r',alpha,s,'b',al1,s,'c',al2,s,'m',alpha_one,s,'k:')
xlabel('alpha');
ylabel('s');
set(gca,'xlim',[-1 5]);
set(gca,'ylim',[-2,Inf]);
legend('real(equation +)','lookup interpolant','root1(W=0)','root2(W=0)','first formula alpha = f(s)');

keyboard