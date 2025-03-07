% fit d'une courbe par une bezier
function ppb = fitbezier(tm,vm)


% creation des matrices
tm = tm(:);
vm = vm(:); 
T = cat(2,ones(size(tm)),tm,tm .^ 2 ,tm .^ 3);
B = [1,0,0,0;-3,3,0,0;3,-6,3,0;-1,3,-3,1];

ppb = (T * B) \ vm;