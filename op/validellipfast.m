% permet de valider la compilation de ellipfast
% flag = 1 si ok
function flag_out = validellipfast

m = linspace(0+eps,1-eps);
[Kf,Ef] = ellipke(m);
Eff      = ellipfast(m,2);
Kff      = ellipfast(m,1);

err = sqrt(sum((Kf -Kff) .^ 2 +(Ef -Eff) .^ 2) ./ length(m) ./ 2);

if ~isfinite(err)
    flag_out = 0;
elseif err > 1e-7
    flag_out = 0;
else
    flag_out = 1;
end

