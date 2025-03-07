function [Rsepa,Zsepa] = protect_sepa_z0(Rsepa,Zsepa,mode)

% recast if needed
Rsepa = double(Rsepa);
Zsepa = double(Zsepa);

% just test
z0_minmax = (max(Zsepa,[],2) + min(Zsepa,[],2)) ./ 2;
a        = (max(Rsepa,[],2) - min(Rsepa,[],2)) ./ 2;
Rmax     = max(Rsepa,[],2) * ones(1,size(Rsepa,2));
dd       = abs(Rsepa - Rmax);
z0_rmax  = sum(Zsepa .* (dd == (min(dd,[],2)* ones(1,size(Rsepa,2)))),2) ./ ...
           sum(dd == (min(dd,[],2)* ones(1,size(Rsepa,2))),2);
       
if any((abs(z0_minmax ./ a) > 5e-2) & (abs(z0_rmax ./ a) > 5e-2))
    fprintf('LCFS given by points is not corrected from vertical shift'); 
    printcom = true;
else
    printcom = false;
end

switch mode
    case 'minmax'
        Zsepa = Zsepa - z0_minmax * ones(1,size(Zsepa,2));
        if printcom
            fprintf(': corrected using center of LCFS\n');
        end
    case 'rmax'
        Zsepa = Zsepa - z0_rmax * ones(1,size(Zsepa,2));
        if printcom
            fprintf(': corrected using Z of maximum R of LCFS\n');
        end
    otherwise
        % leave uncorrected, just display a warning
        if printcom
            fprintf('!\n');
        end
end