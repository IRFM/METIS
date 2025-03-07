% script a inclure dans zero1t
if option.plotconv > 0
   % recherche de imag et nan
   noms = fieldnames(zs);
   stf  = 0;
   for kzw = 1:length(noms)
       nomc = noms{kzw};
       varc = getfield(zs,nomc);
       if any(~isfinite(varc(:)));
            fprintf('NaN in %s\n',nomc);
            stf = 1;
       end
       if any(imag(varc(:)));
            fprintf('Imag in %s\n',nomc);
            stf = 1;
       end
   end
   if stf == 1
            disp('stop in zerot1 : search for NaN or Imag');
            dbstack
            keyboard
   end
end
