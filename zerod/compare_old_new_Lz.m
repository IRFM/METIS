% script de comparaison ancien/nouveau pour les "radiative cooling rate"
[a,z,post]=z0coefpost;
[a_old,z_old,post_old]=z0coefpost(1);
for k=1:26
  figure
  semilogx(post(k).te,post(k).lz,'r',post_old(k).te,post_old(k).lz,'b');
  xlabel('Te (keV)')
  ylabel('L_z');
  title(sprintf('A =%g (%g), Z = %g (%g)',a(k),a_old(k),z(k),z_old(k)));
end

