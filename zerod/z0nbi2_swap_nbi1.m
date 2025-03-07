% script pour echange nbi2 avec nbi 1 si pas de puissance dans nbi1
if all(real(z0dinput.cons.pnbi) < 1e3)
      z0dinput.cons.pnbi     = imag(z0dinput.cons.pnbi);
      z0dinput.cons.ftnbi    = imag(z0dinput.cons.ftnbi);
      z0dinput.option.nb_nbi     = 1;
      z0dinput.option.einj       = z0dinput.option.einj2;
      z0dinput.option.rtang      = z0dinput.option.rtang2;		
      z0dinput.option.zext       = z0dinput.option.zext2;

end 
