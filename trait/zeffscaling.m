function zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz)
% zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz);
% nbar : 1019 m-2
% ptot : MW ( equivalent ploss)
% ip   : MA
% a    : petit rayon (m)
% R0   : grand rayon (m)
% itor : courant toroidal (kA)
% gaz  : 2 -> helium, 1 -> Deuterium
%
if gaz == 2
  alpha = [-0.6907    0.1467    0.3994   -2.0644    0.9234   -0.1183];

  c0 = 3.0160;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);
  a4 = alpha(4);
  a5 = alpha(5);
  a6 = alpha(6);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* ip.^a3 .* a.^a4 .* R0.^a5 .* itor.^a6;
end
if gaz == 1

  alpha = [-0.3931    0.1116   3.6174];
  c0 = 185.2499;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* (a./R0).^a3;


end

