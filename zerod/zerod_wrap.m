function [z0d, info, prof] = zerod_wrap(option, cons, geo, exp0d, fastrun, save_file)


  if nargin < 5
    fastrun = false;
  end

  if nargin >= 6 && ~isempty(save_file)
    save(save_file);
  end

  if fastrun
    [r1, r2, r3] = zerodfast(option, cons, geo, exp0d);
  else
    [r1, r2, r3] = zerod(option, cons, geo, exp0d);
  end

  if nargin >= 6 && ~isempty(save_file)
    save(save_file);
  end

  z0d = r1;
  info = [];
  prof = r3;

end
