% RUN_CERTIFICATION_CASE  Rerun one METIS certification case in GNU Octave and
% compare it with the MATLAB reference stored in certification/metis/<case>.mat.
%
%   octave --no-gui --quiet octave/run_certification_case.m <case> [<result.tsv>]
%
% The case is rerun the way certification/metis_test.m does it (zerod, or
% zerodfast when the reference was made in fast mode). Pass/fail uses the
% largest relative difference between Octave and MATLAB in the flattop means
% (Ip >= 98% of max) of six key quantities: vloop, pfus, li, iboot, te0 and betan.
% The case passes when that difference is within its tolerance in
% octave/certification_tolerances.csv (default 1e-6). Exit status 1 on failure.

more off;
warning ("off", "all");

args = argv ();
if numel (args) < 1
  error ("usage: run_certification_case.m <case> [<result.tsv>]");
end
casename = regexprep (args{1}, '\.mat$', '');
resfile = "";
if numel (args) >= 2
  resfile = args{2};
end

here = fileparts (mfilename ("fullpath"));
addpath (here);
root = metis_octave_setup ();

% ---- tolerance -------------------------------------------------------------
tol = 1e-6;
fid = fopen (fullfile (here, "certification_tolerances.csv"), "r");
if fid > 0
  fgetl (fid);                                   % header
  while true
    l = fgetl (fid);
    if ! ischar (l), break; end
    c = strsplit (l, ",");
    if numel (c) >= 2 && strcmp (strtrim (c{1}), casename)
      tol = str2double (c{2});
    end
  end
  fclose (fid);
end

% ---- run ---------------------------------------------------------------------
key = {"vloop", "pfus", "li", "iboot", "te0", "betan"};
status = "FAIL"; elapsed = NaN; keydev = NaN(1, numel (key)); worst = NaN; worstname = "";
msg = "";
try
  ref = load (fullfile (root, "certification", "metis", [casename ".mat"]));
  ref = ref.post;
  inp = ref.z0dinput;
  inp.option.machine = inp.machine;
  t0 = tic ();
  if numel (ref.zerod.temps) == numel (ref.profil0d.temps)
    zs = zerod (inp.option, inp.cons, inp.geo, inp.exp0d);
  else
    zs = zerodfast (inp.option, inp.cons, inp.geo, inp.exp0d);
  end
  elapsed = toc (t0);

  r = ref.zerod;
  ft = r.ip > 0.98 * max (r.ip);
  for k = 1:numel (key)
    a = mean (r.(key{k})(ft));
    b = mean (zs.(key{k})(ft));
    if abs (a) > 0
      keydev(k) = abs (a - b) / abs (a);
    else
      keydev(k) = abs (b);
    end
  end

  % informative only: worst field, normalised to its range (time derivatives excluded)
  f = fieldnames (r);
  deriv = {"dw", "dpfus", "diboot", "dini", "dwdt", "dwthdt", "dwbpdt", "dwmagtordt", "drmdt"};
  for k = 1:numel (f)
    if any (strcmp (f{k}, deriv)) || ! isfield (zs, f{k}), continue; end
    x = r.(f{k}); y = zs.(f{k});
    if ! isnumeric (x) || ! isequal (size (x), size (y)) || numel (x) < 2, continue; end
    x = real (x(:)); y = real (y(:)); ok = isfinite (x) & isfinite (y);
    sc = max (abs (x(ok)));
    if isempty (sc) || sc == 0, continue; end
    d = max (abs (x(ok) - y(ok))) / sc;
    if isnan (worst) || d > worst
      worst = d; worstname = f{k};
    end
  end

  if any (! isfinite (zs.pfus(:))) || any (! isfinite (zs.vloop(:)))
    msg = "non-finite values in pfus or vloop";
  elseif max (keydev) <= tol
    status = "PASS";
  else
    msg = sprintf ("max key deviation %.2e exceeds tolerance %.1e", max (keydev), tol);
  end
catch err
  msg = strrep (err.message, "\n", " ");
  if ! isempty (err.stack)
    msg = sprintf ("%s (%s:%d)", msg, err.stack(1).name, err.stack(1).line);
  end
end

% ---- report ------------------------------------------------------------------
printf ("\n%s: %s\n", casename, status);
printf ("  run time      %.0f s\n", elapsed);
for k = 1:numel (key)
  printf ("  %-6s        %.2e\n", key{k}, keydev(k));
end
printf ("  tolerance     %.1e\n", tol);
printf ("  worst field   %s (%.2e, informative)\n", worstname, worst);
if ! isempty (msg), printf ("  note          %s\n", msg); end

if ! isempty (resfile)
  fid = fopen (resfile, "w");
  fprintf (fid, "%s\t%s\t%.0f\t%.2e\t%.1e\t%s\t%.2e\t%s\n", casename, status, elapsed, ...
           max (keydev), tol, worstname, worst, msg);
  fclose (fid);
end

exit (! strcmp (status, "PASS"));
