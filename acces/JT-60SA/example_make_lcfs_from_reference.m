% driver for make_lcfs_from_reference
% assume that a METIS file is loaded

sepa_ref.time      = [0     10      20      30      40      50]
sepa_ref.rxup      = [0     0       0       0.4     0.4     0];    % upper triangularity (minor radius unit)
sepa_ref.zxup      = [1     1       1.687   1.687   1.687   1];    % upper altitude X point (minor radius unit), equivalent to upper elongation
sepa_ref.apup      = [0     0       0       0       0       0];    % upper separatrix angle (R,X)  (LFS, degrees)
sepa_ref.amup      = [0     0       0       0       0       0];    % upper separatrix angle (-R,X) (HFS, degrees)
sepa_ref.ra        = [5.8   6.2     6.2     6.2     6.2     6.2];  % major radius R0 (m)
sepa_ref.za        = [0     0.65    0.65    0.65    0.65    0.65]; % altitude of the magnetic axis (m)
sepa_ref.a         = [1     2       2       2       2       1.6];  % minor radius (m)
sepa_ref.rxdo      = [0     0       0.568   0.568   0.568   0];     % lower triangularity (minor radius unit)
sepa_ref.zxdo      = [1     1       2.0     2.0     2.0     2.0];       % lower altitude X point (minor radius unit), equivalent to lower elongation
sepa_ref.apdo      = [0     0       22.46   22.46   22.46   22.46];   % lower separatrix angle (R,X)  (LFS, degrees)
sepa_ref.amdo      = [0     0       67.92   67.92   67.92   67.92];   % lower separatrix angle (-R,X)  (HFS, degrees)

sepa_ref.update_b0 = 'off';      %
sepa_ref.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
sepa_ref.mode      = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
sepa_ref.filename  = '';
sepa_ref.ton       = -Inf;
sepa_ref.toff      = Inf;


z0dinput = make_lcfs_from_references(z0dinput,sepa_ref);
