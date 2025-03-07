function out = extract_st_info(post,tmin,tmax,prev)

if nargin < 3
  tmax = [];
end

% index 
if ~isempty(tmax) && isfinite(tmax)
  ind_valid = find((post.zerod.temps >= tmin) & (post.zerod.temps <= tmax));
else
  ind_valid = find(post.zerod.temps >= tmin);
end
% compute eta_ec
out.eta_ec = mean(post.zerod.ieccd(ind_valid) ./  post.zerod.pecrh(ind_valid));

% search for first and last ST
tloc       = post.zerod.temps(ind_valid);
indice_inv = post.zerod.indice_inv(ind_valid);
indice_st  = find(indice_inv > 0);
% anti aliasing
if ~isempty(indice_st)
    indice_st_ok = indice_st(1);
    for k=2:length(indice_st)
        if indice_st(k-1) == (indice_st(k) - 1)
            % not good
        elseif indice_st(k-1) == (indice_st(k) - 2)
            % not good
        else
            indice_st_ok(end+1) = indice_st(k);
        end
    end
    indice_st  = indice_st_ok;
    out.tmin       = tloc(min(indice_st));
    out.tmax       = tloc(max(indice_st));
    out.nbst       = length(indice_st);
    out.tau_st     = (out.tmax - out.tmin ) ./  out.nbst;
    out.freq_st    = out.nbst ./ (out.tmax - out.tmin );
else
    out.tmin       = NaN;
    out.tmax       = NaN;
    out.nbst       = NaN;
    out.tau_st     = NaN;
    out.freq_st    = NaN;    
end
% serach for rho
rmx      =  post.profil0d.rmx(ind_valid - 1,:);
q        =  post.profil0d.qjli(ind_valid - 1,:);
if ~isempty(indice_st)

    rho_st_v =  NaN * ones(size(indice_st));
    rho_q1_v =  NaN * ones(size(indice_st));
    for k=1:length(indice_st)
        qloc        = q(indice_st(k),:);
        rloc        = rmx(indice_st(k),:) ./ max(rmx(indice_st(k),:));
        rmore       = linspace(0,max(rloc),1001);
        qmore       = pchip(rloc,qloc,rmore);
        ind1        = max(find(qmore <= 1));
        if isempty(ind1)
            dd = abs(qmore - 1);
            ind1 = find(dd == min(dd),1);
        end
        rho_q1_v(k) = rmore(ind1);
        %
        imore       = linspace(1,21,1001);
        rho_st_v(k) = pchip(imore,rmore,indice_inv(indice_st(k)));
    end
    out.rho_st      = mean(rho_st_v);
    out.rho_q1      = mean(rho_q1_v);
    out.std_rho_st  = std(rho_st_v);
    out.std_rho_q1  = std(rho_q1_v);
else
    out.rho_st      = NaN;
    out.rho_q1      = NaN;
    out.std_rho_st  = NaN;
    out.std_rho_q1  = NaN;    
end

% extract power
out.pecrh =  mean(post.zerod.pecrh(ind_valid))/1e6;
% extract maximum
pec        = sum(post.profil0d.pecrh(ind_valid,:),1);
rho        = rmx ./ (max(rmx,[],2) * ones(size(post.profil0d.xli)));
pec_rho    = sum(post.profil0d.pecrh(ind_valid,:) .* rho,1);
out.rho_ec = trapz(post.profil0d.xli,pec_rho,2) ./ max(trapz(post.profil0d.xli,pec,2),eps);
pec_drho2  = sum(post.profil0d.pecrh(ind_valid,:) .* (rho - out.rho_ec) .^ 2,1);
out.delta_rho_ec = 0.076 / 0.0228 * sqrt(trapz(post.profil0d.xli,pec_drho2,2) ./ max(trapz(post.profil0d.xli,pec,2),eps));

% concatenation
if (nargin == 4) && ~isempty(prev)
  noms = fieldnames(prev);
  for k =1:length(noms)
    prev.(noms{k})(end +1) = out.(noms{k});
  end
  out = prev;
end
