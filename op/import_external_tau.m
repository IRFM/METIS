function val_out = import_external_tau(time_in,val_in,val_origin,time_out)

val_out          = interp1_ex(time_in,val_in,time_out,'linear','extrap');
indnok           = find(~isfinite(val_out));
val_out(indnok)  = val_origin(indnok);
val_out          = max(0,val_out);

