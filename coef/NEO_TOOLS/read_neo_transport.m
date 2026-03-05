function struct=read_neo_transport(dirname, opt_rot)
%
%% function struct=read_neo_transport(dirname, opt_rot)
%% opt_rot = 1 -> rotation model 2, reads also rotation output
%

file=load([dirname '/out.neo.transport']);

if ~isempty(file);
nbspecies = (length(file(1,:))-5)./8;
struct.rova  =file(:,1);
struct.phi_norm = file(:,2);
struct.BScur = file(:,3);
struct.vtor_th0 = file(:,4);
struct.upara = file(:,5);
for ij = 1:nbspecies;
    jst = 5+(ij-1).*8;
    struct.(['pflux' num2str(ij)])=file(:,jst+1);
    struct.(['qflux' num2str(ij)])=file(:,jst+2);
    struct.(['vflux' num2str(ij)])=file(:,jst+3);
    struct.(['upara' num2str(ij)])=file(:,jst+4);
    struct.(['kpara' num2str(ij)])=file(:,jst+5);
    struct.(['kparad_' num2str(ij)])=file(:,jst+6);
    struct.(['vtheta' num2str(ij)])=file(:,jst+7);
    struct.(['vphi' num2str(ij)])=file(:,jst+8);    
end;

if opt_rot
file=load([dirname '/out.neo.rotation']);

struct.phi_diff = file(:,2);
for ij = 1:nbspecies;
jst = 2+(2*ij-1);
    struct.(['n_ratio_' num2str(ij)])=file(:,jst);
    struct.(['V_conv_' num2str(ij)])=file(:,jst+1);
end   
struct.phi_theta = file(:,jst+2:end);
else
struct.phi_diff = NaN;
for ij = 1:nbspecies;
struct.(['n_ratio_' num2str(ij)])=NaN;
struct.(['V_conv_' num2str(ij)])=NaN;
end;
struct.phi_theta = NaN;
end;
%

%
%
    else
%
      nbspecies = 6;
struct.rova  = NaN;
struct.phi_norm = NaN;
struct.BScur = NaN;
struct.vtor_th0 = NaN;
struct.upara = NaN;
for ij = 1:nbspecies;
struct.(['pflux' num2str(ij)])=NaN;
struct.(['qflux' num2str(ij)])=NaN;
struct.(['vflux' num2str(ij)])=NaN;
struct.(['upara' num2str(ij)])=NaN;
struct.(['kpara' num2str(ij)])=NaN;
struct.(['kparad_' num2str(ij)])=NaN;
struct.(['vtheta' num2str(ij)])=NaN;
struct.(['vphi' num2str(ij)])=NaN;    
end;
%
if opt_rot;
struct.phi_diff = NaN;
for ij = 1:nbspecies;
struct.(['n_ratio_' num2str(ij)])=NaN;
struct.(['V_conv_' num2str(ij)])=NaN;
end;
struct.phi_theta = NaN;
end;
%
end;   

