function flags = instance_flags

flags.DONT_APPLY_OVERLAY              = py.libmuscle.DONT_APPLY_OVERLAY;
flags.USES_CHECKPOINT_API             = py.libmuscle.USES_CHECKPOINT_API;
flags.KEEPS_NO_STATE_FOR_NEXT_USE     = py.libmuscle.KEEPS_NO_STATE_FOR_NEXT_USE;
flags.STATE_NOT_REQUIRED_FOR_NEXT_USE = py.libmuscle.STATE_NOT_REQUIRED_FOR_NEXT_USE;
for k=1:(8+4+2+1)
    flags.(sprintf('combined_%d',k))  = py.libmuscle.InstanceFlags(uint8(k)); 
end

% if these flags need to be combine, the syntax is :
%     py.libmuscle.InstanceFlags(uint8(8+4+2+1))
