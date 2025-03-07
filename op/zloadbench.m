% fonction de chargement des bench
function [out,ref,jdiff,pe,pepion] = zloadbench
out ={};
jdiff={};
pe={};
pepion ={};
ref =[];
try
   	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn05_amorti06');    
end
try
   	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn05_amorti05');    
end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn05_amorti04');   
end
% try
% 	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn05_amorti03');   
% end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn00_amorti10');    
end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn00_amorti06');    
end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn00_amorti05');    
end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn00_amorti04');   
end
try
	out{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/choc_cn00_amorti03');   
end
try
	ref = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/ref_cn05_amorti05');   
end
try
   	jdiff{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/jdiff_cn05_amorti05');    
end
try
   	jdiff{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/jdiff_cn00_amorti05');    
end
try
   	pe{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/pe_cn05_amorti05');    
end
try
   	pe{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/pe_cn00_amorti05');    
end
try
        pepion{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/pepion_cn05_amorti05');    
end
try
   	pepion{end+1} = zexinfoconv('/usr/drfc/cgc/cgc_data/zineb/bench/pepion_cn00_amorti05');    
end

save /usr/drfc/cgc/cgc_data/zineb/bench/benchdata out ref jdiff pe pepion


