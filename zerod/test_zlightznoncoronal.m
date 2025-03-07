te = logspace(-1,4,1001);
Z = 18;
[Lz,Z_ave] = zlightznoncoronal(te,1e19*ones(size(te)),1,Z);
[ua,uz,post]=z0coefpost;
indc =  find(uz == Z,1);
lzc   = post(indc).lz;
tzc   = post(indc).te;

figure;
loglog(te/1e3,Lz*1e13,tzc,lzc);
