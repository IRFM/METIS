function index = closest(vector,cvalue)

%index = find(vector>cvalue-sigma&vector<cvalue+sigma);

sigma = max(diff(vector))*2;
index = find(vector>cvalue-sigma&vector<cvalue+sigma);

kcount1 = 0;
while length(index>1) & kcount1 < 1000
  kcount1 = kcount1 + 1;
  sigma = sigma / 2;
  index = find(vector>cvalue-sigma&vector<cvalue+sigma);
end

kcount2 = 0;
while isempty(index) & kcount2 < 1000
  kcount2 = kcount2 + 1;
  sigma = sigma+sigma*0.1;
  index = find(vector>cvalue-sigma&vector<cvalue+sigma);
end

% Too much attempts = problem in the function
if kcount1>=1000 | kcount2 >= 1000
  disp(' ')
  disp('Problem in function closest.m')
  disp('Are you sure the arguments are correct?')
  disp(' ')
  index = 9999;
end

% The output must be one index only
if length(index) > 1
  dist1 = abs(cvalue-vector(index(1)));
  dist2 = abs(cvalue-vector(index(2)));
  if dist1 < dist2
    index = index(1);
  else
    index = index(2);
  end
end

% If there are more than two indices found, there is a problem
if length(index) > 2
  disp('Should not arrive: function must be corrected')
end

