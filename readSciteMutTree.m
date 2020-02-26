% reads scite output but skips the root. mutation 1 is the root
% used for mutation rate estimation tests
function [stree, AM] = readSciteMutTree(filepath, m)
fid = fopen(filepath,'r');
AM = zeros(m+1,m+1);

line = fgets(fid);
line = fgets(fid);

for i = 1:m
    line = fgets(fid); 
    data = sscanf(line,'%i %s %i;');
    u = data(1);
    v = data(4);
    AM(u,v) = 1;
end
ind = 1:m;
AM = AM(ind,ind);
for i=1:length(AM)
    stree(i).rate = 0;
    stree(i).children = [];
    stree(i).time = 0;
    for j=1:length(AM)
        if AM(i,j) == 1
         stree(i).children = [stree(i).children j];
         stree(j).parent = i;
        end
        
    end
end

fclose(fid);