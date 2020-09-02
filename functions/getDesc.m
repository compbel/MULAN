function descList = getDesc(AM)
m = length(AM);
descList = cell(1,m);
G = digraph(AM);
for i = 1:m
    descList{i} = bfsearch(G,i);
end