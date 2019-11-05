function AM = stree2AM(stree)
AM = zeros(length(stree), length(stree));
for i=1:length(stree)
    for j=stree(i).children
        AM(i,j) = 1;
    end
end