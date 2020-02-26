function M = streeToMutMatr(stree)

cells = size(stree,2);
mutations = cells;
M = zeros(cells, mutations);
for i=1:cells
    cur = i;
    M(i, cur) = 1;
    while cur ~= 0
        M(i, cur) = 1;
        cur = stree(cur).parent;
    end
end