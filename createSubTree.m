% creates a sub tree for a given vertex, subtree(i) corresponds to
% stree(q(i)) - allows to find corresponded vertexes in original tree
function [subtree, q] = createSubTree(stree, vertex)
q = [vertex];
idx = 1;
subtree(1).parent = 0;
subtree(1).time = 0;
subtree(1).rate = stree(vertex).rate;
subtree(1).children = [];
while idx <= length(q)
    v = q(idx);
    newidx = length(q)+1;
    
    for k=1:length(stree(v).children)
        subtree(newidx).parent = idx;
        subtree(newidx).children = [];
        subtree(newidx).rate = stree(stree(v).children(k)).rate;
        %             subtree(newidx).time = newidx*tStep;
        subtree(idx).children = [subtree(idx).children newidx];
        newidx= newidx+ 1;
        q = [q stree(v).children(k)];
    end
    idx = idx+ 1;
end