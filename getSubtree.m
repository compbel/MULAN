function [substree,subMutOrders,subrates,idsSubtree] = getSubtree(stree,mutOrders,rates,subset)
m = length(stree);
substree = stree(subset);
subMutOrders = mutOrders(subset);
subrates = rates(subset);
idsSubtree = zeros(1,m);
idsSubtree(subset) = 1:length(subset);
for u = subset'
    substree(idsSubtree(u)).children = idsSubtree(substree(idsSubtree(u)).children);
    subMutOrders{idsSubtree(u)} = idsSubtree(subMutOrders{idsSubtree(u)});
    if substree(idsSubtree(u)).parent ~= 0
        if  idsSubtree(substree(idsSubtree(u)).parent) == 0
            substree(idsSubtree(u)).parent = [];
        else
            substree(idsSubtree(u)).parent = idsSubtree(substree(idsSubtree(u)).parent);
        end
    end
end
