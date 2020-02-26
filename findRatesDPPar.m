function [mutOrders,times, likelihood] = findRatesDPPar(stree,rates, minRate,maxRate,eps, AM)
m = length(AM);
G = digraph(AM);
indeg = sum(AM,1);
root = find(indeg == 0);
treeOrder = flip(dfsearch(G,root))';
mutOrders = cell(1,m);

for v = treeOrder
    if isempty(stree(v).children)
        mutOrders{v} = [v];
        continue;
    end
    
%     if length(stree(v).children) == 1
%         mutOrders{v} = [v stree(v).children];
%         continue;
%     end
    
    desc = bfsearch(G,v);
    [substree,subMutOrders,subrates,ids] = getSubtree(stree,mutOrders,rates,desc);
    subAM = AM(desc,desc);
    
    P = perms(flip(stree(v).children));
    likelP = zeros(1,length(P));
    timesP = cell(1,length(P));
    
%     idsP = cell(1,length(P));
%     for j = 1:length(P)
%         idsP{j} = ids;
%     end
    
    parfor j=1:length(P)
%        ids1 = idsP{j}; 
       subMutOrdersP = subMutOrders;
       subMutOrdersP{ids(v)} = [ids(v) ids(P(j,:))];
%        [tP,lP] = findTimeFixedOrderRate1(substree,subAM,subMutOrdersP,subrates,minRate,maxRate);
%        [tP,lP] = findTimeFixedOrderRate2(substree,subAM,subMutOrdersP,subrates,minRate,maxRate);
%        [tP,lP] = findTimeFixedOrderRate3(substree,subAM,subMutOrdersP,subrates,minRate,maxRate);
%        [tP,lP] = findTimeFixedOrderRate4(substree,subAM,subMutOrdersP,subrates,minRate,maxRate);
        [tP,lP] = findTimeFixedOrderRate5(substree,subAM,subMutOrdersP,subrates,minRate,maxRate);
       likelP(j) = -lP; 
       timesP{j} = tP;
    end
    [maxL,imL] = max(round(likelP,5));
    mutOrders{v} = [v P(imL,:)];
    if v == root
        times = zeros(1,m);
        for j = 1:m
            times(j) = timesP{imL}(ids(j));
        end
        likelihood = maxL;
    end
end