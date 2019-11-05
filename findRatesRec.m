function [mutOrders,times, best] = findRatesRec(stree,rates, minRate,maxRate,eps, AM)
tStep = 10;


mutRates = [];
if isempty(stree(1).children)
    times = [0];
    mutRates = [minRate];
    mutOrders = cell(1,1);
    mutOrders{1} = [1];
    best = 1e10;
    return;
end

mOrders = cell(1,length(stree(1).children));
% run recursevely for each child
for i = 1:length(stree(1).children)
    %     create a subtree
    clear subtree;
    q = [stree(1).children(i)];
    idx = 1;
    subtree(1).parent = 0;
    subtree(1).time = 0;
%     subtree(1).rate = stree(stree(1).children(i)).rate;
    subtree(1).children = [];
    subrates = [rates(stree(1).children(i))];
    while idx <= length(q)
        v = q(idx);
        newidx = length(q)+1;
        
        for k=1:length(stree(v).children)
            subtree(newidx).parent = idx;
            subtree(newidx).children = [];
            subrates(newidx) = rates(stree(v).children(k));
            %             subtree(newidx).time = newidx*tStep;
            subtree(idx).children = [subtree(idx).children newidx];
            newidx= newidx+ 1;
            q = [q stree(v).children(k)];
        end
        idx = idx+ 1;
    end
    %      call recursevely
    [moSub, mrSub, tSub] = findRatesRec(subtree,subrates, minRate,maxRate,eps);
    %   restore mutOrders from recursion output
    %   i node in subtree corresponds to q(i) from stree
    for k=1:length(q)
        mOrders{q(k)} = q(moSub{k});
    end
end

% run findRates
P = perms(stree(1).children);
best = 1e10;

for j=1:length(P)
    mOrders{1} = [1 P(j,:)];
    [stree.time] = deal(0);
    % assign time according to order bfs
    %     clear queue;
    %     queue = [1];
    %     idx = 1;
    %     while idx <= length(queue)
    %         cur = queue(idx);
    %         stree(cur).time = idx*tStep;
    %         queue = [queue stree(cur).children];
    %         idx = idx +1;
    %     end
    AMs = stree2AM(stree);
    for o = 1:length(mOrders)
        for k =2:length(mOrders{o})
            AMs(mOrders{o}(k-1), mOrders{o}(k)) = 1;
        end
    end
    G = digraph(AMs);
    
    or = toposort(G);
    step = 0;
    for i=or
        stree(i).time = step*tStep;
        step = step +1;
    end
    [mutRatesT,timesT, L] = findOrder4(stree, rates,AMs,mOrders,minRate,maxRate,eps);
    if L < best
        mutRates = mutRatesT;
        mutOrders = mOrders;
        times = timesT;
        best = L;
    end
    L;
end
% return found order and times