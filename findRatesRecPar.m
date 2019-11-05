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
lperms = zeros(1,length(P));
ordersP = cell(1,length(P));
timesP = cell(1,length(P));
best = 1e10;

parfor j=1:length(P)
    mOrdersP = mOrders;
    streeP = stree;
    mOrdersP{1} = [1 P(j,:)];
    [streeP.time] = deal(0);
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
    AMs = stree2AM(streeP);
    for o = 1:length(mOrdersP)
        for k =2:length(mOrdersP{o})
            AMs(mOrdersP{o}(k-1), mOrdersP{o}(k)) = 1;
        end
    end
    G = digraph(AMs);
    
    or = toposort(G);
    step = 0;
    for i=or
        streeP(i).time = step*tStep;
        step = step +1;
    end
    [mutRatesT,timesT, L] = findOrder4(streeP, rates,AMs,mOrdersP,minRate,maxRate,eps);
    lperms(j) = L;
    ordersP{j} = mOrdersP;
    timesP{j} = timesT;
end
[best,iml] = min(lperms);
times = timesP{iml};
mutOrders = ordersP{iml};

% return found order and times