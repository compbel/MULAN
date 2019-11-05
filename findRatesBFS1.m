function [mutOrders,times, best] = findRatesBFS1(stree,minRate,maxRate,maxTime,eps, AM)
tStep = 10;
G = digraph(AM);
bfs = bfsearch(G,1);
bfs = flip(bfs)';
mutOrders = cell(1,length(AM));
hList = [5,10,20];
myfittype = fittype('a*x - b*log(x) + c', 'dependent',{'y'},'independent',{'x'}, 'coefficients',{'a','b','c'});
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
for v=bfs
    % edge cases
    if isempty(stree(v).children)
        mutOrders{v} = v;
        continue;
    end
    lambda = [];
    if length(stree(v).children) == 1
        mutOrders{v} = [v stree(v).children(1)];
        continue
    end
    
    for c=stree(v).children
        lh = zeros(1,length(hList)); % list of likelihoods for different h
        %create subTree, theta, orders for each child
        [subTree, q] = createSubTree(stree, c);
        n = length(subTree);
        mutOSub = cell(1, n);
        % map stree id -> subTree id
        idMap = containers.Map(q, 1:length(q));
        subTheta = zeros(1, n);
        Theta = 0;
        for i=1:n
            for k=mutOrders{q(i)}
                mutOSub{i} = [mutOSub{i} idMap(k)];
            end
            subTheta(i) = subTree(i).rate;
            Theta = Theta + subTree(i).rate;
        end
        initTime = zeros(1,n);
        step = 0;
        % assign time according to order
        for o = 1:length(mutOSub)
            for k =1:length(mutOSub{o})
                if initTime(mutOSub{o}(k)) == 0
                    initTime(mutOSub{o}(k)) = step*tStep;
                    step = step +1;
                end
            end
        end
        
        %solve for subtree for different h
        lb = zeros(1,n);
        ub = ones(1,n);
        ub(1) = 0;
        for i=1:length(hList)
            l = @(t)likelihoodRatesOrdersH1(t,subTheta,mutOSub,hList(i));
            
            [timeH,lh(i)] = fmincon(l, initTime,[],[],[],[],lb,ub,[],options);
        end
        %find a,b coefficients
        lh = -lh;
        [myfit,gof] = fit(hList',lh',myfittype);
        l = hList(1)*(Theta - myfit.a)+myfit.b+n;
        lambda = [lambda l];
    end
    % find order based of found lambdas
    [B,I] = sort(lambda);
    % save order according to sorted lambda
    mutOrders{v} = [v stree(v).children(I)];
end

% assign time according to order
step = 0;
[stree.time] = deal(0);
for o = 1:length(mutOrders)
    for k =1:length(mutOrders{o})
        if stree(mutOrders{o}(k)).time == 0
            stree(mutOrders{o}(k)).time = step*tStep;
            step = step +1;
        end
    end
end
[mutRatesT,times, best] = findOrder3(stree,[],mutOrders,minRate,maxRate,maxTime,eps);
