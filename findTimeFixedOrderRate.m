function [times,likel] = findTimeFixedOrderRate(stree,AM,mutOrders,rates,minRate,maxRate)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
leafs = find(deg == 0);

A = zeros(m-1,m);
rhs = zeros(m-1,1);
iEq = 0;
for i = 1:m
    if ~isempty(stree(i).children)
        for j = 2:length(mutOrders{i})
            iEq = iEq + 1;
            A(iEq,mutOrders{i}(j-1)) = 1;
            A(iEq,mutOrders{i}(j)) = -1;
        end
    end
end

lb = zeros(1,m);
ub = (m-1)/minRate*ones(1,m);

% lb = deg/maxRate;
% ub = deg/minRate;
% ub(leafs) = (m-1)/minRate;
ub(1) = 0;

theta0 = 0.5*(maxRate + minRate);
currTime = zeros(1,m);
AMOrders = zeros(m,m);
    for i = 1:m
        for j = 2:length(mutOrders{i})
            AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
        end
    end
G = digraph(AMOrders);
% G = digraph(AM);
% v = bfsearch(G,1);
% for i = 1:m
%     currTime(v(i)) = (i-1)/theta0;
% end

tStep = 10;
or = toposort(G);
step = 0;
for i=or
    currTime(i) = step*tStep;
    step = step +1;
end

% options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false);
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
l = @(t)likelihoodRatesOrders2(t,rates,mutOrders);
[times,likel] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
    