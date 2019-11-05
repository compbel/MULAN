function [mutRates,times, newL] = findOrder4(stree,rates, AM,mutOrders,minRate,maxRate,maxTime,eps)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
toCont = true;
leafs = find(deg == 0);
nLeafs = length(leafs);
 
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
% ub = (max(deg)/minRate)*ones(1,m);
ub = (m-1)/minRate*ones(1,m);
ub(1) = 0;
% ub = maxTime*ones(1,m);
 
% theta0 = 0.5*(maxRate + minRate);
% currTheta = theta0*ones(1,m);
% currTime = zeros(1,m);
% G = digraph(AM);
% v = bfsearch(G,1);
% for i = 1:m
%     currTime(v(i)) = (i-1)/theta0;
% end
% currTime(m+1) = max(currTime(leafs));
 
currTheta = rates; %extractfield(stree,'rate');
currTime = extractfield(stree,'time');
% currTime(1) = 0.00001;
 
 
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
l = @(t)likelihoodRatesOrders2(t,currTheta,mutOrders);
[newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
     
 
mutRates = currTheta;
times = newTime;