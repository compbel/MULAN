function [mutRates,times, newL] = findOrder4(stree,rates, AM,mutOrders,minRate,maxRate,eps)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
toCont = true;
leafs = find(deg == 0);
intern = find(deg > 0);
nLeafs = length(leafs);
nIntern = length(intern);
 
A = zeros(m-1+nLeafs,m+1);
rhs = zeros(m-1+nLeafs,1);
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
for i = leafs
    iEq = iEq + 1;
    A(iEq,i) = 1;
    A(iEq,m+1) = -1;
end
% Abd = zeros(nIntern,m+1);
% rhsbd = zeros(nIntern,1);
% iEq = 0;
% for i = 13
%     iEq = iEq + 1;
%     Abd(iEq,i) = 1;
%     Abd(iEq,m+1) = -1;
%     rhsbd(iEq) = -deg(i)/maxRate;
%     iEq = iEq + 1;
%     Abd(iEq,i) = -1;
%     Abd(iEq,m+1) = 1;
%     rhsbd(iEq) = deg(i)/minRate;
% end
% 
% A = [A; Abd];
% rhs = [rhs; rhsbd];
 
lb = zeros(1,m+1);
ub = (m-1)/minRate*ones(1,m+1);
ub(1) = 0;
 
theta0 = 0.5*(maxRate + minRate);
currTheta = rates; %extractfield(stree,'rate');
currTime = extractfield(stree,'time');
 
 
% AMOrders = zeros(m,m);
% EMOrders = zeros(m-1,m+1);
% iEdge = 1;
%     for i = 1:m
%         for j = 2:length(mutOrders{i})
%             AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
%             EMOrders(iEdge,mutOrders{i}(j)) = 1;
%             EMOrders(iEdge,mutOrders{i}(j-1)) = -1;
%             iEdge = iEdge + 1;
%         end
%     end
% % G = digraph(AMOrders);
% % plot(G);
% G = digraph(AM);
% v = bfsearch(G,1);
% for i = 1:m
%     currTime(v(i)) = (i-1)/theta0;
% end
 
currTime(m+1) = max(currTime(leafs));
 
% options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',15000);
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',15000);
l = @(t)likelihoodRatesOrders1(t,currTheta,mutOrders);
[newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
     
 
mutRates = currTheta;
times = newTime;