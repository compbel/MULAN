function [mutRates,times] = findRatesEMorder(stree,AM,mutOrders,minRate,maxRate,maxTime,eps)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
toCont = true;
nLeafs = sum(deg == 0,2);

% A = zeros(m-1,m);
% rhs = zeros(m-1,1);
% iEq = 0;
% for i = 1:m
%     for j = stree(i).children
%         iEq = iEq + 1;
%         A(iEq,j) = -1;
%         A(iEq,i) = 1;
%     end
% end

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
% ub = maxTime*ones(1,m);

theta0 = 0.5*(maxRate + minRate);
currTheta = theta0*ones(1,m);
G = digraph(AM);
v = bfsearch(G,1);
for i = 1:m
    currTime(v(i)) = (i-1)/theta0;
end


% currTheta = extractfield(stree,'rate');
% currTime = extractfield(stree,'time');
% currTime(1) = 0.00001;


options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false);
while toCont
    l = @(t)likelihoodRatesOrders(t,currTheta,mutOrders);
    [newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
    newTheta = deg./(max(newTime)-newTime);
    newTheta(find(isnan(newTheta))) = 0;
    if norm(newTheta-currTheta) <= eps
        toCont = false;
    end
    currTheta = newTheta;
    currTime = newTime;
end

mutRates = currTheta;
times = currTime;