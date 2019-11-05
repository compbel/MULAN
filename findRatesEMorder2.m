function [mutRates,times] = findRatesEMorder2(stree,AM,mutOrders,minRate,maxRate,maxTime,eps)
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
ub = (m-1)/minRate*ones(1,m);

% lb = deg/maxRate;
% ub = deg/minRate;
% ub(leafs) = (m-1)/minRate;
ub(1) = 0;

theta0 = 0.5*(maxRate + minRate);
currTheta = theta0*ones(1,m);
currTime = zeros(1,m);


AMOrders = zeros(m,m);
    for i = 1:m
        for j = 2:length(mutOrders{i})
            AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
        end
    end
% G = digraph(AMOrders);
G = digraph(AM);
v = bfsearch(G,1);
for i = 1:m
    currTime(v(i)) = (i-1)/theta0;
end



% currTime(m+1) = max(currTime(leafs));

% currTheta = extractfield(stree,'rate');
% currTime = extractfield(stree,'time');
% currTime(1) = 0.00001;


options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false);
while toCont
    l = @(t)likelihoodRatesOrders2(t,currTheta,mutOrders);
    [newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
    newTheta = deg./(max(newTime(1:m))-newTime(1:m));   
    ind = find(isnan(newTheta));
    newTheta(ind) = 0;
    ind = setdiff(1:m,ind);
%     if (min(newTheta(ind)) < minRate) || (max(newTheta(ind)) > maxRate)
%         %     newTheta = min(newTheta,maxRate);
% %     newTheta = max(newTheta,minRate);
%         newTheta(ind) = minRate + ((maxRate-minRate)/(sum(newTheta(ind),2)))*newTheta(ind);
%     end

    if norm(newTheta(ind)-currTheta(ind)) <= eps
        toCont = false;
    end
    currTheta = newTheta;
    currTime = newTime;
end

mutRates = currTheta;
times = currTime;