function mutRates = findRatesMHUniform(stree,AM,mutOrders,minRate,maxRate,nIter,sigma)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
toCont = true;
leafs = find(deg == 0);
nLeafs = length(leafs);

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
% for i = 1:m
%     iEq = iEq + 1;
%     A(iEq,i) = 1;
%     A(iEq,m+1) = -1;
%     rhs(iEq) = -deg(i)/maxRate;
%     iEq = iEq + 1;
%     A(iEq,i) = -1;
%     A(iEq,m+1) = 1;
%     rhs(iEq) = deg(i)/minRate;
% end

lb = zeros(1,m+1);
ub = (m-1)/minRate*ones(1,m+1);
ub(1) = 0;

theta0 = 0.5*(maxRate + minRate);
currTheta = theta0*ones(1,m);
currTime = zeros(1,m+1);


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

currTime(m+1) = max(currTime(leafs));

% currTheta = extractfield(stree,'rate');
% currTime = extractfield(stree,'time');
% currTime(1) = 0.00001;


options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',15000);
allLikel = zeros(nIter,1);
allRates = zeros(nIter,m);
l = @(t)likelihoodRatesOrders1(t,currTheta,mutOrders);
[currTime,currL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
for iter = 1:nIter
    if mod(iter,10) == 0
        iter
    end
    newTheta = normrnd(currTheta,sigma);
    if (sum(newTheta < minRate,2) == 0) && (sum(newTheta > maxRate,2) == 0)
        l = @(t)likelihoodRatesOrders1(t,newTheta,mutOrders);
        [newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
%         accprob = exp(-newL + currL);
        accprob = currL/newL;
        p = rand;
        if p <= accprob
            currTheta = newTheta;
            currL = newL;
            currTime = newTime;
        end
    end
    allLikel(iter) = currL;
    allRates(iter,:) = currTheta;
end

figure
plot(allLikel)

[minL,iminL] = min(allLikel);
mutRates = allRates(iminL,:);
