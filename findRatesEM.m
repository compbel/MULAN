function [mutRates,times] = findRatesEM(stree,AM,minRate,maxRate,maxTime,eps)
m = size(stree,2);
deg = zeros(1,m);
for i = 1:m
    deg(i) = length(stree(i).children);
end
toCont = true;

A = zeros(m-1,m);
rhs = zeros(m-1,1);
iEq = 0;
for i = 1:m
    for j = stree(i).children
        iEq = iEq + 1;
        A(iEq,j) = 1;
        A(iEq,i) = -1;
    end
end

minRateV = minRate*ones(1,m);
theta0 = 0.5*(maxRate + minRate);
currTheta = theta0*ones(1,m);
currTime = zeros(1,m);
G = digraph(AM);
v = flip(bfsearch(G,1));
for i = 1:m
    currTime(v(i)) = i/theta0;
end

lb = zeros(1,m);
ub = maxTime*ones(1,m);

options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',false);
while toCont
    l = @(t)likelihoodRates(t,currTheta,deg);
    newTime = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
    newTheta = min([deg./newTime; minRateV]);
    if norm(newTheta-currTheta) <= eps
        toCont = false;
    end
    currTheta = newTheta;
    currTime = newTime;
end

mutRates = currTheta;
times = currTime;