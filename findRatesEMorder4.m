function [mutRates,times,mutOrdersOpt,likel,stopFlag] = findRatesEMorder2(stree,AM,minRate,maxRate,maxTime,maxIter,eps,eps1)
stopFlag = 0;
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
descList = getDesc(AM);


lb = zeros(1,m+1);
ub = (m-1)/minRate*ones(1,m+1);
ub(1) = 0;

theta0 = 0.5*(maxRate + minRate);
currTheta = theta0*ones(1,m);
% [currTime,EMOrders,A,rhs] = getInitTimeConsMatr(mutOrders,stree,minRate,maxRate,maxTime,m,leafs);
currL = -Inf;

% currTheta = extractfield(stree,'rate');
% trueTime = [extractfield(stree,'time') max(extractfield(stree,'time'))];
% currTime(1) = 0.00001;


% options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',15000);
nIter = 1;
ls = zeros(1,maxIter);
while toCont
%     l = @(t)likelihoodRatesOrders1(t,currTheta,mutOrders);
%     [newTime,newL] = fmincon(l,currTime,A,rhs,[],[],lb,ub,[],options);
    
    [mutOrdersOpt,newTime, newL] = findRatesDPPar(stree,currTheta, minRate,maxRate,eps, AM);
    
%     mu = [currTheta, -sum(currTheta,2)];
%     cvx_begin
%         variable t(m+1)
%         minimize -mu*t - sum(log(EMOrders*t))
%         A*t <= rhs
% %         Abd*t <= rhsbd
%         lb' <= t <= ub'
%     cvx_end
%     newL = likelihoodRatesOrders1(t,currTheta,mutOrders);
%     newTime = t';
    
%     newTheta = deg./(newTime(m+1)-newTime(1:m));   
%     newTheta = deg./(max(newTime(1:m))-newTime(1:m));
    newTheta = zeros(1,m);
    for i = intern
        newTheta(i) = deg(i)/(max(newTime(stree(i).children)) - newTime(i));
    end
    
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
%     if abs(currL - newL) < eps1
%         toCont = false;
%     end
    currTheta = newTheta;
    currTime = newTime;
    currL = newL;
    
%     %modify order
%     mutOrdersNew = mutOrders;
%     bigTheta = zeros(1,m);
%     for i = 1:m
%         bigTheta(i) = sum(currTheta(descList{i}),2);
%     end
%     for i = 1:m
%         if ~isempty(stree(i).children)
%             vert1 = intersect(stree(i).children,intern);
%             vert2 = intersect(stree(i).children,leafs);
%             if ~isempty(vert1)
%                 aux = [currTheta(vert1)' vert1'];
%                 aux = sortrows(aux);
% %                 aux = [bigTheta(vert1)' vert1'];
% %                 aux = sortrows(aux,-1);
%                 mutOrdersNew{i} = [i aux(:,2)' vert2];
%             end
%         end
%     end
%     [currTimeNew,EMOrders,Anew,rhsnew] = getInitTimeConsMatr(mutOrdersNew,stree,minRate,maxRate,m,leafs);
%     currTime = currTimeNew;
%     A = Anew;
%     rhs = rhsnew;
%     mutOrders = mutOrdersNew;


    
    nIter = nIter + 1
    ls(nIter) = currL;
    if nIter > maxIter
        stopFlag = 1;
        break;
    end
end

% currTheta = zeros(1,m);
%     for i = intern
%         currTheta(i) = deg(i)/(max(currTime(stree(i).children)) - currTime(i));
%     end

mutRates = currTheta;
times = currTime(1:m);
likel = currL;