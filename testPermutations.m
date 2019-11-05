clear;
m = 10;
minTheta = 1/20;
maxTheta = 1/10;
maxTime = 1000;
eps = 0.001;
[AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(m,minTheta,maxTheta,maxTime);

m = size(AM,1);
G = digraph(AM);
nodeLabels = cell(1,m);
for i = 1:m
    nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
end
plot(G,'NodeLabel',nodeLabels);

mutOrders = cell(1,m);
for i = 1:m
    mutOrders{i} = [i stree(i).children];
end

% sorted order according to mutation rates
sortOrder = cell(1,m);
for i = 1:m
    tmp = [];
    for j = stree(i).children
        tmp = [tmp stree(j).rate];
    end
    [B,I] = sort(tmp);
%     I = flip(I);
    sortOrder{i} = [i stree(i).children(I)];
end

[mo,ti, L(iter)] = findRatesMCMC2(stree,AM,minTheta,maxTheta,maxTime,eps, minTheta, 10);

[mutOrdersI,timesI, LI] = findRatesRec(stree,minTheta,maxTheta,maxTime,eps);
flag = true;
for i=1:m
    if sortOrder{i} ~= mutOrdersI{i}
        disp(i);
        flag = false;
    end
end
if flag
    disp("OK!")
end

% [mutRates,times] = findRatesEMorder2(stree,AM,mutOrders,minTheta,maxTheta,maxTime,eps);

% [mutRates3,times3] = findRatesEMorder3(stree,AM,mutOrders,minTheta,maxTheta,maxTime,eps, minTheta, 1000);

% ratesTrue
% mutRates

% timesTrue
% times

% -likelihoodRatesOrders2(timesTrue,ratesTrue,mutOrders)
% -likelihoodRatesOrders2(times,mutRates,mutOrders)
% -likelihoodRatesOrders2(timesI,mutRatesI,mutOrdersI)
