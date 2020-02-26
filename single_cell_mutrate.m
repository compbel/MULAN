clear;
m = 20;
nTests = 5;
% minTheta = 1/50;
% maxTheta = 1/10;
minTheta = 1/200;
maxTheta = 1/100;
maxTimeSimul = 10000000;
eps = 0.001;
maxIter = 100;

for test = 1:nTests
    test
    [AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(m,minTheta,maxTheta,maxTimeSimul);

    m = size(AM,1);
%     G = digraph(AM);
%     nodeLabels = cell(1,m);
%     for i = 1:m
%         nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
%     end
%     plot(G,'NodeLabel',nodeLabels);
    deg = sum(AM,2)';

    mutOrders = cell(1,m);
    for i = 1:m
        mutOrders{i} = [i stree(i).children];
    end


    -likelihoodRatesOrders2(timesTrue,ratesTrue,mutOrders)
    -likelihoodRatesOrders1([timesTrue max(timesTrue)],ratesTrue,mutOrders)
    
    maxTime = max(timesTrue);

    % [mutOrders,times, best] = findRatesBFS1(stree,minTheta,maxTheta,maxTime,eps, AM);


    [ratesInfer,timesInfer,stopFlag] = findRatesEMorder1(stree,AM,mutOrders,minTheta,maxTheta,maxTime,maxIter,eps);
    % [ratesInfer,timesInfer] = findRatesEMorder2(stree,AM,mutOrders,minTheta,maxTheta,maxTime,eps);

    ind = intersect(find(ratesInfer ~= 0),find(ratesInfer < 1));

    ratesTrue(ind)
    ratesInfer(ind)

    % timesTrue
    % timesInfer
    % 
%     -likelihoodRatesOrders2(timesTrue,ratesTrue,mutOrders)
%     -likelihoodRatesOrders2(timesInfer,ratesInfer,mutOrders)


    accurRate(test) = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesInfer(ind))
    % corrRateP = corr(ratesTrue(ind)', ratesInfer(ind)','Type','Pearson')
    corrRateS(test) = corr(ratesTrue(ind)', ratesInfer(ind)','Type','Spearman')
    corrRateDeg(test) = corr(deg(ind)', ratesInfer(ind)','Type','Spearman')
    stopType(test) = stopFlag;
%     sqrt(sum((ratesTrue(ind)-ratesInfer(ind)).^2,2))/mean(ratesTrue(ind))
    % ind = setdiff(ind,1);

    % accurTime = 1 - mean(abs(timesTrue(ind)-timesInfer(ind))./timesInfer(ind))
end
%% 

% sigma = (maxTheta-minTheta)/20;
% nIterMH = 5000;
% ratesInfer = findRatesMHUniform(stree,AM,mutOrders,minTheta,maxTheta,nIterMH,sigma)
% ind = find(ratesInfer ~= 0);
% accurRate = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesInfer(ind))
% corrRateS = corr(ratesTrue(ind)', ratesInfer(ind)','Type','Spearman')

