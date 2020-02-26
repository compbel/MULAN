clear;
ms = [70 90 110 130 150];
% ms = [10];
nTests = 10;
% minTheta = 1/50;
% maxTheta = 1/10;
minTheta = 1/200;
maxTheta = 1/100;
maxTimeSimul = 10000000;
% eps = 0.01;
eps = 0.001;
eps1 = 0.001;
maxIter = 10;

accurRate = zeros(1,nTests);
accurTime = zeros(1,nTests);
accurRate1 = zeros(1,nTests);
accurTime1 = zeros(1,nTests);
accurOrder1 = zeros(1,nTests);
runTimes1 = zeros(1,nTests);
corrRateS = zeros(1,nTests);
corrRateDeg = zeros(1,nTests);
randRes = zeros(1,nTests);
randRes1 = zeros(1,nTests);

for muts = ms
    for test = 1:nTests
        test
        
        md = 10;
        while md >= 7
            [AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(muts,minTheta,maxTheta,maxTimeSimul,'uniform');

            m = size(AM,1);
            G = digraph(AM);
            nodeLabels = cell(1,m);
            for i = 1:m
                nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
            end
            plot(G,'NodeLabel',nodeLabels);
            deg = sum(AM,2)';
            intern = find(deg > 0);
            md = max(deg);
        end

        mutOrders = cell(1,m);
        for i = 1:m
            mutOrders{i} = [i stree(i).children];
        end

        ratesTrue1 =  deg./(max(timesTrue(1:m))-timesTrue(1:m));
        for i = intern
            ratesTrue2(i) = deg(i)/(max(timesTrue(stree(i).children)) - timesTrue(i));
        end

        1 - errperf(ratesTrue(intern),ratesTrue1(intern),'mare')
        1 - errperf(ratesTrue(intern),ratesTrue2(intern),'mare')

        maxTime = 5*max(timesTrue);

        % [mutOrders,times, best] = findRatesBFS1(stree,minTheta,maxTheta,maxTime,eps, AM);


    %     tic
    %     stree1 = stree;
    %     for i = 1:length(stree)
    %         stree(i).children = flip(stree(i).children);
    %     end
    %     [ratesInfer1,timesInfer1,stopFlag1] = findRatesEMorder1(stree1,AM,mutOrders,minTheta,maxTheta,maxTime,maxIter,eps);
    %     tm1 = toc;

        tic
    %     [ratesInfer,timesInfer,ordersInfer,likelInfer,stopFlag] = findRatesEMorder2(stree,AM,minTheta,maxTheta,maxTime,maxIter,eps,eps1);
        [ratesInfer1,timesInfer1,ordersInfer1,likelInfer1,stopFlag1] = findRatesEMorder4(stree,AM,minTheta,maxTheta,maxTime,maxIter,eps,eps1);
        runTimes1(test) = toc;
    %     [ratesInfer,timesInfer,likelInfer] = findRatesEMfixedOrder(stree,AM,mutOrders,minTheta,maxTheta,maxTime,eps);

    %     ind = intersect(find(ratesInfer ~= 0),find(ratesInfer < 1));
    %     ind1 = intersect(find(ratesInfer1 ~= 0),find(ratesInfer1 < 1));



    %     accurRate(test) = 1 - errperf(ratesTrue(intern),ratesInfer(intern),'mare');
         accurRate1(test) = 1 - errperf(ratesTrue(intern),ratesInfer1(intern),'mare');
    %     accurRate1(test) = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesTrue(ind));
    %     accurTime(test) = 1 - errperf(timesTrue(intern(2:end)),timesInfer(intern(2:end)),'mare');
        accurTime1(test) = 1 - errperf(timesTrue(intern(2:end)),timesInfer1(intern(2:end)),'mare');
    %     accurTime1(test) = 1 - mean(abs(timesTrue(ind(2:end))-timesInfer(ind(2:end)))./timesTrue(ind(2:end)));
        corrRateS(test) = corr(ratesTrue(intern)', ratesInfer1(intern)','Type','Spearman');
        corrRateDeg(test) = corr(deg(intern)', ratesInfer1(intern)','Type','Spearman');
    %     stopType(test) = stopFlag;



        nRandRes = 100;
        rdRes = zeros(1,nRandRes);
        rdRes1 = zeros(1,nRandRes);
        for i = 1:nRandRes
            x = minTheta + rand(1,length(intern))*(maxTheta - minTheta);
            rdRes(i) =  1 - errperf(ratesTrue(intern),x,'mare');
            rdRes1(i) = errperf(ratesTrue(intern),x,'rmse')/errperf(ratesTrue(intern),ratesInfer1(intern),'rmse');
        end
        [b,max_i] = maxk(rdRes,10);
        [c,min_i] = mink(rdRes,10);
        rdRes([max_i min_i]) = [];
        randRes(test) = mean(rdRes);
        [b,max_i] = maxk(rdRes1,10);
        [c,min_i] = mink(rdRes1,10);
        rdRes1([max_i min_i]) = [];
        randRes1(test) = mean(rdRes1);

        intern2 = find(deg >= 2);
        ordCorrs = zeros(1,length(intern2));
        for i = 1:length(intern2)
            ordCorrs(i) = corr(ordersInfer1{intern2(i)}(2:end)',mutOrders{intern2(i)}(2:end)','Type','Kendall');
        end
        accurOrder1(test) = mean(ordCorrs);


        % compare with rand, with simple heuristic; order comparison -
        % spearman, parents in order tree
    end
    ind = find(accurRate1 ~= 0);
    mean(accurRate1(ind))
    mean(accurTime1(ind))
    mean(accurOrder1(ind))
    mean(randRes(ind))
    mean(randRes1(ind))
    save(['res' int2str(muts) 'exh_2_r.mat'],'accurOrder1','accurRate1','accurTime1','corrRateDeg','corrRateS','randRes','randRes1','runTimes1');
end
%% 
% sigma = (maxTheta-minTheta)/20;
% nIterMH = 5000;
% ratesInfer = findRatesMHUniform(stree,AM,mutOrders,minTheta,maxTheta,nIterMH,sigma)
% ind = find(ratesInfer ~= 0);
% accurRate = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesInfer(ind))
% corrRateS = corr(ratesTrue(ind)', ratesInfer(ind)','Type','Spearman')

