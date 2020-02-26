% same as single_cell_mutate_test but the tree is taken from the SCITE
% output. Measures the performance of out algorithm on reconstructed trees


clear;
ms = [30];
% ms = [10];
nTests = 15;
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

% defines whether to use true trees(generated) rather than scite output
trueTest = false;


for muts = ms
    % loads everything from the saved mat file, the tree is getting from
    % scite output folder
    load(['test_scite_high_mut_rate_' int2str(muts) '.mat'], ['true_data'])
    for test = 1:nTests
        test

        trueAM = true_data(test).AM;
        trueStree = true_data(test).stree;
        timesTrue = true_data(test).timesTrue;
        ratesTrue = true_data(test).ratesTrue;
        likelTrue = true_data(test).likelTrue;
        likelEventsTrue = true_data(test).likelEventsTrue;
        aux = true_data(test).aux;
        
        degTrue = sum(trueAM,2)';
        internTrue = find(degTrue > 0);
        
        [stree, AM] = readSciteMutTree(['scite_data_30_output/mutrate_' int2str(muts) '_' int2str(test) '_ml0.gv'], muts);
        if trueTest
           stree = trueStree;
           AM = trueAM;
        end
        m = size(AM,1);
        
        G = digraph(AM);
        indeg = sum(AM,1);
        root = find(indeg == 0);
        
        
        if length(root) ~= 1
            sprintf('Several children in root! %d %d', test, muts)
            continue
        end
        
        G = digraph(trueAM);
        nodeLabels = cell(1,m);
        for i = 1:m
            nodeLabels{i} = [ int2str(i)];
        end
        %plot(G,'NodeLabel',nodeLabels);
        
        
        if isequal(AM,trueAM)
            sprintf('Equal! %d %d', test, muts)
            
        end
        
        deg = sum(AM,2)';
        intern = find(deg > 0);
        
        internIntersect = intersect(intern, internTrue);
        
   
        mutOrders = cell(1,m);
        for i = 1:m
            mutOrders{i} = [i trueStree(i).children];
        end

        % the problem is that intern in true tree and inferred are
        % different
        ratesTrue1 =  deg./(max(timesTrue(1:m))-timesTrue(1:m));
        for i = internIntersect
            ratesTrue2(i) = deg(i)/(max(timesTrue(trueStree(i).children)) - timesTrue(i));
        end

        1 - errperf(ratesTrue(internIntersect),ratesTrue1(internIntersect),'mare');
        1 - errperf(ratesTrue(internIntersect),ratesTrue2(internIntersect),'mare');

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
         accurRate1(test) = 1 - errperf(ratesTrue(internIntersect),ratesInfer1(internIntersect),'mare');
    %     accurRate1(test) = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesTrue(ind));
    %     accurTime(test) = 1 - errperf(timesTrue(intern(2:end)),timesInfer(intern(2:end)),'mare');
        accurTime1(test) = 1 - errperf(timesTrue(internIntersect(2:end)),timesInfer1(internIntersect(2:end)),'mare');
    %     accurTime1(test) = 1 - mean(abs(timesTrue(ind(2:end))-timesInfer(ind(2:end)))./timesTrue(ind(2:end)));
        corrRateS(test) = corr(ratesTrue(internIntersect)', ratesInfer1(internIntersect)','Type','Spearman');
        corrRateDeg(test) = corr(deg(internIntersect)', ratesInfer1(internIntersect)','Type','Spearman');
    %     stopType(test) = stopFlag;



        nRandRes = 100;
        rdRes = zeros(1,nRandRes);
        rdRes1 = zeros(1,nRandRes);
        for i = 1:nRandRes
            x = minTheta + rand(1,length(internIntersect))*(maxTheta - minTheta);
            rdRes(i) =  1 - errperf(ratesTrue(internIntersect),x,'mare');
            rdRes1(i) = errperf(ratesTrue(internIntersect),x,'rmse')/errperf(ratesTrue(internIntersect),ratesInfer1(internIntersect),'rmse');
        end
        [b,max_i] = maxk(rdRes,10);
        [c,min_i] = mink(rdRes,10);
        rdRes([max_i min_i]) = [];
        randRes(test) = mean(rdRes);
        [b,max_i] = maxk(rdRes1,10);
        [c,min_i] = mink(rdRes1,10);
        rdRes1([max_i min_i]) = [];
        randRes1(test) = mean(rdRes1);

        intern2True = find(degTrue >= 2);
        intern2 = find(deg >= 2);
        
        ordCorrs = zeros(1,length(intern2));
        accurOrder1(test) = mean(ordCorrs);
        outputData(test).stree = stree;
        outputData(test).ratesInfer1 = ratesInfer1;
        outputData(test).timesInfer1 = timesInfer1;
        outputData(test).ordersInfer1 = ordersInfer1;
        outputData(test).likelInfer1 = likelInfer1;
        outputData(test).internIntersect = internIntersect;
        outputData(test).AM = AM;
        

        % compare with rand, with simple heuristic; order comparison -
        % spearman, parents in order tree
    end
    ind = find(accurRate1 ~= 0);
    mean(accurRate1(ind))
    mean(accurTime1(ind))
    mean(accurOrder1(ind))
    mean(randRes(ind))
    mean(randRes1(ind))
    trueSuffix = '';
    if trueTest
        trueSuffix = '_true_';
    end
    
    save(['res' int2str(muts) trueSuffix 'scite_high_mut_rate_70.mat'],'accurOrder1','accurRate1','accurTime1','corrRateDeg','corrRateS','randRes','randRes1','runTimes1','outputData');
end
%% 
% sigma = (maxTheta-minTheta)/20;
% nIterMH = 5000;
% ratesInfer = findRatesMHUniform(stree,AM,mutOrders,minTheta,maxTheta,nIterMH,sigma)
% ind = find(ratesInfer ~= 0);
% accurRate = 1 - mean(abs(ratesTrue(ind)-ratesInfer(ind))./ratesInfer(ind))
% corrRateS = corr(ratesTrue(ind)', ratesInfer(ind)','Type','Spearman')

