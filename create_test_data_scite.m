% creates 
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



for muts = ms
    for test = 1:nTests
        test
        
        md = 10;
        s = 1000;
        while md >= 7 || s ~= muts
            [AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(muts,minTheta,maxTheta,maxTimeSimul,'uniform');
            true_data(test).AM = AM;
            true_data(test).timesTrue = timesTrue;
            true_data(test).ratesTrue = ratesTrue;
            true_data(test).stree = stree;
            true_data(test).likelTrue = likelTrue;
            true_data(test).likelEventsTrue = likelEventsTrue;
            true_data(test).aux = aux;
            deg = sum(AM,2)';
            s = size(AM,1);
            intern = find(deg > 0);
            md = max(deg);
        end
            m = size(AM,1);
            G = digraph(AM);
            nodeLabels = cell(1,m);
            for i = 1:m
                nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
            end
            plot(G,'NodeLabel',nodeLabels);
            k = 10;
    end
    save(['test_scite_high_mut_rate_' int2str(muts) '.mat']);
    clear true_data;
end