% finds likelihood of all trees for gawad study based on scite output
clear;
% ms = [10];
nTests = 10;
% minTheta = 1/50;
% maxTheta = 1/10;
minTheta = 1/200;
maxTheta = 1/100;
maxTime = 100000;
% eps = 0.01;
eps = 0.001;
eps1 = 0.001;
maxIter = 10;

pat_size = [20 16 49 78 105 10];
runTimes1 = zeros(1,sum(pat_size));
test = 1;
result = [];
for pat_id=4:5
    s = pat_size(pat_id);
    for mut=0:s
        sprintf('id %d, mut %d',pat_id, mut)
        m_suff = int2str(mut);
        total_mut = s+1;
        if mut==0
            m_suff = '';
            total_mut = s;
        end
        
        
        [stree, AM] = readSciteMutTree(['gawad/pat_' int2str(pat_id) '_cluster_dat' m_suff '_map0.gv' ],total_mut);
        G = digraph(AM);
        m = size(AM,1);
        nodeLabels = cell(1,m);
        for i = 1:m
            nodeLabels{i} = [ int2str(i)];
        end
%         plot(G,'NodeLabel',nodeLabels);
                tic
        [ratesInfer1,timesInfer1,ordersInfer1,likelInfer1,stopFlag1] = findRatesEMorder4(stree,AM,minTheta,maxTheta,maxTime,maxIter,eps,eps1);
        result(test).runTimes = toc;
        result(test).pat_id = pat_id;
        result(test).mut = mut;
        result(test).ratesInfer1 = ratesInfer1;
        result(test).timesInfer1 = timesInfer1;
        result(test).ordersInfer1 = ordersInfer1;
        result(test).likelInfer1 = likelInfer1;
        result(test).stopFlag1 = stopFlag1;
        test = test+1;
    end
    
end
save('gawad_mut_rates_results45.mat','result');