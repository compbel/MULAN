% [stree, AM] = readScite('scite/f_noRep/dataHou18_map0.gv',18);
% m = size(AM,1);
% G = digraph(AM);
% nodeLabels = cell(1,m);
% for i = 1:m
%     nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
% end
% plot(G,'NodeLabel',nodeLabels);
L = zeros(1, 19);

% 1-(1-10^-6)^18
minTheta = 1e-5;
maxTheta = 1e-4;
maxTime = 100000;
eps = 0.001;
mutRate = minTheta;
% for i=1:length(stree)
%     stree(i).rate = mutRate;
% end
% [mutOrders(19),times(19), L(19)] = findRatesMCMC(stree,AM,minTheta,maxTheta,maxTime,eps, minTheta, 3000);
% [mutOrders,times, L(19)] = findRatesRec(stree,minTheta,maxTheta,maxTime,eps, AM);
treesToCalculate = [1,2,16,18,19];
mutOrders = cell(1,1);
times = [];
secRate = [];
secMut = [];
for iter=19
    disp(iter);
    if iter == 19
        [stree, AM] = readScite('scite/f_noRep/dataHou18_map0.gv',18);
    else
        [stree, AM] = readScite(['scite/f_rep' int2str(iter) '/dataHou18_map0.gv' ],19);
    end
    mutRate = minTheta;
    for i=1:length(stree)
        stree(i).rate = mutRate;
    end
    %     [mutOrders,times, L(iter)] = findRatesRec(stree,minTheta,maxTheta,maxTime,eps, AM);
    [mo,ti, L(iter), secRate(iter), sm] = findRatesMCMC2(stree,AM,minTheta,maxTheta,maxTime,eps, minTheta, 2000);
    mutOrders = [mutOrders mo];
    if iter == 19
        ti = [ti 0];
    end
    times = [times ; ti];
    secMut = [secMut; sm];
    
end
L
% save('scite_results_2rates.mat')
