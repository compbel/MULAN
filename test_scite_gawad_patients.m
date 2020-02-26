
L = zeros(1, 6);
Lsingle = zeros(1,6);

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
mutOrders = cell(1,1);
times = [];
secRate = [];
secMut = [];
pat_size = [20 16 49 78 105 10];
vals = [1 2 6];
parfor k=1:3
    iter = vals(k);
    if iter == 5 || iter == 4 || iter == 3
        continue;
    end
    disp(iter);

    [stree, AM] = readScite(['gawad/pat_' int2str(iter) '_cluster_dat_map0.gv' ],pat_size(iter));

    mutRate = minTheta;
    rates = zeros(1,length(stree));
    for i=1:length(stree)
        stree(i).rate = mutRate;
        rates(i) = mutRate;
    end
    [mutOrders,times, Lsingle(k)] = findRatesRec(stree,rates, minTheta,maxTheta,eps, AM);
    [mo,ti, L(k), secRate(k), sm] = findRatesMCMC2(stree,AM,minTheta,maxTheta,maxTime,eps, minTheta, 2000);
    mutOrders = [mutOrders mo];
    %times = [times ; ti];
    secMut = [secMut; sm];
    
end
Lsingle
L
save('scite_results_gawad_patients.mat')
