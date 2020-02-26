nPerm = 1;
T = 200;
range = 1:T;

fun = zeros(nPerm,length(range));
true_fun = zeros(nPerm,length(range));
funH1 = zeros(nPerm,length(range));
violRateOrder = zeros(nPerm,1);
m = 20;
minTheta = 1/25;
maxTheta = 1/8;
maxTime = 1000;
[AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(m,minTheta,maxTheta,maxTime,'binary');

m = size(AM,1);
G = digraph(AM);
nodeLabels = cell(1,m);
for i = 1:m
    nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
end
figure
plot(G,'NodeLabel',nodeLabels);
M = streeToMutMatr(stree);
writeMutAdDNA(M, 'test_Nex_bid_theta_diff_same_pos.nex',[], minTheta);