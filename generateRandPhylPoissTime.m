function [AM,timesMut,rates,stree,likelihood,likelihoodEvents] = generateRandPhylPoissTime(m,minTheta,maxTheta,maxTime)
% stree: [nextChild haplotype parent label frequency timet fitness oldChild rate]
time = 0;
likelihood = 0;
likelihoodEvents = zeros(2,m);
waitTime = zeros(1,m);
timesMut = zeros(1,m);
rates = zeros(1,m);
nCurrMut = 1;
stree(1).children = [];
stree(1).parent = 0;
stree(1).time = 0;
rates(1) = minTheta + rand*(maxTheta - minTheta);
stree(1).rate = rates(1);


while (time <= maxTime) && (nCurrMut < m)
    time = time + 1;
    nCurrMutNew = nCurrMut;
    for i = 1:nCurrMut
        waitTime(i) = waitTime(i)+1;
        prob = exp(-rates(i)*waitTime(i));
        p = rand;
%         if p > prob
        if waitTime(i) > 1/rates(i)
            likelihood = likelihood + log(1-prob);
            likelihoodEvents(1,nCurrMutNew) = log(1-prob);
            waitTime(i) = 0;
            nCurrMutNew = nCurrMutNew + 1;
            stree(i).children = [stree(i).children nCurrMutNew];
            stree(nCurrMutNew).parent = i;
            rates(nCurrMutNew) = minTheta + rand*(maxTheta - minTheta);
            stree(nCurrMutNew).rate = rates(nCurrMutNew);
            stree(nCurrMutNew).time = time;
            timesMut(nCurrMutNew) = time;
        end
    end
    nCurrMut = nCurrMutNew;
end
for i = 1:m
    likelihood = likelihood - rates(i)*waitTime(i);
end
AM = zeros(nCurrMut,nCurrMut);
for i = 1:nCurrMut
    if stree(i).parent ~= 0
        AM(stree(i).parent,i) = 1;
    end
end