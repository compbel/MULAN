function [AM,timesMut,rates,stree,likelihood,likelihoodEvents,aux] = generateRandPhylPoissTime1(m,minTheta,maxTheta,maxTime,prior)
% stree: [nextChild haplotype parent label frequency timet fitness oldChild rate]

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% maxTheta = minTheta;
mu = (log(maxTheta) + log(minTheta))/2;
sigma = (log(maxTheta) - log(minTheta))/6;
h = (maxTheta-minTheta)/(m-1);
time = 0;
likelihood = 0;
likelihoodEvents = zeros(2,m);
aux = zeros(2,m);
waitTime = zeros(1,m);
waitTimeExp = zeros(1,m);
timesMut = zeros(1,m);
rates = zeros(1,m);
nCurrMut = 1;
stree(1).children = [];
stree(1).parent = 0;
stree(1).time = 0;
if strcmp(prior,'uniform')
    rates(1) = minTheta + rand*(maxTheta - minTheta);
end
if strcmp(prior,'lognormal')
    rates(1) = lognrnd(mu,sigma);
end
if strcmp(prior,'unidiscr')
    rates(1) = minTheta;
end
if strcmp(prior,'binary')
    rates(1) = minTheta;
end
stree(1).rate = rates(1);
% waitTimeExp(1) = exprnd(1)/rates(1);
waitTimeExp(1) = 1/rates(1);

while (time <= maxTime) && (nCurrMut < m)
    time = time + 1;
    nCurrMutNew = nCurrMut;
    for i = 1:nCurrMut
        waitTime(i) = waitTime(i)+1;
        if waitTime(i) >= waitTimeExp(i)
            nCurrMutNew = nCurrMutNew + 1;
            likelihood = likelihood + log(rates(i)*waitTime(i)) - rates(i)*waitTime(i);
            likelihoodEvents(1,nCurrMutNew) = log(rates(i)*waitTime(i)) - rates(i)*waitTime(i);
            aux(1,nCurrMutNew) = waitTime(i);
            aux(2,nCurrMutNew) = rates(i);
            waitTime(i) = 0;
            stree(i).children = [stree(i).children nCurrMutNew];
            stree(nCurrMutNew).parent = i;
            if strcmp(prior,'uniform')
                rates(nCurrMutNew) = minTheta + rand*(maxTheta - minTheta);
            end
            if strcmp(prior,'lognormal')
                rates(nCurrMutNew) = lognrnd(mu,sigma);
            end
            if strcmp(prior,'unidiscr')
                rates(nCurrMutNew) = minTheta + (nCurrMutNew-1)*h;
            end
            if strcmp(prior,'binary')
                if rand < 0.5
                    rates(nCurrMutNew) = minTheta;
                else
                    rates(nCurrMutNew) = maxTheta;
                end
            end 
            %
%             z1 = minTheta;
%             z2 = maxTheta;
%             minTheta = z2;
%             maxTheta = 2*z2-z1;
            
            stree(nCurrMutNew).rate = rates(nCurrMutNew);
            stree(nCurrMutNew).time = time;
            timesMut(nCurrMutNew) = time;
%             waitTimeExp(i) = exprnd(1)/rates(i);
%             waitTimeExp(nCurrMutNew) = exprnd(1)/rates(nCurrMutNew);
            waitTimeExp(i) = 1/rates(i);
            waitTimeExp(nCurrMutNew) = 1/rates(nCurrMutNew);
%             waitTimeExp(i) = -log(rand)/rates(i);
%             waitTimeExp(nCurrMutNew) = -log(rand)/rates(nCurrMutNew);
%             waitTimeExp(i) = mean(-log(rand(1,200))/rates(i));
%             waitTimeExp(nCurrMutNew) = mean(-log(rand(1,200))/rates(nCurrMutNew));
        end
    end
    nCurrMut = nCurrMutNew;
end
for i = 1:m
    if waitTime(i) ~= 0
        likelihood = likelihood  - rates(i)*waitTime(i);
        likelihoodEvents(2,nCurrMutNew) = -rates(i)*waitTime(i);
    end
end
AM = zeros(nCurrMut,nCurrMut);
for i = 1:nCurrMut
    if stree(i).parent ~= 0
        AM(stree(i).parent,i) = 1;
    end
end