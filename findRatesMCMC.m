function [mutRates,times, bestL] = findRatesMCMC(stree,AM,minRate,maxRate,maxTime,eps, defaultRate, totalIter)
m = size(stree,2);
normalNodes = datasample(1:m,round(m/2), 'Replace', false);
highMutRateNodes = setdiff(1:m, normalNodes);
secondRate = normrnd((maxRate+minRate)/2,(maxRate+minRate)/4);
while secondRate < minRate || secondRate > maxRate
    secondRate = normrnd((maxRate+minRate)/2,(maxRate+minRate)/4);
end
iter = 0;
bestL = 10e10;
bestNormal = [];
bestSecondRate = 0;
bestTime = [];
bestOrder = [];
curL = 10e10;

% variables to revert move
toMove = 1;
moveToNormal = true;
oldRate = 0;
while iter < totalIter
    iter = iter +1;
    if mod(iter,10) == 0
        fprintf('Current iteration %d \n', [ iter]);
    end
    move = rand;
    % change rate
    if move < 0.5
        oldRate = secondRate;
        secondRate = normrnd(secondRate,(maxRate+minRate)/8);
        while secondRate < minRate || secondRate > maxRate
            secondRate = normrnd(oldRate,(maxRate+minRate)/8);
        end
    else
        % add to normalNodes
        if (rand < 0.5 && length(normalNodes) < m) || isempty(normalNodes)
            moveToNormal = true;
            toMove = datasample(highMutRateNodes,1);
            normalNodes = [normalNodes toMove];
            highMutRateNodes(highMutRateNodes == toMove) = [];
        else % add to highRate
            moveToNormal = false;
            toMove = datasample(normalNodes,1);
            highMutRateNodes = [highMutRateNodes toMove];
            normalNodes(normalNodes == toMove) = [];
        end
    end
    %assign rate
    for i=normalNodes
        stree(i).rate = defaultRate;
    end
    for i=highMutRateNodes
        stree(i).rate = secondRate;
    end
    
    [mutOrdersI,newTime, newL] = findRatesRec(stree,minRate,maxRate,maxTime,eps, AM);
    %save move
    if rand < exp(curL - newL)
        if newL < bestL
            bestL = newL;
            bestNormal = normalNodes;
            bestSecondRate = secondRate;
            bestTime = newTime;
            bestOrder = mutOrdersI;
            fprintf('New L %f iter %d \n', [bestL iter]);
        end
        curL = newL; 
    else %revert move
        if move < 0.5
            secondRate = oldRate;
        else
            if moveToNormal
                normalNodes(normalNodes == toMove) = [];
                highMutRateNodes = [highMutRateNodes toMove];
            else
                highMutRateNodes(highMutRateNodes == toMove) = [];
                normalNodes = [normalNodes toMove];
            end
        end
    end
end

mutRates = bestOrder;
times = bestTime;