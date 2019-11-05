function [mutRates,times, bestL,secondBestRate, bestSecMutations] = findRatesMCMC2(stree,AM,minRate,maxRate,maxTime,eps, defaultRate, totalIter)
m = size(stree,2);
switchNodes = datasample(1:m,2, 'Replace', false);
secondRate = normrnd((maxRate+minRate)/2,(maxRate+minRate)/4);
while secondRate < minRate || secondRate > maxRate
    secondRate = normrnd((maxRate+minRate)/2,(maxRate+minRate)/4);
end
iter = 0;
bestL = 10e10;
bestSwitch = [];
bestSecondRate = 0;
bestTime = [];
bestOrder = [];
curL = 10e10;

% variables to revert move
toMove = 1;
moveToNormal = true;
oldRate = 0;
oldNodes = [];
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
        %change nodes
        oldNodes = switchNodes;
        switchNodes = datasample(1:m,2, 'Replace', false);
    end
    %assign rate
    % first vertex, second 1 - default rate, 0 - second rate
    rates = zeros(1,m);
    q = [1 ~ismember(1, switchNodes)];
    idx = 1;
    while idx <= length(q)
        v = q(idx,:);
        if v(2) == 1
            stree(v(1)).rate = defaultRate;
            rates(v(1)) = defaultRate;
        else
            stree(v(1)).rate = secondRate;
            rates(v(1)) = secondRate;
        end
        
        for c=stree(v(1)).children
            flip = ismember(c, switchNodes );
            if flip
                q = [q ; c ~v(2)];
            else
                q = [q; c v(2)];
            end
        end
        idx = idx+1;
    end
%     for i=normalNodes
%         stree(i).rate = defaultRate;
%     end
%     for i=highMutRateNodes
%         stree(i).rate = secondRate;
%     end
    
    [mutOrdersI,newTime, newL] = findRatesRec(stree,rates, minRate,maxRate,eps, AM);
    %save move
    if rand < exp(curL - newL)
        if newL < bestL
            bestL = newL;
            bestSecMutations = switchNodes;
            secondBestRate = secondRate;
            bestTime = newTime;
            bestOrder = mutOrdersI;
            fprintf('New L %f iter %d \n', [bestL iter]);
        end
        curL = newL; 
    else %revert move
        if move < 0.5
            secondRate = oldRate;
        else
            switchNodes = oldNodes;
        end
    end
end

mutRates = bestOrder;
times = bestTime;