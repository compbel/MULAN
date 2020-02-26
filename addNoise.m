function Mnew = addNoise(M,alpha,beta)
n = size(M,1);
m = size(M,2);
Mnew = M;
for i = 1:n
    for j = 1:m
        if sum(Mnew(:,j),1) == 1
            continue;
        end
        if Mnew(i,j) == 1
            p = rand;
            if p <= alpha
                Mnew(i,j) = 0;
            end
        end
        if Mnew(i,j) == 0
            p = rand;
            if p <= beta
                Mnew(i,j) = 1;
            end
        end
    end
end