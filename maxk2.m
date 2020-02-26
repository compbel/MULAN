function [y,i]=maxk2(A,k)
[A,I] = sort(A, 'descend');
y = A(1:k);
i = I(1:k);
end