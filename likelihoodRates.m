function [f,gradf] = likelihoodRates(t,theta,deg)
f = 0;
m = length(t);
for i = 1:m
    f = f + theta(i)*t(i)-deg(i)*log(theta(i)*t(i));
end

if nargout > 1 % gradient required
    gradf = zeros(1,m);
    for i = 1:m
        gradf(i) = theta(i) - deg(i)/t(i);
    end
end
