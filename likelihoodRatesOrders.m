function f = likelihoodRatesOrders(t,theta,mutOrders)
f = 0;
m = length(t);
for i = 1:m
    for j = 2:length(mutOrders{i})
        u = mutOrders{i}(j-1);
        v = mutOrders{i}(j);
        f = f - log(1-exp(-theta(i)*(t(v) - t(u))));
    end
end
z = 5;

% if nargout > 1 % gradient required
%     gradf = zeros(1,m);
%     for i = 1:m
%         gradf(i) = theta(i) - deg(i)/t(i);
%     end
% end
