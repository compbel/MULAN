function f = likelihoodRatesOrders2(t,theta,mutOrders)
f = 0;
m = length(t);
t1 = max(t) - t;
for i = 1:m
    for j = 2:length(mutOrders{i})
        u = mutOrders{i}(j-1);
        v = mutOrders{i}(j);
        timeInt = t1(u) - t1(v);
        f = f - log(theta(i)*timeInt);
    end
    f = f + theta(i)*t1(i);
end

z = 5;

% if nargout > 1 % gradient required
%     gradf = zeros(1,m);
%     for i = 1:m
%         gradf(i) = theta(i) - deg(i)/t(i);
%     end
% end
