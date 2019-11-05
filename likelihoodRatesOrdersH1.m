function f = likelihoodRatesOrdersH1(t,theta,mutOrders,h)
f = 0;
m = length(t);
for i = 1:m
    for j = 2:length(mutOrders{i})
        u = mutOrders{i}(j-1);
        v = mutOrders{i}(j);
        timeInt = t(v) - t(u);
        f = f - log(timeInt);
    end
    f = f - h*theta(i)*t(i);
end
f = f + h*sum(theta,2) - m*log(h);
