function [currTime,EMOrders] = getInitTime(mutOrders,minRate,maxRate,m,leafs)
    theta0 = 0.5*(maxRate + minRate);
    currTime = zeros(1,m+1);
    AMOrders = zeros(m,m);
    EMOrders = zeros(m-1,m+1);
    iEdge = 1;
        for i = 1:m
            for j = 2:length(mutOrders{i})
                AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
                EMOrders(iEdge,mutOrders{i}(j)) = 1;
                EMOrders(iEdge,mutOrders{i}(j-1)) = -1;
                iEdge = iEdge + 1;
            end
        end
    G = digraph(AMOrders);
    % plot(G);
    % G = digraph(AM);
    v = bfsearch(G,1);
    for i = 1:m
        currTime(v(i)) = (i-1)/theta0;
    end

    currTime(m+1) = max(currTime(leafs));