% clear
% minTheta = 1/10;
% maxTheta = 1/2;
% for i = 1:3
%     r(i) = minTheta + rand*(maxTheta - minTheta);
% end
% R = r(1) + r(2) + r(3);
% 
% f = @(x)(-(r(2)*x(1) + r(3)*x(2) + log(x(1)) + log(x(2)-x(1)) - R*x(3)));
% A = [1 -1 0; 0 1 -1];
% rhs = [0 0];
% initx = [0.25 0.5 0.75];
% lb = [0 0 0];
% ub = [100 100 100];
% options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
% [sol,val] = fmincon(f,initx,A,rhs,[],[],lb,ub,[],options);
% opt = -val
% 
% T = 100;
% range = 1:T;
% fun1 = zeros(1,length(range));
% fun2 = zeros(1,length(range));
% stat = zeros(1,length(range));
% A1 = [1 -1; 0 1];
% A2 = [-1 1; 1 0];
% lb = [0 0];
% initx1 = [0.5 0.75];
% initx2 = [0.75 0.5];
% options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
% for i = 1:length(range)
%     i
%     h = range(i);
%     ub = [h h];
%     rhs = [0 h];
%     f1 = @(x)(-(r(2)*x(1) + r(3)*x(2) + log(x(1)) + log(x(2)-x(1)) - R*h));
%     f2 = @(x)(-(r(2)*x(1) + r(3)*x(2) + log(x(2)) + log(x(1)-x(2)) - R*h));
% %     f1 = @(x)(-(r(2)*x(1) + r(3)*x(2) + log(x(1)) + log(x(2)-x(1))));
%     [sol1,val1] = fmincon(f1,initx1,A1,rhs,[],[],lb,ub,[],options);
%     [sol2,val2] = fmincon(f2,initx2,A2,rhs,[],[],lb,ub,[],options);
%     fun1(i) = val1;
%     fun2(i) = val2;
%     
% %     x = sol1;
% %     stat(i) = r(2)*x(1) + r(3)*x(2) + log(x(1)) + log(x(2)-x(1));
% end
% % figure
% % plot(range,-fun1,range,-fun2)
% 
% opth = -fun1;
% figure
% plot(range,opth)
% 
% % 
% % figure
% % plot(range,-(stat - R*range))
%% 
clear;
nPerm = 1;
T = 200;
range = 1:T;

fun = zeros(nPerm,length(range));
true_fun = zeros(nPerm,length(range));
funH1 = zeros(nPerm,length(range));
violRateOrder = zeros(nPerm,1);
m = 20;
minTheta = 1/20;
maxTheta = 1/10;
maxTime = 1000;
[AM,timesTrue,ratesTrue,stree,likelTrue,likelEventsTrue,aux] = generateRandPhylPoissTime1(m,minTheta,maxTheta,maxTime);
m = size(AM,1);
G = digraph(AM);
nodeLabels = cell(1,m);
for i = 1:m
    nodeLabels{i} = ['(' int2str(i) ',' num2str(round(stree(i).rate,3)) ',' int2str(stree(i).time) ')'];
end
figure
plot(G,'NodeLabel',nodeLabels);

for ord = 1:nPerm
    mutOrders = cell(1,m);
    currTheta = extractfield(stree,'rate');
    sumTheta = sum(currTheta,2)- currTheta(1);
    for i = 1:m
        if ord == 1
            mutOrders{i} = [i stree(i).children];
        else
            pm = randperm(length(stree(i).children));
            mutOrders{i} = [i stree(i).children(pm)];
        end
    end
    
    
    AMOrders = zeros(m,m);
    for i = 1:m
        for j = 2:length(mutOrders{i})
            AMOrders(mutOrders{i}(j-1),mutOrders{i}(j)) = 1;
            if (j > 2) && (currTheta(mutOrders{i}(j-1)) < currTheta(mutOrders{i}(j)))
                violRateOrder(ord) = violRateOrder(ord) + 1;
            end
        end
    end
    G = digraph(AMOrders);
%     figure
%     plot(G)   
    v = bfsearch(G,1);

    A = zeros(m-1,m);
    rhs = zeros(m-1,1);
    iEq = 0;
    for i = 1:m
        if ~isempty(stree(i).children)
            for j = 2:length(mutOrders{i})
                iEq = iEq + 1;
                A(iEq,mutOrders{i}(j-1)) = 1;
                A(iEq,mutOrders{i}(j)) = -1;
            end
        end
    end
    lb = zeros(1,m);

    % currTime = extractfield(stree,'time');
    initTime1 = linspace(0,1,m);
    initTime = initTime1;
    for i = 1:m
        initTime(v(i)) = initTime1(i);
    end
%     initTime(1) = 0.0001;
    initTime(v(end)) = 1 - 0.0001;

    sols = zeros(length(range),m);
    solDiffs = zeros(length(range),m);
    options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false);
    for i = 1:length(range)
        i
        h = range(i);
%         ub = h*ones(1,m);
        ub = ones(1,m);
        ub(1) = 0;
        l = @(t)likelihoodRatesOrdersH1(t,currTheta,mutOrders,h);
%         l = @(t)likelihoodRatesOrders2(t,currTheta,mutOrders);
%         [time,newL] = fmincon(l,initTime,A,rhs,[],[],lb,ub,[],options);
        [time,newL] = fmincon(l,initTime,[],[],[],[],lb,ub,[],options);
        fun(ord,i) = -newL;
        funH1(ord,i) = -likelihoodRatesOrdersH1(time,currTheta,mutOrders,1);
        true_fun(ord,i) = -likelihoodRatesOrdersH1(normalize(extractfield(stree,'time'),'range'),currTheta,mutOrders,h);
        sols(i,:) = time;
        z = 1;
        for j = 1:m
            for u = 2:length(mutOrders{j})
                solDiffs(i,z) = (time(mutOrders{j}(u)) - time(mutOrders{j}(u-1)));
                z = z + 1;
            end
        end
    end
end   

figure
plot(range,fun(1,:))

% figure
% plot(range,fun(1,:),range,fun(2,:),range,fun(3,:),range,fun(4,:),range,fun(5,:))


myfittype = fittype('a*x - b*log(x) + c',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c'})
[myfit1,gof] = fit(range',fun(1,:)',myfittype)
% [myfit2,gof] = fit(range',fun(2,:)',myfittype)
%     plot(myfit,vertSizes',mean(hdAll,2))
approx = myfit1.a*range - myfit1.b*log(range) + myfit1.c;
% figure('NumberTitle', 'off', 'Name', 'Fit approx')
plot(range,fun(1,:),range,approx,'--')

mean(fun(1,:) - approx) 
% figure
% plot(range,sols(:,5))

% 
% figure
% plot(range,sols(:,2))

% figure('NumberTitle', 'off', 'Name', 'Theta approx')
% plot(range,fun(1,:), range, sumTheta*range)
% 
% figure('NumberTitle', 'off', 'Name', 'True time approx')
% plot(range,true_fun(1,:), range, sumTheta*range)

myfit1.a/sumTheta
myfit1.a/sum(currTheta,2)

mean(fun(1,:) - sumTheta*range) 

mean(true_fun(1,:) - sumTheta*range) 

aux = zeros(1,length(range));
aux1 = zeros(1,length(range));
aux3 = zeros(1,length(range));
for i = 1:length(range)
    h = range(i);
    r = currTheta(2);
    t2 = (h*r-3 + sqrt((h*r-3)^2+4*h*r))/(2*h*r);
    aux(i) = t2;
    aux1(i) = h*r*t2 + log(t2) + 2*log(1-t2) + (currTheta(3) + currTheta(4))*h;
%     aux1(i) = (log(t2) + 2*log(1-t2))/h;
end
figure
plot(range,aux,range,sols(:,2))
figure
plot(range,aux1)

figure
plot(range,fun(1,:)./range)

figure
plot(range,range - 3*log(range))

figure
plot(range,funH1(1,:))

aux = [sols(2:end,:); zeros(1,m)] - sols;
aux1 = [fun(1,2:end) 0] - fun(1,:);

[maxL,indL] = max(fun(1,:));
figure
plot(range,maxL - fun(1,:))



%% 
h = 100;
theta1 = 1/10;
theta2 = 1/20;
[x,y] = meshgrid(0:1:h);  
z = theta1*x + theta2*y + log(x) + log(abs(y-x)); 
surf(x,y,-z)
%% 
mus = zeros(size(sols,1),size(sols,2));
for i = 1:length(range)
        i
        h = range(i);
        for j = 1:size(solDiffs,2)
            mus(i,j) = 1/solDiffs(i) + h*sum(currTheta((j+1):end),2);
        end
end
%% 
mus = zeros(1,length(range));
mus1 = zeros(1,length(range));
for q = 1:length(range)
    h = range(q);
    t2 = sols(q,2);
    t3 = sols(q,3);
    r2 = currTheta(2);
    r3 = currTheta(3);
    r4 = currTheta(4);
    h*r2 - 1/(1-t2) +1/t2 - 1/(t3-t2)
    mus(q) = h*r3 + 1/(t3-t2);
end

figure
plot(range,mus,range,(r2+r3+r4)*range+1)
%% 


mus = zeros(1,length(range));
mus1 = zeros(1,length(range));
for q = 1:length(range)
    h = range(q);
    t2 = sols(q,2);
    t3 = sols(q,3);
    t4 = sols(q,4);
    r2 = currTheta(2);
    r3 = currTheta(3);
    r4 = currTheta(4);
    mus(q) = h*r4 + 1/(t4-t3);
end

figure
plot(range,mus,range,(r2+r3+r4)*range+1)

myfittype = fittype('a*x - b*log(x) + c',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b','c'})
[myfit,gof] = fit(range',fun',myfittype)


h*r2 + 1/t2 - 1/(t3-t2)
h*r3 + 1/(t3-t2) - 1/(t4-t3)

