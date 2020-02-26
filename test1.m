clear;
theta = 0.01;
tmax = 500;
t = 1:tmax;
p = (theta*t).*exp(-theta*t);
p1 = 50*normpdf(t,1/theta,0.5/theta);
figure
plot(t,p)

y = randsample(tmax,1000,true,p); 
mean(y)
%% 

n = 100;
x = minTheta + rand(1,n)*(maxTheta - minTheta);
y = minTheta + rand(1,n)*(maxTheta - minTheta);

z = 1-abs(x-y)./y;
mean(z)
histogram(z)

%% 
x = 1:n;
y = (log(x)).^2
plot(y)
%% 
mu = (log(maxTheta) + log(minTheta))/2;
sigma = (log(maxTheta) - log(minTheta))/6;
x = lognrnd(mu,sigma,1,100);
exp(mu + sigma^2/2)
histogram(x)
