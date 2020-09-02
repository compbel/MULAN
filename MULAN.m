disp('Starting MULAN');
addpath('functions');

% filePath = 'hou/f_noRep/dataHou18_map0.gv';
% nMutations=18;

if (~exist('minTheta', 'var'))
    minTheta = 1e-5;
end
if (~exist('maxTheta', 'var'))
    maxTheta = 1e-4;
end
maxTime = 100000;
if (~exist('eps', 'var'))
    eps = 0.001;
end

if (~exist('outputPrefix', 'var'))
    output = 'out';
end


eps1 = 0.001;
maxIter = 10;

[stree, AM] = readScite(filePath,nMutations);

[ratesInfer,timesInfer,ordersInfer,likelInfer,stopFlag] = findRatesEMorder4(stree,AM,minTheta,maxTheta,maxTime,maxIter,eps,eps1);

fileID = fopen([output '.txt'],'w');
fprintf(fileID,'likelihood:\n');
fprintf(fileID,'%6.4f\n ', likelInfer);
fprintf(fileID,'rates:\n');
fprintf(fileID,'%6.6f ', ratesInfer);
fprintf(fileID,'\n');
fprintf(fileID,'times:\n');
fprintf(fileID,'%6.4f ', timesInfer);
fprintf(fileID,'\n');
fprintf(fileID,'order:\n');
for i=1:length(ordersInfer)
    fprintf(fileID,'%4d ', ordersInfer{i});
    fprintf(fileID,'\n');
end
fclose(fileID);
save([output '.mat'],'stree','ratesInfer','timesInfer','ordersInfer','likelInfer','stopFlag');
