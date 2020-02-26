clear;


muts=[70 90 110 130 150];
Data_rate = [];
Data_time = [];
Data_order = [];
Data_runtime = [];
subscript = cell(1,1);
iter = 0;
for i = 1:length(muts)
    outfile_res = ['res' num2str(muts(i)) 'exh_2.mat'];
    load(outfile_res);
%     mean(accurRate)
%     median(accurRate)
    for j = 1:length(accurRate1)
        iter = iter + 1;
        Data_rate = [Data_rate; accurRate1(j)];
        Data_time = [Data_time; accurTime1(j)];
        Data_order = [Data_order; accurOrder1(j)];
        name = ['n=' int2str(muts(i))];
        subscript{iter} = name;
%         iter = iter + 1;
%         subscript{iter} = name;
    end
%     runTimes1 = sort(runTimes1);
    Data_runtime = [Data_runtime; mean(runTimes1)];
end

% figure
% subplot(1,2,1);
% vs = violinplot(Data_err, subscript);
% ylabel('Accuracy');
% xlim([0.5, 3.5]);
% subplot(1,2,2);
% vs = violinplot(Data_corr, subscript);
% ylabel('Spearman correlation');
% xlim([0.5, 3.5]);

figure
vs = violinplot(Data_rate, subscript);
xlim([0.5, 5.5]);
ylim([0.6, 1]);
xlabel('Mutations');
ylabel('MAPA');

figure
vs = violinplot(Data_time, subscript);
xlim([0.5, 5.5]);
ylim([0.6, 1]);
xlabel('Mutations');
ylabel('MAPA');

figure
vs = violinplot(Data_order, subscript);
xlim([0.5, 5.5]);
ylim([0.6, 1]);
xlabel('Times');
ylabel('Kendall tau');

figure
% area(muts,Data_runtime,'LineWidth',2)
plot(muts,Data_runtime,'LineWidth',2)
xlabel('Mutations');
ylabel('Time(seconds)');
corr((muts.^2)',Data_runtime)