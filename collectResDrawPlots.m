clear;


muts=[20 25 30];
Data_err = [];
subscript = cell(1,1);
iter = 0;
for i = 1:length(muts)
    outfile_res = ['res' num2str(muts(i)) 'exh.mat'];
    load(outfile_res);
    for j = 1:length(accurRate)
        iter = iter + 1;
        Data_err = [Data_err; accurRate(j); accurRate(j)];
        name = ['n=' int2str(muts(i))]
        subscript{iter} = name;
        iter = iter + 1;
        subscript{iter} = name;
    end
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
vs = violinplot(Data_err, subscript);
xlim([0.5, 3.5]);
ylim([0, 1]);
xlabel('Mutations');
ylabel('MRA');
