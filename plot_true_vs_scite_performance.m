% gets the results for 70,,90 mut experiments for scite output and true
% data for mutation estimation rates and plots them on violin graphs
clear;
muts=[70];
Data_rate = [];
Data_rate_true = [];
nTests = 10;
accurRate1 = zeros(1,nTests);
accurRateTrue1 = zeros(1,nTests);
subscript = cell(1,1);
iter = 0;
tests = [1 2 3 4 6 7 8];
% alpha = 0.2
% for i = 1:length(muts)
%     load(['res' num2str(muts(i)) '_true_scite_high_mut_rate_input_exh_2_r.mat']);
%     true_data_output = outputData;
%     load(['res' num2str(muts(i)) 'scite_high_mut_rate_input_exh_2_r.mat']);
%     load(['test_scite_high_mut_rate_' num2str(muts(i)) '.mat'], 'true_data');
%     clear accurRate1;
%     clear accurRateTrue1;
%     for test=tests
%         ratesTrue = true_data(test).ratesTrue;
%         internIntersect = outputData(test).internIntersect;
%         accurRate1(test) = 1 - errperf(ratesTrue(internIntersect),outputData(test).ratesInfer1(internIntersect),'mare');
%         accurRateTrue1(test) = 1 - errperf(ratesTrue(internIntersect),true_data_output(test).ratesInfer1(internIntersect),'mare');
%         Data_rate = [Data_rate; accurRate1(test)];
%         Data_rate = [Data_rate; accurRateTrue1(test)];
%         iter = iter+1;
%         name = ['n=' int2str(muts(i)) ' alpha=0.2'];
%         subscript{iter} = name;
%         iter = iter+1;
%         name = ['n=' int2str(muts(i)) ' true alpha=0.2'];
%         subscript{iter} = name;
%     end
%     mean(nonzeros(accurRate1))
%     mean(nonzeros(accurRateTrue1))
%     
% end
%alpha = 01.
for i = 1:length(muts)
    load(['res' num2str(muts(i)) '_true_scite_input_exh_2_r.mat']);
    true_data_output = outputData;
    load(['res' num2str(muts(i)) 'scite_input_exh_2_r.mat']);
    load(['test_scite_' num2str(muts(i)) '.mat'], 'true_data');
    clear accurRate1;
    clear accurRateTrue1;
    for test=1:nTests
        ratesTrue = true_data(test).ratesTrue;
        internIntersect = outputData(test).internIntersect;
        accurRate1(test) = 1 - errperf(ratesTrue(internIntersect),outputData(test).ratesInfer1(internIntersect),'mare');
        accurRateTrue1(test) = 1 - errperf(ratesTrue(internIntersect),true_data_output(test).ratesInfer1(internIntersect),'mare');
        Data_rate = [Data_rate; accurRate1(test)];
        Data_rate = [Data_rate; accurRateTrue1(test)];
        iter = iter+1;
        name = ['scite tree'];
        subscript{iter} = name;
        iter = iter+1;
        name = ['clean tree'];
        subscript{iter} = name;
    end
    mean(nonzeros(accurRate1))
    mean(nonzeros(accurRateTrue1))
    
end

figure
vs = violinplot(Data_rate, subscript);
xlim([0.5, 2.5]);
ylim([0.6, 1]);
xlabel('Mutations');
ylabel('MAPA');

