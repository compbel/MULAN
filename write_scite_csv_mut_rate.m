clear;
ms = [30];

alpha = 0.2;
beta = 0.00001;
out_folder = 'scite_data_30/';
multiply = 1;
for muts = ms
    load(['test_scite_high_mut_rate_' int2str(muts) '.mat'])
    for test = 1:nTests
        stree = true_data(test).stree;
        M = streeToMutMatr(stree);
        tmp = M;
        for i=2:multiply
            tmp = [tmp;M];
        end
        M = tmp;
        M = M';
        M_noise = addNoise(M, alpha, beta);
        fid = fopen([out_folder 'mutrate_' int2str(muts) '_' int2str(test) '.csv'], 'w');

        for i=1:size(M_noise,1)
            for j=1:size(M_noise,2)
                fprintf(fid, [int2str(M_noise(i,j)) ' ']);
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end