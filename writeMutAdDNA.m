function writeMutAdDNA(M, file, rates, minTheta)
c = size(M,1);
m = size(M,2);
fid = fopen(file, 'w');
fprintf(fid, '#NEXUS\n');

fprintf(fid,'Begin data;\n');
fprintf(fid,'	Dimensions ntax=%d nchar=%d;\n', c, 3*m+6);
fprintf(fid,'	Format datatype=nucleotide gap=-;\n');
fprintf(fid,'	Matrix\n');
start_codon = 'ATG';
major = 'CCC';
minor = 'CCA';
pos1_major = 'TTA';
pos1_minor = 'CTA';
for i=1:c
    fprintf(fid,'c_%d\tATG', i);
    for j=1:m
        if isempty(rates) || rates(j) == minTheta
            if M(i,j) == 0
                fprintf(fid, major);
            else
                fprintf(fid, minor);
            end
        else
            if M(i,j) == 0
                fprintf(fid, pos1_major);
            else
                fprintf(fid, pos1_minor);
            end
        end
    end
    fprintf(fid, 'TAA\n');
end
fprintf(fid,';\nEnd;');