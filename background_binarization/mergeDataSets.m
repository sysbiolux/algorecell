function [data,genes]=mergeDataSets(data1,genes1,data2,genes2,quantileNorm,tpmNorm);

[genes,i1,i2]=intersect(genes1,genes2);
data=[data1(i1,:), data2(i2,:)];

if quantileNorm
    disp('Quantile normalization ...')
    data=quantilenorm(data);
end

if tpmNorm %including length normalization
    disp('Convert reads to TPM')
    [geneLengths,geneLengthNames]=xlsread('GeneMaxLengths_Julia2019.xlsx');
    % normalize reads by gene length
    [genesWithLength,ia,ib]=intersect(genes,geneLengthNames);
    disp('Length found for genes:')
    numel(genesWithLength)
    lengths=geneLengths(ib);
    [genesWithLength(1:20), num2cell(lengths(1:20))]
    %     [genesWithLength(end-5:end), num2cell(lengths(end-5:end))]
    dataLength=data(ia,:)./repmat(lengths,1,size(data,2));
    % tpm normalization
    dataTPM=1e6*dataLength./repmat(sum(dataLength,1),size(dataLength,1),1);
    [genesWithLength(1:20), num2cell(dataTPM(1:20,1))]
    data=dataTPM;
    genes=genesWithLength;
    %     data=1e6*data./repmat(sum(data,1),size(data,1),1);
end
end