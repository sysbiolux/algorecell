% load TPM datasets (query and background)
quantileNorm=1
tpmNorm=1 %including length norm

disp('Loading data ...')
load('dataTPM_query_v2.mat') %Here only corrected raw counts
if quantileNorm
    data1=dataReads;
    genes1=genes;
else
    data1=dataTPM;
    genes1=genesWithLength;
end
load('dataTPM_Background.mat')
if quantileNorm
    temp1=(sum(dataReads)>1e7) & (sum(dataReads)<1e8); %read counts between 10 Mio and 100 Mio reads
    sum(temp1)
    disp('Pre quantile normalization of background data ...')
    NormData = quantilenorm(dataReads);
    temp2=median(NormData)>=1; %median read count >= 1 to filter out unusual distributions
    sum(temp2)
    selSamples=temp1.*temp2;
    sum(selSamples)
    data2=dataReads(:,find(selSamples));
    genes2=genes;
else
    data2=seldataTPM;
    genes2=genesWithLength;
end

[data,genes]=mergeDataSets(data1,genes1,data2,genes2,quantileNorm,tpmNorm);
disp('Data loaded and merged ...')

% save datagenes_final data genes qsamples TFs selSamples
%% re-load preprocessed data & genes (also includes manually added qsamples & TFs)
load('datagenes_final')
disp('Loading: Done!')

%% check tissue coverage in selected samples
disp('check tissue coverage in selected samples')
load SamplesTissues2666.mat
selnames=samples(find(selSamples));
selTissues=tissues(find(selSamples))
[C,ia,ic]=unique(selTissues); numel(C)
% [C,ia,ic]=unique(tissues); numel(C) %all tissues
%
figure, hist(ic,numel(C));
title('Selected samples per tissue')
[n,centers]=hist(ic,numel(C));
text(centers,n+1,C,'Color','r')

%% check tpm distributions of selected samples
disp('check tpm distributions of selected samples')
samplesToPlot=1:size(data,2);
boxesPerPlot=100
Nrsplot=ceil(numel(samplesToPlot)/boxesPerPlot)
for counter=1:Nrsplot
    %     figure
    figure('units','normalized','outerposition',[0 0 1 1])
    temp=samplesToPlot(((counter-1)*boxesPerPlot+1):min(numel(samplesToPlot),((counter-1)*boxesPerPlot+boxesPerPlot)));
    boxplot(log10(data(:,temp)+1))
    title('TPM distributions of selected samples (pseudo, log10)')
end

%% check tpm distributions of selected genes
disp('check tpm distributions of selected genes')
% highly expressed mouse housekeepers: https://www.nature.com/articles/s41598-017-04520-z/tables/2
% genesToPlotNames=sort({'HK1','HK2','HK3','AKT1','ARF1','EIF1','UBL5','VCP','RUNX2','OSR1','AHR'});
% extract TF data
disp('Mouse TFs:')
temp=tdfread('Mus_musculus_TF.txt');
TFs=upper(cellstr(temp.Symbol));
TFs=unique([TFs', 'LPL','CEBPA','SP7','ADIPOQ','FABP4','ALPL','HEY1']');
numel(TFs)
% genesToPlotNames=sort(TFs(1:end))

higha=0 % 75-100 // 0-70 % 75-100 // 0-100
lowa=75
higho=75
lowo=100
TFsonly=1
minmax=0 % 0=mean, 1=minmax, 2=all/individual (to be implemented)

uca=prctile(data(:,34:end)',higha)';
lca=prctile(data(:,34:end)',lowa)';
uco=prctile(data(:,34:end)',higho)';
lco=prctile(data(:,34:end)',lowo)';
contr=data(:,[1,12,23]);
adip=data(:,[2:6,13:17,24:28]);
osteo=data(:,[7:11,18:22,29:33]);

disp('Selected Genes:');

if minmax==1 %min/max
    temp1=max(adip,[],2)>uca;
    sum(temp1)
    temp2=min(adip,[],2)<lca;
    sum(temp2)
    temp3=max(osteo,[],2)>uco;
    sum(temp3)
    temp4=min(osteo,[],2)<lco;
    sum(temp4)
elseif minmax==2 %all/individual (to be implemented)
    
else %mean
    temp1=mean(adip,2)>uca;
    sum(temp1)
    temp2=mean(adip,2)<lca;
    sum(temp2)
    temp3=mean(osteo,2)>uco;
    sum(temp3)
    temp4=mean(osteo,2)<lco;
    sum(temp4)
end

selGenes=temp1.*temp2.*temp3.*temp4;
sum(selGenes)
genesToPlotNames=sort(genes(find(selGenes)))

% candidate marker genes:
% genesToPlotNames=sort({'PPARG','LPL','CEBPA','RUNX2','SP7','BGLAP','BMP4','ERBB2','ADIPOQ','FABP4','ALPL','AHR','SPARC','SPP1','ALPI','ALPPL2','AKP3','GLIS1','MEF2C','PLAGL1','RUNX3','DLX5','DLX6','HEY1','ATF3','KLF4','RBPJ','IRF1'})
% marker genes:
genesToPlotNames=sort({'LPL','CEBPA','SP7','ADIPOQ','FABP4','ALPL','HEY1'})
genesToPlotNames=sort({'MYC','HEY1','TRP63','GATA3','RORA','SMAD3','TCF7L2','PAX8'}) %predictedToBeUPregulated
genesToPlotNames=sort({'SP1','CEBPA','ESR1','NR2F2','XBP1','E2F3'}) %predictedToBeDOWNregulated
% genesToPlotNames=union(permDOWN,tempDOWN) %tempUP %load: nodeUPDOWN_perturbations_Loic_0120.mat
% genesToPlotNames=suspects %tempUP %load: suspects.mat
% genesToPlotNames=sort({'HEY1','ESR1','SP1'})
% genesToPlotNames=sort({'LPL','SP7','FABP4','HEY1','COL4A6'})
% genesToPlotNames=sort({'SOX2','NANOG','POU5F1'})
% genesToPlotNames={'AHR', 'AR', 'ARNT', 'ARNTL', 'ATF2', 'ATF4', 'CEBPA', 'CEBPB', 'CLOCK', 'CTCF', 'CUX1', 'DBP', 'DDIT3', 'E2F1', 'E2F3', 'E2F4', 'E2F8', 'EGR1', 'EGR3', 'ESR1', 'ESR2', 'ETS1', 'FOS', 'FOXA2', 'FOXH1', 'FOXM1', 'FOXO1', 'FOXO3', 'FOXP3', 'GATA2', 'GATA3', 'GATA6', 'GFI1', 'GLI1', 'HEY1', 'HMGB2', 'HNF1A', 'HNF4A', 'ID1', 'ID3', 'JUN', 'JUNB', 'KLF2', 'KLF4', 'LEF1', 'MAFG', 'MEF2B', 'MEIS1', 'MYC', 'NANOG', 'NCOR2', 'NFATC3', 'NFE2', 'NFE2L2', 'NFIA', 'NFKB1', 'NFYA', 'NFYB', 'NFYC', 'NKX2-1', 'NR0B2', 'NR1D1', 'NR1I2', 'NR2F2', 'NR3C1', 'NR4A1', 'ONECUT1', 'PAX8', 'PBX1', 'PDX1', 'POU5F1', 'REL', 'RELA', 'RELB', 'RORA', 'RUNX1', 'RXRA', 'SALL4', 'SMAD1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD5', 'SMAD7', 'SOX2', 'SOX4', 'SOX9', 'SP1', 'SP3', 'SREBF1', 'STAT3', 'STAT5A', 'STAT5B', 'TCF7L2', 'TFAP2C', 'TGIF1', 'TRP53', 'TRP63', 'TTF1', 'USF1', 'USF2', 'XBP1', 'YY1', 'ZFP281', 'ZSCAN10'}
% genesToPlotNames={'E2F7','FOXM1','FOXO1','GATA2','GATA6','LEF1','NANOG','NFE2L2','NFYB','RELA','RELB','RORA','RUNX1','SMAD3','TCF7L2','TEAD2','TGIF1','USF1','USF2','ZFP57'}
% genesToPlotNames=a %%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! uncomment for plotting of genes selcted below

% [~,genesToPlotNames]=xlsread('Pathway_FattyAcid.xlsx')
% [~,genesToPlotNames]=xlsread('Pathway_TriAcylGlycerolSynthesis.xlsx')
% [~,genesToPlotNames]=xlsread('Pathway_Glycolysis.xlsx')

if TFsonly
    genesToPlotNames=intersect(TFs,genesToPlotNames)
end
genesToPlot=find(ismember(genes,genesToPlotNames))
% plotNames=genesWithLength(genesToPlot)
plotNames=genes(genesToPlot)
boxesPerPlot=25
Nrsplot=ceil(numel(genesToPlot)/boxesPerPlot)
for counter=1:Nrsplot
    %     figure
    figure('units','normalized','outerposition',[0 0 1 1])
    temp=((counter-1)*boxesPerPlot+1):min(numel(genesToPlot),((counter-1)*boxesPerPlot+boxesPerPlot));
    boxplot(log10(data(genesToPlot(temp),34:end)+1)',plotNames(temp))
    ylim([0 4])
    hold on
    spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
    for i = 1:numel(genesToPlot(temp)) %control
        allData=log10(data(genesToPlot(temp(i)),[1,12,23])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'k*','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %adipo
        allData=log10(data(genesToPlot(temp(i)),[2:3,13:14,24:25])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'ro','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %osteo
        allData=log10(data(genesToPlot(temp(i)),[7:8,18:19,29:30])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'go','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %adipo late
        allData=log10(data(genesToPlot(temp(i)),[4:6,15:17,26:28])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'ro','linewidth', 2)
    end
    for i = 1:numel(genesToPlot(temp)) %osteo late
        allData=log10(data(genesToPlot(temp(i)),[9:11,20:22,31:33])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'go','linewidth', 2)
    end
end

%% determine TPM cutoff
for counter=1:5
    figure
    hist(log10(data(:,counter)),30)
end
disp('TPM OFF cutoff: 10^0 = 1')
1
disp('TPM ON cutoff: 10^0 = 1')
1

%% discretizing query data per sample: median/upper quartile based
disp('discretizing query data per sample ...')
lowao=50 % 0 if in lower half of background data
highao=75 % 1 if in upper quartile of background data
TFsonly=1

lcao=prctile(data(:,34:end)',lowao)';
ucao=prctile(data(:,34:end)',highao)';

ddata=0.5*ones(size(data,1),numel(qsamples)); %nan
dscore=zeros(size(ddata));
for counter=1:numel(qsamples)
    disp('Discretizing sample ...')
    counter
    temp0=(data(:,counter)>=ucao);
    ddata(temp0,counter)=1;
    temp1=repmat(data(:,counter),1,size(data(:,34:end),2));
    temp2=sum((temp1-data(:,34:end))>=0,2);
    temp3=temp2*100/size(data(:,34:end),2);
    temp4=max(0,(temp3-highao)/(100-highao));
    dscore(temp0,counter)=temp4(temp0);
    
    temp0=(data(:,counter)<lcao);
    ddata(temp0,counter)=0;
    temp1=repmat(data(:,counter),1,size(data(:,34:end),2));
    temp2=sum((temp1-data(:,34:end))>=0,2);
    temp3=temp2*100/size(data(:,34:end),2);
    temp4=max(0,(lowao-temp3)/(lowao-0));
    dscore(temp0,counter)=temp4(temp0);
    
% % %     temp0=((data(:,counter)<1).*(ddata(:,counter)==1));
% % %     ddata(find(temp0),counter)=0.5; %TPM ON cutoff
% % %     dscore(find(temp0),counter)=1; %score=1
    temp0=(data(:,counter)<1);
    ddata(temp0,counter)=0; %TPM cutoff
    dscore(temp0,counter)=1; %score=1 (not expressed)
end
dgenes=genes;
if TFsonly
    [C,ia,ib]=intersect(genes,TFs);
    dataTFs=data(ia,:);
    ddataTFs=ddata(ia,:);
    dscoreTFs=dscore(ia,:);
    dgenesTFs=genes(ia);
end
disp('Done!')
disp('ON x3 / OFF x3 over C/A1-15/O1-15:')
if TFsonly
    temp=[sum(ddataTFs==1)',sum(ddataTFs==0)']';
else
    temp=[sum(ddata==1)',sum(ddata==0)']';
end
[temp(1,1:11);temp(1,12:22);temp(1,23:33);temp(2,1:11);temp(2,12:22);temp(2,23:33)]

%%
if TFsonly
    dataDif=dataTFs;
    ddataDif=ddataTFs;
    dscoreDif=dscoreTFs;
    dgenesDif=dgenesTFs;
else
    dataDif=data;
    ddataDif=ddata;
    dscoreDif=dscore;
    dgenesDif=dgenes;
end

% find genes always being 1 (quality control via plotting)
temp=nansum(ddataDif,2);
a=dgenesDif(find(temp==33))

%% check background distributions of selected genes (quartile cutoff ok?)
genesToPlot=find(ismember(dgenesDif,a(1:min(numel(a),30))))
for counter=1:numel(genesToPlot)
    figure
    %     figure('units','normalized','outerposition',[0 0 1 1])
    title('TPM distributions of selected genes over all selected samples (1. TPM 2. pseudo, log10)')
    subplot(2,1,1)
    temp=dataDif(genesToPlot(counter),34:end)';
    hist(temp,30)
    temp2=prctile(temp,lowao);
    temp3=axis;
    hold on, plot([temp2 temp2],[temp3(3) temp3(4)],'g','LineWidth',2)
    temp2=prctile(temp,highao);
    temp3=axis;
    hold on, plot([temp2 temp2],[temp3(3) temp3(4)],'r','LineWidth',2)
    xlabel(dgenesDif(genesToPlot(counter)))
    
    subplot(2,1,2)
    temp=log10(dataDif(genesToPlot(counter),34:end)+1)';
    hist(temp,30)
    temp2=prctile(temp,lowao);
    temp3=axis;
    hold on, plot([temp2 temp2],[temp3(3) temp3(4)],'g','LineWidth',2)
    temp2=prctile(temp,highao);
    temp3=axis;
    hold on, plot([temp2 temp2],[temp3(3) temp3(4)],'r','LineWidth',2)
    xlabel(dgenesDif(genesToPlot(counter)))
end

%% TPM distributions of selected genes over all selected samples (1. TPM 2. pseudo, log10)
temp=nansum(ddataDif,2);
a=dgenesDif(find(temp==33))
% plot selected genes a in plotting section above, by uncommenting ...=a

%% find diff expressed genes via kmeans, minimum 3 fc of centroid locations,
% at least 3  data points per kmeans cluster
% and fullfilling a ttest2 between the two clusters
diffGenes=dgenesDif;

addiffGenes=diffGenes;
for counter2=1:10 %re-run multiple times & remove border cases
    addiffGenesTemp={};
    for counter=1:numel(diffGenes)
        temp=dataDif(find(ismember(dgenesDif,diffGenes(counter))),[1:23,24,25:33]);
        temp=log10(temp+1);
        [idx,C]=kmeans(temp',2);
        temp2=[min(temp(idx==1)), max(temp(idx==1)), min(temp(idx==2)), max(temp(idx==2))];
        %         if abs(C(1)-C(2))>=log10(3) & (min(sum(idx==1),sum(idx==2))>=3) %still having 3-fold diff between centroid locations
        if ttest2(temp(idx==1),temp(idx==2)) & abs(C(1)-C(2))>=log10(3) & (min(sum(idx==1),sum(idx==2))>=3) %still having 3-fold diff between centroid locations
            diffGenes(counter)
            [median(temp2), sum(idx==1)*(C(1)<C(2))+sum(idx==2)*(C(2)<C(1)), sum(idx==1)*(C(1)>C(2))+sum(idx==2)*(C(2)>C(1))]
            addiffGenesTemp=[addiffGenesTemp; diffGenes(counter)];
        else
            %     diffGenes(counter)
            %         [median(temp2)]
            %         C
        end
    end
    addiffGenes=intersect(addiffGenes,addiffGenesTemp);
end
addiffGenes
numel(addiffGenes)
a=addiffGenes;
% consolidated list over multiple runs:
addiffGenes={'ADIPOQ' 'AHR' 'ALPL' 'ARNT2' 'ATF5' 'ATOH8' 'CEBPA' 'DBP' 'DDIT3' 'DMRTC2' 'EGR1' 'EGR2' 'EGR3' 'ESR1' 'ETV4' 'FABP4' 'FOS' 'GM20422' 'GM3055' 'GM44973' 'HEY1'  'HMGA2' 'HOXC13' 'IRF7' 'LPL' 'NR4A1' 'PLAGL1' 'SCX' 'SOX8' 'SP7' 'TLX1' 'TRP63' 'TSC22D3' 'VDR' 'ZFP658' 'ZHX3'}%     37
% 36

%% diff expressed genes: upper cluster=1; lower cluster=0
for counter=1:numel(addiffGenes)
    addiffGenes(counter)
    temp0=find(ismember(dgenesDif,addiffGenes(counter)))
    temp=dataDif(temp0,[1:23,24,25:33])
    temp=log10(temp+1)
    [idx,C]=kmeans(temp',2);
    if C(1)>C(2)
        C=[C(2);C(1)];
        idx=idx-1;
        idx(idx==0)=2;
    end
    temp2=[min(temp(idx==1)), max(temp(idx==1)), min(temp(idx==2)), max(temp(idx==2))];
    temp3=median(temp2)
    
    ddataDif(temp0,idx==1)=0; %lower cluster=0
    dscoreDif(temp0,idx==1)=(temp3-temp(idx==1))/(temp3-min(temp));
    ddataDif(temp0,idx==2)=1; %upper cluster=1
    dscoreDif(temp0,idx==2)=(temp(idx==2)-temp3)/(max(temp)-temp3);
    
    ddataDif(temp0,:)
    dscoreDif(temp0,:)
end
save queryDataDiscretized_MedianUpperQuartile_minTPM1 qsamples ddataDif dscoreDif dgenesDif lowao highao

%% some statistics
figure, hist(ddataDif), title('Frequency of discretized values per sample')
figure, plot(sort(sum(ddataDif'==1))), title('Genes called ON in all sample')
figure, plot(sort(sum(ddataDif'==0.5))), title('Genes not discretized in all samples')
figure, plot(sort(sum(ddataDif'==0),'descend')), title('Genes called OFF in all samples')
figure, hist(dscoreDif), title('Frequency of discretization scores per sample')

%%
toCompare=[2:11:33, 7:11:33; 3:11:33, 8:11:33; 4:11:33, 9:11:33; 5:11:33, 10:11:33; 6:11:33, 11:11:33]
genesToPlotNames={}; a=[];
for counter=1:size(toCompare,1) %time-points
    disp(' ')
    disp(['Time point:' num2str(counter)])
    temp1=ddataDif(:,toCompare(counter,1:3));
    temp2=ddataDif(:,toCompare(counter,4:6));
    temp3=(sum(temp1==1,2)>=2);
    disp(['Genes Adipo ON:' num2str(sum(temp3))])
    %         a=dgenesDif(find(temp3))
    genesToPlotNames=unique([genesToPlotNames; a]);
    temp3=(sum(temp1==1,2)>=2).*(sum((temp2>=0.5),2)<=1);
    disp(['Genes Adipo ON & Osteo OFF:' num2str(sum(temp3))])
            a=dgenesDif(find(temp3))
    genesToPlotNames=unique([genesToPlotNames; a]);
    temp3=(sum((temp2==1),2)>=2);
    disp(['Genes Osteo ON:' num2str(sum(temp3))])
%     a=dgenesDif(find(temp3))
    genesToPlotNames=unique([genesToPlotNames; a]);
    temp3=(sum(temp1>=0.5,2)<=1).*(sum((temp2==1),2)>=2);
    disp(['Genes Adipo OFF & Osteo ON:' num2str(sum(temp3))])
%             a=dgenesDif(find(temp3))
    genesToPlotNames=unique([genesToPlotNames; a]);
end
genesToPlotNames
numel(genesToPlotNames)

%% check tpm distributions of selected genes
disp('check tpm distributions of selected genes')
% genesToPlotNames=a
genesToPlot=find(ismember(genes,genesToPlotNames))
plotNames=genes(genesToPlot)
boxesPerPlot=25
Nrsplot=ceil(numel(genesToPlot)/boxesPerPlot)
for counter=1:Nrsplot
    %     figure
    figure('units','normalized','outerposition',[0 0 1 1])
    temp=((counter-1)*boxesPerPlot+1):min(numel(genesToPlot),((counter-1)*boxesPerPlot+boxesPerPlot));
    boxplot(log10(data(genesToPlot(temp),34:end)+1)',plotNames(temp))
    ylim([0 4])
    hold on
    spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
    for i = 1:numel(genesToPlot(temp)) %control
        allData=log10(data(genesToPlot(temp(i)),[1,12,23])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'k*','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %adipo
        allData=log10(data(genesToPlot(temp(i)),[2:3,13:14,24:25])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'ro','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %osteo
        allData=log10(data(genesToPlot(temp(i)),[7:8,18:19,29:30])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'go','linewidth', 1)
    end
    for i = 1:numel(genesToPlot(temp)) %adipo late
        allData=log10(data(genesToPlot(temp(i)),[4:6,15:17,26:28])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'ro','linewidth', 2)
    end
    for i = 1:numel(genesToPlot(temp)) %osteo late
        allData=log10(data(genesToPlot(temp(i)),[9:11,20:22,31:33])+1)';
        xCenter = 1:numel(genesToPlot);
        plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'go','linewidth', 2)
    end
end
