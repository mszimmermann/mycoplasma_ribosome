%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add interior point distance calculator for calculation of distances
% between ribosomes in each tomogram
addpath('InterPointDistanceMatrix')

% folder with data files
fileFolder = 'data_ribosomes';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first read polysome data (for distance of 7nm)
fileName = [fileFolder, filesep, 'polysome_tomoNum_pid_riboclass_sequence.txt'];

% read polysome table from file
polysomeTable = readtable(fileName);
% convert data table to matrix
polysomeData = table2array(polysomeTable);
% read type names from variable names of the table
typenames = polysomeTable.Properties.VariableNames;

% first two columns are polysome IDs. separate data and IDs in two matrices
polysomeDataIDs = polysomeData(:,1:2);
polysomeData = polysomeData(:,3:end);

% get the unique names of the ribosome types (clusters)
cluster_names = unique(polysomeData(:));
cluster_names(isnan(cluster_names))=[];
maxcluster = length(cluster_names);
% create a matrix of zeros to count pairs of ribosome states that occur
% after each other
ribosome_pairs_experimental = zeros(maxcluster, maxcluster);

% for each row in exampleData, record pairs of ribosome types (clusters)
for i = 1:size(polysomeData,1)
    % go only until first nan
    for j=1:(nnz(~isnan(polysomeData(i,:)))-1) 
        ribosome_pairs_experimental(cluster_names == polysomeData(i,j),...
                       cluster_names == polysomeData(i,j+1)) = ...
            ribosome_pairs_experimental(cluster_names == polysomeData(i,j),...
                       cluster_names == polysomeData(i,j+1)) + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all ribosome data across tomograms with coordinates 
% - to compare total distribution of ribosome classes to polysome-engaged ribosomes
% - to perform distance threshold analysis for polysome definitions 
fileName = [fileFolder, filesep, 'motl_annoted_addpolysome_allcombined.txt'];

polysome_total_data = table2array(readtable(fileName, 'delimiter',','));

polysome_total_distance = zeros(5*size(polysome_total_data,2), 7);
poly_dist_idx = 1;
tom_unique = unique(polysome_total_data(5,:)); 
% perform polysome analysis per tomogram
for tom_i = 1:length(tom_unique)
    % get ribosomes from the current tomogram
    curpoly = polysome_total_data(:, polysome_total_data(5,:)==tom_unique(tom_i));
    
    % get coordinates of RNA entry (rows 11-13)
    A = curpoly(11:13,:)';
    % get coordinates of RNA exit (rows 14-16)
    B = curpoly(14:16,:)';
    % calculate pairwise distances between RNA entry and exit sites
    d = ipdm(A,B);
    
    % loop through all ribosome pairs in the tomogram
    % and record tomogram ID, left and right ribosome ID, 
    % left and right ribosome states, and distance 
    curdist = zeros(size(curpoly,2)*(size(curpoly,2)-1)/2,1);
    curtomos = polysome_total_data(4, polysome_total_data(5,:)==tom_unique(tom_i));
    curstates = polysome_total_data(20, polysome_total_data(5,:)==tom_unique(tom_i));
    curstart = zeros(size(curpoly,2)*(size(curpoly,2)-1)/2,1);
    curend = zeros(size(curpoly,2)*(size(curpoly,2)-1)/2,1);
    curstate1 = zeros(size(curpoly,2)*(size(curpoly,2)-1)/2,1);
    curstate2 = zeros(size(curpoly,2)*(size(curpoly,2)-1)/2,1);
    idx=1;
    for i = 1:size(d,1)
        for j=i+1:size(d,2)
            % for each pair take the smallest distance between one RNA
            % entry and the other RNA exit sites to depermine order
            if d(i,j)<d(j,i)
                curdist(idx) = d(i,j);
                curstate1(idx) = curstates(j);
                curstate2(idx) = curstates(i);
                curstart(idx) = curtomos(j);
                curend(idx) = curtomos(i);
            else
                curdist(idx) = d(j,i);
                curstart(idx) = curtomos(i);
                curend(idx) = curtomos(j);
                curstate1(idx) = curstates(i);
                curstate2(idx) = curstates(j);
            end
            idx = idx+1;
        end
    end
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,1) = tom_unique(tom_i);
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,2) = curstart;
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,3) = curend;
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,4) = curdist;
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,6) = curstate1;
    polysome_total_distance(poly_dist_idx:poly_dist_idx+length(curdist)-1,7) = curstate2;
    poly_dist_idx = poly_dist_idx+length(curdist);
end
            
% convert pixels to nm
polysome_total_distance(:,5) = polysome_total_distance(:,4)*0.17;

% plot distribution of lengths
figure
subplot(1,2,1)
histogram(polysome_total_distance(:,5))
xlabel('Distance, nm')
ylabel('Number of ribosome pairs')
axis square

subplot(1,2,2)
histogram(polysome_total_distance(polysome_total_distance(:,5)<20,5))
xlabel('Distance, nm')
ylabel('Number of ribosome pairs')
suptitle('Distribution of distances between ribosomes')
axis square

print(gcf, '-painters', '-dpng', '-r600', ...
    'figure_plot_distribution_of_ribosome_distances')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate percentage of ribosome classes per cell (total, monosomes and
% engaged in polysomes)
% first check which ribosomes are not engaged in polysomes
% keep only ribosome-engaged in distance
% keep distance threshold at 7nm
distance_threshold = 7;
polysome_total_distance_filtered = polysome_total_distance(...
    polysome_total_distance(:,5)<=distance_threshold,:);


% get number of unique tomgrams (5 row in total polysome data)
all_tomograms = unique(polysome_total_data(5,:));
% calculate ribosome class percentage in total, poly- and mono-ribosomes
all_ribosomes_percentage = zeros(maxcluster,length(all_tomograms));
polysome_ribosomes_percentage = zeros(maxcluster,length(all_tomograms));
monosome_ribosomes_percentage = zeros(maxcluster,length(all_tomograms));

for i=1:maxcluster
    for j=1:length(all_tomograms)
        % last row (20) in all data is polysome class, count number per tomogram
        all_ribosomes_percentage(i,j) = nnz(polysome_total_data(20,polysome_total_data(5,:)==all_tomograms(j))==(i-1))/...
                                        nnz(polysome_total_data(5,:)==all_tomograms(j));
        % get IDs of polysome-engaged ribosomes in this tomogram
        polysome_indeces = unique(polysome_total_distance_filtered(...
            polysome_total_distance_filtered(:,1)==all_tomograms(j),2:3));
        % get all current ribosomes
        all_current_ribosomes = polysome_total_data(:,polysome_total_data(5,:)==all_tomograms(j));
        % get indices of ribosomes not engaged in poysomes
        [~,monoidx] = setdiff(all_current_ribosomes(4,:), polysome_indeces);
        all_current_monosomes = all_current_ribosomes(:,monoidx);
        % calculate percentage of current class
        monosome_ribosomes_percentage(i,j) = nnz(all_current_monosomes(20,:)==(i-1))/...
            length(all_current_monosomes(20,:));
        % check if there are polysomes in the current tomogram
        if nnz(polysomeDataIDs(:,1)==all_tomograms(j))
            polysome_ribosomes_percentage(i,j) = nnz(polysomeData(polysomeDataIDs(:,1)==all_tomograms(j),:)==(i-1))/...
                nnz(~isnan(polysomeData(polysomeDataIDs(:,1)==all_tomograms(j),:)));
        end
    end
end

% replace zeros with nans
all_ribosomes_percentage(all_ribosomes_percentage==0) = nan;
monosome_ribosomes_percentage(monosome_ribosomes_percentage==0) = nan;
polysome_ribosomes_percentage(polysome_ribosomes_percentage==0) = nan;

% calculate meand and medians ignoring tomograms with 0 counts
all_percentage_mean = nanmean(all_ribosomes_percentage,2);
all_percentage_std = nanstd(all_ribosomes_percentage,[],2);

mono_percentage_mean = nanmean(monosome_ribosomes_percentage,2);
mono_percentage_std = nanstd(monosome_ribosomes_percentage,[],2);

poly_percentage_mean = nanmean(polysome_ribosomes_percentage,2);
poly_percentage_std = nanstd(polysome_ribosomes_percentage,[],2);

%calculate WMW p-values of medians between poly and all
poly_percentage_p = zeros(size(polysome_ribosomes_percentage,1),1);
poly_percentage_FC = zeros(size(polysome_ribosomes_percentage,1),1);
for i=1:size(polysome_ribosomes_percentage,1)
    poly_percentage_p(i) = ranksum(all_ribosomes_percentage(i,:),...
        polysome_ribosomes_percentage(i,:));
    poly_percentage_FC(i) = log2(nanmedian(polysome_ribosomes_percentage(i,:))./...
        nanmedian(all_ribosomes_percentage(i,:)));
end
poly_percentage_fdr = mafdr(poly_percentage_p,'bhfdr',1);

%calculate WMW p-values of medians between poly and mono
poly_mono_percentage_p = zeros(size(polysome_ribosomes_percentage,1),1);
poly_mono_percentage_FC = zeros(size(polysome_ribosomes_percentage,1),1);
for i=1:size(polysome_ribosomes_percentage,1)
    poly_mono_percentage_p(i) = ranksum(monosome_ribosomes_percentage(i,:),...
        polysome_ribosomes_percentage(i,:));
    poly_mono_percentage_FC(i) = log2(nanmedian(polysome_ribosomes_percentage(i,:))./...
        nanmedian(monosome_ribosomes_percentage(i,:)));
end
poly_mono_percentage_fdr = mafdr(poly_mono_percentage_p,'bhfdr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotdata = [all_percentage_mean, mono_percentage_mean, poly_percentage_mean];
plotdatastd = [all_percentage_std, mono_percentage_std, poly_percentage_std];
plotdatafdr = poly_percentage_fdr;
plotdatamonofdr = poly_mono_percentage_fdr;
plotdatamedianfc = poly_percentage_FC;
plotdatamonomedianfc = poly_mono_percentage_FC;

% define order of plotting the states
plotorder = [9 5 8 6 4 2 0 3 1 7] + 1;

plotdata = plotdata(plotorder,:);
plotdatastd = plotdatastd(plotorder,:);
plotdatafdr = plotdatafdr(plotorder);
plotdatamonofdr = plotdatamonofdr(plotorder);
plotdatamedianfc = plotdatamedianfc(plotorder);
plotdatamonomedianfc = plotdatamonomedianfc(plotorder);


figure
barh(plotdata)
set(gca, 'YTick', 1:length(plotorder))
set(gca, 'YTickLabel', plotorder-1)
hold on 
errorbar(plotdata(:,1),[1:length(plotdata)]-0.24,...
    plotdatastd(:,1), '.k', 'horizontal');
errorbar(plotdata(:,2),[1:length(plotdata)],...
    plotdatastd(:,2), '.k', 'horizontal');
errorbar(plotdata(:,3),[1:length(plotdata)]+0.24,...
    plotdatastd(:,3), '.k', 'horizontal');
for i=1:size(plotdata,1)
    if abs(plotdatamedianfc(i))>=log2(1.5)
        text(max(plotdata(i,:))+max(plotdatastd(i,:))*1.3, i-0.24, ...
            sprintf('%.3f', plotdatafdr(i)))
    end
    if abs(plotdatamonomedianfc(i))>=log2(1.5)
        text(max(plotdata(i,:))+max(plotdatastd(i,:))*1.3, i+0.24, ...
            sprintf('%.3f', plotdatamonofdr(i)))
    end
end
xlim([0 0.6])
ylim([0.5 length(plotorder)+0.5])
axis ij
legend({'All', 'Monosome', 'Polysomes'}, 'location', 'SouthEast')
axis square
xlabel('Fraction')
title({'Differences in ribosome class fractions'...
       'Highlighted differences >50%'})
print(gcf, '-painters', '-dpng', '-r600', ...
    'figure_plot_bar_main_ribosome_percentages_per_tomogram_all_and_mono_vs_polysome_wmw')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%   'figure_plot_bar_main_ribosome_percentages_per_tomogram_all_and_mono_vs_polysome_wmw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate theoretical distributions of ribosome pairs based on 
% abundance of different classes engaged in polysomes
polysome_pairs_theoretical = zeros(maxcluster);
for i=1:maxcluster
    for j=1:maxcluster
        polysome_pairs_theoretical(i,j) = (nnz(polysome_total_data(20,:)==(i-1))/nnz(polysome_total_data(20,:))) *...
             (nnz(polysome_total_data(20,:)==(j-1))/nnz(polysome_total_data(20,:)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform random shuffling of polysomes to calculate distribution of pair
% frequencies
% the three column is the type of triplet member 1, 2, 3)
ribosome_pair_ids = zeros(maxcluster^2,2);
idx=1;
for i=1:maxcluster
    for j=1:maxcluster
        ribosome_pair_ids(idx,:) = [i-1 j-1];
        idx=idx+1;
    end
end
            
% define number of shuffles. The more shuffles, the longer the
% calculations. Tip: change to 10000 for better resolution of p-values. 
shufflen=1000;
shuffled_matrixes = zeros(length(ribosome_pair_ids), shufflen);
% calculate number of polysomes of certain length
shuffled_polysome_length = zeros(shufflen, size(polysomeData,2));

for shi = 1:shufflen
    %tic
    exampleData_shuffled = polysomeData(:);
    exampleData_shuffled = exampleData_shuffled(randperm(numel(polysomeData)));
    exampleData_shuffled = reshape(exampleData_shuffled, size(polysomeData));
    for i=1:size(exampleData_shuffled,1)
        curvals = exampleData_shuffled(i,:);
        exampleData_shuffled(i,:) = [curvals(~isnan(curvals)),...
            nan(1,nnz(isnan(curvals)))];
    end
    exampleData_shuffled(sum(~isnan(exampleData_shuffled),2)<2,:)=[];
    
    % create a matrix of zeros to count pairs of ribosome states that occur
    % after each other
    % get the unique names of the clusters
    ribosome_pairs_shuffled = zeros(maxcluster, maxcluster);

    % for each row in exampleData, record pairs of ribosome clusters
    for i = 1:size(exampleData_shuffled,1)
        % go only until first nan
        for j = 1:(nnz(~isnan(exampleData_shuffled(i,:)))-1) 
            curidx = (ribosome_pair_ids(:,1) == exampleData_shuffled(i,j)) &...
                (ribosome_pair_ids(:,2) == exampleData_shuffled(i,j+1));
            if nnz(curidx)
                shuffled_matrixes(curidx, shi) = ...
                    shuffled_matrixes(curidx, shi) + 1;
                ribosome_pairs_shuffled(exampleData_shuffled(i,j)+1, exampleData_shuffled(i,j+1)+1) = ...
                    ribosome_pairs_shuffled(exampleData_shuffled(i,j)+1, exampleData_shuffled(i,j+1)+1)+1;
            else
                fprintf('Not found %d %d\n', i, j);
            end
        end
        curlen = nnz(~isnan(exampleData_shuffled(i,:)));
        shuffled_polysome_length(shi, curlen) = shuffled_polysome_length(shi, curlen)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate p-values and fdr and plot distributions

% save p-values and fold changes between experimental and shuffled
% distributions
pmat_shuffled = zeros(length(plotorder));
fcmat_shuffled = zeros(length(plotorder));
% calcuate p-values and fold changes for pairs regardless of pair order
pmat_shuffled_noorder = zeros(length(plotorder));
fcmat_shuffled_noorder = zeros(length(plotorder));
% save frequency of shuffled ribosomes
ribosome_pairs_shuffled_freq = zeros(maxcluster);
for i=1:length(plotorder)
    rowi = plotorder(i);
    for j=1:length(plotorder)
        coli = plotorder(j);
    
        curidx = (ribosome_pair_ids(:,1) == (rowi-1)) &...
                (ribosome_pair_ids(:,2) == (coli-1));
            
        % get opposite pair for no-order analysis
        curidx_opposite = (ribosome_pair_ids(:,1) == (coli-1)) &...
                (ribosome_pair_ids(:,2) == (rowi-1));
        
        % get distribution of frequencies for the current ribosome pair
        shuffled_dist = zeros(shufflen,1);
        shuffled_dist_noorder = zeros(shufflen,1);
        for shi=1:shufflen
            
            shuffled_dist(shi) = shuffled_matrixes(curidx, shi)./...
                sum(shuffled_matrixes(:,shi));
            
            if coli~=rowi
                shuffled_dist_noorder(shi) = (shuffled_matrixes(curidx, shi)+...
                    shuffled_matrixes(curidx_opposite, shi))./...
                    sum(shuffled_matrixes(:,shi));
            else
                shuffled_dist_noorder(shi) = shuffled_matrixes(curidx, shi)./...
                    sum(shuffled_matrixes(:,shi));
            end
        end

        % get experiemntal value for pair frequency
        pairs_exp_all = ribosome_pairs_experimental(rowi,coli)./...
            sum(ribosome_pairs_experimental(:));
        
        if coli~=rowi
            pairs_exp_all_noorder = (ribosome_pairs_experimental(rowi,coli) + ...
                ribosome_pairs_experimental(coli,rowi))./...
                sum(ribosome_pairs_experimental(:));
        else
            pairs_exp_all_noorder = ribosome_pairs_experimental(rowi,coli)./...
                sum(ribosome_pairs_experimental(:));
        end
        
        % calculate mean of shuffled frequencies
        ribosome_pairs_shuffled_freq(rowi, coli) = ...
            mean(shuffled_dist);
        
        % to calculate permutation p-value, add the experimental value to
        % the vector of shuffled freqencies
        shuffled_dist(shufflen+1) = pairs_exp_all;
        % calculate number of equal or more extreme frequency values
        curP = nnz(shuffled_dist<=pairs_exp_all)/(shufflen+1);
        curP_shuffled = min(curP, 1-curP);
        if curP_shuffled == 0
            curP_shuffled = 1/(shufflen+1);
        end
        % save p-value and log2 fc to matrix
        pmat_shuffled(rowi, coli) = curP_shuffled;
        fcmat_shuffled(rowi, coli) = log2(pairs_exp_all/mean(shuffled_dist));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % no-order pair p-values
        % to calculate permutation p-value, add the experimental value to
        % the vector of shuffled freqencies
        shuffled_dist_noorder(shufflen+1) = pairs_exp_all_noorder;
        % calculate number of equal or more extreme frequency values
        curP = nnz(shuffled_dist_noorder<=pairs_exp_all_noorder)/(shufflen+1);
        curP_shuffled = min(curP, 1-curP);
        % put minimal value instead of 0
        if curP_shuffled == 0
            curP_shuffled = 1/(shufflen+1);
        end
        pmat_shuffled_noorder(rowi, coli) = curP_shuffled;
        fcmat_shuffled_noorder(rowi, coli) = log2(pairs_exp_all_noorder/mean(shuffled_dist_noorder));
       
        
    end
end
% adjust for multiple hypotheses testing
pmat_shuffled_fdr = pmat_shuffled(:);
pmat_shuffled_fdr = mafdr(pmat_shuffled_fdr, 'bhfdr', 'true');
pmat_shuffled_fdr = reshape(pmat_shuffled_fdr, size(pmat_shuffled));
% adjust for multiple hypotheses testing pairs without order
pmat_shuffled_noorder = tril(pmat_shuffled_noorder);
pmat_shuffled_noorder_fdr = pmat_shuffled_noorder;
pmat_shuffled_noorder_fdr(pmat_shuffled_noorder_fdr==0)=nan;
pmat_shuffled_noorder_fdr = pmat_shuffled_noorder_fdr(:);

pmat_shuffled_noorder_fdr = mafdr(pmat_shuffled_noorder_fdr, 'bhfdr', 'true');
pmat_shuffled_noorder_fdr = reshape(pmat_shuffled_noorder_fdr, size(pmat_shuffled));


% plot shuffled results to figure
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spidx=1; % index of subplots
for i=1:length(plotorder)
    rowi = plotorder(i);
    for j=1:length(plotorder)
        coli = plotorder(j);
    
        curidx = (ribosome_pair_ids(:,1) == (rowi-1)) &...
                (ribosome_pair_ids(:,2) == (coli-1));
        
        % get distribution of frequencies for the current ribosome pair
        shuffled_dist = zeros(shufflen,1);
        for shi=1:shufflen
            
            shuffled_dist(shi) = shuffled_matrixes(curidx, shi)./...
                sum(shuffled_matrixes(:,shi));
        end

        % multiply by 1000 for plotting
        shuffled_dist = shuffled_dist*1000;
        
        subplot(length(plotorder),length(plotorder),spidx)
        h = histogram(shuffled_dist,10);
        hold on

        % get experiemntal value for pair frequency
        pairs_exp_all = ribosome_pairs_experimental(rowi,coli)./...
            sum(ribosome_pairs_experimental(:));
        
        % multiply by 1000 for plotting
        pairs_exp_all= pairs_exp_all*1000;
        
        plot(pairs_exp_all*[1 1],...
            [1 max(h.Values)], 'r', 'LineWidth', 2)
        
        % get theoretical value for pair frequency
               % multiply by 1000 for plotting
        pairs_theor_all= polysome_pairs_theoretical(rowi,coli)*1000;
        
        plot(pairs_theor_all*[1 1],...
            [1 max(h.Values)], 'g', 'LineWidth', 2)
              
        title(sprintf('%d %d', cluster_names(rowi), cluster_names(coli)))
    
        curP_shuffled = pmat_shuffled_fdr(rowi, coli);
       
        % highlight significant differences in yellow
        if curP_shuffled<=0.01
            set(gca,'Color','yellow')
        end
        
        % set the same y limits
        ylim([0 0.4*shufflen])
        % keep y ticks only for forst column
        if j>1
            set(gca, 'YTick', [])
        end
        % add x and y labels for one plot
        if (i==round(length(plotorder)/2)) && (j==1)
            ylabel('Number of simulations')
        end
        if (i==length(plotorder)) && (j==round(length(plotorder)/2))
            xlabel('Fraction of pairs, x 10e-3')
        end
        % increase subplot counter
        spidx=spidx+1;
          
    end
end
suptitle('Experimental (red), theoretical (green) and simulated (blue) frequencies, highlighted in yellow pFDR<=0.01')
orient landscape
set(gcf, 'InvertHardCopy', 'off');
print(gcf, '-painters', '-dpng', '-r600', ...
    'figure_plot_exp_vs_suffled_distributions')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%   'figure_plot_exp_vs_suffled_distributions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot heatmap of experimental vs shuffled distributions
ribosome_pairs_exp_freq = ribosome_pairs_experimental./sum(ribosome_pairs_experimental(:));

plotdata = (ribosome_pairs_exp_freq./...
    ribosome_pairs_shuffled_freq);
plotdatafdr = pmat_shuffled_fdr;

plotdata=plotdata( plotorder, plotorder );
plotdatafdr=plotdatafdr( plotorder, plotorder );

figure
imagesc(plotdata)
set(gca, 'XTick', 1:length(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:length(plotorder))
set(gca, 'YTickLabel', plotorder-1) 
h = colorbar;
set(get(h,'label'),'string','Experimental/Shuffled');
caxis([0 2])

title('Experimental / shuffled frequencies')
axis square

xlabel('Following ribosome')
ylabel('Leading ribosome')
% reverse y axis
%set(gca,'YDir','normal')

% plot significant changes
hold on
jindent=0.1;
for i =1:length(plotdata)
    for j=1:length(plotdata)
        if plotdatafdr(i,j)<=0.01
            text(j-jindent,i+jindent,'*');
        end
    end
end

orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     'figure_imagesc_experimental_fs_shuffled')
print(gcf, '-painters', '-dpng', '-r600', ...
    'figure_imagesc_experimental_fs_shuffled')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot experimental, theoretical, shuffled, and ratio thereof
% plot pair count matrix
fig = figure('units','normalized','outerposition',[0 0 1 1]);

plot_log_flag = 1;

% experimental
plotdata = ribosome_pairs_experimental;
% calculate fraction of pairs
plotdata = plotdata/sum(plotdata(:));
plotdata = plotdata(plotorder, plotorder);

if plot_log_flag
    plotdata = log10(plotdata);
end

subplot(2,3,1)
imagesc(plotdata)
set(gca, 'XTick', 1:max(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:max(plotorder))
set(gca, 'YTickLabel', plotorder-1) 

h = colorbar;
set(get(h,'label'),'string','Fraction');
if plot_log_flag
    set(get(h,'label'),'string','Fraction, log10');
    caxis([-4 -0.75])
end
xlabel('Following ribosome')
ylabel('Leading ribosome')
title('Experimental pairs')
axis square

% theoretical
plotdata = polysome_pairs_theoretical;
plotdata = plotdata(plotorder, plotorder);

if plot_log_flag
    plotdata = log10(plotdata);
end

% plot pair count matrix
subplot(2,3,2)
imagesc(plotdata)
set(gca, 'XTick', 1:max(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:max(plotorder))
set(gca, 'YTickLabel', plotorder-1) 

h = colorbar;
set(get(h,'label'),'string','Fraction');
if plot_log_flag
    set(get(h,'label'),'string','Fraction, log10');
    caxis([-4 -0.75])
end

xlabel('Following ribosome')
ylabel('Leading ribosome')
title('Theoretical pairs')
axis square

% shuffled
plotdata = ribosome_pairs_shuffled_freq;
plotdata = plotdata(plotorder, plotorder);

if plot_log_flag
    plotdata = log10(plotdata);
end

% plot pair count matrix
subplot(2,3,3)
imagesc(plotdata)
set(gca, 'XTick', 1:max(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:max(plotorder))
set(gca, 'YTickLabel', plotorder-1) 

h = colorbar;
set(get(h,'label'),'string','Fraction');
if plot_log_flag
    set(get(h,'label'),'string','Fraction, log10');
    caxis([-4 -0.75])
end

xlabel('Following ribosome')
ylabel('Leading ribosome')
title('Shuffled pairs')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ratio Exp / Theoretical
plotdata = (ribosome_pairs_exp_freq./...
    polysome_pairs_theoretical);

plotdata=plotdata( plotorder, plotorder );

subplot(2,3,5)
imagesc(plotdata)
set(gca, 'XTick', 1:length(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:length(plotorder))
set(gca, 'YTickLabel', plotorder-1) 
h = colorbar;
set(get(h,'label'),'string','Experimental/Theoretical');
caxis([0.25 4])


title('Experimental / theoretical frequencies')
axis square

xlabel('Following ribosome')
ylabel('Leading ribosome')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental / shuffled
plotdata = (ribosome_pairs_exp_freq./...
    ribosome_pairs_shuffled_freq);
plotdatafdr = pmat_shuffled_fdr;

plotdata=plotdata( plotorder, plotorder );
plotdatafdr=plotdatafdr( plotorder, plotorder );

subplot(2,3,6)
imagesc(plotdata)
set(gca, 'XTick', 1:length(plotorder))
set(gca, 'XTickLabel', plotorder-1) 
set(gca, 'YTick', 1:length(plotorder))
set(gca, 'YTickLabel', plotorder-1) 
h = colorbar;
set(get(h,'label'),'string','Experimental/Shuffled');
caxis([0.5 2])

title('Experimental / shuffled frequencies')
axis square

xlabel('Following ribosome')
ylabel('Leading ribosome')
% reverse y axis
%set(gca,'YDir','normal')

% plot significant changes
hold on
jindent=0.1;
for i =1:length(plotdata)
    for j=1:length(plotdata)
        if plotdatafdr(i,j)<=0.01
            text(j-jindent,i+jindent,'*');
        end
    end
end

suptitle('Experimental, theoretical and shuffled occurance of ribosome type pairs')

orient landscape

if plot_log_flag
     % print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    %     'figure_imagesc_fraction_of_pairs_experimental_theoreticalALL_shuffled_caxis0_5_2')
    print(gcf, '-painters', '-dpng', '-r600', ...
        'figure_imagesc_fraction_of_pairs_experimental_theoreticalALL_shuffled_caxis0_5_2_log10')
    
else
    % print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    %     'figure_imagesc_fraction_of_pairs_experimental_theoreticalALL_shuffled_caxis0_5_2')
    print(gcf, '-painters', '-dpng', '-r600', ...
        'figure_imagesc_fraction_of_pairs_experimental_theoreticalALL_shuffled_caxis0_5_2')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform analysis of polysomes at different thresholds of ribosome
% definition
    
% define range of threshold distances
dist_range = [3 4 5 6 7 8 9 10 100];
% get polysomes at different thresholds
poly_idx = 1;
curpolysome = nan * ones(60000, 300);

for i=1:length(dist_range)
    curthreshold = dist_range(i);
    
    curpoly = polysome_total_distance(polysome_total_distance(:,5)<=curthreshold,:);
    curpoly_names = [strcat(arrayfun(@(x) num2str(x), curpoly(:,1), 'unif', 0), '_', ...
                                    arrayfun(@(x) num2str(x), curpoly(:,2), 'unif', 0)),...
                     strcat(arrayfun(@(x) num2str(x), curpoly(:,1), 'unif', 0), '_', ...
                                    arrayfun(@(x) num2str(x), curpoly(:,3), 'unif', 0))];
    curpoly_names_taken = zeros(length(curpoly_names),1);
    
    poly_start = setdiff(curpoly_names(:,1), curpoly_names(:,2));
    [~, poly_start_idx] = intersect(curpoly_names(:,1), poly_start, 'stable'); 

    for j=1:length(poly_start_idx)
        poly_next = curpoly_names(poly_start_idx(j),2);
        curpolysome(poly_idx,1) = curthreshold;
        curpolysome(poly_idx,2) = curpoly(poly_start_idx(j),1);
        curpolysome(poly_idx,3) = j;
        curpolysome(poly_idx,4) = curpoly(poly_start_idx(j),6);
        curpolysome(poly_idx,5) = curpoly(poly_start_idx(j),7);
        curpoly_names_taken(poly_start_idx(j)) = 1;
        cur_j = 6;
        while ismember(poly_next, curpoly_names(:,1))
            idx = (ismember(curpoly_names(:,1), poly_next) & curpoly_names_taken==0);
            if nnz(idx)==0
                break
            end
            if nnz(idx)>1
                idx = find(idx);
                idx = idx(1);
            end
            poly_next = curpoly_names(idx,2);
            curpolysome(poly_idx,cur_j) = curpoly(idx,7);
            curpoly_names_taken(idx) = 1;
            cur_j = cur_j+1;
        end
        poly_idx = poly_idx+1;
    end
end
curpolysome(poly_idx:end,:) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate representation of pairs in different thresholds

% create a matrix of zeros to count pairs of ribosome states that occur
% after each other

ribosome_pairs_threshold = zeros(size(ribosome_pair_ids,1),length(dist_range)-1);

for dist_i = 1:length(dist_range)-1
    % first three columns are threshold, tomogram and polysome IDs
    curpolysomeData = curpolysome(curpolysome(:,1)==dist_range(dist_i),:);
    curpolysomeData = curpolysomeData(:,4:end);
    % for each row in exampleData, record pairs of ribosome clusters and add
    % one to the corresponding vector of pair counts ribosome_triples
    for i = 1:size(curpolysomeData,1)
        % go only until first nan
        for j=1:(nnz(~isnan(curpolysomeData(i,:)))-1)
            curidx = (ribosome_pair_ids(:,1) == curpolysomeData(i,j)) &...
                (ribosome_pair_ids(:,2) == curpolysomeData(i,j+1));
            if nnz(curidx)
                if dist_range(dist_i)<100
                    ribosome_pairs_threshold(curidx, dist_i) = ribosome_pairs_threshold(curidx, dist_i) + 1;
                else
                    % for all polysomes, calculate only odd pairs
                    %if mod(j,2)==1
                        ribosome_pairs_threshold(curidx, dist_i) = ribosome_pairs_threshold(curidx, dist_i) + 1;
                   % end
                end     
            else
                fprintf('Not found %d, %d\n', curpolysomeData(i,j), curpolysomeData(i,j+1))
            end
        end
    end
end
clear curpolysomeData

ribosome_pairs_threshold_fraction = ribosome_pairs_threshold;
for i=1:size(ribosome_pairs_threshold,2)
    ribosome_pairs_threshold_fraction(:,i) = ribosome_pairs_threshold_fraction(:,i)/sum(ribosome_pairs_threshold(:,i));
end

% calculate theoretical numbers from the full ribosome dataset
ribosomes_theor = polysome_total_data(20,:);
ribosome_pair_theor = zeros(length(ribosome_pair_ids),1);
for i=1:length(ribosome_pair_ids)
    curR1 = ribosome_pair_ids(i,1);
    curR2 = ribosome_pair_ids(i,2);
    curfract1 = nnz(ribosomes_theor==curR1)/length(ribosomes_theor);
    curfract2 = nnz(ribosomes_theor==curR2)/length(ribosomes_theor);
    % frequency of pair equals product of ribosome frequencies
    ribosome_pair_theor(i) = curfract1*curfract2;
end

% calculate difference between experimental and theoretical frequency at
% differen thresholds
ribosome_pairs_threshold_fraction_difference = ribosome_pairs_threshold_fraction;
for i=1:size(ribosome_pairs_threshold_fraction,1)
    ribosome_pairs_threshold_fraction_difference(i,:) = ribosome_pairs_threshold_fraction(i,:)-...
        ribosome_pair_theor(i);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate number of engaged ribosomes per class
ribosome_class_engaged = zeros(maxcluster, length(dist_range));
for dist_i = 1:length(dist_range)
    % first three columns are threshold, tomogram and polysome IDs
    polysomeData = curpolysome(curpolysome(:,1)==dist_range(dist_i),:);
    polysomeData = polysomeData(:,4:end);
    for j=1:maxcluster
        ribosome_class_engaged(j, dist_i) = nnz(polysomeData==j-1);
    end
end

% get right (following) and left (leading) neighbors engagement
ribosome_class_engaged_left = zeros(maxcluster, length(dist_range));
ribosome_class_engaged_right = zeros(maxcluster, length(dist_range));

for dist_i=1:length(dist_range)
    curthreshold = dist_range(dist_i);
    
    curpoly = polysome_total_distance(polysome_total_distance(:,5)<=curthreshold,:);
    for j=1:maxcluster
        ribosome_class_engaged_left(j, dist_i) = nnz(curpoly(:,6)==j-1);
        ribosome_class_engaged_right(j, dist_i) = nnz(curpoly(:,7)==j-1);
    end
end

% for each ribosome class, calculate total fraction of ribosomes engaged in
% polysomes, and leading and following neighbor separately
ribosome_class_engaged_fraction = ribosome_class_engaged;
ribosome_class_engaged_fraction_left = ribosome_class_engaged_left;
ribosome_class_engaged_fraction_right = ribosome_class_engaged_right;

for j=1:size(ribosome_class_engaged,2)
    ribosome_class_engaged_fraction(:,j) = ribosome_class_engaged_fraction(:,j)/sum(ribosome_class_engaged(:,j));
    ribosome_class_engaged_fraction_left(:,j) = ribosome_class_engaged_fraction_left(:,j)/sum(ribosome_class_engaged_left(:,j));
    ribosome_class_engaged_fraction_right(:,j) = ribosome_class_engaged_fraction_right(:,j)/sum(ribosome_class_engaged_right(:,j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figure
% plot experimental minus theoretical frequency
% plot leading divided by following ribosome class fractions

fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1)
hold on
for i=1:size(ribosome_pairs_threshold_fraction_difference, 1)
    % plot only main classes (0..9)
    if (ribosome_pair_ids(i,1)<10) && ( (ribosome_pair_ids(i,2)<10))
        plot(ribosome_pairs_threshold_fraction_difference(i,:), 'Color', [.5 .5 .5]);
        % label pairs if they are different from theoretical
        if max(abs(ribosome_pairs_threshold_fraction_difference(i,1:3)))>0.015
            text(1,ribosome_pairs_threshold_fraction_difference(i,1),...
                 sprintf('%d %d',ribosome_pair_ids(i,:)))
        end
    end
end

axis square
set(gca, 'XTick', 1:length(dist_range)-1)
set(gca, 'XTickLabels', dist_range(1:end-1))
xlim([1 length(dist_range)-1])
ylim([-0.08 0.08])
xlabel('Distance threshold, nm')
ylabel('Experimental - theoretical fraction')
title('Polysome pair fraction compared to theoretical')
% plot leading divided by following ribosome class fractions
% Elongation factor bound colors (bluish green)
% elongation factor binding is required to exit state color (orange)
state_colors_EF = [[.8 .4 0];... %orange
    [0 .6 .5];...  % bluish green,...
    [.8 .4 0];...
    [0 .6 .5];...
    [.5 .5 .5];...
    [.8 .4 0];...
    [0 .6 .5];...
    [0 .6 .5];...
    [0 .6 .5];...
    [.8 .4 0]];

% plot left to right ratio vs distance
subplot(1,2,2)
hold on
axis square
xlabel('Distance threshold, nm')
title('Ribosome engagement in polysome pairs')
ylabel('Leading to following ratio')
% plot only main classes (0..9)
for i=1:10
       
    curplot = ribosome_class_engaged_fraction_left(i,:)./ribosome_class_engaged_fraction_right(i,:);
   
    plot( [1:length(curplot)], curplot,...
        'color',state_colors_EF(i,:), 'LineWidth', 2)
    
    ylim([0 2])
    text(1, curplot(1),sprintf('%d', i-1))   
end
plot([1,length(curplot)], [1 1], 'k--')
set(gca, 'XTick', 1:length(dist_range)-1)
set(gca, 'XTickLabels', dist_range(1:end-1))
xlim([1 length(dist_range)-1])
legend({'EF binding required to change state',...
    'EF bound to ribosome'})
orient landscape
print(gcf, '-painters', '-dpng', '-r600', ...
    'figure_plot_ribosome_pairs_distance_and_LR_colored_by_EF')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     'figure_plot_ribosome_pairs_distance_and_LR_colored_by_EF')
