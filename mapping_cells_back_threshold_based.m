%% Preparation of parameter matrices using a threshold (threshold based filtering happens after z-scoring)

% Define the folder names and file names
root_path = ''; % insert file path to a root folder containing subfolders of all tips
folder_names = {'Sample1', 'Sample4', 'Sample5', 'Sample9', 'Sample14', 'Sample19', 'Sample20', 'Sample21', 'Sample22','Sample23','Sample24','Sample25','Sample26','Sample27','Sample28','Sample29'};
stats_matrix_file_name = 'stats_matrix_all.mat';

% Set the threshold based on z-score
threshold = 1;

% Loop through the folders
for i = 1:length(folder_names)
    % Load the stats_matrix_all from the current folder
    stats_matrix_file_path = fullfile(folder_names{i}, stats_matrix_file_name);
    load(stats_matrix_file_path, 'stats_matrix_all');
    
    % Z-score all columns of stats_matrix_all
    stats_matrix_all_zscored = zscore(stats_matrix_all);

    % Find indices of rows in the z-scored matrix that exceed the threshold in the fourth column
    indices = find(stats_matrix_all_zscored(:,3) >= threshold);

    % Extract the rows of the z-scored matrix corresponding to these indices
    stats_matrix_ellipthresh_zscored = stats_matrix_all_zscored(indices,:);

    % Save the thresholded matrix in z-scored form to the current folder
    save_path = fullfile(folder_names{i}, 'stats_matrix_ellipthresh_zscored.mat');
    save(save_path, 'stats_matrix_ellipthresh_zscored');
end

%% Visualize the joint heatmap of filtered parameter matrices

%HEATMAP ANALYSIS OF THRESHOLDED MATRICES + MAPPING BACK
root_path = ''; % insert file path to a root folder containing subfolders of all tips
folders   = {'Sample1', 'Sample4', 'Sample5', 'Sample9', 'Sample14', 'Sample19', 'Sample20', 'Sample21', 'Sample22','Sample23','Sample24','Sample25','Sample26','Sample27','Sample28','Sample29'}; % folders here (you need save combined stats_matrix_all for each group, ex. ampullae, bifurcated, etc. into separate folders)
counter = 0;
for xp_num = 1:length(folders)
    xp(xp_num).stat_all = load(fullfile(root_path,folders{xp_num},'stats'));

    xp(xp_num).Centroid = xp(xp_num).stat_all.stats.Centroid; %centroid from the stats matrix
       
    xp(xp_num).M  = load(fullfile(root_path,folders{xp_num},'stats_matrix_ellipthresh_zscored'));
    xp(xp_num).M  = xp(xp_num).M.stats_matrix_ellipthresh_zscored;
    xp(xp_num).Mt = xp(xp_num).M.';
    xp(xp_num).M_length = size(xp(xp_num).M,1);
    xp(xp_num).CellIdentities    = find(xp(xp_num).stat_all.stats.Volume);
    xp(xp_num).Big_ID = counter+(1:xp(xp_num).M_length)';
    counter = counter+xp(xp_num).M_length;
end
%% 

%Here we concatenate our thresholded matrices, modify them to fit the
%clustergram and make a clustergram
M = cat(1,xp.M);
Mt = M';
% heatmap for nc-genes, fixed columns
params     = {'Volume','Volume/Surface ratio','Ellipticity','Elongation','Longest Axis'};

heatm   = clustergram(Mt,'RowLabels',params',...
    'RowPDist','cosine','ColumnPDist','cosine','linkage','average',...
    'DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
%% Mapping cells back from the heatmap

% in rest of the script ONLY change "xp_ref_num"
for xp_ref_num = 1:length(folders);
 % type here the embryo number 1)E11-1  2)E11  3)3c2  4)4a1  5)4b2
counter    = 6;  % this pointer is for the color
indices = 1042; % group number from heatmap
tiffname = 'green'; %write the name of pseudocolored cell group that we will save as tiff
%=========================================================================
cmp = jet(length(indices));
cmp(1,:)  = [0.1,0.9,1];   % light blue    
cmp(2,:)  = [0.8,0.6,0.9]; % lilac
cmp(3,:)  = [0.6,0.9,0.4]; % green   
cmp(4,:)  = [0.8,0.2,0.4]; % red pink   
cmp(5,:)  = [0.9,1,0];     % yellow
cmp(6,:)  = [0.1,0.4,1];   % dark blue 
cmp(7,:)  = [0,0.8,0.8];   % cyan 
cmp(8,:)  = [0.9,0,1];     % blueish purple 
cmp(9,:)  = [1,0.6,0];     % burnt-ish orange   
cmp(10,:) = [0.9,0.4,0.5]; % baby pink      
%figure
%imshow(image1(:,:,40),[])
sub_cluster_BIG = struct('cells_of_interest',[],'cells_of_interest_small',[]);
CellIdentities = cat(1,xp.Big_ID); 
bds = [0;cumsum(cat(1,xp.M_length))];

for num = indices 
    
    group_of_interest = clusterGroup(heatm, num, 'col'); 
    Col_Labels        = group_of_interest.ColumnLabels; % here we have extracted some column labels from the clustergram; 
    
    Double_Labels     = cell(0); % pre-allocating
    % improvement: Determine way how to get all clusters out of a clustergram. 
    for i = 1 : length(Col_Labels)
        Double_Labels{i} = str2double(Col_Labels{i});
    end
    cells_of_Interest_Big = CellIdentities(cell2mat(Double_Labels));
    sub_cluster_BIG(counter).cells_of_interest = cells_of_Interest_Big;

    for xp_num = 1 : length(xp) % Dont change this value!!
        % identify the cells that come from expertiment xp_num
        in = intersect(find(cells_of_Interest_Big>bds(xp_num)),find(cells_of_Interest_Big<bds(xp_num+1)));
        % find the corresponding identities in that experiment 
        template = xp(xp_num).Big_ID; 
        pattern = cells_of_Interest_Big(in);
        distance = pdist2(template,pattern);
        in_xp_num = find(min(distance,[],2)==0);
        small_id = xp(xp_num).CellIdentities(in_xp_num);
        sub_cluster_BIG(counter).cells_of_interest_small{xp_num} = small_id; 
    end

    counter = counter+1;
end
% here the end result is image of membrain staining with centroids for
% selected subgroups

sub_cluster_BIG_pruned = sub_cluster_BIG; % pre-allocating
for i = 1 : length(sub_cluster_BIG)
    template = sub_cluster_BIG(i).cells_of_interest;
    for j = 1 : length(sub_cluster_BIG)
        pattern = sub_cluster_BIG(j).cells_of_interest;
        if isempty(setdiff(pattern,template)) &&(i~=j)
            sub_cluster_BIG_pruned(i).cells_of_interest = setdiff(template,pattern); 
        end
    end

    cells_of_Interest_Big = sub_cluster_BIG_pruned(i).cells_of_interest; 

    for xp_num = 1 : length(xp)
        % identify the cells that come from expertiment xp_num
        in = intersect(find(cells_of_Interest_Big>bds(xp_num)),find(cells_of_Interest_Big<bds(xp_num+1)));
        % find the corresponding identities in that experiment 
        template = xp(xp_num).Big_ID; 
        pattern = cells_of_Interest_Big(in);
        distance = pdist2(template,pattern);
        in_xp_num = find(min(distance,[],2)==0);
        small_id = xp(xp_num).CellIdentities(in_xp_num);
        sub_cluster_BIG_pruned(i).cells_of_interest_small{xp_num} = small_id; 
    end
end
% Make culster 3D stack visualization with cell outlines
% first define new sub group of cells then start with loading the 
% label watershed matrix for that embryo you want to image

Label = load(fullfile(root_path,folders{xp_ref_num},'Final_Label.mat')); %dont forget to change to Label2!!
Label = Label.Final_Label;
Label_sub = 0*Label; % pre-allocating
for c = 1 : length(sub_cluster_BIG_pruned)
    cells_of_Interest = sub_cluster_BIG_pruned(c).cells_of_interest_small{xp_ref_num};
    for i = 1 : length(cells_of_Interest) 
        Label_sub(Label==cells_of_Interest(i)) = c;
    end
end
for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(sub_cluster_BIG_pruned)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=c) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(cmp(c,1).*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(cmp(c,2).*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(cmp(c,3).*255*double(tmp))+temp(:,:,3);
    end
 imwrite(temp,['-',num2str(xp_ref_num),'1042ellip.tif'],'tiff','Compression','none','WriteMode','append');

end
end