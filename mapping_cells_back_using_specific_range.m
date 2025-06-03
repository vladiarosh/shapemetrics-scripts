%% *Mapping cells back for a parameter of choice above the threshold of choice (based on raw values)*

param_name = {'Volume', 'vol-surfArea ratio', 'ellipticity', 'elongation', 'long_axis_length', 'minor_axis', 'intermediate_axis'};
parent_folder = ''; % insert file path to a root folder containing subfolders of all tips
subfolders = {'Initial bud 1', 'Initial bud 2','Ampulla 1', 'Ampulla 2', 'Ampulla 3', 'Asymmetrical ampulla 1', 'Asymmetrical ampulla 2', 'Asymmetrical ampulla 3', 'Bifurcated tip 1', 'Bifurcated tip 2', 'Bifurcated tip 3',  'DKO tip 1', 'DKO tip 2', 'DKO tip 3'};

num_cells_of_interest = zeros(length(subfolders), length(param_name));

for folder_idx = 1:length(subfolders)
    current_folder = fullfile(parent_folder, subfolders{folder_idx});
    
    % Load stats_matrix_all and Final_Label from the current folder
    load(fullfile(current_folder, 'stats_matrix_all'));  % Replace with your actual file names
    load(fullfile(current_folder, 'Final_Label'));  % Replace with your actual file names

    for parameter = 1
        threshold = 100;
        stat_vector = stats_matrix_all(:, parameter);
        index = 1:length(stat_vector);
        cells_of_interest = index(stat_vector < threshold);
        num_cells_of_interest(folder_idx, parameter) = length(cells_of_interest);

        Label_sub = 0 * Final_Label;

        for i = 1:length(cells_of_interest)
            Label_sub(Final_Label == cells_of_interest(i)) = cells_of_interest(i);
        
        end
        for z = 1:size(Label_sub, 3)
            temp = zeros(size(Label_sub, 1), size(Label_sub, 2), 3, 'uint8');
            for c = 1:length(cells_of_interest)
                tmp = Label_sub(:, :, z);
                tmp(tmp ~= cells_of_interest(c)) = 0;
                tmp = tmp > 0;

                cmp = zeros(10, 3);
                cmp(1, :) = [0.1, 0.9, 1];
                cmp(2, :) = [0.8, 0.6, 0.9];
                cmp(3, :) = [0.6, 0.9, 0.4];
                cmp(4, :) = [0.8, 0.2, 0.4];
                cmp(5, :) = [0.9, 1, 0];
                cmp(6, :) = [0.1, 0.4, 1];
                cmp(7, :) = [0, 0.8, 0.8];
                cmp(8, :) = [0.9, 0, 1];
                cmp(9, :) = [1, 0.6, 0];
                cmp(10, :) = [0.9, 0.4, 0.5];

                temp(:, :, 1) = uint8(cmp(6, 1) * 255 * double(tmp)) + temp(:, :, 1);
                temp(:, :, 2) = uint8(cmp(6, 2) * 255 * double(tmp)) + temp(:, :, 2);
                temp(:, :, 3) = uint8(cmp(6, 3) * 255 * double(tmp)) + temp(:, :, 3);
            end
    image_filename = [param_name{parameter}, 'below100cubicmicrons.tif'];
    imwrite(temp, fullfile(current_folder, image_filename), 'tif', 'Compression', 'packbits', 'WriteMode', 'append');
            
        end
    end
end

% Save the number of cells_of_interest as a separate file
save('num_cells_of_interest.mat', 'num_cells_of_interest');
%% 
%% Visualize single parameters (all 5 at the same time) on a range of iterest
% parameter 1: volume parameter 2: vol-surfArea ratio parameter 3: ellipticity 
% parameter 4: elongation parameter 5: the length of the longest axis

param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};
for parameter=1

rangeA = [200 250]; % change range minimum (left) and maximum (rigth) value
%Important note - in this version we use raw values, not zscored-values. If
%you want to use zscored-matrices, zscore the stats_matrix_all file and use
%zscored_matrix as stat_vector


    stat_vector=stats_matrix_all(:,parameter);
    index=1:length(stat_vector);
    
    cells_of_interest=index(stat_vector'>rangeA(1) & stat_vector' < rangeA(2)); %use > or < depending on larger or smaller stat required

Label_sub = 0*Final_Label;  


  
    for i = 1 : length(cells_of_interest) 
        Label_sub(Final_Label == cells_of_interest(i)) = cells_of_interest(i);
    end
    %load Label_sub in volume viewer to explore 


for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(cells_of_interest)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest(c)) = 0;
        tmp = tmp>0;
cmp(1,:)  = [0.1,0.9,1];   % light blue    
cmp(2,:)  = [0.8,0.6,0.9]; % lilac
cmp(3,:)  = [0.6,0.9,0.4]; % green   
cmp(4,:)  = [0.8,0.2,0.4]; % red pink   
cmp(5,:)  = [0.9,1,0];     % yellow
cmp(6,:)  = [0.1,0.4,1];   % dark blue 
cmp(7,:)  = [0,0.8,0.8];   % cyan 
cmp(8,:)  = [0.9,0,1];     % purple 
cmp(9,:)  = [1,0.6,0];     % orange   
cmp(10,:) = [0.9,0.4,0.5]; % baby pink 
        temp(:,:,1) = uint8(cmp(6,1).*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(cmp(6,2).*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(cmp(6,3).*255*double(tmp))+temp(:,:,3);
    end
    imwrite(temp,[param_name{parameter},'_200-250.tif'],'tif','Compression','packbits','WriteMode','append'); %change '_small.tif' to '_large.tif' when viewing the larger populations
end

end
%% Visualize single parameter on a set of ranges

%MULTIPLE RANGE ANALYSIS FOR RAW DATA
param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};

% define the different ranges of volume
volume_ranges = [50 80; 80 110; 110 140; 140 170; 170 200; 200 250];

for range_index = 1:size(volume_ranges, 1)
rangeA = volume_ranges(range_index, :);
% rest of the code remains the same
stat_vector=stats_matrix_all(:,1);
index=1:length(stat_vector);

cells_of_interest=index(stat_vector'>rangeA(1) & stat_vector' < rangeA(2));

Label_sub = 0*Final_Label;  

for i = 1 : length(cells_of_interest) 
    Label_sub(Final_Label == cells_of_interest(i)) = cells_of_interest(i);
end

for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(cells_of_interest)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest(c)) = 0;
        tmp = tmp>0;
        cmp(1,:)  = [0.1,0.9,1];   % light blue    
        cmp(2,:)  = [0.8,0.6,0.9]; % lilac
        cmp(3,:)  = [0.6,0.9,0.4]; % green   
        cmp(4,:)  = [0.8,0.2,0.4]; % red pink   
        cmp(5,:)  = [0.9,1,0];     % yellow
        cmp(6,:)  = [0.1,0.4,1];   % dark blue 
        cmp(7,:)  = [0,0.8,0.8];   % cyan 
        cmp(8,:)  = [0.9,0,1];     % purple 
        cmp(9,:)  = [1,0.6,0];     % orange   
        cmp(10,:) = [0.9,0.4,0.5]; % baby pink 
        temp(:,:,1) = uint8(cmp(6,1).*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(cmp(6,2).*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(cmp(6,3).*255*double(tmp))+temp(:,:,3);
    end
    imwrite(temp,[param_name{1},'_',num2str(rangeA(1)),'-',num2str(rangeA(2)),'.tif'],'tif','Compression','packbits','WriteMode','append');
end
end