%% Top- and bottom-10% of cells for a parameter of choice

   %TAILS ANALYSIS
%% Visualize single parameters (all 5 at the same time) and save 'mapped back' cell images as tifs)
% parameter 1: volume              
% parameter 2: vol-surfArea ratio              
% parameter 3: ellipticity          
% parameter 4: elongation
% parameter 5: the length of the longest axis
param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};
for parameter=1
zscored_matrix = zscore(stats_matrix_all); 
bot_treshold_in_z_score=quantile(zscored_matrix(:,parameter),0.1);
top_treshold_in_z_score=quantile(zscored_matrix(:,parameter),0.9);

    stat_vector=zscored_matrix(:,parameter);
    index=1:length(stat_vector);
   
    cells_of_interest_top_10=index(stat_vector'>top_treshold_in_z_score) %use > or < depending on larger or smaller stat required
    cells_of_interest_bot_10=index(stat_vector'<bot_treshold_in_z_score)
    Label_sub = 0*Final_Label;  
 
    for i = 1 : length(cells_of_interest_top_10)
        Label_sub(Final_Label == cells_of_interest_top_10(i)) = cells_of_interest_top_10(i);
    end
    for i = 1 : length(cells_of_interest_bot_10)
        Label_sub(Final_Label == cells_of_interest_bot_10(i)) = cells_of_interest_bot_10(i);
    end
    %load Label_sub in volume viewer to explore
for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(cells_of_interest_top_10)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest_top_10(c)) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(1*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(0*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(0*255*double(tmp))+temp(:,:,3);
    end
    for c = 1 : length(cells_of_interest_bot_10)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest_bot_10(c)) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(0*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(0*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(1*255*double(tmp))+temp(:,:,3);
    end
    imwrite(temp,[param_name{parameter},'_tails.tif'],'tif','Compression','packbits','WriteMode','append'); %change '_small.tif' to '_large.tif' when viewing the larger populations
end
end
 
%% Top-10% of cells for a parameter of choice

   %TAILS ANALYSIS
%% Visualize single parameters (all 5 at the same time) and save 'mapped back' cell images as tifs)
% parameter 1: volume              
% parameter 2: vol-surfArea ratio              
% parameter 3: ellipticity          
% parameter 4: elongation
% parameter 5: the length of the longest axis
param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};
for parameter=1
zscored_matrix = zscore(stats_matrix_all); 
bot_treshold_in_z_score=quantile(zscored_matrix(:,parameter),0.1);
top_treshold_in_z_score=quantile(zscored_matrix(:,parameter),0.9);

    stat_vector=zscored_matrix(:,parameter);
    index=1:length(stat_vector);
   
    cells_of_interest_top_10=index(stat_vector'>top_treshold_in_z_score) %use > or < depending on larger or smaller stat required
    
    Label_sub = 0*Final_Label;  
 
    for i = 1 : length(cells_of_interest_top_10)
        Label_sub(Final_Label == cells_of_interest_top_10(i)) = cells_of_interest_top_10(i);
    end
   
    end
    %load Label_sub in volume viewer to explore
for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(cells_of_interest_top_10)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest_top_10(c)) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(1*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(0*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(0*255*double(tmp))+temp(:,:,3);
    end
  
    imwrite(temp,[param_name{parameter},'_toptail.tif'],'tif','Compression','packbits','WriteMode','append'); %change '_small.tif' to '_large.tif' when viewing the larger populations
end
%% *Top-10% of cells for a parameter of choice for all UB tips at once*

param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};
parent_folder = ''; % insert file path to a root folder containing subfolders of all tips

% List of subfolders within the parent folder
subfolders = {'Ampulla 1', 'Ampulla 2' ,'Ampulla 3', 'Asymmetrical ampulla 1','Asymmetrical ampulla 2','Asymmetrical ampulla 3','Bifurcated tip 1','Bifurcated tip 2','Bifurcated tip 3','Initial bud 1','Initial bud 2','DKO tip 1','DKO tip 2','DKO tip 3'};  % Add more folders as needed
cells_of_interest_top_10_table = {};
for folder_idx = 1:length(subfolders)
    current_folder = fullfile(parent_folder, subfolders{folder_idx});
    
    % Load stats_matrix_all and Final_Label from the current folder
    load(fullfile(current_folder, 'stats_matrix_all'));  % Replace with your actual file names
    load(fullfile(current_folder, 'Final_Label'));  % Replace with your actual file names
    cells_of_interest_top_10_row = {};
    for parameter = 4
        zscored_matrix = zscore(stats_matrix_all);
        top_treshold_in_z_score = quantile(zscored_matrix(:, parameter), 0.9);
        stat_vector = zscored_matrix(:, parameter);
        index = 1:length(stat_vector);
        
        cells_of_interest_top_10 = index(stat_vector' > top_treshold_in_z_score); 
       
        Label_sub = 0 * Final_Label;
        cells_of_interest_top_10_row{parameter} = cells_of_interest_top_10;
      
        for i = 1 : length(cells_of_interest_top_10)
            Label_sub(Final_Label == cells_of_interest_top_10(i)) = cells_of_interest_top_10(i);
           
        end
        cells_of_interest_top_10_table{folder_idx} = cells_of_interest_top_10_row;
        % Rest of your code for processing and saving images
        for z = 1 : size(Label_sub, 3)
            temp = zeros(size(Label_sub, 1), size(Label_sub, 2), 3, 'uint8');
            for c = 1 : length(cells_of_interest_top_10)
                tmp = Label_sub(:,:,z);
                tmp(tmp ~= cells_of_interest_top_10(c)) = 0;
                tmp = tmp > 0;
                temp(:,:,1) = uint8(0 * double(tmp)) + temp(:,:,1);
                temp(:,:,2) = uint8(0 * double(tmp)) + temp(:,:,2);
                temp(:,:,3) = uint8(255 * double(tmp)) + temp(:,:,3);
            end
            image_filename = [param_name{parameter}, '_toptail_cb_adj_2.tif'];
            imwrite(temp, fullfile(current_folder, image_filename), 'tif', 'Compression', 'packbits', 'WriteMode', 'append');

           
            
        end
    end
   
end

%% *Bottom-10% of cells for a parameter of choice for all UB tips at once*

param_name={'Volume','vol-surfArea ratio','ellipticity','elongation','long_axis_length'};
parent_folder = ''; % insert file path to a root folder containing subfolders of all tips

% List of subfolders within the parent folder
subfolders = {'Ampulla 1', 'Ampulla 2', 'Ampulla 3', 'Asymmetrical ampulla 1','Asymmetrical ampulla 2','Asymmetrical ampulla 3','Bifurcated tip 1','Bifurcated tip 2','Bifurcated tip 3','Initial bud 1','Initial bud 2','DKO tip 1','DKO tip 2','DKO tip 3'};  % Add more folders as needed

for folder_idx = 1:length(subfolders)
    current_folder = fullfile(parent_folder, subfolders{folder_idx});
    
    % Load stats_matrix_all and Final_Label from the current folder
    load(fullfile(current_folder, 'stats_matrix_all'));  % Replace with your actual file names
    load(fullfile(current_folder, 'Final_Label'));  % Replace with your actual file names
    
    for parameter = 1
        zscored_matrix = zscore(stats_matrix_all);
        bot_treshold_in_z_score=quantile(zscored_matrix(:,parameter),0.1);
        stat_vector = zscored_matrix(:, parameter);
        index = 1:length(stat_vector);
        
        cells_of_interest_bot_10=index(stat_vector'<bot_treshold_in_z_score);
        Label_sub = 0 * Final_Label;
        
        for i = 1 : length(cells_of_interest_bot_10)
        Label_sub(Final_Label == cells_of_interest_bot_10(i)) = cells_of_interest_bot_10(i);
        end
        
        % Rest of your code for processing and saving images
        for z = 1 : size(Label_sub, 3)
            temp = zeros(size(Label_sub, 1), size(Label_sub, 2), 3, 'uint8');
            for c = 1 : length(cells_of_interest_bot_10)
        tmp = Label_sub(:,:,z);
        tmp(tmp~=cells_of_interest_bot_10(c)) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(2*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(2*127*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(2*80*double(tmp))+temp(:,:,3);
            end

            % Save the image in the current folder
            image_filename = [param_name{parameter}, '_bottomtail.tif'];
            imwrite(temp, fullfile(current_folder, image_filename), 'tif', 'Compression', 'packbits', 'WriteMode', 'append');
        end
    end
end