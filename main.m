
% USER INPUTS ============================================================

image_type = 'sum';% 'max';

projection_method = 'linear';

data_path = fullfile(pwd,'test_data','ap_sagittal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO WORK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add paths ===============================================================

src_path = fullfile(pwd,'src');
addpath(src_path)


% Setup inputs ============================================================
image_path = fullfile(data_path,'images',image_type);
image_files = dir(fullfile(image_path,'*.tif'));
image_data_path = strcat({image_files.folder},'\',{image_files.name});

labels_path = fullfile(data_path,'labels');
label_files = dir(fullfile(labels_path,'*.mat'));
label_data_path = strcat({label_files.folder},'\',{label_files.name});

figure_save_path = fullfile(data_path,'projections',image_type,'figures','raw');
obj_save_path = fullfile(data_path,'projections',image_type,'objects');
data_save_path = fullfile(data_path,'projections',image_type,'data');

if ~exist(figure_save_path,'dir')
    mkdir(figure_save_path);
end

if ~exist(obj_save_path,'dir')
    mkdir(obj_save_path);
end

if ~exist(data_save_path,'dir')
    mkdir(data_save_path);
end
%% Make projections
% Step through images------------------------------------------------------

for i = 1:numel(image_data_path)
    
    disp(['Data set ',num2str(i),' of ', num2str(numel(image_data_path))])
    
    % Make image matrix----------------------------------------------------
    imageInfo = imfinfo(image_data_path{i});
    n_channels = numel(imageInfo);
    im = zeros(imageInfo(1).Height,imageInfo(1).Width,n_channels);
    
    for j = 1:n_channels
        im(:,:,j) = imread(image_data_path{i},j);
    end

    % Apply labels --------------------------------------------------------
    label = load(label_data_path{i});
    
    % make roi mask--------------------------------------------------------    
    hull_points =  label.gTruth.LabelData.roi{:};
    image_size = [imageInfo(1).Height,imageInfo(1).Width];
    
    [roi_mask, roi_mask_pixels] = makeRoiMask(hull_points,image_size);
    
    % make intensity matrix------------------------------------------------
    idx = sub2ind(size(im(:,:,1)),...
                       roi_mask_pixels(:,2),...
                       roi_mask_pixels(:,1));
    
    intensity = zeros(numel(idx),n_channels);
    
    for j = 1:n_channels
        x = im(:,:,j);
        intensity(:,j) = x(idx);
    end


    % make manifolds-------------------------------------------------------
    
    manifold = makeManifold(label,projection_method); 

    obj = projectMyPixel(roi_mask_pixels,intensity,manifold,projection_method);
       

    % Save data and objects -----------------------------------------------
    image_name = extractBefore(image_files(i).name,'.tif');

    save_name = fullfile(obj_save_path,['obj_',image_name,'.mat']);
    save(save_name, 'obj')

    csv_save_name = fullfile(data_save_path,['data_',image_name,'.txt']);
    writeData2csv(obj, csv_save_name)

    % Make figures --------------------------------------------------------
    stage = extractBetween(image_name,'stage=','_');
    sample = extractBetween(image_name,'sample=','_');
    projection_name = extractBetween(image_name,'projection=','_');
    channel_names = strcat('c0',[extractAfter(image_name,'_c0')]);
    
    marker_names = split(channel_names,'_');

    marker_colors = [0.3 0.3 0.3;
                     0.0 0.7 0.0;
                     0.0 0.0 0.7;
                     0.7 0.0 0.0];

    save_folder = string(fullfile(figure_save_path,...
        extractBefore(image_name,'_projection')));

    if ~exist(save_folder,'dir')
        mkdir(save_folder);
    end

%--------------------------------------------------------------------------
% PLOT RAW DATA -----------------------------------------------------------
    plotRawFluorescence(obj,save_folder,marker_names);

%--------------------------------------------------------------------------
% PLOT PROJECTIONS --------------------------------------------------------

    plotData(obj, projection_name, marker_names, marker_colors, save_folder)

end

%% Make Percentile Traces
trace_save_path = fullfile(data_path,'raw_percentile_traces',image_type);

if ~exist(trace_save_path,'dir')
    mkdir(trace_save_path);
end

obj_files = dir(fullfile(obj_save_path,'*.mat'));

for i = 1:numel(obj_files)
    
    load(fullfile(obj_files(i).folder,obj_files(i).name))
    
    traces = percentileTraces(obj,obj_files(i),trace_save_path);    
    
end



%% Static functions========================================================
function writeData2csv(obj, save_name)

    data = [obj.coordinates, obj.projection_absolute, obj.projection_fractional, obj.intensity];
    
    writematrix(data,save_name,'Delimiter',",")

 
end


function [mask, mask_pixels] = makeRoiMask(hull_points,image_size)%--------

    % Make a hull around the segmented ROI
    shp = alphaShape(hull_points,inf);
        
    % Make pixel coordinate pairs 
    [idx, jdx] = meshgrid([1:image_size(2)],[1:image_size(1)]);
    pixel_coordinates = [idx(:), jdx(:)];

    % get the pixels within the hull
    mask_pixel_idx = shp.inShape(pixel_coordinates);
    
    % make mask matrix
    mask_pixels= pixel_coordinates(mask_pixel_idx,:);
    
    mask = reshape(mask_pixel_idx,image_size(1),image_size(2))==1;
        
end%-----------------------------------------------------------------------


function manifold = makeManifold(segmentation, method)%----------------------------
      
    field_names = fieldnames(segmentation.gTruth.LabelData);
                
    switch method

        case 'linear'
            
            manifold = segmentation.gTruth.LabelData.(field_names{2}){:};
            
        case 'radial'
            
            ring = segmentation.gTruth.LabelData.(field_names{2}){:};
            focus = segmentation.gTruth.LabelData.(field_names{3}){:};

            manifold{1} =  ring ;
            manifold{2} = focus(1:2) + focus(3:4)/2;

        case 'conical'
            
            arc = segmentation.gTruth.LabelData.(field_names{2}){:}{:};% Arc
            apex = segmentation.gTruth.LabelData.(field_names{3}){:}; % Apex
            
            manifold{1} =  arc ;
            manifold{2} = apex(1:2) + apex(3:4)/2;
                        
    end
    
end%-----------------------------------------------------------------------

function f = plotRawFluorescence(obj,save_folder,marker_names )        


    f = obj.intensity;
    f_logged = log10(f+1);     
    
    x = obj.coordinates(:,1);
    y = obj.coordinates(:,2);
    
    % ROI with signal -----------------------------------------------------
    
    fig_1 = figure(1);
    
    for i = 1:size(f_logged,2)
        
        clf
        
        scatter3(x,y,f(:,i),1,f_logged(:,i))
        view(2)
        
        caxis([min(f_logged(:,i)), max(f_logged(:,i))])
        
        colorbar

        title(marker_names{i})

        axis equal
        axis ij

        save_name = ['roi_' marker_names{i}];
        saveFigure(fig_1,save_folder,save_name)
                
    end

    close(fig_1)

    % Fluorescence distribution -------------------------------------------

    fig_1 = figure(1);
    clf
    fig_1.Position = [0   0  1000  1000];
    plotmatrix(f)
    axis equal
    title({'Raw Fluorescence',strjoin([strcat(marker_names(1:end-1)',','),marker_names(end)])})

    save_name = 'fluorescence_raw';
    saveFigure(fig_1,save_folder,save_name)
    close(fig_1)

    fig_2 = figure(2);
    clf
    fig_2.Position = [0   0  1000  1000];
    plotmatrix(f_logged)
    axis equal
    title({'log10(Raw Fluorescence + 1)',strjoin([strcat(marker_names(1:end-1)',','),marker_names(end)])})

    save_name = 'fluorescence_logged';
    saveFigure(fig_2,save_folder,save_name)
    close(fig_2)

end

function plotData(obj,projection_name,marker_names,marker_colors,save_path)
 
    fluorescence = log10(obj.intensity+1);
    
    projection_type = {'projection_absolute', 'projection_fractional'};

    for i = 1:2

        fig = figure(1);
        fig.Position = [0   0  1000  1000];
        clf
        hold on   

        for j = 1:size(fluorescence,2)    
        
            subplot(size(fluorescence,2)  ,1,j)

            scatter(obj.(projection_type{i}),fluorescence(:,j),3,...
                    'MarkerFaceColor',marker_colors(j,:),...
                    'MarkerFaceAlpha',0.1, 'MarkerEdgeColor','none')                       
            
            title(marker_names{j}) 

            ylim([min(fluorescence(:,j)), max(fluorescence(:,j))])
            

        end

        save_name = strcat(projection_name, '_', projection_type{i});
        saveFigure(fig,save_path,save_name)      

    end

    close(fig)

end

function saveFigure(fig,save_folder,save_name)
            
            fullSaveName = fullfile(save_folder,save_name);
             
            savefig(fig,strcat(fullSaveName,'.fig'))
            saveas(fig,strcat(fullSaveName,'.png'),'png')
            saveas(fig,strcat(fullSaveName,'.svg'),'svg')
            
end
