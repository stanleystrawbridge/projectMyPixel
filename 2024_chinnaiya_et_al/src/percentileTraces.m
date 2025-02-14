classdef percentileTraces


    properties %===========================================================
        projection
        file
        save_path

        n_bins = 500
        centiles = [0 10 25 50 75 90]
        bins_absolute
        bins_relative
        values

        % Sampel Metadata
        stage
        sample_number
        marker_names
        number_of_markers



    end %==================================================================


    methods %==============================================================


        function traces = percentileTraces(projection,file,save_path) %---- 

            % initializ inputs
            traces.projection = projection;
            traces.file = file;
            
            traces = setIntput(traces,save_path);

            traces = traces.makeTraces;

            traces.plotTraces;

            traces.overlayTraces;

            traces.overlayTracesWithIQR;

            save(traces.save_path,'traces')

        end %--------------------------------------------------------------

        
        function traces = setIntput(traces,save_path) %--------------------

            % Sampel Metadata
            metadata = traces.file.name;

            traces.stage = extractBetween(metadata,'stage=','_');
            traces.sample_number = extractBetween(metadata,'sample=','_');
 
            folder = extractBetween(metadata,'stage=','_projection');
            traces.save_path = fullfile(save_path,folder{:});
            if ~exist(traces.save_path,'dir')
                mkdir(traces.save_path);
            end
            
            channel_names = strcat('c0',[extractAfter(metadata,'_c0')]);
            traces.marker_names = split(channel_names,'_');
            traces.number_of_markers = numel(traces.marker_names);

            if size(traces.projection.intensity,1) < traces.n_bins
                traces.n_bins = size(traces.projection.intensity,2);
            end

            traces.values = zeros(traces.n_bins,...
                                  numel(traces.centiles),...
                                  traces.number_of_markers);

        end %--------------------------------------------------------------


        function traces = makeTraces(traces) %-----------------------------
            
            [min_d,max_d] = bounds(traces.projection.projection_absolute);
            
            bins = linspace(min_d,max_d+1,traces.n_bins+1)';

            traces.bins_absolute = bins(1:end-1);
            traces.bins_relative =  traces.bins_absolute - min(traces.bins_absolute);
            traces.bins_relative =  traces.bins_absolute./max(traces.bins_absolute);

            for i = 1:traces.n_bins

                for j = 1:traces.number_of_markers
                
                    start = bins(i);
                    stop =  bins(i+1);  
                    
                    idx = traces.projection.projection_absolute >=start & ...
                    traces.projection.projection_absolute < stop;
                        
                    traces.values(i,:,j) = prctile(traces.projection.intensity(idx,j),traces.centiles);
    
                end

            end

        end %--------------------------------------------------------------


        function plotTraces(traces) %--------------------------------------

            smoothing = 20;

            n_centiles = numel(traces.centiles);
            c = [linspace(0,0.3922,n_centiles); 
                linspace(0,0.5843,n_centiles);
                linspace(0,0.9294,n_centiles)]';

            line_widths = [2.5*ones(1,3), 5, 2.5*ones(1,2)];

            for i = 1:traces.number_of_markers

                fig1 = figure(1);clf;hold on;
                fig2 = figure(2);clf;hold on;

                for j = 1:n_centiles

                    figure(1)
                    plot(traces.bins_absolute, smooth(traces.values(:,j,i),...
                        smoothing),LineWidth=line_widths(j),Color=c(j,:))
                    
                    figure(2)
                    plot(traces.bins_relative, smooth(traces.values(:,j,i),...
                        smoothing),LineWidth=line_widths(j),Color=c(j,:))

                end

                if  isempty(extractBefore(traces.marker_names{i},'.'))
                    marker = traces.marker_names{i};
                else
                    marker = extractBefore(traces.marker_names{i},'.');
                end                          

                title(strcat(traces.stage{:},'/',traces.sample_number ,'/',marker))

                figure(1)

                legend(['min',string(traces.centiles(2:end))],...
                    'Location','southoutside','Orientation','horizontal')
                xlabel('Distance (\mum)')
                ylabel('Fluorescence')                
                title(strcat(traces.stage{:},'/',traces.sample_number ,'/',marker))
                set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')

                axis tight
                
                save_name = [traces.marker_names{i}, '_absolute_distance'];
                traces.saveFigure(fig1,save_name);

                figure(2)
                legend(['min',string(traces.centiles(2:end))],...
                    'Location','southoutside','Orientation','horizontal')
                xlabel('Distance')
                ylabel('Fluorescence')

                set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')
                
                axis tight
                
                save_name = [traces.marker_names{i}, '_relative_distance'];             
                traces.saveFigure(fig2,save_name);

            end

        end %--------------------------------------------------------------
       
        function overlayTraces(traces) %--------------------------------------
            
            smoothing = 20;

            x = permute(traces.values(:,4,:), [1 3 2]);
            
            for i = 1:traces.number_of_markers
                x(:,i) = smooth(x(:,i),smoothing);
            end
            % 
            % x = x - repmat(min(x),size(x,1),1);
            % % x = x ./ repmat(max(x),size(x,1),1);
            
            fig1 = figure(1);
            clf
            hold on

            fig2 = figure(2);
            clf
            hold on

            c = {'g','r','b'};
            labels = {};

            for i = 2:traces.number_of_markers
                
                figure(1);

                plot(traces.bins_relative,x(:,i),c{i-1},'LineWidth',3)
                shg

                axis([0 max(traces.bins_relative) 0 1])

                figure(2);

                plot(traces.bins_absolute,x(:,i),c{i-1},'LineWidth',3)

                axis([0 max(traces.bins_absolute) 0 1])

                shg
                          
                if  isempty(extractBefore(traces.marker_names{i},'.'))
                    labels = [labels;traces.marker_names{i}];
                else
                    label = {extractBefore(traces.marker_names{i},'.')};
                    labels = [labels;label];
                end

            end

            figure(1);
            grid on
            xlabel('Distance')
            ylabel('Fluorescence')
            title(strcat(traces.stage{:},'/',traces.sample_number))
            set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')
            legend(labels,'Location','southoutside',Orientation='horizontal' )
            
            figure(2);
            grid on
            xlabel('Distance  (\mum)')
            ylabel('Fluorescence')
            title(strcat(traces.stage{:},'/',traces.sample_number))
            set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')
            legend(labels,'Location','southoutside',Orientation='horizontal' )

            save_name = 'overlay_relative_distance';             
            traces.saveFigure(fig1,save_name);
                        
            save_name = 'overlay_absolute_distance';             
            traces.saveFigure(fig2,save_name);

        end %--------------------------------------------------------------


        function overlayTracesWithIQR(traces) %--------------------------------------
            
            smoothing = 20;
            
            x_25 = permute(traces.values(:,3,:), [1 3 2]);
            x_50 = permute(traces.values(:,4,:), [1 3 2]);
            x_75 = permute(traces.values(:,5,:), [1 3 2]);
            
            for i = 1:traces.number_of_markers

                x_25(:,i) = smooth(x_25(:,i),smoothing);
                x_50(:,i) = smooth(x_50(:,i),smoothing);
                x_75(:,i) = smooth(x_75(:,i),smoothing);

            end

            x_25 = x_25 - repmat(min(x_50),size(x_25,1),1);
            x_75 = x_75 - repmat(min(x_50),size(x_75,1),1);
            x_50 = x_50 - repmat(min(x_50),size(x_50,1),1);

    
            x_25 = x_25 ./ repmat(max(x_50),size(x_25,1),1);
            x_75 = x_75 ./ repmat(max(x_50),size(x_75,1),1);
            x_50 = x_50 ./ repmat(max(x_50),size(x_50,1),1);
                   
            fig1 = figure(1);
            clf
            hold on

            fig2 = figure(2);
            clf
            hold on

            c = {'g','r','b'};
            labels = {};

            for i = 2:traces.number_of_markers
                
                figure(1);
                hold on
                x_grid = [traces.bins_relative; flip(traces.bins_relative)];
                y = [x_25(:,i); flip(x_75(:,i))];
                
                fill(x_grid,y,c{i-1},'EdgeColor','none','FaceAlpha',0.2)

                plot(traces.bins_relative,x_50(:,i),c{i-1},'LineWidth',5)

                axis([0 max(traces.bins_relative) 0 1])

                shg

                figure(2);
                hold on
                x_grid = [traces.bins_absolute; flip(traces.bins_absolute)];
                y = [x_25(:,i); flip(x_75(:,i))];
                
                fill(x_grid,y,c{i-1},'EdgeColor','none','FaceAlpha',0.2)

                plot(traces.bins_absolute,x_50(:,i),c{i-1},'LineWidth',5)

                axis([0 max(traces.bins_absolute) 0 1])

                shg
                          
                if  isempty(extractBefore(traces.marker_names{i},'.'))
                    labels = [labels;traces.marker_names{i};traces.marker_names{i}];
                else
                    label = {extractBefore(traces.marker_names{i},'.')};
                    labels = [labels;label; label];
                end

            end

            figure(1);
            grid on
            xlabel('Distance')
            ylabel('Fluorescence')
            ylim([0 1])
            title(strcat(traces.stage{:},'/',traces.sample_number))
            set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')
            legend(labels,'Location','southoutside',Orientation='horizontal' )
            set(gcf,'Position' ,[0 0 950 950])


            figure(2);
            grid on
            xlabel('Distance  (\mum)')
            ylabel('Fluorescence')
            ylim([0 1])
            title(strcat(traces.stage{:},'/',traces.sample_number))
            set(gca,'FontName','Arial','FontSize',15,'FontWeight','bold')
            legend(labels,'Location','southoutside',Orientation='horizontal' )
            set(gcf,'Position' ,[0 0 950 950])

            save_name = 'overlayIQR_relative_distance';             
            traces.saveFigure(fig1,save_name);
            close(fig1)
            
            save_name = 'overlayIQR_absolute_distance';             
            traces.saveFigure(fig2,save_name);
            close(fig2)
            
        end %--------------------------------------------------------------


        function saveFigure(traces,fig,save_name) %-----------------------
            
            full_save_name = fullfile(traces.save_path,save_name);
                        
            savefig(fig,strcat(full_save_name,'.fig'))
            saveas(fig,strcat(full_save_name,'.png'),'png')
            saveas(fig,strcat(full_save_name,'.svg'),'svg')

        end %--------------------------------------------------------------


    end %==================================================================
   
end