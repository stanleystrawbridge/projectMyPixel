classdef projectMyPixel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties%============================================================
        
        % Inputs
        coordinates
        intensity
        manifold        
        method                              

        n_coordinates
        
        % Outputs
        projection_absolute
        projection_fractional

    end%===================================================================

    
    methods%===============================================================        
        
        function obj = projectMyPixel(varargin) %--------------------------

            obj = obj.importData(varargin);            

            obj = obj.linearMethod;
            
        end %--------------------------------------------------------------

        
        %------------------------------------------------------------------
        % INITIALIZATION---------------------------------------------------
               
        function obj = importData(obj,args) %------------------------------
            
            obj.coordinates = args{1};
            obj.intensity = args{2};
            obj.manifold = args{3};
            obj.method = args{4};

            obj.n_coordinates = size(obj.coordinates,1);

        end %--------------------------------------------------------------                                  


        %------------------------------------------------------------------
        % PROJECTION METHODS-----------------------------------------------
        
        function obj = linearMethod(obj) %---------------------------------
                                    
            line = obj.manifold;

            % Fix           
                       
            line_distance = cumsum(sqrt(sum(diff(line).^2,2)));
            
            n_line_points = size(line,1);
            
            obj.projection_absolute = zeros(size(obj.coordinates,1),1);
            obj.projection_fractional = obj.projection_absolute;
                        
            % Calculate the distance of each query point to each
            % linesegment
            for i = 1:obj.n_coordinates
            
                current_point = obj.coordinates(i,:);
                
                idx = projectMyPixel.minDistance2Point(line,current_point);                                                     
                
                % find the second closest line point which is adjecent to
                % the closest line point                
                if idx==1
                    
                    segment = line(1:2,:);
                    segment_number = 1;

                    % project that point onto the line segment
                    proj_p = projectMyPixel.point2segment(segment,current_point); 
                    
                elseif idx == n_line_points
                    
                    segment = line(end-1:end,:);
                    segment_number = n_line_points-1;
                   
                    % project that point onto the line segment
                    proj_p = projectMyPixel.point2segment(segment,current_point); 
                    
                else
                
                    adjacent_points = line([idx-1,idx+1],:);
                
                    % project that point onto the line segment
                    segment_1 = line([idx(1)-1,idx(1)],:);
                    proj_p_1 = projectMyPixel.point2segment(segment_1,current_point); 
                    
                    segment_2= line([idx(1),idx(1)+1],:);
                    proj_p_2 = projectMyPixel.point2segment(segment_2,current_point); 


                    d = sqrt(sum(([proj_p_1;proj_p_2] - repmat(current_point,2,1)).^2,2));

                    [~,jdx] = min(d);                                   

                    if jdx == 1
                        proj_p = proj_p_1;
                        segment = [adjacent_points(jdx,:);
                              line(idx,:)];
                        segment_number = idx-1;
                    else
                        proj_p = proj_p_2;
                        segment = [line(idx,:);
                                   adjacent_points(jdx,:)];
                        segment_number = idx;
                    end
                    
                end
                                              
                % Calculate distance along segment
                distance_along_segment = sqrt(sum((proj_p - segment(1,:)).^2));
                
                if segment_number>1
                    distance_along_line = line_distance(segment_number-1) + ...
                                          distance_along_segment;
                else
                    distance_along_line = distance_along_segment;
                end

                % Save distance along line
                obj.projection_absolute(i) = distance_along_line;                            
                
            end 

            obj.projection_absolute =   obj.projection_absolute ...
                                      - min(obj.projection_absolute);
            
            % Save fractional distance along line
                obj.projection_fractional = ...
                   obj.projection_absolute./max(obj.projection_absolute);


        end %--------------------------------------------------------------
        

    end %==================================================================
        
     
    methods (Static) %=====================================================
                     
        
        function idx = minDistance2Point(reference_points,query_point) %---

            n_ref_points = size(reference_points,1);
            Q = repmat(query_point,n_ref_points,1);

            d = reference_points - Q;
            d = sqrt(sum(d.^2,2));

            [~,idx] = sort(d);  
            idx = idx(1);

        end%---------------------------------------------------------------


        function  proj_p = point2segment(segment,p)  %---------------------

            dims = size(p,2);
            
            if dims<3           
                p = [p, 0];            
                segment = [ segment, [0;0]];                                  
            end

            a = segment(1,:);
            b = segment(2,:);

            S = sqrt(sum((a-b).^2));
            A = sqrt(sum((a-p).^2));
            B = sqrt(sum((b-p).^2));
            
            ba = b - a; 
            pa = p - a;

            dProjP = sqrt(sum(cross(ba,pa).^2))./ sqrt(sum(ba.^2));
            proj_p = a + dot(pa,ba)./dot(ba,ba).*ba;
            
            % Determing if the projected point is within the line segment
            onOrOff = A.^2 + B.^2 - 2.*(dProjP.^2) < S.^2;
            
            % if it is not on it, take the closest point.
            if ~onOrOff

                d = sqrt(sum((segment - repmat(p,2,1)).^2,2));
                
                [~,jdx] = min(d); 

                proj_p = segment(jdx,:);

            end
            
            if dims < size(proj_p,2)
                proj_p = proj_p(1:2);
            end
            
        end %---------------------------------------------------------------

    
        function [radius_length, distance_along_radius] = line2segment(focus,boundary,point)
            
            n_segments = size(boundary,1)-1;

            reference_segment = [focus; point];

            intersection_points = zeros(n_segments,2);
            isPositiveDirection = zeros(n_segments,1);
            isInRingSegment = zeros(n_segments,1);

            for i = 1:n_segments

                boundary_segment = boundary(i:i+1,:);               
                
                % calculate intersection point
                x_0 = projectMyPixel.lineLineIntersection([boundary_segment;reference_segment]);

                intersection_points(i,:) = x_0;

                % determine if intersection point is in the positive
                % direction from focus to point                               
                isPositiveDirection(i) = ...
                    all(((point - focus) > 1e-5) == ((x_0 - focus) > 1e-5),2);
                
                % determine if intersection point is within ring segment
                ring_segment_length = sqrt(sum(diff(boundary_segment).^2));
                subsegment_length_sum = sum(sqrt(sum((boundary_segment - [x_0;x_0]).^2,2)));

                

                isInRingSegment(i) =  ...
                    abs(ring_segment_length-subsegment_length_sum) < 1e-5; % less than a small value, numerically 0
            
            end

            % Calculate distances
            radius_lengths = sqrt(sum((intersection_points...
                            - repmat(focus,n_segments,1)).^2,2));

            distance_along_radius = sqrt(sum((point - focus).^2));

            % get the correct distance
            keep_idx = isPositiveDirection & isInRingSegment;
            if sum(keep_idx) == 1

                radius_length = radius_lengths(keep_idx);
                
                
            elseif sum(isPositiveDirection) >0 

                %Take minimumem distance 
                
                radius_lengths(~isPositiveDirection) = inf;

                radius_length = min(radius_lengths);

            else

                % If/when a user comes back to this, you will need to
                % select the correct intersection point based on vetor
                % magnitude.
                error('no unique intersection point, Get Stan to fix it!')

            end

        end %--------------------------------------------------------------


        function x_0 = lineLineIntersection(data) %------------------------
 
            x = data(:,1);
            y = data(:,2);
            % wikipedia 
            % https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
            D = (x(1)-x(2))*(y(3)-y(4)) - (x(3)-x(4))*(y(1)-y(2));

            N_x = (x(1)*y(2)-x(2)*y(1))*(x(3)-x(4)) ...
                 - (x(3)*y(4)-x(4)*y(3))*(x(1)-x(2));
            N_y = (x(1)*y(2)-x(2)*y(1))*(y(3)-y(4)) ...
                 - (x(3)*y(4)-x(4)*y(3))*(y(1)-y(2));

            x_0 =[N_x/D, N_y/D];
        end %--------------------------------------------------------------

    end %==================================================================
    
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

