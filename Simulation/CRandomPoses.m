classdef CRandomPoses
    %CRandomPoses Class for generation of random poses
    %   Detailed explanation goes here
    
    properties
        min_d
        max_d
        min_ang
        ang_z_max
        ang_x_max
        ang_y_max
        device
    end
    
    methods
        function this = CRandomPoses( S )
            % S is an input struct with fields
            % 	min_d
            %   max_d
            %   min_ang
            %   ang_z_max
            %   ang_x_max
            %   ang_y_max
            %   device
            %
            % It's possible to use the output of configuration file
            % as input:
            % S = readConfigFile( config_file, config_section );
            
            % Copy fields
            this.min_d = S.min_d;
            this.max_d = S.max_d;
            this.min_ang = S.min_ang;
            this.ang_z_max = S.ang_z_max;
            this.ang_x_max = S.ang_x_max;
            this.ang_y_max = S.ang_y_max;
            this.device = S.device;    
        end
        
        function [R, t] = gen( this, N )
            % this.gen(N) returns N poses fulfilling
            % the conditions imposed in constructor
            extractStructFields( this );
            
            max_ang = 90 - min_ang;
            R_Cam_LRF = [ 0 -1 0 ; 0 0 -1 ; 1 0 0 ];
            if strcmp(device,'Lidar')
                % Rotation matrix to correct axes
                R_rel = R_Cam_LRF;
            else
                R_rel = eye(3);
            end
            
            %% World to Cam translation
            % Directions from W to LIDAR
            % Generate N valid values
            cont = 0;
            v_w_c = zeros(3,N);
            t = zeros(3,N);
            while cont < N
                v_w_c_ = rand( 3, 1 );
                v_w_c_ = snormalize( v_w_c_ );
                % Distance from W to LIDAR
                d_w_c_  = min_d + (max_d-min_d) * rand( 1 );
                ang_ = acosd( v_w_c_ );
                if ~any( ang_ > max_ang )
                    cont = cont + 1;
                    v_w_c(:,cont) = v_w_c_;
                    t(:,cont) = v_w_c_ * d_w_c_;
                end
            end
            
            %% World to Cam rotation
            % Set firstly camera pointing to World origin with X axis in XY world plane
            % direction
            % TODO: Impose origin-pointing axis (instead of device)
            R_z = -v_w_c;
            R_x = -skew([0 0 1]') * R_z;
            R_x = snormalize( R_x );
            R_y = cross(R_z,R_x,1);
            % Add random rotation on Cam X and Z axis
            rand_ang_z = ang_z_max * (-1 + 2*rand(1,N)); % ang_z_max typically pi
            rand_ang_x = ang_x_max * (-1 + 2*rand(1,N)); % ang_x_max typically FOV
            rand_ang_y = ang_y_max * (-1 + 2*rand(1,N));
            R = reshape( [R_x;R_y;R_z], 3,3,N );
            for i=1:N
                % Check premultiplication and good order
                %     R(:,:,i) = R(:,:,i) * RotationX(ang_x(i)) * RotationZ(ang_z(i));
                %     R(:,:,i) = R(:,:,i) * RotationX(rand_ang_x1(i)) * RotationZ(rand_ang_z(i)) * RotationX(rand_ang_x(i));
                %     R(:,:,i) = R(:,:,i) * RotationY(rand_ang_x1(i)) * RotationZ(rand_ang_z(i)) * RotationX(rand_ang_x(i));
                R(:,:,i) = R(:,:,i) * ...
                           RotationZ(rand_ang_z(i)) * ...
                           RotationX(rand_ang_x(i)) * ...
                           RotationY(rand_ang_y(i)) * ...
                           R_rel;
            end
            
            if N > 1
                % Transform output to cell arrays:
                t = mat2cell( t, 3, ones(1,N) );
                R = mat2cell( R, 3, 3, ones(1,N) );
                R = permute( R, [1 3 2] );
            end
        end
    end
    
end

