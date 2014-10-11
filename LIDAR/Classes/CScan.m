classdef CScan < handle
    %CScan Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Minimal data (as stored in Rawlog)
        r   % Range data (in m)
        ts  % timesample
        
        % Scan data
%         xy
        mask
        
        % Some extra options
        crop
        
        % Timesample
        delta_ts % Difference with associated data (camera)
        
        % Storage information
        path
        metafile
    end
    
    methods
        function obj = CScan( in1, ts )
            % obj = CScan( S )
            % obj = CScan( r, ts )
            if nargin~=0
                switch nargin
                    case 1
                        S = in1;
                    case 2
                        r = in1; %#ok<*PROP>
                end
                if exist('r','var')
                    % Array constructor
                    if iscell( r )
                        obj(1,numel(ts)) = CScan;
                        [obj.r] = deal( r{:} );
                        [obj.ts] = deal( ts{:} );
                    else
                        obj.r = r;
                        obj.ts = ts;
                    end
                    
                elseif exist( 'S', 'var' )
                    % Array constructor
                    obj(1,numel(S)) = CScan;
                    C = fieldnames( S );
                    C = { C{:} };
                    for field = C
                        field = field{:}; %#ok<FXSET>
                        if isprop(obj,field)
                            [obj.(field)] = deal(S.(field));
                        end
                    end
                else
                    error('Not valid input for this class');
                end
            end
        end
        
        function filterRangeMedian( obj, n )
            % filterRangeMedian( obj, n )
            % Applies median filter to neighbour distances
            % in every group of 2n+1 elements
            head = repmat( obj.r(1),   1, 2*n );
            tail = repmat( obj.r(end), 1, 2*n );
            augmented = [ head, obj.r, tail ];
            dist = augmented - circshift( augmented, [0 1] ); % Distance between neighbours
            rr = repmat( augmented, n, 1 );
            for i=1:n
                rr(i,:) = circshift( rr(i,:), [0 i-n-1] );
            end
            rr = median( rr, 1 );
            obj.r = rr(2*n+1 : end-2*n);
        end
        
        function filterRangeNeighDist( obj, d )
            % filterRangeMedian( obj, n )
            % Applies median filter to neighbour distances
            % in every group of 2n+1 elements

            % r change wrt point in left and right positions, respectively
            dist_l = obj.r - circshift( obj.r, [0 +1] );
            dist_r = obj.r - circshift( obj.r, [0 -1] );
            mask_l = abs(dist_l) > d;
            mask_r = abs(dist_r) > d;
            
            mask = mask_l & mask_r;
            
            obj.r(mask) = 0;
        end
        
        function r = filterRangeBilateral( obj, n, s_sigma, r_dist )
            S = fspecial('gaussian',[1,2*n+1], s_sigma); % Spacial filter (gaussian)
            
            head = repmat( obj.r(1),   1, n );
            tail = repmat( obj.r(end), 1, n );
            augmented = [ head, obj.r, tail ];
            newval = zeros( size( augmented ) );
            
            for i = (1:length(obj.r)) + n
                vals = augmented( i-n : i+n );
                R = abs(augmented(i) - vals) < r_dist;
                filt = S .* R;
                newval(i) = dot(filt,vals) / sum(filt);
            end
            
            if nargout == 1
                r = newval( n+1 : end-n );
            else
                obj.r = newval( n+1 : end-n );
            end
        end
        
        function k = polarCurvature( obj )
            vd  = [1/280 	-4/105 	1/5 	-4/5 	0 	4/5 	-1/5 	4/105 	-1/280]';
            vdd = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]';
            
            n = 4;
            head = repmat( obj.r(1),   1, n );
            tail = repmat( obj.r(end), 1, n );
            augmented = [ head, obj.r, tail ];
            rr = repmat( augmented, 2*n+1, 1 );
            for i=1:n
                rr(i,:) = circshift( rr(i,:), [0 i-n-1] );
            end
            
            r  = obj.r;
            D  = vd'  * rr / obj.resolution_r;
            DD = vdd' * rr / obj.resolution_r^2;
            D  =  D(n+1 : end-n); % First derivative
            DD = DD(n+1 : end-n); % Second derivative
            
            % Formula given in Curvature of graph
            % http://en.wikipedia.org/wiki/Curvature
            k = abs( r.^2 + 2*D.^2 - r.*DD ) ./ ...
                ( r.^2 + D.^2 ).^1.5;
        end
        
        function loadScan( obj )
            % Blender case
            obj.xy = loadBlenderPCD( obj.path );
        end 
    end
    
    methods (Static)
        function scan = average( varargin )
            % CScan.average(Scan_1, Scan_2, ..., Scan_n)
            % Averaged sum of several scans
            r = cell2mat( varargin' );
            r = mean( r, 1 );
            
            scan = CScan( r, [] );
        end
        
    end
    
end

