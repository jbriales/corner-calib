classdef CGT < handle
    % CGT Class with persistent variables to store GT and access it from
    % anywhere (without passing as an argument)

properties ( GetAccess = 'public', SetAccess = 'private' )
    % Empty
end

methods ( Static )
    function R = R( R )
        persistent mem_R;
        if nargin == 1
            mem_R = R;
        elseif nargin == 0
            R = mem_R;
        else
            disp('Something went wrong');
        end
    end
    
    function t = t( t )
        persistent mem_t;
        if nargin == 1
            mem_t = t;
        elseif nargin == 0
            t = mem_t;
        else
            disp('Something went wrong');
        end
    end
    
    function triplets = triplets( triplets )
        % triplets = triplets( triplets )
        % Storage of all feasible inlier configuration triplets
        % Intended for Rotation computation
        persistent mem_triplets;
        if nargin == 1
            mem_triplets = triplets;
        elseif nargin == 0
            triplets = mem_triplets;
        else
            disp('Something went wrong');
        end
    end
    
    function pairs = pairs( pairs )
        % pairs = pairs( pairs )
        % Storage of all feasible inlier configuration pairs
        % Intended for Translation computation
        persistent mem_pairs;
        if nargin == 1
            mem_pairs = pairs;
        elseif nargin == 0
            pairs = mem_pairs;
        else
            disp('Something went wrong');
        end
    end
    
    function triplets_ = triplets_( triplets_ )
        % All combinations (triplets) feasible for current iteration segment tags,
        % both inliers and outliers included
        persistent mem_triplets_;
        if nargin == 1
            mem_triplets_ = triplets_;
        elseif nargin == 0
            triplets_ = mem_triplets_;
        else
            disp('Something went wrong');
        end
    end
    
    function pairs_ = pairs_( pairs_ )
        % All combinations (pairs) feasible for current iteration segment tags,
        % both inliers and outliers included
        persistent mem_pairs_;
        if nargin == 1
            mem_pairs_ = pairs_;
        elseif nargin == 0
            pairs_ = mem_pairs_;
        else
            disp('Something went wrong');
        end
    end
    
    function segs3D = segs3D( segs3D )
        persistent mem_segs3D;
        if nargin == 1
            mem_segs3D = segs3D;
        elseif nargin == 0
            segs3D = mem_segs3D;
        else
            disp('Something went wrong');
        end
    end
    
    function out = any( name, in )
        eval(['persistent ',name]); % persistent name
        if nargin == 2
            eval([name,' = in;']);
        elseif nargin == 1
            out = eval(name);
        else
            disp('Something went wrong');
        end
    end
    
    function out = other( id, val )
        persistent others;
        if nargin == 2
            others{id} = val;
        elseif nargin == 1
            out = others{id};
        else
            disp('Something went wrong');
        end
    end
    
    % Useful methods to check and compare values to GT
    function [missed, idx_j]  = missedInliers( triplets )
        assert(~isempty(CGT.triplets),'CGT.triplets is empty');
        assert(size(triplets,1)==3,'triplets input must be 3xN');
        % Intersect all possible CGT triplets with currently feasible ones
        triplets_gt = intersect( CGT.triplets', CGT.triplets_','rows' )';
        % Get difference of triplets
        [missed, idx_j] = setdiff(triplets_gt',triplets','rows');
        missed = missed';
    end
       
    function [slipped, idx_j] = slippedOutliers( triplets )
        assert(~isempty(CGT.triplets),'CGT.triplets is empty');
        assert(size(triplets,1)==3,'triplets input must be 3xN');
        % Intersect all possible CGT triplets with currently feasible ones
        triplets_gt = intersect( CGT.triplets', CGT.triplets_','rows' )';
        % Get difference of triplets
        [slipped, idx_j] = setdiff(triplets',triplets_gt','rows');
        slipped = slipped'; % Transpose to get column
    end
    
    function [deg, idx_j] = degenerateInliers( triplets )
        % TODO
    end
    
        function [missed, idx_j]  = missedInliersT( pairs )
        assert(~isempty(CGT.pairs),'CGT.pairs is empty');
        assert(size(pairs,1)==2,'triplets input must be 3xN');
        % Intersect all possible CGT pairs with currently feasible ones
        pairs_gt = intersect( CGT.pairs', CGT.pairs_','rows' )';
        % Get difference of pairs
        [missed, idx_j] = setdiff(pairs_gt',pairs','rows');
        missed = missed';
    end
       
    function [slipped, idx_j] = slippedOutliersT( pairs )
        assert(~isempty(CGT.pairs),'CGT.pairs is empty');
        assert(size(pairs,1)==2,'pairs input must be 2xN');
        % Intersect all possible CGT pairs with currently feasible ones
        pairs_gt = intersect( CGT.pairs', CGT.pairs_','rows' )';
        % Get difference of pairs
        [slipped, idx_j] = setdiff(pairs',pairs_gt','rows');
        slipped = slipped'; % Transpose to get column
    end
end

end