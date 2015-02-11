classdef CNetSegments < handle
    %CNetSegments Set of segments with common points
    %   This class is a group of 2D segments enriched with connectivity
    %   constraints and related methods
    
    properties
        pts % Array of points in net
        C   % Matrix of connection (indeces of each pair of points)
    end
    properties (SetAccess=private)
        segs
    end
    properties (Dependent)
        Nsegs
        Npts
    end
    
    methods
        function this = CNetSegments( pts, C )
            assert( size(pts,1) == 2, 'Points must be 2xN' );
            assert( size(C,1) == 2, 'Connection matrix must be 2xN' );
            
            % Load inputs
            this.pts = pts;
            this.C   = C;
            
            this.segs = CSegment2D.empty(0,this.Nsegs);
            this.updateSegs;
        end
        
        % Overloading functions
        function varargout = subsref( this, S )
            switch S(1).type
                case '()'
                    if numel(S) == 1
                        % Extract referenced segments
                        varargout{:} = this.segs( S.subs{:} );
                    else
%                         S(1) = [];
%                         varargout(1:numel(S(1).subs{:}) = ...
%                         var = ...
%                             { this.segs( S(1).subs{:} ).(S(2).subs) };
%                         nargout = 2;
                        [varargout{ 1 : numel(S(1).subs{:}) }] = ...
                            deal( this.segs( S(1).subs{:} ).(S(2).subs) );
%                           var = deal( ...
%                             this.segs( S(1).subs{:} ).(S(2).subs) );
%                         [varargout{1:nargout}] = ...
%                             builtin('subsref',this.segs,S);
                    end
                otherwise
                    % Act as usual
%                     builtin('subsref',this,S)
                    [varargout{1:nargout}] = builtin('subsref',this,S);
%                     out = subsref( this, S );
            end
            
        end
        function this = subsasgn( this, S, in )
            switch S(1).type
                case '()'
                    if numel(S) > 1
                        warning('Check case implementation');
                        keyboard;
                    end
                    % Modify referenced segment and its connected
                    k = S.subs{:};
                    ij = this.C(:,k); i=ij(1); j=ij(2);
                    this.pts(:,i) = in.p1;
                    this.pts(:,j) = in.p2;
                    this.updateSegs;
                    
                otherwise
                    % Act as usual
%                     builtin('subsref',this,S)
                    this = builtin('subsasgn',this,S,in);
%                     out = subsref( this, S );
            end
            
        end
        
        function Nsegs = get.Nsegs( this )
%             Nsegs = numel( this.segs );
            Nsegs = size(this.C,2);
        end
        function Npts  = get.Npts( this )
            Npts = size( this.pts, 2 );
        end
        function tags  = tags( this, mask )
%             tags = cellfun(@(x)num2str(x), num2cell(find(mask)),...
%                         'UniformOutput', false);
            tags = cellfun(@(x)num2str(x), {this.segs(mask).tag},...
                'UniformOutput', false);
        end
        
        function updateSegs( this )
            for k=1:this.Nsegs
                this.segs(k) = CSegment2D( ...
                    this.pts(:,this.C(1,k)),...
                    this.pts(:,this.C(2,k)) );
                % Add tag corresponding to index
                this.segs(k).tag = k;
            end
        end
        
        function transform( this, type, T )
            switch type
                case 'euclidean2D'
                    R = T(:,1:2);
                    t = T(:,3);
                    this.pts = R * this.pts + repmat(t,1,this.Npts);
                otherwise
                    warning('Unknown transformation %s\n',type);
            end
            this.updateSegs;
        end
        
        function [segs1,segs2] = getNeighbours( this, k )
            % [segs1,segs2] = getNeighbours( this, k )
            % Input: index k of interest segment
            % Return two arrays of segments:
            % - segs1 are the segments which have in common point i
            % - segs2 are the segments which have in common point j
            % Input segment is excluded from output lists

            [J1,J2] = this.getNeighboursIdxs( k );
            segs1 = this.segs( J1 );
            segs2 = this.segs( J2 );
            
        end
        function [J1,J2] = getNeighboursIdxs( this, k )
            % [segs1,segs2] = getNeighbours( this, k )
            % Input: index k of interest segment
            % Return two arrays of segments:
            % - segs1 are the segments which have in common point i
            % - segs2 are the segments which have in common point j
            % Input segment is excluded from output lists

            ij = this.C(:,k); i=ij(1); j=ij(2);
            [~,J1] = find( this.C == i );
            [~,J2] = find( this.C == j );
            
            % Remove input segment from output list
            J1 = setdiff( J1, k );
            J2 = setdiff( J2, k );
        end
        
        function [J1,J2] = getPointContainers( this, k )
            % [J1,J2] = getPointContainers( this, k )
            % k is the index of searched point
            % J1 are those segments whose p1 point is k
            % J2 are those segments whose p2 point is k
            [~,J1] = find( this.C(1,:) == k );
            [~,J2] = find( this.C(2,:) == k );
        end
        
        function h = plot( this, format, tag )
            this.segs.plot( format, tag );
        end
    end
    
end