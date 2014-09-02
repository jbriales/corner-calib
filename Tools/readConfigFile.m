function varargout = readConfigFile( varargin )
% S = readConfigFile( filename )
% Returns a struct where fields are different lines in config file
% 
% ConfigFile Example:
%   dir = [1 0]
%   M = [1 1 ; 2 3]
%   str = 'this is a string'
%   FOVd = 270.2


for i=1:length(varargin)
    filename = varargin{i};
    FID = fopen( filename, 'r' );
    
    if FID < 0
        error('File %s was not found\n', filename)
    else
        % Initialize struct
        S = struct;
        
        while(~feof(FID))
            % Read new line
            line = fgetl( FID );
            
            % Remove comments
            pos.comment = [ strfind( line, '#' ), strfind( line, '//' ) ];
            if ~isempty( pos.comment )
                line(pos.comment:end) = [];
            end
            
            % Find = symbol
            pos.eq = strfind( line, '=' );
            if isempty( pos.eq )
                continue
            end
            
            % Find key word
            pos.key = pos.eq - 1;
            while line(pos.key)==' '
                pos.key = pos.key - 1;
            end
            key = line(1:pos.key);
            
            % Find value
            pos.val = pos.eq + 1;
            str = strtrim( line(pos.val:end) );
            
            % Read matrices
            mat = str;
            while mat(end)==';'
                mat = fgetl( FID );
                mat = strtrim( mat );
                str = [str, mat];
            end
            %         if mat(end)~=']'
            %             error('Matrix has to finish with ]')
            %         end
            
            % Store value in struct
            eval( strcat( 'S.(key)=', str,';' ) )
            
            % Delete previous parameters
            clear pos
            
        end
    end
    varargout(i) = {S};
end
    
    