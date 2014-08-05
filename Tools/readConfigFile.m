function S = readConfigFile( filename )
% S = readConfigFile( filename )
% Returns a struct where fields are different lines in config file
% 
% ConfigFile Example:
%   dir = 1 0
%   Npts = 1081
%   range = 60
%   FOVd = 270.2


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
        pos.comment = strfind( line, '//' );
        if ~isempty( pos.comment )
            line(pos.comment:end) = [];
        end
        
        % Find = symbol
        pos.eq = strfind( line, '=' );
        
        % Find key word
        pos.key = pos.eq - 1;
        while line(pos.key)==' '
            pos.key = pos.key - 1;
        end
        key = line(1:pos.key);
        
        % Find value
        pos.val = pos.eq + 1;
        val = textscan(line(pos.val:end), '%f');
        
        % Store value in struct
        S.(key) = val{:}'; % Using dynamic fieldname
        
        % Deletes previous parameters
        clear pos
       
    end
end
    
    
    