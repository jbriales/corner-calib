function saveConfigFile( filename, varargin )
% saveConfigFile( S )
% Save a struct fields in different lines in config file
% 
% ConfigFile Example:
%   dir = [1 0]
%   M = [1 1 ; 2 3]
%   str = 'this is a string'
%   FOVd = 270.2

FID = fopen( filename, 'w' );

if FID < 0
    error('File %s was not created\n', filename)
else
    for i=1:length(varargin)
        S = varargin{i};

        if isstruct(S)
            vars = fieldnames( S );
            for k=1:length(vars)
                key = vars{k};
                var = S.(key);
                write( key, var );
            end
        else
            key = inputname( 1 + i );
            var = S;
            write( key, var );
        end
    end
end
    
% Nested function to write variable and its name
    function write( key, var )
        if isinteger( var )
            format = '%d\t';
        else
            format = '%.3e\t';
        end
        fprintf(FID, '%s = ', key);
        if isnumeric( var )
            fprintf(FID, '[');
            if isempty( var )
                fprintf(FID,' ');
            elseif isvector( var )
                fprintf(FID, format, var);
            else % Is matrix
                fprintf(FID, format, var(1,:));
                for row=2:size(var,1)
                    fprintf(FID, ';\n\t');
                    fprintf(FID, format, var(row,:));
                end
            end
            fprintf(FID, ']\n');
        elseif ischar( var )
            fprintf(FID, ' ''%s''\n', var);
        end
        fprintf(FID,'\n');
    end
end

