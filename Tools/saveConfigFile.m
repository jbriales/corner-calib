function saveConfigFile( filename, S )
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
    vars = fieldnames( S );
    for k=1:length(vars)
        key = vars{k};
        fprintf(FID, '%s = ', key);
        var = S.(key);
        if isnumeric( var )
            fprintf(FID, '[');
            if isempty( var )
                fprintf(FID,' ');
            elseif isvector( var )
                fprintf(FID, '%.3e ', var);
            else % Is matrix
                fprintf(FID, '%.3e ', var(1,:));
                for row=2:size(var,1)
                    fprintf(FID, ';');
                    fprintf(FID, '%.3e ', var(row,:));
                end
            end
            fprintf(FID, ']\n');
        elseif ischar( var )
            fprintf(FID, ' ''%s''\n', var);
        end
    end
end
    
end