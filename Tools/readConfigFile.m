function S = readConfigFile( filename, section, field )
% S = readConfigFile( filename )
% Returns a struct where fields are different lines in config file
% S = readConfigFile( filename, section )
% Returns a struct with fields at the beginning of .ini and inside certain
% section specified by '[section_name]'
% S = readConfigFile( filename, section, field )
% Returns a variable with name 'field' in section '[section]'
% If no section is given, use '' in section
% 
% ConfigFile Example:
%   dir = [1 0]
%   M = [1 1 ; 2 3]
%   str = 'this is a string'
%   FOVd = 270.2

% for i=1:length(varargin)
%     filename = varargin{i};
    FID = fopen( filename, 'r' );
    
    % Look for section parameter
    
    otherSection = false;
%     if i+1 <= length(varargin)
    if exist('section','var')
        if ~isempty(section)
%         next = varargin{i+1};
%         if isSectionHeader(next) % This is a section specifier
%             section = next;
%             i = i+1; % Fix index for next loop
            inSection = false;
%         end
        end
    else
        section = [];
    end
    if ~exist('field','var')
        field = [];
    end
    
    if FID < 0
        error('File %s was not found\n', filename)
    else
        % Initialize struct
        S = [];
        
        while(~feof(FID))
            % Read new line
            line = fgetl( FID );
            
            % Remove comments
            line = removeComment( line );
            if isempty(line) % Only comment line
                continue
            end
            
            % Check section if necessary and control pipeline
            if ~isempty(section)
                if isSectionHeader( line )
                    if strcmp(section,line)
                        inSection = true;
                        otherSection = false;
                    else
                        otherSection = true;
                        if inSection % Have already reached interest sect
                            break
                        end
                    end
                end
            end
            if otherSection % If in other section, skip line
                continue
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
            
            % In case of specific variable, check its name and skip if not
            if ~isempty(field)
                if ~strcmp(field,key)
                    continue
                end
            end
            
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
            str = removeComment( str );
            %         if mat(end)~=']'
            %             error('Matrix has to finish with ]')
            %         end
            
            % Store value in struct
            if isempty(field) % All fields into struct
                eval( strcat( 'S.(key)=', str,';' ) )
            else % Specific variable is output
                eval( strcat( 'S=', str,';' ) );
                return
            end
            
            % Delete previous parameters
            clear pos
            
        end
    end
    fclose(FID);
    
    if isempty(S)
        warning('Output is empty!');
    end
%     varargout(i) = {S}; %#ok<AGROW>
% end
end

function line = removeComment( line )
pos_comment = [ strfind( line, '#' ), strfind( line, '//' ) ];
if ~isempty( pos_comment )
    line(pos_comment(1):end) = [];
end
line = strtrim( line );
end

function out = isSectionHeader( line )
    out = (line(1)=='[' && line(end)==']');
end