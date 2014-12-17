% Get all currently existing variables
input_vars = who;
% Remove this from list
input_vars( strcmp(input_vars,'this') ) = [];

% Assign each variable to object property if:
% - Input variable exists
% - Object property named as input variable exists
thisprops = fieldnames(this);
n = numel(input_vars);
for idx = 1:n
    var = input_vars{idx};
    if any( strcmp(var, thisprops) )
        this.(var) = eval(var);
    end
end