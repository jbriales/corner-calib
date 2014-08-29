function extractStructFields( S )
% extractStructFields( S )
% Copy struct S field values to variables with the same name
vars = fieldnames( S );
for k=1:length(vars)
    field = vars{k};
    assignin('caller',field,S.(field));
    % 'base' for base WS
    % 'caller' for WS in parent function
end

end