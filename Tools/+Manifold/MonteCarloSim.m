function [out, YY, inc] = MonteCarloSim( F, X0, varargin )
% out = simulate( F, X0, options )
% Output:
%   mu is the mean value (Montecarlo)
%   A is the covariance matrix (Montecarlo)
% 
% Input:
%   F - function handle of the form F(X) where X is a manifold
%   X0 - initial data point in a certain manifold, with cov data
% 
% Options:
%   Ref - reference output manifold to compare with
%   N - number of samples for simulation
%   verbose - set true to show various information
% 
% Notes:
%   - N should be of the order of 1e4 at least, though it can become slow
%   - N or the order of 1e3 works faster (for online time) but results
%       fluctuates in the range of 10% error
%   - This function is based on the use of manifold classes with internal
%   operations defined

%% Preallocation of variables

allfields = {'N', 1e4
             'Ref', []
             'Par', false
             'verbose', true
            };
opt = setVarargin(allfields, varargin{:});
% Copy struct values to variables with field names
vars = fieldnames( opt );
for i=1:length(vars)
    field = vars{i};
    eval([ field,'= opt.',field,';' ]);
end
clear i allfields field vars

if verbose
    fprintf('Montecarlo analysis:\n')
    fprintf('--------------------\n')
    fprintf('Nsamples = %d\n',N)
end

%% Generation of random samples
% if isa(X0,'Manifold.Dyn')
%     keyboard
%     if all( cellfun(@(X)isa(X,'Manifold.S2'), X0.vars ) )
%         X = mvnrnd( X0.X', X0.A_X, N )';
%         CX = mat2cell(X,[3 3 3],ones(1,N));
%         clear X
%         for ii=1:N
%             X(ii) = Manifold.Dyn( ...
%                 Manifold.S2(CX{1,ii}),...
%                 Manifold.S2(CX{2,ii}),...
%                 Manifold.S2(CX{3,ii}) );
%         end
%     else
%         inc_eps = mvnrnd( zeros(1,X0.dim), X0.A_x, N )';
%         [~,X] = X0 + inc_eps; %#ok<RHSFN>
%         clear inc_eps
%     end
% else
    inc_eps = mvnrnd( zeros(1,X0.dim), X0.A_x, N )';
    [~,X] = X0 + inc_eps; %#ok<RHSFN>
    clear inc_eps
% end

%% Loop solving of problem
Y = Ref.empty(0,N);
% Computation of results is done in Parallel or For form
if ~opt.Par
    for i=1:N
        Y(i) = F( X(i) ); % Out is the output object
    end
else
    Y = F(X);
end

%% Statistics computation
mu = Ref.mean( Y );
for i=1:length(Y)
    inc_eps{i} = Y(i) - mu;
end
inc_eps = cell2mat( inc_eps );
A  = inc_eps * inc_eps' / (N-1);
mu.setMinimalCov( A );
out = mu;

%% Uncertainty propagation (if AF)
if ~isempty(Ref)
    %% Verbose information output
    if opt.verbose
        A_err = norm(Ref.A_x-out.A_x,'fro') / norm(out.A_x,'fro') * 100;
        %             [D_KL, rel_norm] = divergence_KL(mu, A, muF, AF);
        fprintf('-----------------------------------------\n')
        fprintf('The mean of F(X) computed by Montecarlo was:\n')
        disp(out.X)
        fprintf('The grountruth mean predicted as F(X0) was:\n')
        disp(Ref.X)
        fprintf('The covariance detected by Montecarlo in R resulting from corrupted data:\n')
        disp(out.A_x)
        fprintf('The covariance detected by propagation in R resulting from corrupted data:\n')
        disp(Ref.A_x)
        fprintf('The relative Fro error is: %f/100\n\n', A_err);
        fprintf('-----------------------------------------\n\n')
    end
end

if nargout >= 2
    % Output also all computed output
    YY  = Y;
    inc = inc_eps;
end
end