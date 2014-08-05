function [ mu, A, muF, AF ] = simulate( F, X0, A0, varargin )
% [ mu, A ] = simulate( F, X0, A0, options )
% [ mu, A, muF, AF ] = simulate( F, X0, A0, varargin )
% Output:
%   mu is the mean value (Montecarlo)
%   A is the covariance matrix (Montecarlo)
%   muF is the mean value (calculated from X0)
%   AF is the covariance matrix (propagated with JF)
% 
% Input:
%   F - function handle of the form F(X), parameters are listed in a
%   column vector
%   X0 - initial data point (column vector)
%   A0 - nX x nX covariance matrix in minimal space of X0 parameter
% 
% Options:
%   N - number of samples for simulation
%   manSum - handle to function that returns array of samples given a initial point
%       and an array with increments in manifold (see below) - default Rn
%   manSubtraction - handle to function that returns array of increments in manifold
%       given a reference point and array of samples - default Rn
%   manMean - handle to function that returns mean parameter
%       given an array of samples - default Rn
%   verbose - set true to show various information
% 
% Further explanations:
% The input functions are of the form:
%   X = manSum(X0,inc_eps)
%   inc_eps = manSubtraction(X,X0)
%   mu = manMean( X )
% where:
%   X - NX x Nsamples array with colum vectos of disturbed parameters
%   X0 - NX column vector with parameters in the space
%   inc_eps - nX x Nsamples array with increments in manifold tangent space
%   mu - NX column vector with average value in corresponding space
% These functions can be generated with auxiliar function
% MonteCarlo.spaceOps (see help)
% Currently supported spaces and manifolds are:
%   - R
%   - S1
%   - SO(3)
% 
% Examples of use:
% 1.- Estimate covariance of rotation computed from R2 and 3 S1 direction vectors
% [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( {'R','R','S1','S1','S1'}, {'SO(3)'} );
% [mu,A] = MonteCarlo.simulate( F_R, X0, A_CO, 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
%                               'AF', A_R, 'N', 1e3 );
% 2.- Estimate covariance of S1 direction vector computed from R2n
% [manSum,manSubtraction,manMean] = MonteCarlo.spaceOps( repmat({'R'},1,size(A_pts,1)), {'S1'} );
% [mu,A] = MonteCarlo.simulate( F_L, X0, A_pts, 'manSum',manSum, 'manSubtraction',manSubtraction, 'manMean',manMean,...
%                               'AF', A_l, 'N', 1e4 );
% 
% Notes:
%   - N should be of the order of 1e4 at least, though it can become slow
%   - N or the order of 1e3 works faster (for online time) but results
%       fluctuates in the range of 10% error
%   - Input functions can be manually programmed for non considered spaces

%% Preallocation of variables
defaultManSum = @(X0,inc)repmat(X0,1,size(inc,2))+inc;
defaultManSubtraction = @(X,X0)X-repmat(X0,1,size(X,2));
defaultManMean = @(X)mean(X,2);

allfields = {'N', 1e4
             'manSum', defaultManSum
             'manSubtraction', defaultManSubtraction
             'manMean', defaultManMean
             'AF', []
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

% Manifold and space dimensions
nX = size( A0, 1 );
NX = size( X0, 1 );

%% Generation of random samples
inc_eps = mvnrnd( zeros(1,nX), A0, N )';
X = manSum( X0, inc_eps );
clear inc_eps

%% Loop solving of problem
Y  = zeros( numel(F(X0)) , N);
% Computation of results is done in Parallel or For form
if ~opt.Par
    for i=1:N
        out = F( X(:,i) ); % Out is the vector of studied parameters
        Y(:,i) = out(:);   % Vectorized result
    end
else
    Y = F(X);
end

%% Statistics computation
mu = manMean( Y );
inc_eps = manSubtraction( Y, mu );
A  = inc_eps * inc_eps' / (N-1);

%% Uncertainty propagation (if AF)
if ~isempty(AF)
    muF = F(X0);
    %% Verbose information output
    if opt.verbose
        if numel(A) == numel(AF)
            A_err = norm(A-AF,'fro') / norm(A,'fro') * 100;
        else
            error('Dimensions do not fit')
        end
        %             [D_KL, rel_norm] = divergence_KL(mu, A, muF, AF);
        fprintf('-----------------------------------------\n')
        fprintf('The mean of F(X) computed by Montecarlo was:\n')
        disp(mu)
        fprintf('The grountruth mean predicted as F(X0) was:\n')
        disp(muF)
        fprintf('The covariance detected by Montecarlo in R resulting from corrupted data:\n')
        disp(A)
        fprintf('The covariance detected by propagation in R resulting from corrupted data:\n')
        disp(AF)
        fprintf('The relative error is: %f/100\n\n', A_err);
        fprintf('-----------------------------------------\n\n')
    end
end

end