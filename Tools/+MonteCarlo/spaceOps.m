function [fSum, fSubtraction, fMean] = spaceOps( cellInSpaces, cellOutSpaces )
% [fSum, fSubtraction, fMean] = spaceOps( cellInSpaces, cellOutSpaces )
% Generate a function handle to
%   - manSum
%   - manSubtraction
%   - manMean
% that Montecarlo.simulate needs as input for non-Euclidean spaces
% 
% The functions will be a function of the form:
%   X = manSum(X0,inc_eps)
%   inc_eps = manSubtraction(X,X0)
%   mu = manMean( X )
% where:
%   X - NX x Nsamples array with colum vectos of disturbed parameters
%   X0 - NX column vector with parameters in the space
%   inc_eps - nX x Nsamples array with increments in manifold tangent space
%   mu - NX column vector with average value in corresponding space
% 
% Example of use:
% [...] = MonteCarlo.spaceOps( {'R','R','S1','SO(3)'}, {'SO(3)','SO(3)' )
% Function handle for an input space consisting of R2, S1 and SO(3)
% subspaces concatenated and output space of SO(3) and SO(3)
% 
% [...] = MonteCarlo.spaceSubtraction( repmat({'R'},1,3), {'SO(3)'} )
% Function handle for a space consisting of purely Euclidean R3 and output
% SO(3)
% 
% After defining this function handle this is passed to Montecarlo.simulate
% MonteCarlo.simulate( F, X0, 'manSum', fSum, ...
%                             'manSubtraction',fSubtraction, ...
%                             'manMean', fMean, ... );

%% For input space
cellSpaces = cellInSpaces;
[sum_funs, in_spans_spa, in_spans_man] = deal( cell(1,length(cellSpaces)) );
idx_spa = 0;
idx_man = 0;

% Sum operation
for k=1:length(cellSpaces)
    man = cellSpaces{k};
    switch man.topo
        case 'R'
            Nspa = man.dim;
            Nman = man.dim;
            sum_funs{k} = @MonteCarlo.manSumRn;
        case 'SO'
            if man.dim == 3
                Nspa = 9;
                Nman = 3;
                sum_funs{k} = @MonteCarlo.manSumRot;
            else
                error('%s(%d) not implemented yet',man.topo,man.dim);
            end
        case 'S'
            if man.dim == 1
                Nspa = 2;
                Nman = 1;
                sum_funs{k} = @MonteCarlo.manSumS1;
            else
                error('%s(%d) not implemented yet',man.topo,man.dim);
            end
        otherwise
            error('The specified space does not exist')
    end
    in_spans_spa{k} = idx_spa + (1:Nspa);
    in_spans_man{k} = idx_man + (1:Nman);
    idx_spa = idx_spa + Nspa;
    idx_man = idx_man + Nman;
end
NX = idx_spa;
nX = idx_man;

function X = manSum(X0,inc_eps)
    N  = size(inc_eps,2);
    X = zeros( NX, N );
    
    for i=1:length(sum_funs)
        idx_spa = in_spans_spa{i};
        idx_man = in_spans_man{i};
        fun = sum_funs{i};
        X(idx_spa,:) = fun( X0(idx_spa), inc_eps(idx_man,:) );
    end
end
fSum = @manSum;

%% For output space
cellSpaces = cellOutSpaces;
[sub_funs, mean_funs, out_spans_spa, out_spans_man] = deal( cell(1,length(cellSpaces)) );
idx_spa = 0;
idx_man = 0;

for k=1:length(cellSpaces)
    man = cellSpaces{k};
    switch man.topo
        case 'R'
            Nspa = man.dim;
            Nman = man.dim;
            sub_funs{k} = @MonteCarlo.manSubtractionRn;
            mean_funs{k} = @MonteCarlo.manMeanRn;
        case 'SO'
            if man.dim == 3
                Nspa = 9;
                Nman = 3;
                sub_funs{k} = @MonteCarlo.manSubtractionRot;
                mean_funs{k} = @MonteCarlo.manMeanRot;
            else
                error('%s(%d) not implemented yet',man.topo,man.dim);
            end
        case 'S'
            if man.dim == 1
                Nspa = 2;
                Nman = 1;
                sub_funs{k} = @MonteCarlo.manSubtractionS1;
                mean_funs{k} = @MonteCarlo.manMeanS1;
            else
                error('%s(%d) not implemented yet',man.topo,man.dim);
            end
        otherwise
            error('The specified space does not exist')
    end
    out_spans_spa{k} = idx_spa + (1:Nspa);
    out_spans_man{k} = idx_man + (1:Nman);
    idx_spa = idx_spa + Nspa;
    idx_man = idx_man + Nman;
end
NY = idx_spa-1;
nY = idx_man-1;

function inc_eps = manSubtraction(Y,Y0)
    N  = size(Y,2);
    inc_eps = zeros( nY, N );
    
    for i=1:length(sub_funs)
        idx_spa = out_spans_spa{i};
        idx_man = out_spans_man{i};
        fun = sub_funs{i};
        inc_eps(idx_man,:) = fun( Y(idx_spa,:), Y0(idx_spa) );
    end
end
fSubtraction = @manSubtraction;

function Y0 = manMean(Y)
    Y0 = zeros( NY, 1 );
    
    for i=1:length(sub_funs)
        idx_spa = out_spans_spa{i};
        fun = mean_funs{i};
        Y0(idx_spa,:) = fun( Y(idx_spa,:) );
    end
end
fMean = @manMean;

end