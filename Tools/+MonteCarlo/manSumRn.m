function X = manSumRn( X0, inc )

N = size( inc, 2 );
X = repmat(X0,1,N) + inc;
end