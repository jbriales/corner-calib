function X = manSumS1( X0, inc_alpha )

N = size( inc_alpha, 2 );
n  = zeros(2,N);
n0 = X0;
for i=1:N
    ca = cos(inc_alpha(i));
    sa = sin(inc_alpha(i));
    R = [ca -sa ; sa ca];
    n(:,i) = R * n0;
end
X = n;
end