function inc_alpha = manSubtractionS1( n, n0 )

N = size( n, 2 );
ort_n0 = [ -n0(2) +n0(1) ]';
ort_n0 = repmat( ort_n0, 1, N );
inc_alpha = asin( dot( ort_n0, n, 1 ) );

end
