function plotEstimatedImg( R_c_w )

title('Approximate image with current estimated rotation in Camera')
hold on
rgb = 'rgb';
for i=1:3
    p  = makeinhomogeneous(R_c_w(:,i));
    op = [ zeros(2,1) p ];
    plot( op(1,:), op(2,:), rgb(i) )
end
axis equal

end