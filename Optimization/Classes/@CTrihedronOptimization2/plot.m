function plot( obj )

% TODO
warning('This function is TODO yet')

% Plotting results and weights
hold on
bar_width = 0.25;

% title(sprintf('Optimization with WEIGHTED=%d',weighted))

all_label( all_label==0 ) = [];
rgb = 'rgb';
subplot(411), hold on
% for i=1:length(all_label)
%     bar(i,err(i), rgb(all_label(i)))
% end
for i=1:3
    bar( find(all_label==i), err(all_label==i), bar_width, rgb(i) )
end

subplot(412), hold on
errW = W*err;
for i=1:3
    bar( find(all_label==i), errW(all_label==i), bar_width, rgb(i) )
end

subplot(413), hold on
% for i=1:length(all_label)
%     bar(i,W(i,i), rgb(all_label(i)))
% end
for i=1:3
    v = diag(W);
    bar( find(all_label==i), v(all_label==i), bar_width, rgb(i) )
end
subplot(414), hold on
axis ij, colormap(gray)
handle = bar3(abs(W)); % abs to set 0 as black
colorbar
for k=1:length(handle)
    zdata = get(handle(k),'ZData');
    set(handle(k),'CData',zdata,'FaceColor','interp');
end

end