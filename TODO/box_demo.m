function box_demo

r = 100; % graph x-axis range

% filter
filter = ones(1,19);
%filter = rand(1,19); % random filter!
filter = filter / sum(filter); % normalise

% filter inserted into the graph x=axis
a = zeros(1,r);
a(41:(41+length(filter)-1)) = filter;

subplot(3,2,1)
plot(a);
title('filter');

for i=1:5
    a = conv(a, filter, 'same');
    a = a(1:r);

    % fit a gaussian
    n = sum(a); % number of elements
    m = sum(a .* [1:r]) / n; % mean
    s = sqrt(sum(a .* ([1:r] - m).^2) / n); % standard deviation

    x = [1:r];
    scale = max(a);
    data = scale*exp(-(x-m).^2/(2*s*s)); 

    % plot the graph
    subplot(3,2,i+1);
    plot(x,a,x,data);
    txt = sprintf('iteration %d', i);
    title(txt);
end