function box_demo2
colormap(gray(256));

w = 128; % image width
h = 128; % image height
r = 32; % radius of centre circle

filter = ones(21,21);
filter = filter / sum(filter(:));

img = zeros(h,w);

% draw a circle
for y=1:h
    for x=1:w
        dx = x - w/2;
        dy = y - h/2;
        r2 = dx*dx + dy*dy;

        if(r2 < r*r)
            img(y,x) = 255;
        end
    end
end

subplot(1,2,1);
image(img);
title('original');
axis ([1, w, 1, h], 'square');

for i=1:4
    img = conv2(img, filter, 'same');
end

subplot(1,2,2);
image(img);
title('after 4 iterations of 21x21 box filter');