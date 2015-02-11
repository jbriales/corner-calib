function h = plotHomLineWin( l, c )
% plotHomLineWin( l, c )
% Plot an homogeneous line (3x1) in the current window (inside current axis
% configuration)

if ~exist('c','var')
    c = 'k-';
end

if numel(l) ~= 3
    error('PLOTHOMLINEWIN: %d elements in homogeneous line, should be 3', numel(l))
end
l = l(:); % Make column vector

border = [ 1 1 0 0
           0 0 1 1
            -axis  ];
        
if 1
    p = makeinhomogeneous( skew(l) * border );
    in = insideAxis( p );
    
    p = p(:,in);
    
    h = plot(p(1,:),p(2,:),c);
else % Change size_img
    n = l(1:2);
    v = [ -n(2) n(1) ]';
    p = -l(3)*n;
    
    p = zeros(2,2);
    p(:,1) = getBorderIntersection( p, +v, size_img );
    p(:,2) = getBorderIntersection( p, -v, size_img );
end

end

function inside = insideAxis( p )

lim = axis;
lim = lim + 0.0001 * [ -1 +1 -1 +1 ];

inside = p(1,:) >= lim(1) & p(1,:) <= lim(2) & p(2,:) >= lim(3) & p(2,:) <= lim(4);

end

function q = getBorderIntersection( p, v, size_img )
% Get intersection with image border of line going from p in direction v

n = [ -v(2) v(1) ]';
l = [ n' -n'*p ]';
ax = [ 1 size_img(2) 1 size_img(1) ];
borders = [ 1 1 0 0
    0 0 1 1
    -ax   ]; % Homogeneous lines for image borders
q = makeinhomogeneous( skew(l) * borders );
pq = q - repmat(p,1,4);
% Find closest q of the 4 results in positive v direction
q  =  q(:, v'*pq > 0);
pq = pq(:, v'*pq > 0 );  % Remove back intersections
[~,Iq] = min( v'*pq ); % Take closest intersection
q = q(:,Iq);
clear Iq
end
        