function stack = stackImages( imgs )

N4 = length( imgs );

[N1, N2, N3] = size( imgs(1).I );

% MK = zeros(N1,N2,N3,N4);
stack = repmat(uint8(0),[N1 N2 N3 N4]);

tiled = [imgs.I];

for k=1:3
    stack(:,:,k,:) = reshape( tiled(:,:,k), N1, N2, N4 );
end
end