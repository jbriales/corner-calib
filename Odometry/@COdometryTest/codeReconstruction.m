function [s, PTS] = codeReconstruction( this,...
    im_pts,...
    Rrel, trel, K, imsize_,...
    reconThres )

    % TODO: Set in object
    imsize = repmat( imsize_, 3,1 )';

    % Compute 3D points for pair 1-2
    P{1} = K * eye(3)  * [ eye(3) zeros(3,1) ]; % Camera matrix ref Cam1
    P{2} = K * Rrel{1}'* [ eye(3) -trel{1}  ];  % Camera matrix for non-ref Cam2
    
    Npts = size(im_pts{1},2);
    X = zeros(4, Npts);
    for i=1:Npts
        x = makeinhomogeneous( [ im_pts{1}(:,i), im_pts{2}(:,i) ] );
%         X(:,i) = vgg_X_from_xP_lin(x,P,imsize); % Linear algorithm
        X(:,i) = vgg_X_from_xP_nonlin(x,P,imsize); % Non-linear algorithm
    end
% %     % Plot reconstruction
% %     mask_finite = abs(X(4,:)) > 0.01;
% %     pts3D{1} = makeinhomogeneous( X(:,mask_finite) );
% %     % Transform points to World coordinate system
% %     pts3D_W = this.Cam.R * pts3D{1} + repmat(this.Cam.t,1,size(pts3D{1},2));
% %     subplot( this.hF_scene )
% %     hold on, plot3(pts3D_W(1,:),pts3D_W(2,:),pts3D_W(3,:),'or')
% %     subplot( this.hF_image )

    % Compute 3D points for pair 2-3
    P{1} = K * eye(3)  * [ eye(3) zeros(3,1) ]; % Camera matrix ref Cam1
    P{2} = K * Rrel{2}'* [ eye(3) -trel{2}  ];  % Camera matrix for non-ref Cam2
    
    Y = zeros(4, Npts);
    for i=1:Npts
        y = makeinhomogeneous( [ im_pts{2}(:,i), im_pts{3}(:,i) ] );
        Y(:,i) = vgg_X_from_xP_lin(y,P,imsize);
    end
    
    % Remove very far points (more uncertain) wrt average values
%     Q = quantile(sqrt(sum(makeinhomogeneous(X).^2,1)),[.25 .75]);
    % Use both list of points
    Q = quantile( [sqrt(sum(makeinhomogeneous(X).^2,1)),...
                   sqrt(sum(makeinhomogeneous(Y).^2,1))],...
                   [.25 .75] );
               % Parameter 5 is hard coded, BE CAREFUL
    mask_finite = (sqrt(sum(makeinhomogeneous(X).^2,1)) < Q(2)+5*(Q(2)-Q(1))) & ...
                  (sqrt(sum(makeinhomogeneous(Y).^2,1)) < Q(2)+5*(Q(2)-Q(1)));
%     mask_finite = (sqrt(sum(makeinhomogeneous(X).^2,1)) < 100) | ...
%                   (sqrt(sum(makeinhomogeneous(Y).^2,1)) < 100);
%     mask_finite = (abs(X(4,:)) > 0.01) & (abs(Y(4,:)) > 0.01);
%     pts3D{1} = X(:,mask_finite);
%     pts3D{2} = Y(:,mask_finite);
    pts3D{1} = makeinhomogeneous( X(:,mask_finite) );
    pts3D{2} = makeinhomogeneous( Y(:,mask_finite) );
    % Correct Z+ sign to avoid points behind camera
    for k=1:2
        pts3D{k} = repmat(sign(pts3D{k}(3,:)),3,1).* pts3D{k};
    end
%     pts3D{1} = makeinhomogeneous( X );
%     pts3D{2} = makeinhomogeneous( Y );
    
    % Find scale factor s = ||trel{2}|| / ||trel{1}||
    A = pts3D{1} - repmat(trel{1},1,size(pts3D{1},2));
    B = Rrel{1} * pts3D{2};
    s = sum( dot(A,B,1) ) / sum( dot(B,B,1) );
    
    % Return triangulated points for second pair (2-3), wrt Camera 2
    PTS = pts3D{2};
end