function [Rrel, cell_vanishing_directions, trel, triplets] = codeElqursh( this, segs1, segs2, matches, img_size, tracker1, tracker2 )
            
    % Format specified in HELP
%     segs1_ = [segs1.p1 ; segs1.p2];
%     segs2_ = [segs2.p1 ; segs2.p2];
    
    % Reshape segments to cells? What format for Elqursh!!?? PFFFF!!
    N1 = numel(segs1);
    segs1_ = cell(1,numel(segs1));
    for i=1:N1
        segs1_{i} = [ segs1(i).p1, segs1(i).p2 ];
    end
    N2 = numel(segs2);
    segs2_ = cell(1,numel(segs2));
    for i=1:N2
        segs2_{i} = [ segs2(i).p1, segs2(i).p2 ];
    end
    
    if size(matches,2)~=2
        warning('matches for ali_relpose should be Mx2')
        keyboard
        if size(matches,1)==2
            matches = matches';
        else
            error('Bad input');
        end
    end
    
    size_im = img_size;
    K = this.Cam.K;
    
    params = struct('min_rot_inliers',1,...
                    'min_tran_inliers',1,...
                    'filter_triplets',true,...
                    'min_candidate_triplets',1,...
                    'rot_thresh',1); % Default 2 is too big (get outliers!?)
            
	if 0
    [Rrel,trel,tripletpairs,inliers,err,...
        p1,p2,inliersT,seg_pairs1,seg_pairs2,...
        cell_vanishing_directions] = ...
        ali_relpose(segs1_,segs2_,matches,size_im,K,params); %#ok<ASGLU>
    else
        [Rrel,trel,tripletpairs,inliers,err,...
        p1,p2,inliersT,seg_pairs1,seg_pairs2,...
        cell_vanishing_directions] = ...
        ali_relpose(segs1_,segs2_,matches,size_im,K,params); %#ok<ASGLU>
    end
    triplets = tripletpairs(1:3,inliers);
    
    % Visually check results
    if 0
        hF = figure('Name','Visually check Elqursh result');
        hI = imshow( tracker1.img ); hold on;
        set(hF,'Position', [771   202   808   570] );
        freezeColors;
        for i=1:numel(inliers)
            %         tags = tracker1.segs.tags( tracker1.maskSegs );
            title( num2str(err(i)) );
            idxs = tripletpairs(1:3,inliers(i));
            hsegs = segs1(idxs).plot('r');
            delete(hsegs);
        end
        hI = imshowpair( tracker1.img, tracker2.img, 'montage' );
    end
    
end