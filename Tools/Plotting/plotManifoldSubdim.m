        function h = plotManifoldSubdim( gv, FC, Fx, x0 )
            Ngv = length(gv);
            [w_i,w_j] = meshgrid( gv );
            
            % Linearize code
            inc_eps = [ w_i(:), w_j(:), zeros(Ngv^2,1) ]';
            labels = { 'x', 'y', 'z'
                        'y', 'z', 'x' };
            titles = { 'X-Y', 'Y-Z', 'Z-X' };
            for k=1:3
                inc_eps_ = circshift( inc_eps, k-1 ); % Shift inc_eps rows
                x = Fx(x0,inc_eps_);
                cost = FC(x);
                cost = reshape(cost,Ngv,Ngv);
                
                cost_GT = FC( x0 );
                subplot(1,3,k); hold on;
                xlabel(labels{1,k});
                ylabel(labels{2,k});
                title( titles{k} );
                if 1 % Surf plot
                    surf(w_i,w_j, cost);
                    contour3(w_i,w_j, cost, 'k');
                    plot3(0,0,cost_GT, '.y', 'LineWidth', 3);
                    shading interp;
                    view([90 90]);
                    axis([gv(1) gv(end) gv(1) gv(end)])
                    %                         axis([gv(1) gv(end) gv(1) gv(end)...
                    %                               min(Err{k_sub}(:)) median(Err{k_sub}(:))]);
                else
                    contourf(w_i,w_j, cost, 20);
                    shading interp;
                    view([90 90]);
                    plot(0,0,'.y', 'LineWidth', 3);
                end
            end
            h = []; % TODO
        end