function dispTransImg(arr)
    % 2 or 3 dim means colocated
    % 4 or 5 dims means normal DOT data
    
    %% check if it is ok
    num_dims = ndims(arr);
    if num_dims < 3 || num_dims > 5
        error('Number of array dimensions must be -- 3 or 5.\n');
    end
    
    %% create figure, display integrated image, init middle pixel
    f = figure('Position',[384 396 1301 559]);
    ax_trans = subplot(1,2,1); 
    
    imagesc(sum(arr,3)); axis off; axis square; colorbar;
    r = round(size(arr,1)); c = round(size(arr,2));
        
    
    %% initialize the middle pixel as the transient displayed
    figure(f);
    subplot(1,2,2);
    plot(squeeze(arr(r,c,:))); xlabel('Time (ps)');
    title(sprintf('Transient @ pix (%d,%d)',r,c));
    
    %% hovering functionality
    set(f,'WindowButtonMotionFcn', @hoverCallback);
    function hoverCallback(src, evt)
        figure(f);
        subplot(1,2,2);
        mousePoint = get(ax_trans, 'CurrentPoint');
        mouseX = round(mousePoint(1,1));
        mouseY = round(mousePoint(1,2));

        if mouseX >= 1 && mouseX <= size(arr,2) && mouseY >= 1 && mouseY <= size(arr,1)
            plot(squeeze(arr(mouseY,mouseX,:))); xlabel('Time (ps)');
            title(sprintf('Transient @ pix (%d,%d)',mouseY,mouseX));
        end
    end
    
    %% initialize image for scanning

    
    %% click functionality for scanning through the source position
    set(f,'WindowButtonMotionFcn', @hoverCallback);
    function clickCallback(src, evt)
        figure(f);
        subplot(1,2,2);
        mousePoint = get(ax_trans, 'CurrentPoint');
        mouseX = round(mousePoint(1,1));
        mouseY = round(mousePoint(1,2));

        if mouseX >= 1 && mouseX <= size(arr,2) && mouseY >= 1 && mouseY <= size(arr,1)
            plot(squeeze(arr(mouseY,mouseX,:))); xlabel('Time (ps)');
            title(sprintf('Transient @ pix (%d,%d)',mouseY,mouseX));
        end
    end
end