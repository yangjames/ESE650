function generate_train_paths
    %% establish figure
    %clear all
    close all
    map_rgb = imread('aerial_color.jpg');
    f = figure('Name','Gui Classifier','Visible','off','Position',[0, 0, 1800, 1000]);
    ax = axes('Units','pixels');
    h = imshow(map_rgb);
    hold on
    path_plot = plot(0,0,'r-*');
    
    % Create buttons
    select_path = uicontrol('Style', 'pushbutton', 'String', 'Select Region',...
        'Position', [200 20 85 20],...
        'Callback', @get_path);
    
    save_data = uicontrol('Style', 'edit',...
        'Position', [370 20 250 20],...
        'Callback', @save_to_file,...
        'HorizontalAlignment','left');

    f.Visible = 'on';                    
    
    data_idx = 1;
    paths = {};
    function get_path(source, callbackdata)
        %{d
        [~,x,y] = roipoly;
        path = [x(1:end-1) y(1:end-1)];
        %path = ginput();
        x_start = round(path(1:end-1,1));
        y_start = round(path(1:end-1,2));
        x_end = path(2:end,1);
        y_end = path(2:end,2);
        x_b = [];
        y_b = [];
        %{d
        for i = 1:size(path,1)-1
            [x_b_i,y_b_i] = getMapCellsFromRay(x_start(i),y_start(i),x_end(i),y_end(i));
            [~,min_idx] = min(abs(x_start(i)-x_b_i));
            if min_idx ~= 1
                x_b_i = flip(x_b_i);
                y_b_i = flip(y_b_i);
            end
            x_b = [x_b;x_b_i(1:end)];
            y_b = [y_b; y_b_i(1:end)];
        end
        %}
        
        %[x_b,y_b] = getMapCellsFromRay(int32(x_start'),int32(y_start'),x_end',y_end');
        set(path_plot,'xdata',x_b,'ydata',y_b)
        paths{data_idx} = [x_b y_b];
        data_idx = data_idx + 1;
    end

    function save_to_file(source,callbackdata)
        if ~isempty(source.String)
            save([source.String '.mat'],'paths');
        end
        set(save_data,'String','');
    end
end