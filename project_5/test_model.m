function test_model
    %% establish figure
    clear all
    close all
    data = load('test_car_model_3.mat');
    map_rgb = imresize(imread('aerial_color.jpg'),data.scale);
    
    f = figure('Name','Gui Classifier','Visible','off','Position',[0, 0, 1800, 1000]);
    ax = axes('Units','pixels');
    h = imshow(map_rgb);
    hold on
    path_plot = plot(0,0,'r-');
    
    % Create buttons
    select_path = uicontrol('Style', 'pushbutton', 'String', 'Select Region',...
        'Position', [200 20 85 20],...
        'Callback', @get_path);
    f.Visible = 'on';                    
    
    function get_path(source, callbackdata)
        path = ginput();
        start = round(path(1,:));
        goal = round(path(end,:));
        
        ctg = dijkstra_matrix(data.cost_map,goal(2),goal(1));
        [ip1, jp1] = dijkstra_path(ctg, data.cost_map, start(2), start(1));
        set(path_plot,'xdata',jp1,'ydata',ip1)
    end
end