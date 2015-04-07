function map_labeler
    
    %% set classification values and pixel coordinate values
    class_dict = {'car_traversible','human_traversible'};
    class_idx = 1;
    coords = cell(length(class_dict),1);
    
    %% establish figure
    map_rgb = imread('aerial_color.jpg');
    f = figure('Name','Gui Classifier','Visible','off','Position',[50, 50, 1800, 1000]);
    ax = axes('Units','pixels');
    h = imshow(map_rgb);
    
    % Create push button
    select_region = uicontrol('Style', 'pushbutton', 'String', 'Select Region',...
        'Position', [200 20 85 20],...
        'Callback', @get_roi);
    
    class_select = uicontrol('Style', 'popup',...
       'String', class_dict,...
       'Position', [290 20 80 20],...
       'Callback', @set_class);
   
    save_current = uicontrol('Style', 'edit',...
        'Position', [370 20 250 20],...
        'Callback', @save_to_file,...
        'HorizontalAlignment','left');
        
    f.Visible = 'on';
    
    function get_roi(source, callbackdata)
        mask = roipoly;
        indices = find(mask);
        coords{class_index} = [coords{class_index}; indices];
    end

    function set_class(source,callbackdata)
        class_idx = source.Value;
    end

    function save_to_file(source,callbackdata)
        if ~isempty(source.String)
            im_size = size(map_rgb);
            save([source.String '.mat'],'coords','colors','class_dict','im_size');
        end
        set(save_current,'String','');
    end
end