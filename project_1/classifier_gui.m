function classifier_gui
    %% get file names
    a = dir;
    file_names = {a.name};
    valid_files = ~cellfun(@isempty,regexp(file_names,'^\d+\.\d+\.png\s*$'));
    file_names = file_names(valid_files);
    im_index = 1;
    
    %% set classification values and pixel coordinate values
    color_dict = {'red','green','blue','brown','grey','black','yellow'};
    color_idx = 1;
    
    coords = cell(length(file_names),1);
    colors = cell(length(file_names),1);
    
    %% establish figure
    f = figure('Visible','off');
    ax = axes('Units','pixels');
    
    rgb_im = imread(file_names{im_index});
    h = imshow(rgb_im);
    title(['File ' num2str(im_index) ' of ' num2str(length(file_names))])
    
    % Create push button
    next = uicontrol('Style', 'pushbutton', 'String', 'Next',...
        'Position', [80 20 50 20],...
        'Callback', @get_next_image);
    
    previous = uicontrol('Style', 'pushbutton', 'String', 'Previous',...
        'Position', [20 20 50 20],...
        'Callback', @get_previous_image);
    
    select_region = uicontrol('Style', 'pushbutton', 'String', 'Select Region',...
        'Position', [150 20 80 20],...
        'Callback', @get_roi);
    
    color_select = uicontrol('Style', 'popup',...
       'String', color_dict,...
       'Position', [20 340 100 50],...
       'Callback', @set_color);
   
    f.Visible = 'on';
    
    function get_next_image(source, callbackdata)
        if im_index < length(file_names)
            im_index = im_index+1;
            rgb_im = imread(file_names{im_index});
            set(h,'CData',rgb_im);
            title(['File ' num2str(im_index) ' of ' num2str(length(file_names))])
            %drawnow
        end
    end

    function get_previous_image(source,callbackdata)
        if im_index > 1
            im_index = im_index-1;
            rgb_im = imread(file_names{im_index});
            set(h,'CData',rgb_im);
            title(['File ' num2str(im_index) ' of ' num2str(length(file_names))])
            %drawnow
        end
    end

    function get_roi(source, callbackdata)
        indices = find(roipoly);
        coords{im_index} = [coords{im_index}; indices]
        current_color = color_dict{color_idx};
        colors{im_index} = vertcat(colors{im_index}, repmat({current_color},length(indices),1));
        colors{im_index}(:);
    end

    function set_color(source,callbackdata)
        color_idx = source.Value;
    end
end