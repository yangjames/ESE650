function classifier_gui
    %% get file names
    file_path = '../dataset/sequences/00/image_2/';
    a = dir(file_path);
    file_names = {a.name};
    valid_files = ~cellfun(@isempty,regexp(file_names,'^\d+\.png\s*$'));
    file_names = file_names(valid_files);
    for i = 1:length(file_names)
        distances{i} = file_names{i}(1:end-4);
    end
    distances = cellfun(@str2num,distances)';
    im_index = 1;
    
    %% set classification values and pixel coordinate values
    color_dict = {'road','not road'};
    %color_dict = {'red','green','blue','brown','grey','black','yellow'};
    color_idx = 1;
    
    coords = cell(length(file_names),1);
    colors = cell(length(file_names),1);
    
    %% establish figure
    f = figure('Name','Gui Classifier','Visible','off','Position',[100, 100, 1080, 810]);
    ax = axes('Units','pixels');
    
    rgb_im = imread([file_path file_names{im_index}]);
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
        'Position', [200 20 85 20],...
        'Callback', @get_roi);
    
    color_select = uicontrol('Style', 'popup',...
       'String', color_dict,...
       'Position', [290 20 80 20],...
       'Callback', @set_color);
   
    save_current = uicontrol('Style', 'edit',...
        'Position', [370 20 250 20],...
        'Callback', @save_to_file,...
        'HorizontalAlignment','left');
        
    f.Visible = 'on';
    
    function get_next_image(source, callbackdata)
        if im_index < length(file_names)
            im_index = im_index+1;
            rgb_im = imread([file_path file_names{im_index}]);
            set(h,'CData',rgb_im);
            title(['File ' num2str(im_index) ' of ' num2str(length(file_names))])
            %drawnow
        end
    end

    function get_previous_image(source,callbackdata)
        if im_index > 1
            im_index = im_index-1;
            rgb_im = imread([file_path file_names{im_index}]);
            set(h,'CData',rgb_im);
            title(['File ' num2str(im_index) ' of ' num2str(length(file_names))])
            %drawnow
        end
    end

    function get_roi(source, callbackdata)
        mask = roipoly;
        indices = find(mask);
        coords{im_index} = [coords{im_index}; indices];
        %current_color = color_dict{color_idx};
        colors{im_index} = [colors{im_index}; repmat(color_idx,length(indices),1)];%vertcat(colors{im_index}, repmat({current_color},length(indices),1));
    end

    function set_color(source,callbackdata)
        color_idx = source.Value;
    end

    function save_to_file(source,callbackdata)
        if ~isempty(source.String)
            im_size = size(rgb_im);
            save([source.String '.mat'],'coords','colors','color_dict','file_names','im_size','distances');
        end
        set(save_current,'String','');
    end
end