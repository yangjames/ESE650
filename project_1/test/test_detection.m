clear all
dirstruct = dir('*.png');
xf = cell(length(dirstruct),1);
df = xf;
yf = df;

for i = 1:length(dirstruct),

    % Current test image
    im = imread(dirstruct(i).name);

    % Run detection algorithm
    [x, y, d] = detect_barrel(im);

    xf{i} = x;
    yf{i} = y;
    df{i} = d;

    %{
    % Display results:
    hf = figure(1);
    image(im);
    hold on;
    plot(x, y, 'g+');
    title(sprintf('Barrel distance: %.1f m', d));
    hold off;
    pause;
    %}
end
for i = 1:length(xf)
    fprintf('\niteration: %d\n',i)
    disp(xf{i})
    disp(yf{i})
    disp(df{i})
end