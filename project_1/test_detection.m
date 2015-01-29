clear all
dirstruct = dir('*.png');

for i = 1:length(dirstruct),

    % Current test image
    im = imread(dirstruct(i).name);

    % Run detection algorithm
    [x, y, d] = detect_barrel(im);

    % Display results:
    hf = figure(1);
    image(im);
    hold on;
    plot(x, y, 'g+');
    title(sprintf('Barrel distance: %.1f m', d));
    hold off;
    pause;
end