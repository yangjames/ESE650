clear all
dirstruct = dir('*.png');

for i = 1:length(dirstruct),

    % Current test image
    im = imread(dirstruct(i).name);

    % Your computations here!
    [x, y, d] = detect_barrel(im);

    % Display results:
    hf = figure(1);
    image(im);
    hold on;
    plot(x, y, 'g+');
    title(sprintf('Barrel distance: %.1f m', d));

    % You may also want to plot and display other
    % diagnostic information such as the outlines
    % of connected regions, etc.
    hold off;
    pause;
end