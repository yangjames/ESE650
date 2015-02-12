function coord = pixel_to_world(pixel)
f = 3.7;
w = 4.8;
h = 3.6;
s = 0.015;
coord = [repmat(f,1,size(pixel,2));...
        -(pixel(1,:)-320/2)*s; ...
        -(pixel(2,:)-240/2)*s];
                    