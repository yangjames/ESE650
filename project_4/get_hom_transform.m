function H = get_hom_transform(R,t)
H = [R t(:);0 0 0 1];