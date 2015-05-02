function E = EssentialMatrixFromFundamentalMatrix(F,K)
E_init=K'*F*K;
[U,~,V] = svd(E_init);
E = U*[1 0 0; 0 1 0; 0 0 0]*V';