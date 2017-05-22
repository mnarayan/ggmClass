% p = 15;
% m = 100;
% n = 10;
% % R = ones(1,p); rho = .5; R(2:p) = rho.^[2:p];
% % popSigi = toeplitz(R);
% popSigi = eye(p);
% F = randn(p,p); [V D] = eig(F+F'); [sorteig sortidx] = sort(diag(D),'descend');
% popSig =  V(:,sortidx(1:2))*diag(sorteig(1:2))*V(:,sortidx(1:2))';
% A = zeros(m,p,n);
% %
% % popSigi = popSigi/trace(popSigi);
% % popSig = popSig/trace(popSig);
%
% for cc=1:n
%     A(:,:,cc) = randn(m,p)*sqrtm(popSigi+popSig);
% end

% Check that code runs without error and that hatSigma_i are nearly empty

[hatSigma hatSigma_i] = covariance.common_lowrank_decomposition(A,{});
figure(1)
subplot(1,3,1); imagesc(popSig); colorbar; axis image equal;
subplot(1,3,2); imagesc(popSigi); colorbar; axis image equal;
subplot(1,3,3); imagesc(hatSigma); colorbar; axis image equal;
figure(2)
subplot(2,2,1); imagesc(corr(A(:,:,1))); colorbar; axis image equal;
subplot(2,2,2); imagesc(corr(A(:,:,2))); colorbar; axis image equal; 
subplot(2,2,3); imagesc(hatSigma_i(:,:,1)); colorbar; axis image equal;
title(sprintf('Frob.: %.2f',norm(popSigi(:,:,1)-hatSigma_i(:,:,1))))
subplot(2,2,4); imagesc(hatSigma_i(:,:,2)); colorbar; axis image equal; 
title(sprintf('Frob.: %.2f',norm(popSigi(:,:,1)-hatSigma_i(:,:,2))))
