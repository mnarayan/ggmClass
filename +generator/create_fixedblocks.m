function [blockmat] = create_fixedblocks(n,blocklength)
%CREATE_FIXEDBLOCKS Generate block resamples of fixed length
% Inputs
% 	N - Number of Resamples to Generate
%	n - size of the number of available observations
% 	b - size of the bootstrap

	%y = randsample(population,k)
	b2 = blocklength; %ceil(2*n^(1/3));
	%% moving block 
	%N2 = n-b2+1;
	% submat2 = zeros([N2 b2]);
	% for i = 1:N2
	% 	submat2(i,:) = [i:min(i+b2-1,n)];
	% end	
	
	%% fixed block length
	N2 = floor(n/(b2));
	submat2 = zeros([N2 b2]);
	for ii = 1:N2
		submat2(ii,:) = [max(1,(ii-1)*b2+1):min(ii*b2,n)];
	end	
	blockmat = submat2;

	%blklength = ceil(b/(b2-1));	
	%submat = zeros([N blklength*b2]);
	% if(replacement)
	% 	blockresamples = genbootstraps(N,N2,blklength);
	% else
	% 	blockresamples = gensubsamples(N,N2,blklength);
	% end
	% for ii = 1:N
	% 	tmprs = submat2(blockresamples(ii,:),:);
	% 	submat(ii,:) = sort(reshape(tmprs,[1 numel(tmprs)]));
	% end

end