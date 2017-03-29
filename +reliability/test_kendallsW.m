X = repmat([1:10]',[1 2]); 
W = reliability.kendallsW(X,[])


X(2:3,:) = 2;
W = reliability.kendallsW(X,[])

X(2:3,2) = 3;
W = reliability.kendallsW(X,[])

X(:,2) = flipud(X(:,1)); 
W = reliability.kendallsW(X,[])
