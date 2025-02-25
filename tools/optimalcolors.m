function c = optimalcolors(n)
% build optimal colorscale
% n : how many colors? (int)


c = zeros(4,3);
c(2,1)=1;c(3,2)=1;c(4,3)=1;
c=[c; ~c(2:end,:) ];
if n<=size(c,1)
    c = c(1:n,:);
else
    c = [c; rand([n-size(c,1), 3])];
end