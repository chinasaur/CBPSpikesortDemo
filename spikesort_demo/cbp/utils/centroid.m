function t = centroid(x)

t = (1:length(x))*(x(:)./sum(x));