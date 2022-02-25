function path = pickPath(P)
x = cumsum([0 P(:).'/sum(P(:))]);
x(end) = 1e3*eps + x(end);
path = vec2ind(histcounts(rand,x)');