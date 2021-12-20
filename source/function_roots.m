function roots = function_roots(coeffs, monomials)
    f = fimplicit(dot(coeffs, monomials));
    set(gca, 'YTick', [], 'XTick', []);
    pbaspect([1 1 1]);

    xlow  = f.XRange(1) - 0.1;
    xhigh = f.XRange(2) + 0.1;
    ylow  = f.YRange(1) - 0.1;
    yhigh = f.YRange(2) + 0.1;

    axis([xlow, xhigh, ylow, yhigh]);
    fig = frame2im(getframe(gcf));
    close;
    
    %%
    fig = min(fig, [], 3);
    
    [~, i1] = min(sum(fig, 2));
    fig(i1, :) = 255;
    [~, i2] = min(sum(fig, 2));
    fig(i2, :) = 255;
    
    [~, j1] = min(sum(fig, 1));
    fig(:, j1) = 255;
    [~, j2] = min(sum(fig, 1));
    fig(:, j2) = 255;
    
    if i1 > i2, tmp = i1; i1 = i2; i2 = tmp; end
    if j1 > j2, tmp = j1; j1 = j2; j2 = tmp; end

    fig = fig(i1:i2, j1:j2);
    
    %%
    [n, m] = size(fig);
    fig(fig > 100) = 0;
    
    [y, x] = find(fig);
    roots = [x / n * (xhigh - xlow) + xlow, - (y / m * (yhigh - ylow) + ylow)];
end