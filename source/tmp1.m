 syms x y;
func = x ^ 2 + y ^ 2 - 1;
f = fimplicit(func);
set(gca, 'YTick', [], 'XTick', []);
pbaspect([1 1 1]);

xlow  = f.XRange(1) - 0.1;
xhigh = f.XRange(2) + 0.1;
ylow  = f.YRange(1) - 0.1;
yhigh = f.YRange(2) + 0.1;

axis([xlow, xhigh, ylow, yhigh]);
fig = frame2im(getframe(gcf));
close;
imwrite(fig, "image.png");