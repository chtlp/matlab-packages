% h: handle of the figure (e.g., gcf)
% a: resize x
% b: resize y
function resize_figure(handle, a, b)

A = get(handle, 'Position');
x = A(1); y = A(2); w = A(3); h = A(4);

xc = x + w/2;
yc = y + h/2;

w = w * a;
h = h * b;

x = xc - w/2;
y = yc - h/2;

set(handle, 'Position', [x, y, w, h]);