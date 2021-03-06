function startup

[p, ~, ~] = fileparts(mfilename('fullpath'));

addpath(fullfile(p, 'bgamv120'))
addpath(fullfile(p, 'matlab'))
addpath(fullfile(p, 'mtimesx'))
addpath(fullfile(p, 'baxfun'))
addpath(fullfile(p, 'Multiprod_2009'))
addpath(fullfile(p, 'Vector_algebra_2009'))

addpath(fullfile(p, 'cvx'))
addpath(fullfile(p, 'cvx/structures'))
addpath(fullfile(p, 'cvx/lib'))
addpath(fullfile(p, 'cvx/functions'))
addpath(fullfile(p, 'cvx/commands'))
addpath(fullfile(p, 'cvx/builtins'))