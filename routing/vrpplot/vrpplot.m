%% initialize
clear ; close all; clc

%%  Load points and routes data

points_file = '../data/vrp_421_41_1';
routes_file = 'vrp_routes_file';
points = load(points_file);
routes = load(routes_file);

N = points(1, 1);
V = points(1, 2);
capacity = points(1, 3);
demand = points(2:end, 1);
X = points(2:end, 2);
Y = points(2:end, 3);

%% plot the data
%setenv GNUTERM 'x11';

% plot warehouse
x1 = X(1);
y1 = Y(1);
plot(x1, y1, "r*");
hold on

% plot points
plot(X(2:end), Y(2:end), 'b+');
hold on

% plot routes
[m, n] = size(routes);
colors = ["r", "g", "b", "k", "y", "m"];
c = 1;

for i = 1:n-1
  if (routes(i) == 0 && routes(i+1) == 0) 
    c = c + 1;
    if (c > size(colors)) c = 1; endif
    %if (c == 6) break; endif
  endif
  j = routes(i) + 1;
  jp1 = routes(i+1) + 1;
  x = [X(j), X(jp1)];
  y = [Y(j), Y(jp1)];
  plot(x, y, colors(c));
  hold on
endfor

axis([-40 40 -40 40 -40 40]); 
axis square;


