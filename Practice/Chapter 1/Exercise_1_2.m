%% 
% *Convex hull*
% 
% Prove that conv(C) = {all convex combinations of points in C}.

% Generate some random points
C = rand(5, 2);

% Find the convex hull using convhull
K = convhull(C);

% Generate some random coefficients for a convex combination
alpha = rand(5, 1);
alpha = alpha / sum(alpha); % ensure the coefficients sum to 1

% Compute the convex combination of points in C
x = C.' * alpha;

% Check if x lies within the convex hull K
in_hull = inpolygon(x(1), x(2), C(K, 1), C(K, 2));

% Display the results
disp(['Convex hull of C: ' num2str(K.')]);
disp(['Convex combination of points in C: ' num2str(x.')]);
disp(['Is the convex combination in the convex hull? ' num2str(in_hull)]);