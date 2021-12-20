%%
clear; clc;
addpaths;

pop.F = [2 0 -1; 0 2 -2];
pop.n = 2;
pop.k = 1;
pop.d = 1;

% pop.G{1} = [2 0 0 0 0 -1; 0 2 0 0 0 -1; 0 0 2 0 0 -1; 0 0 0 2 0 -1; 0 0 0 0 2 -1; 0 0 0 0 0  1];
% pop.G{2} = [2 0 0 0 0  1; 0 2 0 0 0  1; 0 0 2 0 0  1; 0 0 0 2 0  1; 0 0 0 0 2  1; 0 0 0 0 0 -1];
pop.G{1} = [2 0 -1; 0 2 -1; 0 0  1];
pop.G{2} = [2 0  1; 0 2  1; 0 0 -1];
pop.G{3} = [0 1  1];
 
pop.I = {1:2};
pop.J = {1:size(pop.G, 2)};

sdp = gendata2(pop,'SBSOS');
sol = csol(sdp,'sdpt3');
psol = postproc(pop,sdp,sol);
param = psol.YY{1};
disp(param);

y = zeros(2);
y(tril(true(2))) = param(4:end);
y = y + triu(y', 1);
disp(y);


% sdp = gendata2(pop, 'SBSOS', 'sedumi');
% sol = csol(sdp, 'sedumi');
% psol = postproc(pop, sdp, sol);
% disp(psol.YY{1});