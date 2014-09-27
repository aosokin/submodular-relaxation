% viterbiPottsMex performs Viterbi algorithm to minimize energy on a chain graph
% 
% energy = viterbiPottsMex(unary, costs)
% [energy, labels] = viterbiPottsMex(unary, costs)
% 
% INPUT
%     unary     -   input sequence, N x K double matrix, where N - number of objects, K -
%               number of labels
%               unary potential for object i to be of label j is unary(i, j);
%     costs     -   pairwise coefficients matrix, (N - 1) x 1 vector;
% 
%   Paiwise potential for variable(x_i) to be of label k and variable (x_{i+1}) to be of label l is costs(i)
% 
% OUTPUT
%     energy    -   energy of the best labeling
%     labels    -   the best labeling
%     
% Anton Osokin (firstname.lastname@gmail.com),  22.05.2013