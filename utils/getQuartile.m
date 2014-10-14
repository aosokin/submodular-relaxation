function [med, lower, upper] = getQuartile(vec)
%getQuartile finds the median and the two quartiles ov vector vec
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

vec = sort(vec);

if mod(length(vec), 2) == 1
    med = vec((length(vec) - 1) / 2 + 1);
else
    med = (vec(length(vec)  / 2 + 1) + vec(length(vec)  / 2 )) / 2;
end

lower = vec( round( length(vec) / 4 ) );
upper = vec( round( length(vec) / 4 * 3 ) );

end
