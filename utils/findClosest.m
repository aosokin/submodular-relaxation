function pos = findClosest( vector, values )
%findClosest finds the closest element in vector for each element in values
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  14.10.2014

vector = vector(:);
values = values(:)';

[~, pos] = min(abs(bsxfun(@minus, vector, values)),[], 1);

end

