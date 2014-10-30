function [ classColor, className ] = createMsrcClassColorTable
% createMsrcClassColorTable is a helper function to work with the MSRC-21 dataset 
%
% See createDataset_MSRC for instructions
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2014

className = cell(0, 1);
classColor = nan(0, 3);
className{end + 1, 1} = 'void';
classColor(end + 1, 1 : 3) = [ 0	0	0];
className{end + 1, 1} = 'building';
classColor(end + 1, 1 : 3) = [ 128	0	0];	 
className{end + 1, 1} = 'grass';
classColor(end + 1, 1 : 3) = [ 0	128	0];	 
className{end + 1, 1} = 'tree';
classColor(end + 1, 1 : 3) = [ 128	128	0	]; 
className{end + 1, 1} = 'cow';
classColor(end + 1, 1 : 3) = [ 0	0	128	 ];
% className{end + 1, 1} = 'horse';   % get gid of horses and mountains as suggested
% classColor(end + 1, 1 : 3) = [ 128	0	128	 ];
className{end + 1, 1} = 'sheep';
classColor(end + 1, 1 : 3) = [ 0	128	128	 ];
className{end + 1, 1} = 'sky';
classColor(end + 1, 1 : 3) = [ 128	128	128	 ];
% className{end + 1, 1} = 'mountain'; % get gid of horses and mountains as suggested
% classColor(end + 1, 1 : 3) = [ 64	0	0	 ];
className{end + 1, 1} = 'aeroplane';
classColor(end + 1, 1 : 3) = [ 192	0	0	 ];
className{end + 1, 1} = 'water';
classColor(end + 1, 1 : 3) = [ 64	128	0	 ];
className{end + 1, 1} = 'face';
classColor(end + 1, 1 : 3) = [ 192	128	0	 ];
className{end + 1, 1} = 'car';
classColor(end + 1, 1 : 3) = [ 64	0	128	 ];
className{end + 1, 1} = 'bicycle';
classColor(end + 1, 1 : 3) = [ 192	0	128	 ];
className{end + 1, 1} = 'flower';
classColor(end + 1, 1 : 3) = [ 64	128	128	 ];
className{end + 1, 1} = 'sign';
classColor(end + 1, 1 : 3) = [ 192	128	128	 ];
className{end + 1, 1} = 'bird';
classColor(end + 1, 1 : 3) = [ 0	64	0	 ];
className{end + 1, 1} = 'book';
classColor(end + 1, 1 : 3) = [ 128	64	0	 ];
className{end + 1, 1} = 'chair';
classColor(end + 1, 1 : 3) = [ 0	192	0	 ];
className{end + 1, 1} = 'road';
classColor(end + 1, 1 : 3) = [ 128	64	128	 ];
className{end + 1, 1} = 'cat';
classColor(end + 1, 1 : 3) = [ 0	192	128	 ];
className{end + 1, 1} = 'dog';
classColor(end + 1, 1 : 3) = [ 128	192	128	 ];
className{end + 1, 1} = 'body';
classColor(end + 1, 1 : 3) = [ 64	64	0	 ];
className{end + 1, 1} = 'boat';
classColor(end + 1, 1 : 3) = [ 192	64	0	 ];

end

