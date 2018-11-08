#!/usr/bin/env bash

export PATH=$PATH:$(dirname $0)
CELLORGANIZER=/usr/local/tools/cellorganizer/2.5
WORKING_DIRECTORY=`pwd`
MATLAB=/usr/local/ufrb/MATLAB/bin/matlab

SEED=$1
NUMBER_OF_IMAGES=1
COMPRESSION=$2

echo "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
cd ./cellorganizer
setup();
cd('$WORKING_DIRECTORY')

options.targetDirectory = pwd;
options.prefix = 'output';
options.debug = false;
options.verbose = true;
options.display = false;

options.seed=$SEED;
try
    state = rng( options.seed );
catch err
    state = rand( 'seed', options.seed ); %#ok<RAND>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.numberOfSynthesizedImages=$NUMBER_OF_IMAGES;
options.compression = '$COMPRESSION';
options.targetDirectory = pwd;
options.debug = false;
options.verbose = true;
options.display = false;

tic
answer = slml2img( {'./cellorganizer/models/2D/nucleolus.mat', ...
  './cellorganizer/models/2D/endosome.mat'}, options );
toc

files = dir( 'output*.tif' );

for i=1:1:length(files)
  file = files(i).name;
  img = tif2img( file );
  img2 = reshape( img, size(img, 1 ), [] );
  img2 = uint8(img2);
  imwrite( img2, ['output' num2str(i) '.png'] );
end

exit;" > script.m

ln -s $CELLORGANIZER $(pwd)/cellorganizer
$MATLAB -nodesktop -nosplash -r "script;"
