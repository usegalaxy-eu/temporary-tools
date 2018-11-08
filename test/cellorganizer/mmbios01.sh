#!/usr/bin/env bash

export PATH=$PATH:$(dirname $0)
CELLORGANIZER=/usr/local/tools/cellorganizer/2.5
WORKING_DIRECTORY=`pwd`
MATLAB=/usr/local/ufrb/MATLAB/bin/matlab

SEED=$1
NUMBER_OF_SIMPLE_GEOMETRIES=$2
IMAGE_SIZE=$3
DOWNSAMPLE_FACTOR=$4

echo "
% mmbios01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
addpath( genpath([pwd filesep 'cellorganizer']));

seed = $SEED;
try
    state = rng(seed);
catch
    state = RandStream.create( 'mt19937ar', 'seed', seed );
    RandStream.setDefaultStream( state );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE SIMPLE GEOMETRIES
if ~exist( './synthetic_images' )
    mkdir( './synthetic_images' );
end

for i=1:1:$NUMBER_OF_SIMPLE_GEOMETRIES
    options.image_size = $IMAGE_SIZE;

    a = randi([1/4*options.image_size 1/2*options.image_size]);
    b = randi([1/4*options.image_size 1/2*options.image_size]);
    c = randi([1/4*options.image_size 1/2*options.image_size]);
    filename = ['./synthetic_images/synthetic' num2str(sprintf('%04d',i)) '_0.tif'];
    if ~exist( filename )
        disp(['Making image ' filename]);
        img = generate_ellipsoid(a,b,c,options);
        img2tif( img, filename, 'lzw' );
    end

    a = randi([1/2*options.image_size 3/4*options.image_size]);
    b = randi([1/2*options.image_size 3/4*options.image_size]);
    c = randi([1/2*options.image_size 3/4*options.image_size]);
    filename = ['./synthetic_images/synthetic' num2str(sprintf('%04d',i)) '_1.tif'];
    if ~exist( filename )
        disp(['Making image ' filename]);
        img = generate_ellipsoid(a,b,c,options);
        img2tif( img, filename, 'lzw' );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAIN MODEL WITH SIMPLE GEOMETRIES
clear options
tic

% location of TIF files and resolution
dna = './synthetic_images/*_0.tif';
cell = './synthetic_images/*_1.tif';
options.model.resolution = [0.049, 0.049, 0.2000];

%training flag is set to framework because want to train a model for
%nuclear and cell shape
options.train.flag = 'framework';

%train model at full resolution
options.downsampling = [$DOWNSAMPLE_FACTOR $DOWNSAMPLE_FACTOR $DOWNSAMPLE_FACTOR];

%combination of supported model class and type
options.nucleus.name = 'ellipsoids';
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'cylindrical_surface';

%combination of supported model class and type
options.cell.name = 'ellipsoids';
options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';

%output filename
options.model.filename = 'model.mat';

%helper options
options.verbose = true;
options.debug = true;

%CellOrganizer call
img2slml( '3D', dna, cell, [], options );
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exit;" > script.m

echo "Running the following script in Matlab"
cat script.m

echo $WORKING_DIRECTORY
ln -s $CELLORGANIZER $(pwd)/cellorganizer
$MATLAB -nodesktop -nosplash -r "script;"
