#!/usr/bin/env bash

export PATH=$PATH:$(dirname $0)
CELLORGANIZER=/usr/local/tools/cellorganizer/2.5
WORKING_DIRECTORY=`pwd`
MATLAB=/usr/local/ufrb/MATLAB/bin/matlab

INPUT=$1
cp -v $INPUT $(pwd)/model.mat

echo "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
cd ./cellorganizer
setup(true);
cd('$WORKING_DIRECTORY')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify model name and load model
model_name = './model.mat';
load(model_name);

% show shape space by calling the function
f = figure('visible','off');
showShapeSpaceFigure(model);
saveas( f, 'show_shape_space.png', 'png' );

exit;" > script.m

echo "Running the following script in Matlab"
cat script.m

echo $WORKING_DIRECTORY
ln -s $CELLORGANIZER $(pwd)/cellorganizer
$MATLAB -nodesktop -nosplash -r "script;"
