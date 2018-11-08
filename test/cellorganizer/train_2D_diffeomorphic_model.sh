#!/usr/bin/env bash

export PATH=$PATH:$(dirname $0)
CELLORGANIZER=/usr/local/tools/cellorganizer/2.5
WORKING_DIRECTORY=`pwd`
MATLAB=/usr/local/ufrb/MATLAB/bin/matlab

DATASET=$1
NUMBER_OF_IMAGES=$2
DOWNSAMPLE_FACTOR=$3

echo "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
cd ./cellorganizer
setup(true)
cd('$WORKING_DIRECTORY')

options.verbose = true;
options.debug = false;
options.display = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Training component: main call to img2slml

% collection, modify these according to your needs
directory = './cellorganizer/images/HeLa/2D/$DATASET/';

number_of_images = $NUMBER_OF_IMAGES;
files = dir([ directory filesep 'orgdna' filesep 'cell*.tif' ]);

if number_of_images > length(files)
  warning('You are trying to select more images than those available in this dataset. Selecting the maximum available images.')
  number_of_images = length(files)
end

files = files(1:number_of_images);
dna = {};
for i=1:1:length(files)
    dna{length(dna)+1} = [ directory filesep 'orgdna' filesep files(i).name];
end

files = dir([ directory filesep 'orgcell' filesep 'cell*.tif' ]);
files = files(1:number_of_images);
cellm = {};
for i=1:1:length(files)
    cellm{length(cellm)+1} = [ directory filesep 'orgcell' filesep files(i).name];
end

files = dir([ directory filesep 'orgprot' filesep 'cell*.tif' ]);
files = files(1:number_of_images);
protein = {};
for i=1:1:length(files)
    protein{length(protein)+1} = [ directory filesep 'orgprot' filesep files(i).name];
end

files = dir([ directory filesep 'crop' filesep 'cell*.tif' ]);
files = files(1:number_of_images);
masks = {};
for i=1:1:length(files)
    masks{length(masks)+1} = [ directory filesep 'crop' filesep files(i).name];
end
options.masks = masks;

options.model.resolution = [ 0.49, 0.49 ];
options.model.filename = 'model.mat';
options.model.id = num2str(now);
options.model.name = 'cellorganizer-on-galaxy';

%set nuclei and cell model name
options.nucleus.name = 'num2str(now)';
options.cell.model = 'num2str(now)';

%set the dimensionality of the model
dimensionality = '2D';

%documentation
options.documentation.description = 'This model has been trained on CellOrganizer on Galaxy+Bridges';

%set model type
options.nucleus.class = 'framework';
options.cell.class = 'framework';
options.nucleus.type = 'diffeomorphic';
options.cell.type = 'diffeomorphic';
options.train.flag = 'framework';

options.downsampling = [$DOWNSAMPLE_FACTOR,$DOWNSAMPLE_FACTOR];

%train the model
success = img2slml( dimensionality, dna, cellm,[], options );

exit;" > script.m

echo "Running the following script in Matlab"
cat script.m

echo $WORKING_DIRECTORY
ln -s $CELLORGANIZER $(pwd)/cellorganizer
$MATLAB -nodesktop -nosplash -r "script;"
