%% Load dataframe

% Columns 1 through 8
%     'psd_segid'    'BBOX_bx'    'BBOX_by'    'BBOX_bz'    'BBOX_ex'    'BBOX_ey'    'BBOX_ez'    'postsyn_sz' 
% Columns 9 through 16 
%     'postsyn_wt'    'postsyn_x'    'postsyn_y'    'postsyn_z'    'presyn_sz'    'presyn_wt'    'presyn_x'    'presyn_y' 
% Columns 17 through 23
%     'presyn_z'    'size'    'postsyn_segid'    'presyn_segid'    'centroid_x'    'centroid_y'    'centroid_z'

% examples use cases

%eg = (df.presyn_seg(df.postsyn_seg==76181));% find presynaptic partners of a cell
%eg = (df.postsyn_seg(df.presyn_seg==76181));% find postsynaptic partners of a cell


% path to synapses edge list 

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end