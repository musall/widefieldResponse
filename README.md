# Analysis of cortical wiefield and hippocampal activity
This repo contains a Matlab-based code collection for  analysis of cortical widefield, behavior, and electrophysiological data in the study by Mitlasoczki et al. 2024
Region-specific spreading depolarization drives aberrant post-ictal behavior

[https://doi.org/10.1101/2024.10.12.618012)](https://doi.org/10.1101/2024.10.12.618012)

To get to the results in the paper, download the dataset from the study with the data in the folder 'F:\widefieldResponse_data\'. Otherwhise, change the datapath to the recordings in the function wr_WidefieldCA1stim.m
You can then use the main function 'wr_WidefieldCA1stim' to create the panels for figure 3 in the manuscript. For the code to work, all dependencies in the dependency folder need to be added to the matlab path or need to shift the current working directory to the dependency folder so all required subfunctions are usable.
