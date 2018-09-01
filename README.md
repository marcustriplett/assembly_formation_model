# Computational model of neural assembly formation

MATLAB code for simulating the development of spontaneously active neural assemblies. See corresponding journal article "Emergence of spontaneous assembly activity in developing neural networks without afferent input." M. A. Triplett et al. (2018). To appear in PLOS Computational Biology.

## Usage
Calling `assembly_formation_network.m` in MATLAB will begin the simulation. There are three variables controlling the basic behaviour and output of the model:
* Set `draw_on = 1` to enable replotting of the weight matrix and spontaneous activity raster at regular intervals (given by `pause_interval`).
* Set `pause_on = 1` to also pause whenever the weight matrix and raster are replotted. Press any key to resume the simulation.
* Set `save_on = 1` to write basic summary statistics used in the journal article to file.

Once the `duration` has elapsed, the code will produce several graphs reporting how summary statistics change over time. If `num_trials > 1` the graphs will show the mean plus standard error. The autosimilarity plot will only be calculated for time points exceeding `cluster_sample_time`.
