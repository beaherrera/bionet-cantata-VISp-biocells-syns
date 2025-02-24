## Compiling NEURON mechanisms
NEURON mechanisms are located components/mechanisms/modfiles. Run the following commands to compile them:

```bash
$ cd ../components/mechanisms
$ nrnivmodl modfiles 
$ cd -
```

## Running:
To run a simulation, install bmtk and run the following:
```bash
$ python run_bionet.py config.imulation.aibs_axon.syns.json
```

to run using multiple cores using MPI:
```bash
$ mpirun -np <N> python run_bionet.py config.simulation.aibs_axon.syns.json
```

If successful, will create a *output* directory containing log, spike trains and recorded cell variables. By default only the "internal" spike times are recorded to a spikes.h5 file - as specified in the config.

Parameters such as run-time, dL and dt can be changed updating the config file ```config.imulation.aibs_axon.syns.json```.

## Output directories:


## The Network:
The network files have already been built and stored as SONATA files in the *network_axon/* directory. The bmtk Builder script used to create the network files is *build_network.py*. To rebuild the network files, run:

```
$ python build_network.py
```
or
```bash
$ mpirun -np <N> python build_network.py
```

This will overwrite the existing files in the network directory. Note that there is some randomness in how the network is built, so expect (slightly) different simulation results everytime the network is rebuilt.

## Simulation Parameters
Parameters to run the simulation, including run-time, inputs, recorded variables, and networks are stored in config.json and can modified with a text editor.

## Plotting the voltage traces
Before running the following command to plot the voltage traces and save them to .png files, unzip the membrane_potential.h5 file inside the output directories:
```bash
$ python plot_traces.py\
    --bmtk_json_path config.imulation.aibs_axon.syns.json\
    --canata_output_dir <canata_output_dir>\
    --node_types network_axon/v1_node_types.csv\
    --save-dir figs\
    --save-as bionet-VISp-biocells-syns
```

## Plotting the spike times and mean firing rate
To plot the spike times and mean firing rate, run the following command:
```bash
$ python compare_spikes.py\
    --bmtk_json_path config.imulation.aibs_axon.syns.json\
    --canata-spikes-csv <canata_output_dir/spikes.csv>\
    --node_types network_axon/v1_node_types.csv\
    --save-dir figs\
    --save-as bionet-VISp-biocells-syns
```