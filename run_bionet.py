"""Simulates an example network of 450 cell receiving two kinds of exernal input as defined in the configuration file"""

import sys, os
from bmtk.simulator import bionet
from bmtk.simulator.bionet.default_setters.cell_models import set_params_peri
from bmtk.analyzer.spike_trains import plot_raster
from neuron import h


@bionet.model_processing
def aibs_axon(hobj, template_name, dynamics_params):
    if dynamics_params is not None:
        # fix_axon_peri(hobj)
        set_params_peri(hobj, dynamics_params)

    return hobj


def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    net = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=net)
    h.secondorder = 0
    sim.run()

    plot_raster(config_file=config_file, group_by="pop_name")

    bionet.nrn.quit_execution()


if __name__ == "__main__":
    # Find the appropriate config.json file
    config_path = None
    if __file__ != sys.argv[-1]:
        config_path = sys.argv[-1]
        if not os.path.exists(config_path):
            raise AttributeError(
                "configuration file {} does not exist.".format(config_path)
            )
    else:
        for cfg_path in ["config.json", "config.simulation.json"]:
            if os.path.exists(cfg_path):
                config_path = cfg_path
                break
        else:
            raise AttributeError("Could not find configuration json file.")

    run(config_path)
