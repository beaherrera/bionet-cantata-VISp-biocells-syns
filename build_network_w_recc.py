import os
import shutil
import math
import numpy as np
import pandas as pd
import argparse
import json
from pprint import pprint

from bmtk.builder.networks import NetworkBuilder
from bmtk.builder.bionet.swc_reader import get_swc
from bmtk.builder.bionet import rand_syn_locations


def gaussianLL(src_props, trg_props, w0, sigma):
    src_tuning = src_props["tuning_angle"]
    tar_tuning = trg_props["tuning_angle"]

    delta_tuning = abs(
        abs(abs(180.0 - abs(float(tar_tuning) - float(src_tuning)) % 360.0) - 90.0)
        - 90.0
    )
    weight = w0 * math.exp(-((delta_tuning / sigma) ** 2))

    return weight


def load_models(models_path, model_processing=None, N_default=None):
    models = []
    models_dict = json.load(open(models_path, "r"))
    for layer, layer_dict in models_dict["locations"].items():
        for model_name, models_dict in layer_dict.items():
            ei = models_dict.get("ei", None)
            for model in models_dict["models"]:
                if not model.get("enabled", True):
                    continue

                tuning_angle = np.linspace(
                    0.0, 360.0, N_default or model.get("N", 1), endpoint=False
                )

                models.append(
                    {
                        "N": N_default or model.get("N", 1),
                        "pop_name": model_name,
                        "ei": ei,
                        "model_type": model["model_type"],
                        "model_template": model["model_template"],
                        "model_processing": model_processing
                        or model.get("model_processing", None),
                        "dynamics_params": model["dynamics_params"],
                        "morphology": model.get("morphology", None),
                        "tuning_angle": tuning_angle,
                    }
                )
    return models


def build_network(axon_type, rng_seed=None):
    if rng_seed is not None:
        np.random.seed(rng_seed)

    if axon_type == "stub":
        model_processing = "aibs_perisomatic"
        output_dir = "network_stub_recc"
    elif axon_type in ["axon", "full"]:
        model_processing = "aibs_axon"
        output_dir = "network_axon_recc"
    else:
        raise ValueError("Unknown axon_type")

    models = load_models(
        "model_props/v1_node_models.biophysical_simplified.json",
        model_processing=model_processing,
        N_default=1,
    )

    net = NetworkBuilder("v1")
    for model in models:
        net.add_nodes(**model)

    # Step 2: We want to connect our network. Just like how we have node-types concept we group our connections into
    # "edge-types" that share rules and properties
    cm = net.add_edges(
        source=net.nodes(ei="i"),
        # {
        #     "ei": "i"
        # },  # For synaptic source cells select all inhibitory cells, ei==i, both biophys and point
        target=net.nodes(ei="i"),
        # {
        #     "ei": "i",
        #     "model_type": "biophysical",
        # },  # For synaptic targets select all inhibitory biophys cells
        connection_rule=5,  # All matching source/target pairs will have
        # syn_weight=0.0002,  # synaptic weight
        # target_sections=[
        #     "somatic",
        #     "basal",
        # ],  # Gives the simulator the target sections and
        # distance_range=[0.0, 1e20],  # distances (from soma) when creating connections
        delay=2.0,
        dynamics_params="GABA_InhToInh.json",
        model_template="exp2syn",
    )
    cm.add_properties("syn_weight", rule=0.0002, dtypes=float)
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": [
                "somatic",
                "basal",
            ],
            "distance_range": [0.0, 1e20],
            "return_coords": False,
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    cm = net.add_edges(
        source=net.nodes(ei="i"),  # {"ei": "i"},
        target=net.nodes(ei="e"),  # {"ei": "e", "model_type": "biophysical"},
        connection_rule=lambda trg, src: 5,
        # syn_weight=0.00018,
        # weight_function="wmax",
        # distance_range=[0.0, 50.0],
        # target_sections=["somatic", "basal", "apical"],
        delay=2.0,
        dynamics_params="GABA_InhToExc.json",
        model_template="exp2syn",
    )
    cm.add_properties("syn_weight", rule=0.00018, dtypes=float)
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": ["somatic", "basal", "apical"],
            "distance_range": [0.0, 50.0],
            "return_coords": False,
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    cm = net.add_edges(
        source=net.nodes(ei="e"),  # {"ei": "e"},
        target=net.nodes(ei="i"),  # {"pop_name": "PV1"},
        connection_rule=5,
        # syn_weight=0.00035,
        # weight_function="wmax",
        # distance_range=[0.0, 1e20],
        # target_sections=["somatic", "basal"],
        delay=2.0,
        dynamics_params="AMPA_ExcToInh.json",
        model_template="exp2syn",
    )
    cm.add_properties("syn_weight", rule=0.00035, dtypes=float)
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": ["somatic", "basal"],
            "distance_range": [0.0, 1e20],
            "return_coords": False,
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    cm = net.add_edges(
        source=net.nodes(ei="e"),  # {"ei": "e"},
        target=net.nodes(ei="e"),  # {"pop_name": "Scnn1a"},
        connection_rule=5,
        # syn_weight=6.4e-05,
        # weight_function="gaussianLL",
        # weight_sigma=50.0,
        # distance_range=[30.0, 150.0],
        # target_sections=["basal", "apical"],
        delay=2.0,
        dynamics_params="AMPA_ExcToExc.json",
        model_template="exp2syn",
    )
    cm.add_properties(
        "syn_weight",
        rule=gaussianLL,
        rule_params={
            "w0": 6.4e-05,
            "sigma": 50.0,
        },
        dtypes=float,
    )
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": ["basal", "apical"],
            "distance_range": [30.0, 150.0],
            "return_coords": False,
            "morphology_dir": "components/mechanisms",
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    net.build()
    net.save(output_dir=output_dir)

    virt = NetworkBuilder("virt")
    virt.add_nodes(N=100, ei="e", model_type="virtual")

    cm = virt.add_edges(
        source=virt.nodes(),
        target=net.nodes(ei="e"),
        connection_rule=lambda *_: np.random.randint(10, 20),
        delay=2.0,
        dynamics_params="AMPA_ExcToExc.json",
        model_template="exp2syn",
    )
    cm.add_properties("syn_weight", rule=0.0001, dtypes=float)
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": ["somatic", "basal", "apical"],
            "distance_range": [0.0, 1.0e20],
            "return_coords": True,
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    cm = virt.add_edges(
        source=virt.nodes(),
        target=net.nodes(ei="i"),
        connection_rule=lambda *_: np.random.randint(2, 10),
        delay=2.0,
        dynamics_params="AMPA_ExcToInh.json",
        model_template="exp2syn",
    )
    cm.add_properties("syn_weight", rule=0.0008, dtypes=float)
    cm.add_properties(
        [
            "afferent_section_id",
            "afferent_section_pos",
            "afferent_swc_id",
            "afferent_swc_pos",
            "afferent_section_xcoords",
            "afferent_section_ycoords",
            "afferent_section_zcoords",
        ],
        rule=rand_syn_locations,
        rule_params={
            "sections": ["somatic", "basal", "apical"],
            "distance_range": [0.0, 1.0e20],
            "return_coords": True,
        },
        dtypes=[int, float, int, float, float, float, float],
    )

    virt.build()
    virt.save(output_dir=output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("axon_types", nargs="*", type=str, default=["stub", "axon"])
    parser.add_argument("--rng-seed", nargs="?", type=int, default=100)
    args = parser.parse_args()

    for axon_type in args.axon_types:
        build_network(axon_type=axon_type, rng_seed=args.rng_seed)
