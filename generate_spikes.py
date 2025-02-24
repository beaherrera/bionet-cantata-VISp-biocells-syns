from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

psg = PoissonSpikeGenerator()
psg.add(
    node_ids="network_axon/virt_nodes.h5",
    firing_rate=20.0,
    times=(0.5, 2.5),
    population="virt",
)
psg.to_sonata("inputs/virt_spikes.h5")
