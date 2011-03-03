def rydberg_to_ev(energy_rydberg):
    return 13.6056923 * energy_rydberg

def hartree_to_ev(energy_rydberg):
    return rydberg_to_ev(energy_rydberg) * 2
