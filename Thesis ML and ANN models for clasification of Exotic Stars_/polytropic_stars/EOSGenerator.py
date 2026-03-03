from itertools import product
from EOSPath import EOSPath
import numpy as np


class EOSGenerator:
    def __init__(self, p_saturation, initial_density, segments, gamma_options):
        self.p_saturation = p_saturation
        self.initial_density = initial_density
        self.segments = segments
        self.gamma_options = gamma_options
        self.eos_models = []

    def generate(self):
        segment_rhos = np.exp(np.linspace(np.log(self.initial_density), np.log(self.initial_density * self.segments), self.segments))
        all_gamma_paths = np.array(list(product(self.gamma_options, repeat=self.segments)))

        for i, gammas in enumerate(all_gamma_paths):
            gamma0 = gammas[0]
            K_base = self.p_saturation / (self.initial_density ** gamma0)
            eos = EOSPath(segment_rhos, gammas, K_base=K_base, name=f"EOS_{i}")
            self.eos_models.append(eos)

    def get_models(self):
        return self.eos_models


