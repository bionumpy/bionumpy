
_bionumpy_random_seed = None

def seed(seed: int):
    global _bionumpy_random_seed
    _bionumpy_random_seed = seed


def default_rng():
    pass