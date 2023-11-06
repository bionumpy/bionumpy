import dataclasses
import hypothesis.strategies as st


@dataclasses.dataclass
class BaseStrategy:
    strategy: callable

    def sample(self, n):
        pass


String = BaseStrategy(st.text)



def oo_strategy(cls):
    return cls
