from dataclasses import dataclass, field

from pygmid import sweep

# @pytest.mark.parametrize("config_file", [Path("tests/pygmid/sweep/config_xt018.cfg")])
# def mytest():
#
# from pygmid.sweep.simulator import SpectreSimulator

# @dataclass
# class testSpectreSimulator:
#     # _simulator: SpectreSimulator = field(default_factory=lambda: SpectreSimulator(['blablabla']), repr=False)
#     def __init__(self):
#         self._simulator = SpectreSimulator(['blablabla'])


sweep.run("config_xt018.cfg")
# testSpectreSimulator()
