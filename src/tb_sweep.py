from pygmid import sweep

# @pytest.mark.parametrize("config_file", [Path("tests/pygmid/sweep/config_xt018.cfg")])
def mytest():
    sweep.run("./tests/pygmid/sweep/config_xt018.cfg")


