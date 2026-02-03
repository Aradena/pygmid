import path
import sys

# directory reach
directory = path.path(__file__).abspath()

# setting path
sys.path.append(directory.parent.parent)

from src import pygmid

if __name__ == "__main__":
    pygmid.sweep.run("config_xt018_1v8.cfg")

