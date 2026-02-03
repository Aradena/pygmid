import sys

# setting path
sys.path.append('../src')

from src import pygmid

if __name__ == "__main__":
    pygmid.sweep.run("config_xt018_1v8.cfg")

