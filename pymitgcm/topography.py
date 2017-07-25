import numpy as np

from . import gcmutils


class Topography:
    def __init__(self, model):
        """ Instantiate a topography object. """

        self.filename = 'topo.bin'
        
        self.model = model
        self.nx, self.ny = model.nx, model.ny


    def save(self, savedir):
        """ Save topography to a binary file. """
        topovar = { self.filename : self.data }
        gcmutils.savegcminput(topovar, savedir)

