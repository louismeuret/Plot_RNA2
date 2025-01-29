# -*- coding: utf-8 -*- 
class MissingReference(Exception):
    '''
    Raised when a method which needs a reference structure is called on a trajectory without one
    '''

class FrameOutOfBounds(Exception):
    '''
    Raised when trying to access a frame which does not exists
    '''

class EmptyEnsemble(Exception):
    '''
    Raised when trying to compute functions on an empty ensemble
    '''

class BiasPropertiesError(Exception):
    '''
    Raised when a property of a BiasProperties object is needed but not found
    '''

class GuessedTimestep(Warning):
    '''
    Raised when timestep had to be guessed from trajectory and might be wrong
    '''