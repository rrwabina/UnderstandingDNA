import numpy as np 
import pandas as pd 

class FSM:
    '''
    Finite State Machine is an abstract machine which consists of a set of states 
    (including the initial state and one or more end states), 
    a set of input events, a set of output events, and a state transition function. 
    '''
    def __init__(self):
        self.handlers = {}
        self.startState = None
        self.endStates = []

    def add_state(self, name, handler, end_state = 0):
        name = name.upper()
        self.handlers[name] = handler
        if end_state:
            self.endStates.append(name)

    def set_start(self, name):
        self.startState = name.upper()

    def run(self, cargo):
        try:
            handler = self.handlers[self.startState]
        except:
            raise InitializationError('Must call .set_start() before .run()')
        if not self.endStates:
            raise InitializationError('At least one state must be an end state')
        while True:
            (newState, cargo) = handler(cargo)
            print('Reached', newState)
            if newState.upper() in self.endStates:
                print('Success!', newState)
                print(cargo)
                break
            else:
                handler = self.handlers[newState.upper()]



