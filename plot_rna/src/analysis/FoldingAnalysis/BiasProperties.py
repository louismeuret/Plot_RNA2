# -*- coding: utf-8 -*- 
class BiasProperties:
    def __init__(self, time = None, bias_potential = None, bias_penality = None, cum_tot_force = None, bias_inst_force = None, z_coord = None,
                s_coord = None, w_coord = None, z_min = None, s_min = None, w_min = None, closeness = None, 
                progress = None):
        self.time = time
        self.bias_potential = bias_potential
        self.bias_penality = bias_penality
        self.cum_tot_force = cum_tot_force
        self.bias_inst_force = bias_inst_force
        self.z_coord = z_coord
        self.s_coord = s_coord
        self.w_coord = w_coord
        self.z_min = z_min
        self.s_min = s_min
        self.w_min = w_min
        self.closeness = closeness
        self.progress = progress
