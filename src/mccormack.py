from typing import Union
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import warnings


class McCormack:
    valid_disturbances = ['density factor', 'flow factor', 'max flow value', 'rho jam']

    def __init__(self, L, nx, dt, T, rho_jam, v0, tau, clip_rho=(None, None), assert_cfl = True):
        self.L = L
        self.nx = nx 
        self.dx = L / nx
        self.x = np.linspace(0, L, nx)

        self.dt = dt 
        self.T = T
        self.t = np.arange(0, T, dt)

        if self.dt > self.dx / v0 and assert_cfl:
            raise ValueError(f"CFL (dt < dx / v0) condition not satisfied with dt: {self.dt}, dx: {self.dx}, v0: {v0}")

        self.v0 = v0 
        self.rho_jam = rho_jam
        self.tau = tau

        self.clip_rho = clip_rho

    def Ve(self, rho):
        if any(self.clip_rho):
            rho = np.clip(rho, self.clip_rho[0], self.clip_rho[1])
        
        return self.v0 * (1 - rho / self.rho_jam)

    def pressure(self, rho):
        if any(self.clip_rho):
            rho = np.clip(rho, self.clip_rho[0], self.clip_rho[1])
        P = (self.v0 - self.Ve(rho)) / (2 * self.tau) 
        return P

    def flux(self, rho, Q):
        if any(self.clip_rho):
            rho = np.clip(rho, self.clip_rho[0], self.clip_rho[1])

        f1 = Q 
        f2 = Q**2 / rho + self.pressure(rho)
        return f1, f2

    def initialize(self, rho0, v0):
        if rho0 is None:
            rho0 = np.ones(self.nx)

        if v0 is None:
            v0 = self.Ve(rho0)

        Q0 = rho0 * v0

        self.rho = [rho0.copy()]
        self.Q = [Q0.copy()]

        return rho0.copy(), Q0.copy()
    
    def source(self, rho, Q):
        if any(self.clip_rho):
            rho = np.clip(rho, self.clip_rho[0], self.clip_rho[1])

        s1 = np.zeros_like(rho)
        s2 = (rho * self.Ve(rho) - Q) / self.tau
        return s1, s2

    def predict(self, rho, Q):
        f1, f2 = self.flux(rho, Q)
        s1, s2 = self.source(rho, Q)

        rho_pred = rho.copy() 
        Q_pred = Q.copy()

        rho_pred[1:-1] -= self.dt / self.dx * (f1[2:] - f1[1:-1]) + self.dt * s1[1:-1]
        Q_pred[1:-1]   -= self.dt / self.dx * (f2[2:] - f2[1:-1]) + self.dt * s2[1:-1]

        # Slight hack for boundary conditions by copying neighboring cell states
        rho_pred[0] = rho_pred[1]
        rho_pred[-1] = rho_pred[-2]

        Q_pred[0] = Q_pred[1]
        Q_pred[-1] = Q_pred[-2]

        return rho_pred, Q_pred

    def correct(self, rho, Q, rho_pred, Q_pred):
        f1_pred, f2_pred = self.flux(rho_pred, Q_pred)
        s1_pred, s2_pred = self.source(rho_pred, Q_pred)

        rho_new = rho.copy() 
        Q_new = Q.copy()

        rho_new[1:-1] = 0.5 * (rho_pred[1:-1] + rho[1:-1] - self.dt / self.dx * (f1_pred[1:-1] - f1_pred[0:-2]) + self.dt * s1_pred[1:-1])
        Q_new[1:-1]   = 0.5 * (Q_pred[1:-1]   + Q[1:-1]   - self.dt / self.dx * (f2_pred[1:-1] - f2_pred[0:-2]) + self.dt * s2_pred[1:-1])


        rho_new[0]  = rho_new[1]
        rho_new[-1] = rho_new[-2]

        Q_new[0]  = Q_new[1]
        Q_new[-1] = Q_new[-2]


        return rho_new, Q_new

    def run(self, rho0: np.array = None, v0: np.array = None, disturbance = {}) -> Union[pd.DataFrame, pd.DataFrame]:
        self.default_rho_jam = self.rho_jam

        rho, Q = self.initialize(rho0, v0)

        self.rho_pred = []
        self.Q_pred = []

        for ti in self.t[1:]:
            rho_pred, Q_pred = self.predict(rho, Q)

            self.rho_pred.append(rho_pred.copy())
            self.Q_pred.append(Q_pred.copy())

            rho, Q = self.correct(rho, Q, rho_pred, Q_pred)

            # Possible to include a disturbance during the simulation by changing density and/or flow for given time-steps at given locations
            if disturbance:
                if disturbance['time'][0] <= ti < disturbance['time'][1]:
                    disturb_idx = np.where((self.x > disturbance['location'][0]) & (self.x < disturbance['location'][1]))

                    if not any([disturbance_type in self.valid_disturbances for disturbance_type in disturbance]):
                        
                        warnings.warn("\n\nNo disturbance defined from " + ",".join(self.valid_disturbances) + "\n\n")
                    
                    if 'density factor' in disturbance:
                        rho[disturb_idx] *= disturbance['density factor']

                    if 'flow factor' in disturbance:
                        Q[disturb_idx] *= disturbance['flow factor']

                    if 'max flow value' in disturbance:
                        print(Q[disturb_idx], disturbance['max flow value'])
                        Q[disturb_idx] = np.clip(Q[disturb_idx], None, disturbance['max flow value'])

                    if 'rho jam' in disturbance:
                        self.rho_jam = disturbance['rho jam']

            if 'rho jam' in disturbance and ti >= disturbance['time'][1]:
                self.rho_jam = self.default_rho_jam
            
            if any(self.clip_rho):
                rho = np.clip(rho, self.clip_rho[0], self.clip_rho[1])

            self.rho.append(rho.copy())
            self.Q.append(Q.copy())

        print("Simulation finished!..")
            
        return pd.DataFrame(self.rho, index = self.t, columns = self.x), pd.DataFrame(self.Q, index = self.t, columns = self.x)
