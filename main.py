from src.mccormack import McCormack
from src.scenarios import gaussian_density_example, with_disturbance, numerical_oscilations, lecture_high_res, lecture_low_res

import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt


if __name__ == "__main__":
    rho0, v0, params = gaussian_density_example()  # Produces negative values in flow and speed
    model = McCormack(**params)
    df_rho, df_Q = model.run(rho0, v0, disturbance = {})

    # Example with reduced jam density in during simulation
    # rho0, v0, params, disturbance = with_disturbance()
    # model = McCormack(**params)
    # df_rho, df_Q = model.run(rho0, v0, disturbance = disturbance)

    # Example with numerical oscilations by time step too large for changes
    # rho0, v0, params, disturbance = numerical_oscilations()
    # model = McCormack(**params)
    # df_rho, df_Q = model.run(rho0, v0, disturbance = disturbance)

    # Example with high resolution
    # rho0, v0, params = lecture_high_res()
    # model = McCormack(**params)
    # df_rho, df_Q = model.run(rho0, v0, disturbance = {})

    # Example with low resolution
    # rho0, v0, params = lecture_low_res()
    # model = McCormack(**params)
    # df_rho, df_Q = model.run(rho0, v0, disturbance = {})

    print(f"min rho: {df_rho.min().min()}  max rho: {df_rho.max().max()}")
    print(f"min Q: {df_Q.min().min()}  max Q: {df_Q.max().max()}")
    
    kwds = dict(
        cmap = 'seismic', 
        cbar = True, 
        cbar_kws = {'orientation': 'horizontal', 'shrink': .8, 'aspect': 30, 'location': 'top'},
        # center = 0 # for negative and positive values -- e.g.
        )  

    fig, ax = plt.subplots(1, 3, figsize=(21, 7))
    sns.heatmap(df_rho.T, ax=ax[0], **kwds)
    ax[0].set_title('Density ($\\rho$)')
    sns.heatmap(df_Q.T, ax=ax[1], **kwds)
    ax[1].set_title('Flow ($Q$)')
    sns.heatmap((df_Q / df_rho).T, ax=ax[2], **kwds)
    ax[2].set_title('Density / Flow ($v = \\frac{\\rho}{Q}$)')

    for axis in ax:
        axis.set_xticks(np.arange(0, len(model.t), 100))
        axis.set_xticklabels(np.round(model.t[::100], 2))

        axis.set_yticks(np.arange(0, len(model.x), 10))
        axis.set_yticklabels(np.round(model.x[::10], 2))
        axis.set_xlabel('t (hrs)')
        axis.set_ylabel('x (km)')
        
    plt.tight_layout()
    plt.show()
