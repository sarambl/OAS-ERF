# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Adding a temperature dependency

# %% [markdown]
# ### Paper reference: Yu et al (2017)
# Yu, F., Luo, G., Nadykto, A. B., & Herb, J. (2017). Impact of temperature dependence on the possible contribution of organics to new particle formation in the atmosphere. Atmospheric Chemistry and Physics, 17(8), 4997–5005. https://doi.org/10.5194/acp-17-4997-2017
#

# %% [markdown]
# \begin{equation}
# f_T = exp \Big[ \frac{\Delta H}{k} \big( \frac{1}{T} - \frac{1}{T_0} \big) \Big]
# \end{equation}

# %% [markdown]
# - Best fit with obsevations according to Yu et al (2017) is by $\Delta H = 35$ kcal/mol.
# - $T_0 = 278.$ 
# - $k= 1.380649 \times 10^{−23}$ J/K
# - $k = 0.001985875$ kcal/mol/K
#

# %% [markdown]
# It is unclear how the authors set the maximum value, but they say in the review that they initially used 10 as a max and that the difference was negligible. 

# %%
import numpy as np

# %%
T_0 = 278
k_kcalmolK =  0.001985875 
dH = 35

def fT(T, max_v=10):
    f =  np.exp(dH/k_kcalmolK*(1/T-1/T_0))
    if type(T) is np.ndarray:
        f[f>max_v] = max_v
    else:
        f = min(f,max_v)
    return f


# %%
fT(270)

# %%
T = np.linspace(270,310)


f = fT(T)

import matplotlib.pyplot as plt
plt.plot(T,f)
plt.yscale('log')
plt.xlim([270,310])
plt.ylim([1e-3,10])

# %%
a = [1,2,4]

# %%
[]

# %%
type(T)

# %%
