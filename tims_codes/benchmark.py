
import pickle
flux = pickle.load(open('fluxes_1200.p', 'rb'))
from matplotlib import pyplot as plt
plt.plot(flux['vel_kms'],flux['total'])
plt.show()