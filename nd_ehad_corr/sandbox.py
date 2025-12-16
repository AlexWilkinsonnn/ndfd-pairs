import matplotlib.pyplot as plt
import numpy as np
from lmfit import Minimizer, create_params, report_fit
from lmfit.lineshapes import gaussian, lorentzian
def residual(pars, x, data):
  model = (gaussian(x, pars['amp_g'], pars['cen_g'], pars['wid_g']) +
  lorentzian(x, pars['amp_l'], pars['cen_l'], pars['wid_l']))
  return model - data


def fit_model(x, y):
  print(x)
  print(y)
  pfit = create_params(amp_g=10, cen_g=5, wid_g=1, amp_l=10,
    peak_split=dict(value=2.5, min=0, max=5),
    cen_l=dict(expr='peak_split+cen_g'),
    wid_l=dict(expr='wid_g'))
    
  mini = Minimizer(residual, pfit, fcn_args=(x, y))
  
  out = mini.leastsq()
  
  best_fit = y + out.residual
  
  report_fit(out.params)
  
  plt.plot(x, y, 'o')
  plt.plot(x, best_fit, '--', label='best fit')
  plt.legend()
  plt.show()
  
  return out.params['cen_g'], out.params['wid_g']

def main():
  np.random.seed(0)
  x = np.linspace(0, 20.0, 20)
  data = (gaussian(x, 21, 6.1, 1.2) + lorentzian(x, 10, 9.6, 1.3) +
    np.random.normal(scale=0.1, size=x.size))
    
  mean, unc = fit_model(x, data)
  
if __name__ == "__main__":
  main()


