from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.close('figure 1')

data = open('test_gold.dat', 'r')
x, y = np.loadtxt(data, delimiter=' ', usecols=(0, 1), unpack=True)

ref = open('ref.dat', 'r')
x_ref, y_ref = np.loadtxt(ref, delimiter=' ', usecols=(0, 1), unpack=True)

def func(s, depth, alpha, eq,a):
    return depth*(1-np.exp(-alpha*(s-eq)))*(1-np.exp(-alpha*(s-eq))) + a

popt, pcov = curve_fit(func, x, y)
std_dev = np.sqrt(np.diag(pcov)) #Standard deviation errors on the parameters
print popt
print std_dev

table = open('latex_result_table.txt','w')
for i in xrange(0,len(x_ref)):
	s = x_ref[i]
	E_calc = func(s,*popt)
	E_ref = y_ref[i]
	delta = np.absolute(E_ref-E_calc)
	err1 = np.absolute(func(s,popt[0]+std_dev[0],popt[1]+std_dev[1],popt[2]+std_dev[2],popt[3]+std_dev[3])-func(s,popt[0]-std_dev[0],popt[1]-std_dev[1],popt[2]-std_dev[2],popt[3]-std_dev[3]))
	err2 = np.absolute(func(s,popt[0]+std_dev[0],popt[1]-std_dev[1],popt[2]-std_dev[2],popt[3]+std_dev[3])-func(s,popt[0]-std_dev[0],popt[1]+std_dev[1],popt[2]+std_dev[2],popt[3]-std_dev[3]))
	err = np.maximum(err1,err2)
	table.write("\\hline\n%f & %f & %f & %f & %f \\\\ \n " % (s, E_calc, E_ref, delta, err))


plot(x,y,label = "Raw data")
plot(x,func(x,*popt), label = "Morse potential fit")
plot(x_ref,y_ref, label = "Litterature values")
legend = plt.legend()
show()

