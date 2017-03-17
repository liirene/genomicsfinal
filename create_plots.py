import numpy as np
import matplotlib.pyplot as plt

# use this as a module and a script

# Store info here w/in Protein_values

x = np.array([0,4,8,12,24,48])
y = np.array([1,2,3,4,5,6])

# Initializing figure values
fig, ax = plt.subplots() # initializes figure and axes object
plt.scatter(x= x, y = y, c= 'red')
            # cmap = colormapinstnace
plt.ylabel('Log2 Fold Change')
plt.xlabel('Hours')
plt.ylim(0, max(y) + 20)
plt.xlim(0, 50.0)

# regression lines

# calculate coefficients and store them: array outputs are:
# polynomial coeffs, highest power first (3,2,1,0)
fit = np.polyfit(x=x, y=y, deg=3)
ax.plot(x, fit[0]*x**3 + fit[1]*x**2 + fit[2]*x + fit[3], color = 'blue')
ax.scatter(x,y)

# save shit to disk
# Make multiple graphs at once?

# plot all proteins involved in one pathway on the same graph?

plt.close('all') # closes all open figures
# plt.close('fig') will close fig object (given that you assigned it)
