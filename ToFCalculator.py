#############################################################################################

# *    Title: TimeDifferences.py
# *    Author: Reilly Schaffer
# *    Date: 03/31/2025
# *    Code version: 1
# *    Availability: https://github.com/phann21/ResearchProgress/blob/main/TimeDifferences.py

##############################################################################################
import math
from scipy.optimize import fsolve
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

# Initial Transducer locations
T1 = [0, 0, 0]
T2 = [-10, 17.321, -2]
T3 = [-10, -17.321, -2]
T4 = [20, 0, -2]

ME = [0,0,-20]

# Rate of propagation
rate = 1481*1000

# System of equations
def system_of_equations(vars):
    t1, t2, t3, t4 = vars
    eq1 = math.sqrt((ME[0] - T1[0])**2 + (ME[1] - T1[1])**2 + (ME[2] - T1[2])**2) - rate * (t1)
    eq2 = math.sqrt((ME[0] - T2[0])**2 + (ME[1] - T2[1])**2 + (ME[2] - T2[2])**2) - rate * (t2)
    eq3 = math.sqrt((ME[0] - T3[0])**2 + (ME[1] - T3[1])**2 + (ME[2] - T3[2])**2) - rate * (t3)
    eq4 = math.sqrt((ME[0] - T4[0])**2 + (ME[1] - T4[1])**2 + (ME[2] - T4[2])**2) - rate * (t4)

    return [eq1, eq2, eq3, eq4]

# Initial guess for the solution
initial_guess = [1,1, 5, 1]

# Solve the system of equations
solution = fsolve(system_of_equations, initial_guess)
t1, t2, t3, t4 = solution[0], solution[1],solution[2],solution[3]
print("Solution: t1 =", t1, "t2 =", t2,"t3 =", t3, "t4 =", t4)
fig = plt.figure()

# syntax for 3-D projection
ax = plt.axes(projection ='3d')

# plotting
ax.scatter3D(T1[0], T1[1], T1[2], color='red')
ax.scatter3D(T2[0], T2[1], T2[2], color='red')
ax.scatter3D(T3[0], T3[1], T3[2], color='red')
ax.scatter3D(T4[0], T4[1], T4[2], color='red')
ax.scatter3D(ME[0], ME[1], ME[2], color='green')

plt.show()