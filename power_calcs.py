'''
Define initial power requirements for the turbine from the fundamental equation.
If one variable is missing, calculate it.
If 2 variables are missing, plot the missing variables.

Equation:

P = effiency * rho * A * V^3

'''
import numpy as np
import matplotlib.pyplot as plt

def P_eq(effiency, rho, A, V):
    return effiency * rho * A * V**3

def effiency_eq(P, rho, A, V):
    return P / (rho * A * V**3)

def rho_eq(P, effiency, A, V):
    return P / (effiency * A * V**3)

def A_eq(P, effiency, rho, V):
    return P / (effiency * rho * V**3)

def V_eq(P, effiency, rho, A):
    return (P / (effiency * rho * A))**(1/3)


def power_calc(P=False, effiency=False, A=False, V=False):
    # if one variable is missing, calculate it
    # if 2 variables are missing return a plot of the missing variables
    # if 3 or more missing, return an error

    # define variable ranges
    res = 10
    P_range = np.linspace(400, 800, res) # W
    effiency_range = np.linspace(0, 1, res) # decimal
    rho = 1000 # kg/m^3 - density of water
    A_range = np.linspace(0, 2, res) # m^2
    V_range = np.linspace(0, 3, res) # m/s

    inputs = [P, effiency, A, V]
    input_names = ['P', 'effiency', 'A', 'V']
    ranges = [P_range, effiency_range, A_range, V_range]
    missing_var = []
    for i in inputs:
        if i == False:
            missing_var.append(False)
        else:
            missing_var.append(True)

    if missing_var.count(False) == 1:
        # calculate the missing variable
        if missing_var[0] == False:
            P = P_eq(effiency, rho, A, V)
            print('Power = ', P)
        elif missing_var[1] == False:
            effiency = effiency_eq(P, rho, A, V)
            print('Effiency = ', effiency)
        elif missing_var[2] == False:
            A = A_eq(P, effiency, rho, V)
            print('Area = ', A)
        elif missing_var[3] == False:
            V = V_eq(P, effiency, rho, A)
            print('Velocity = ', V)

        return P, effiency, rho, A, V
    
    elif missing_var.count(False) == 2:
        # plot the missing variables
        # either P + effiency, P + A, P + V, effiency + A, effiency + V, A + V
        # identify which variables are missing and replace with a range
        vars = [P, effiency, A, V]
        int_vars = [] # vars of interest
        for i in range(len(missing_var)):
            if missing_var[i] == False:
                vars[i] = ranges[i]
                int_vars.append(i)
        
        # Define a dictionary mapping indices to plot titles
        plot_equations = {
            (0, 1): P_eq,
            (0, 2): P_eq,
            (0, 3): P_eq,
            (1, 2): effiency_eq,
            (1, 3): effiency_eq,
            (2, 3): A_eq
        }

        # Iterate over the dictionary
        for indices, equation in plot_equations.items():
            if int_vars == list(indices):
                print(f"Plotting {equation}")
                
                y = equation(vars[0], vars[1], vars[2], vars[3])
                x = vars[int_vars[1]]

                # plot the equation
                plt.figure()
                plt.plot(x, y)  
                plt.title(f"{input_names[int_vars[0]]} vs {input_names[int_vars[1]]}")
                plt.xlabel(input_names[int_vars[1]])
                plt.ylabel(input_names[int_vars[0]])
                plt.show()
                pass

       

        return P, effiency, rho, A, V
    else:
        # return an error
        raise ValueError('Too many missing variables. Please enter 3 or more variables.')
    

if __name__ == '__main__':


    P, efficiency, rho, A, V = power_calc(V=1.5, effiency=0.35)