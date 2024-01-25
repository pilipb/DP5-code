import numpy as np
import matplotlib.pyplot as plt
import argparse

# define equations
def P_eq(efficiency, rho, A, V):
    return efficiency * rho * A * V**3

def efficiency_eq(P, rho, A, V):
    return P / (rho * A * V**3)

def rho_eq(P, efficiency, A, V):
    return P / (efficiency * A * V**3)

def A_eq(P, efficiency, rho, V):
    return P / (efficiency * rho * V**3)

def V_eq(P, efficiency, rho, A):
    return (P / (efficiency * rho * A))**(1/3)


def power_calc(P=False, efficiency=False, A=False, V=False):
    # if one variable is missing, calculate it
    # if 2 variables are missing return a plot of the missing variables
    # if 3 or more missing, return an error

    # define variable ranges
    res = 50
    P_range = np.linspace(400, 800, res) # W
    efficiency_range = np.linspace(0, 1, res) # decimal
    rho = 1000 # kg/m^3 - density of water
    A_range = np.linspace(0, 2, res) # m^2
    V_range = np.linspace(0.16, 2.6, res) # m/s

    inputs = [P, efficiency, A, V]
    input_names = ['P', 'efficiency', 'A', 'V']
    ranges = [P_range, efficiency_range, A_range, V_range]
    missing_var = []
    for i in inputs:
        if i == False:
            missing_var.append(False)
        else:
            missing_var.append(True)

    print('\n-------------------------------------------\n')
    if missing_var.count(False) == 1:
        # calculate the missing variable
        if missing_var[0] == False:
            print('\n---------Calculating Power---------')
            P = P_eq(efficiency, rho, A, V)
            print('Power = ', P)
        elif missing_var[1] == False:
            print('\n---------Calculating efficiency---------')
            efficiency = efficiency_eq(P, rho, A, V)
            print('efficiency = ', efficiency)
        elif missing_var[2] == False:
            print('\n---------Calculating Area---------')
            A = A_eq(P, efficiency, rho, V)
            print('Area = ', A)
        elif missing_var[3] == False:
            print('\n---------Calculating Velocity---------')
            V = V_eq(P, efficiency, rho, A)
            print('Velocity = ', V)

        print('\n-------------------------------------------')
        print('All variables calculated.')
        print('\nPower = ', P)
        print('\nefficiency = ', efficiency)
        print('\nArea = ', A)
        print('\nVelocity = ', V)
        print('\n-------------------------------------------')

        return P, efficiency, rho, A, V
    
    elif missing_var.count(False) == 2:
        print('\n-------------------------------------------')
        # plot the missing variables
        # either P + efficiency, P + A, P + V, efficiency + A, efficiency + V, A + V
        # identify which variables are missing and replace with a range
        vars = [P, efficiency, A, V]
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
            (1, 2): efficiency_eq,
            (1, 3): efficiency_eq,
            (2, 3): A_eq
        }

        # Iterate over the dictionary
        for indices, equation in plot_equations.items():
            if int_vars == list(indices):
                print('\nMissing variables: %s, %s'  %( input_names[int_vars[0]], input_names[int_vars[1]] ))
                # print("Plotting {equation}")
                
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
        print('\n-------------------------------------------')

        return P, efficiency, rho, A, V
    
    elif missing_var.count(False) == 0:
        print('\n-------------------------------------------')
        print('All variables entered.')
        print('\nPower = ', P)
        print('\nefficiency = ', efficiency)
        print('\nArea = ', A)
        print('\nVelocity = ', V)
        print('\n-------------------------------------------')
        return P, efficiency, rho, A, V
    
    else:
        # return an error
        raise ValueError('Too many missing variables. Please enter 3 or more variables.')


    

if __name__ == '__main__':
    '''
    Define initial power requirements for the turbine from the fundamental equation.
    If one variable is missing, calculate it.
    If 2 variables are missing, plot the missing variables.

    Equation:

    P = efficiency * rho * A * V^3

    P = Power (W)
    efficiency = efficiency (decimal)
    rho = density of water (kg/m^3)
    A = Area (m^2)
    V = Velocity (m/s)

    Usage:

    python3 power_calcs.py -P 400 -e 0.35 -A 0.8 -V 1
    (or any combination of 2 or more variables)

    Help:

    python3 power_calcs.py -h

    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-P', '--power', default=False, type=float, help='Power (W)')
    parser.add_argument('-e', '--efficiency', default=False, type=float, help='Efficiency (decimal)')
    parser.add_argument('-A', '--area', default=False, type=float, help='Area (m^2)')
    parser.add_argument('-V', '--velocity', default=False, type=float, help='Velocity (m/s)')
    args = parser.parse_args()

    P = args.power
    efficiency = args.efficiency
    A = args.area
    V = args.velocity

    P, efficiency, rho, A, V = power_calc(V=V, efficiency=efficiency, A=A, P=P)