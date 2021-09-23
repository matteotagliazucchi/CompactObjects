import numpy as np  

# 4th order Runge Kutta step

def rk_fixedstep(x, y, step, func, args):
    """
    4th order Runge Kutta with fixed step. Works for vectorial ODE dy/dx = f(x,y).
    It computes y_{n+1} once we provide x_n, y_n and a step
    
    Parameters:
    ===========
    x: float
        independent variable 
    y: array of float
        dependent variables "y"
    step: float
        step size/increment
    func: callable
        function that calculate the derivatives of dependent variables y w.r.t. x. This arguments basically is f(x,y).
    args: tuple
        additional, and optional, arguments of dxdy
    -------------
    y_next: array of float
        array containing the points y_{n+1} = y_n + (h/6) (k_1 + 2 k_2 + 2k_3 + k_4) 
        where k_i are functions of y_n, x_n and f(x_n,y_n) (Valid for n=0,1,2,3...) and the "x" points. 
    """

    n = len(np.atleast_1d(y))  # number of differential equations
    
    k1 = np.zeros(n)
    k2 = np.zeros(n)
    k3 = np.zeros(n)
    k4 = np.zeros(n)

    y_next = np.zeros(n)

    # inizialization of k_i's

    if args is None :
        k1 = step*func(x, y)
        k2 = step*func(x + step/2, y + k1/2)
        k3 = step*func(x + step/2, y + k2/2)
        k4 = step*func(x + step, y + k3)
    else :
        k1 = step*func(x, y, args)
        k2 = step*func(x + step/2, y + k1/2, args)
        k3 = step*func(x + step/2, y + k2/2, args)
        k4 = step*func(x + step, y + k3, args)

    # filling y_next
   
    y_next = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    return y_next


# Ode solver (hopefully similar to scipy.integrate.odeint)

def rk_odesolver(func, y0, x_interval, args = None, stop_event = None):
    """
    4th order Runge Kutta algorithm without any condition. Works for vectorial ODE dy/dx = f(x,y).

    Parameters:
    ===========
    x_interval : array of float
        array with equal-spaced interval of the independent variable 
    y0: array of floats
        array with initial conditions for the dependent variables 
    func: callable
        function that calculates the derivatives of dependent variables y w.r.t. x . This arguments basically is f(x,y).
    args : tuple
        additional, and optional, arguments of dxdy
    stop_event: boolean
        condition for stopping the loop. 
    ------------
    solutions: multidimensional array of floats
        multi-dim array that contains the points y_{n+1} = y_n + (h/6) (k_1 + 2 k_2 + 2k_3 + k_4) 
        where k_i are functions of y_n, x_n and f(x_n,y_n) (Valid for n=0,1,2,3...).
    """

    n = len(np.atleast_1d(y0))  # number of differential equations
    number_of_iterations = len(x_interval)
    initial_step = np.diff(x_interval)[0] 

    solutions = np.zeros((number_of_iterations, n)) # we store all the y_n points

    if stop_event is not None:
        
        # inizialization to check if we reached the event
        check_event = False
    
        for i in range(number_of_iterations):
        
            # we make a step forward
            solutions[i,:] = rk_fixedstep(
                x_interval[i], y0, initial_step, func, args)
            
            # check if we reached our solution
            if bool(stop_event(x_interval[i],y0)): 
                #print("Event reached. Return all the points until the reached event.")
                check_event = True
                return solutions[0:i,:], check_event

            # if we do not reach the event we redefine the variable to make the loop continue
            else :  
                 y0 = solutions[i, :]
        
        # if we end the loop and we did not reach the event
        #print("Event NOT reached. Return all the points computed so far.")   
        return solutions, check_event 

       
    else: # if we did not set a stop event we return all the points
        #print("Return all the calculated points.") 
        return solutions