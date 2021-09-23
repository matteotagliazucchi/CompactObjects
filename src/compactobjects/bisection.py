import numpy as np

def bisection(func,a,b,N):
    '''Approximate solution of f(x)=0 on interval [a,b] by bisection method.

    Parameters
    ===========
    func : callable
        The function for which we are trying to approximate a solution f(x)=0.
    a,b : float
        The interval in which to search for a solution. The function returns
        None if f(a)*f(b) >= 0 since a solution is not guaranteed.
    N : (positive) integer
        The number of iterations to implement.
    -------
    x_N : float
        The midpoint of the Nth interval computed by the bisection method. The
        initial interval [a_0,b_0] is given by [a,b]. If f(m_n) == 0 for some
        midpoint m_n = (a_n + b_n)/2, then the function returns this solution.
        If all signs of values f(a_n), f(b_n) and f(m_n) are the same at any
        iteration, the bisection method fails and return None.

    Examples:
    >>> f = lambda x: x**2 - x - 1
    >>> bisection(f,1,2,25)
    1.618033990263939
    >>> f = lambda x: (2*x - 1)*(x - 3)
    >>> bisection(f,0,1,10)
    0.5
    '''
    if func(a)*func(b) > 0: 
        #print("No solutions")
        return None
    if func(a)*func(b) == 0:
        #print("Found exact solution.")
        if func(a) == 0 : return a
        else: return b
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = func(m_n)
        if func(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif func(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            #print("Found exact solution.")
            return m_n
        else:
            #print("No solutions with bisection. Return 'None'.")
            return None
    #print("Found approximated solution.")
    return (a_n + b_n)/2