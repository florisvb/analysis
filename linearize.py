# Author: Floris van Breugel
# Date: 15 October 2009

# this function returns the A,B,C,D matrices for the linear state space representation of 
# an arbitrary first order (state) representation of a dynamical system. 
# xdot = A*(x-x0) + B*(u-u0) + xdot(x0,u0)
# y = C*(x-x0) + D*(u-u0) + y(x0,u0)


# xdot, x, y, and u, can be sympy Symbols, sympy Matrices, or lists (containing sympy Symbols)

# x is the list of states, in symbolic form, as used in xdot

# xdot is the list of the state derivatives, in terms of the symbols used in x

# u is the list of inputs    

# y is the desired output, in terms of the symbols used in x and u, this is optional. if not entered, only  A and B (and F0) matrices are returned

# x0 and u0 are optional, they are the point about which you wish to linearize. If not entered the function will linearize about zero for all states/inputs.


import numpy as np
import sympy as sp
from sympy import sin, cos
pi = np.pi
import scipy
#from control.matlab import ss



def linearize(xdot, x, u, y=None, x0=None, u0=None):

    
   

    # we need matrices... if input is a symbol or list of symbols, convert them to matrices
    if type(xdot).__name__!='Matrix':
    
        if type(xdot).__name__=='Symbol':
            xdot = sp.Matrix([xdot])
        elif type(xdot).__name__=='list':
            xdot = sp.Matrix(xdot)
        elif 1:
            raise ValueError, "Please enter xdot as a sympy Matrix, Symbol, or list of Symbols"
            
    if type(x).__name__!='Matrix':
    
        if type(x).__name__=='Symbol':
            x = sp.Matrix([x])
        elif type(x).__name__=='list':
            x = sp.Matrix(x)
        elif 1:
            raise ValueError, "Please enter x as a sympy Matrix, Symbol, or list of Symbols"
    
    if type(u).__name__!='Matrix':
    
        if type(u).__name__=='Symbol':
            u = sp.Matrix([u])
        elif type(u).__name__=='list':
            u = sp.Matrix(u)
        elif 1:
            raise ValueError, "Please enter u as a sympy Matrix, Symbol, or list of Symbols"
            
    if u0 is not None:
        if type(u0).__name__!='Matrix':
        
            if type(u0).__name__=='list':
                u0 = sp.Matrix(u0)
            elif 1:
                raise ValueError, "Please enter u0 as a sympy Matrix, or list of values"
                
    if x0 is not None:        
        if type(x0).__name__!='Matrix':
        
            if type(x0).__name__=='list':
                x0 = sp.Matrix(x0)
            elif 1:
                raise ValueError, "Please enter x0 as a sympy Matrix, or list of values"
        
    if y is not None:
        if type(y).__name__!='Matrix':
        
            if type(y).__name__=='Symbol':
                y = sp.Matrix([y])
            if type(y).__name__=='list':
                y = sp.Matrix(y)
            elif 1:
                raise ValueError, "Please enter y as a sympy Matrix, Symbol, or list of Symbols"

    # number of states
    len_x = max(x.shape)
    len_xdot = max(xdot.shape)
    len_u = max(u.shape)
    if y is not None:
        len_y = max(y.shape)
        
    # check sizes, should probably have more checks...
    if len_x != len_xdot:
        print('Error: incompatible x and xdot')

    # check operating point: if None, make it zeros
    if x0 is None:
        x0 = sp.zeros([len_x,1])
    if u0 is None: 
        u0 = sp.zeros([len_u,1])  

    
        
    # substitution list: for operating point
    sub_arr = []
    for i in range(0,len_x):
        sub_arr.append((x[i],x0[i]))
    for i in range(0,len_u):
        sub_arr.append((u[i],u0[i]))  
    
    ## Form A Matrix ##
    DA = xdot.jacobian(x)   

    # plug in operating point
    A = DA.subs(sub_arr)
    
    ## Form B Matrix ##
    DB = sp.zeros([len_xdot,len_u])
    for i in range(0,len_xdot):
        for j in range(0,len_u):
            DB[i,j] = xdot[i].diff(u[j])
    
    B = DB.subs(sub_arr)

    if y is not None:

        ## Form C Matrix ##
        DC = sp.zeros([len_y,len_x])
        for i in range(0,len_y):
            for j in range(0,len_x):
                DC[i,j] = y[i].diff(x[j])

        C = DC.subs(sub_arr)

        
        ## Form D Matrix ##
        DD = sp.zeros([len_y,len_u])
        for i in range(0,len_y):
            for j in range(0,len_u):
                DD[i,j] = y[i].diff(u[j])
        
        D = DD.subs(sub_arr)
    
        # constants in y
        G0 = y.subs(sub_arr)


    # note if F0 is not zero, we're not linearizing about an eq. point
    F0 = xdot.subs(sub_arr)
    
    
    if y is not None: 
        return [A,B,C,D,F0,G0]
    if y is None: 
        return [A,B,F0]
    

def linsys(xdot, x, u, y, x0=None, u0=None):

    # y is required for linsys, but not linearize
    # apparently 'ss' does not support multiple outputs; linearize does
    
    
    As,Bs,Cs,Ds,F0,G0 = linearize(xdot, x, u, y, x0, u0)
    
    sumF0 = 0
    for i in F0:
        sumF0 += i
    if sumF0 > 0.001:
        print('Warning: The system was not linearized about an equilibrium point!')
        print
        print 'xdot at x0 = ', F0
        
        
    if Cs.shape[0] > 1:
        raise ValueError, "C matrix cannot have more than one row; system must be SISO"
        

    A = scipy.matrix(As).astype(np.float)
    B = scipy.matrix(Bs).astype(np.float)
    C = scipy.matrix(Cs).astype(np.float)
    D = scipy.matrix(Ds).astype(np.float)
    
    
    
    sys = 0 #ss(A,B,C,D)
    
    return sys
    
def atest_1():
    
    theta, thetadot = sp.symbols('theta, thetadot')

    u = sp.Symbol('u')


    u0 = [0]

    x = [theta, thetadot]
    xdot = [thetadot, -sin(theta)]
    x0 = [pi/4.0, 0]
    y = [2*theta, thetadot**(1/2)]
    
    A,B,C,D,F0 = linearize(xdot, x, u, y, x0, u0) 
    
    print
    print
    print 'xdot = A*(x-x0) + B*(u-u0) + xdot(x0,u0)'
    print 'y = C*(x-x0) + D*(u-u0) + y(x0,u0)'
    print
    
    print'A= '
    print A
    print
    print'B= '
    print B
    print
    print'C= '
    print F0
    print
    print'D= '
    print F0
    print
    print'F0= '
    print F0
    
    return 1
    
def test_2():

    x1,x2,x3 = sp.symbols('x1,x2,x3')
    x = sp.Matrix([x1,x2,x3])
    A = sp.Matrix([[-4, -5, -2,],[1,0,0],[0,1,0]])
    B = sp.Matrix([[1,0,0],[1,2,3]]).T
    u1,u2 = sp.symbols('u1,u2')
    u = sp.Matrix([u1,u2])


    xdot = A*x+B*u
    y = [x1+1]
    
    A,B,C,D,F0,G0 = linearize(xdot,x,u,y)
    
    print
    print
    print 'xdot = A*(x-x0) + B*(u-u0) + xdot(x0,u0)'
    print
    
    print'A= '
    print A
    print
    print'B= '
    print B
    print
    print'C= '
    print C
    print
    print'D= '
    print D
    print
    print'F0= '
    print F0
    print
    print'G0= '
    print G0
    
    sys = linsys(xdot,x,u,y)
    print sys
    
    return 1
    
if __name__=='__main__':
    import nose
    nose.runmodule()
    
    # to run tests: nosetests -s linearize.py
