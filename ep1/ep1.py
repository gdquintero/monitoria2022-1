import numpy as np

# Computing LU decompotition of a tridiagonal matrix
def lu(a,b,c,n):
    l = np.zeros(n-1)
    u = np.zeros(n)
    u[0] = b[0]
    
    for i in range(1,n):
        l[i-1] = a[i] / u[i-1]
        u[i] = b[i] - l[i-1] * c[i-1]

    return l, u

# Computing the solution of LUx = d
def solve_lu(l,u,d,n):
    x = np.zeros(n)
    y = np.zeros(n)
    y[0] = d[0]

    for i in range(1,n):
        y[i] = d[i] - l[i-1] * y[i-1]

    x[n-1] = y[n-1] / u[n-1]

    for i in range(n-2,-1,-1):
        x[i] = (y[i] - c[i] * x[i+1]) / u[i]
    return x

# Computing the solution for the cyclic tridiagonal system
def solve_lu_cyclic(a_cyc,b_cyc,c_cyc,d_cyc,n):
    a = np.zeros(n-1)
    c = np.zeros(n-1)
    v = np.zeros(n-1)
    w = np.zeros(n-1)
    x_cyc = np.zeros(n)

    a[1:]  = a_cyc[1:-1]
    c[:-1] = c_cyc[:-2]
    w[0]   = c_cyc[-1]; w[-1] = a_cyc[-1]
    v[0]   = a_cyc[0]; v[-1] = c_cyc[-2]

    l, u = lu(a,b_cyc[:-1],c,n-1)
    y = solve_lu(l,u,d[:-1],n-1)
    z = solve_lu(l,u,v,n-1)

    x_cyc[-1] = (d_cyc[-1] - c_cyc[-1]*y[0] - a_cyc[-1]*y[-1]) / \
                (b_cyc[-1] - c_cyc[-1]*z[0] - a_cyc[-1]*z[-1])

    x_cyc[:-1] = y[:] - x_cyc[-1]*z[:]

    return x_cyc

# Starting Task 1
n = 5

a = np.zeros(n)
b = np.zeros(n)
c = np.zeros(n)
d = np.zeros(n)

for i in range(n):
    ind = i + 1
    if (i < n - 1):
        a[i] = ((2. * ind) - 1.) / (4. * ind)
    else:
        a[i] = ((2 * ind) -1) / (2. * ind)
    
    c[i] = 1. - a[i]
    b[i] = 2.
    d[i] = np.cos((2. * np.pi * (ind**2)) / (n**2))


x = solve_lu_cyclic(a,b,c,d,n)
print("The solution of the cyclic tridiagonal system is: ","\n")
print(x)


