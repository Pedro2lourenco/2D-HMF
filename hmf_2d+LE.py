import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,sqrt,pi,square,log


E_in = .498  # Initial energy of the system

h = 0.5      # Time step
N = 5        # Number of particles

# Phase space variables (positions and momenta)
x = np.zeros(N)
y = np.zeros(N)
px = np.zeros(N)
py = np.zeros(N)

n = 4  # Controls periodic box size

# Coupling constants
c1 = c2 = 1
c = -1
d = 1

lyap_num = 2*N  # Number of Lyapunov exponents (dimension of tangent space)

# Initial perturbation vector (random normalized vector)
w = 2*np.random.randint(4*N)-np.ones(4*N)
w = w/np.dot(w,w)


# =============================
# EQUATIONS OF MOTION
# =============================

def evolveQ(X,Y,PX,PY,delta):
    # Update positions using momenta (Hamiltonian dynamics)
    X += PX*delta
    Y += PY*delta
    
    return(X,Y)

def evolveP(X,Y,PX,PY,delta):
    # Compute mean-field terms (magnetizations and correlations)
    M1x = np.sum(np.cos(X))/N
    M1y = np.sum(np.sin(X))/N
    M2x = np.sum(np.cos(Y))/N
    M2y = np.sum(np.sin(Y))/N
    Ppx = np.sum(np.cos(X+Y))/N
    Ppy = np.sum(np.sin(X+Y))/N
    Pmx = np.sum(np.cos(X-Y))/N
    Pmy = np.sum(np.sin(X-Y))/N
    
    # Forces derived from the Hamiltonian
    Fx = c*((M1y*np.cos(X))-(M1x*np.sin(X))) + \
         0.5*d*((Ppy*np.cos(X+Y))-(Ppx*np.sin(X+Y)) +
                (Pmy*np.cos(X-Y))-(Pmx*np.sin(X-Y)))
    
    Fy = c*((M2y*np.cos(Y))-(M2x*np.sin(Y))) + \
         0.5*d*((Ppy*np.cos(X+Y))-(Ppx*np.sin(X+Y)) -
                (Pmy*np.cos(X-Y))+(Pmx*np.sin(X-Y)))
    
    # Update momenta
    PX += Fx*delta
    PY += Fy*delta
    
    return(PX,PY)


# =============================
# SUZUKI-TROTTER 4th ORDER COEFFICIENTS
# =============================

s = 1/(4-4**(1/3))

delta1 = h*s
delta2 = 0.5*(h*s)
delta3 = 0.5*h*(1-(3*s))
delta4 = (1-(4*s))*h


# =============================
# LINEARIZED (TANGENT) DYNAMICS
# Used for Lyapunov exponent calculation
# =============================

def evolvedQ(dX,dY,dPx,dPy,delta):
    # Linearized position evolution
    dX += dPx*delta 
    dY += dPy*delta 
    
    return(dX,dY)

def evolvedP(X,Y,dX,dY,dPx,dPy,delta):
    
    # Compute base system averages
    M1x = np.sum(np.cos(X))/N
    M1y = np.sum(np.sin(X))/N
    M2x = np.sum(np.cos(Y))/N
    M2y = np.sum(np.sin(Y))/N
    Ppx = np.sum(np.cos(X+Y))/N
    Ppy = np.sum(np.sin(X+Y))/N
    Pmx = np.sum(np.cos(X-Y))/N
    Pmy = np.sum(np.sin(X-Y))/N
    
    # Linearized variations of those quantities
    dM1x = np.sum(dX*-np.sin(X))/N
    dM1y = np.sum(dX*np.cos(X))/N
    dM2x = np.sum(dY*-np.sin(Y))/N
    dM2y = np.sum(dY*np.cos(Y))/N
    dPpx = np.sum(-np.sin(X+Y)*(dX+dY))/N
    dPpy = np.sum(np.cos(X+Y)*(dX+dY))/N
    dPmx = np.sum(-np.sin(X-Y)*(dX-dY))/N
    dPmy = np.sum(np.cos(X-Y)*(dX-dY))/N
    
    # Linearized forces (Jacobian action)
    dFx = c*((dM1y*np.cos(X)-M1y*np.sin(X)*dX) -
             (dM1x*np.sin(X)+M1x*np.cos(X)*dX)) + \
          0.5*d*((dPpy*np.cos(X+Y)-Ppy*np.sin(X+Y)*dX) -
                 (dPpx*np.sin(X+Y)+Ppx*np.cos(X+Y)*dX) +
                 (dPmy*cos(X-Y)-Pmy*sin(X-Y)*dX) -
                 (dPmx*np.sin(X-Y)+Pmx*np.cos(X-Y)*dX))

    dFy = c*((dM2y*np.cos(Y)-M2y*np.sin(Y)*dY) -
             (dM2x*np.sin(Y)+M2x*np.cos(Y)*dY)) + \
          0.5*d*((dPpy*np.cos(X+Y)-Ppy*np.sin(X+Y)*dY) -
                 (dPpx*np.sin(X+Y)+Ppx*np.cos(X+Y)*dY) -
                 (dPmy*cos(X-Y)+Pmy*sin(X-Y)*dY) +
                 (dPmx*np.sin(X-Y)-Pmx*np.cos(X-Y)*dY))
    
    # Update tangent momenta
    dPx += dFx*delta
    dPy += dFy*delta
    
    return(dPx,dPy)


# =============================
# SUZUKI-TROTTER INTEGRATOR
# =============================

def Suzuki4_linear(X,Y,dX,dY,dPx,dPy):
    # 4th order integrator for tangent dynamics
    evolvedQ(dX,dY,dPx,dPy,delta2)
    evolvedP(X,Y,dX,dY,dPx,dPy,delta1)
    evolvedQ(dX,dY,dPx,dPy,delta1)
    evolvedP(X,Y,dX,dY,dPx,dPy,delta1)
    evolvedQ(dX,dY,dPx,dPy,delta3)
    evolvedP(X,Y,dX,dY,dPx,dPy,delta4)
    evolvedQ(dX,dY,dPx,dPy,delta3)
    evolvedP(X,Y,dX,dY,dPx,dPy,delta1)
    evolvedQ(dX,dY,dPx,dPy,delta1)
    evolvedP(X,Y,dX,dY,dPx,dPy,delta1)
    evolvedQ(dX,dY,dPx,dPy,delta2)
    
    return(dX,dY,dPx,dPy)
    

def Suzuki4(X,Y,PX,PY):
    # 4th order symplectic integrator for original system
    evolveQ(X,Y,PX,PY,delta2)
    evolveP(X,Y,PX,PY,delta1)
    evolveQ(X,Y,PX,PY,delta1)
    evolveP(X,Y,PX,PY,delta1)
    evolveQ(X,Y,PX,PY,delta3)
    evolveP(X,Y,PX,PY,delta4)
    evolveQ(X,Y,PX,PY,delta3)
    evolveP(X,Y,PX,PY,delta1)
    evolveQ(X,Y,PX,PY,delta1)
    evolveP(X,Y,PX,PY,delta1)
    evolveQ(X,Y,PX,PY,delta2)
    
    return(X,Y,PX,PY)


# =============================
# ENERGY FUNCTIONS
# =============================

    
def K(PX,PY):
    # Kinetic energy per particle
    return 0.5*np.sum(np.square(PX)+np.square(PY))/N

def V(X,Y):
    # Potential energy (mean-field interactions)
    M1x = np.sum(np.cos(X))/N
    M1y = np.sum(np.sin(X))/N
    M2x = np.sum(np.cos(Y))/N
    M2y = np.sum(np.sin(Y))/N
    Ppx = np.sum(np.cos(X+Y))/N
    Ppy = np.sum(np.sin(X+Y))/N
    Pmx = np.sum(np.cos(X-Y))/N
    Pmy = np.sum(np.sin(X-Y))/N
    
    
    M1 = np.sqrt(np.square(M1x)+np.square(M1y))
    M2 = np.sqrt(np.square(M2x)+np.square(M2y))
    
    Pp = np.sqrt(np.square(Ppx)+np.square(Ppy))
    Pm = np.sqrt(np.square(Pmx)+np.square(Pmy))
    
    #return (2-((M1x*np.cos(X))+(M1y*np.sin(X)))-((M2x*np.cos(Y))+(M2y*np.sin(Y)))-0.5*(2-((Ppx*np.cos(X+Y))+(Ppy*np.sin(X+Y)))-((Pmx*np.cos(X-Y))+(Pmy*np.sin(X-Y)))))
    return 0.5*c*(2-np.square(M1)-np.square(M2)) + 0.25*d*(2-np.square(Pp)-np.square(Pm))
    
def H(K,V):
    return (K+V)


# =============================
# INITIAL CONDITION: WATERBAG
# =============================

def waterbag(X,Y,PX,PY):
    
    # Generates uniform initial distribution in phase space

    deltax = deltay = 2*n*np.pi
    deltapy = 1e-5
    
    M1x = (np.sin(deltax)/deltax)
    M1y = ((1.-np.cos(deltax))/deltax)
    M1_2 = np.square(M1x) + np.square(M1y)
    
    M2x = (np.sin(deltay)/deltay)
    M2y = ((1.-np.cos(deltay))/deltay)
    M2_2 = np.square(M2x) + np.square(M2y)
    
    a = deltax*deltay
    
    P_plusx = (1./a)*(np.sin(deltax)*np.sin(deltay) - (1.-np.cos(deltax))*(1.-np.cos(deltay)))
    P_plusy = (1./a)*((1.-np.cos(deltax))*np.sin(deltay)+(1.-np.cos(deltay))*np.sin(deltax))
    P_plus2 = np.square(P_plusx) + np.square(P_plusy)
    
    P_minusx = (1./a)*(np.sin(deltax)*np.sin(deltay) + (1.-np.cos(deltax))*(1.-np.cos(deltay)))
    P_minusy = (1./a)*((1.-np.cos(deltax))*np.sin(deltay) - (1.-np.cos(deltay))*np.sin(deltax))
    P_minus2 = np.square(P_minusx) + np.square(P_minusy)
    
    U1 = 0.5*c*(2. - M1_2 - M2_2)
    U2 = 0.25*d*(2. - P_plus2 - P_minus2)
    
    alpha = 1./24.*np.square(deltapy) + U1 + U2

    deltapx = np.sqrt(24.*(E_in - alpha))
    
    for i in range(N):
        X = np.random.rand(N)*(deltax)
        Y = np.random.rand(N)*(deltay)
        PX = (2.*np.random.rand(N) - np.ones(N))*(0.5*deltapx)
        PY = (2.*np.random.rand(N) - np.ones(N))*(0.5*deltapy)


    return (X,Y,PX,PY)



# Initialize system
x,y,px,py = waterbag(x,y,px,py)

# =============================
# TIME EVOLUTION + LYAPUNOV
# =============================

time = np.arange(0,100,h)

L = []     # Lyapunov exponent evolution
exp = []   # Final values for each run
lamb = 0   # Accumulated logarithm


for k in range(len(w)):

  print(k)

  for i in range(1,len(time)+1):

    #print(i)

    # Evolve system
    x,y,px,py = Suzuki4(x, y, px, py)
    
    # Evolve tangent vector
    w[0:N],w[N:2*N],w[2*N:3*N],w[3*N:4*N] = Suzuki4_linear(x, y, w[0:N], w[N:2*N], w[2*N:3*N], w[3*N:4*N])

    # Accumulate logarithmic growth (Lyapunov)
    lamb += log(np.linalg.norm(w))
    
    LE = lamb/(float(i)*h) # Finite-time Lyapunov exponent

    L.append(LE)

    # Renormalize tangent vector to avoid overflow
    w = (1/(np.sqrt(np.dot(w,w))))*w
    
    ##periodic boundary conditions
    for jj in range(N):
        if( np.abs(x[jj]) > 2*n*np.pi):
            x[jj] -= np.sign(x[jj])*2.*n*np.pi
        if( np.abs(y[jj]) > 2*n*np.pi):
            y[jj] -= np.sign(y[jj])*2.*n*np.pi

# Catch the last LLE in time series
last = L[-1]
exp.append(last)

