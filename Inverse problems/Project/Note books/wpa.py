from numpy import *

# -------------------------------------
# Wave propagation algorithm
# --------------------------------------

def limiter(r,lim_choice='MC'):
    wlimiter = empty(r.shape)
    if lim_choice == None:
        wlimiter = ones(r.shape)
    elif lim_choice == 'minmod':        
        # wlimitr = dmax1(0.d0, dmin1(1.d0, r))
        wlimiter = maximum(0,minimum(1,r))
    elif lim_choice == 'superbee':
        # wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
        a1 = minimum(1,2*r)
        a2 = minimum(2,r)        
        wlimiter = maximum(0,maximum(a1,a2))
    if lim_choice == 'MC':
        # c = (1.d0 + r)/2.d0
        # wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
        c = (1 + r)/2
        wlimiter = maximum(0,minimum(c,minimum(2,2*r)))
    elif lim_choice == 'vanleer':
        # wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
        wlimiter = (r + abs(r))/(1 + abs(r))
            
    return wlimiter


def evolve_q(mx,ax,bx,dt, dx, q0, F, bc, tv, lim_choice=None, uvel=1):

    # Velocities defined at interfaces
    up = max(array([uvel,0]))
    um = min(array([uvel,0]))    

    M = len(tv) - 1

    t = 0
    Q = zeros((mx,M+1))
    Q[:,0] = q0

    second_order = True

    q = q0.copy()
    alpha = uvel*dt/dx
    for n in range(0,M):
        t = tv[n]

        # Apply boundary conditions
         #q_ext = concatenate(([0,0], q,[0,0]))
        q_ext = bc(q)
        waves_x = q_ext[1:] - q_ext[:-1]

        # Compute fluctuations
        apdq = up*waves_x[1:-1]
        amdq = um*waves_x[1:-1]
        if second_order:    
            # waves_ext = concatenate(([waves_x[0]],waves_x,[waves_x[-1]]))
            wl = waves_x[:-2]*waves_x[1:-1]
            wr = waves_x[2:]*waves_x[1:-1]
            w2 = waves_x[1:-1]*waves_x[1:-1]
            m = w2 > 0
            r = ones(w2.shape)

            if uvel > 0:
                r[m] = wl[m]/w2[m]
            else:
                r[m] = wr[m]/w2[m]

            wlimiter = limiter(r,lim_choice)

            cxx_coeff = 0.5*abs(uvel)*(1 - abs(uvel)*dt/dx)    
            cxx = cxx_coeff*waves_x[1:-1]*wlimiter
            apdq -= cxx
            amdq += cxx

            # Update the solution and include source term F
            q = q - dt/dx*(apdq[:-1] + amdq[1:]) +  dt*F[:,n]
        
        Q[:,n+1] = q
        
    return Q
