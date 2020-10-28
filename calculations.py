from Beamer import *

if __name__ == '__main__':
    with openBeamer('calculations', 'calculations') as doc:
        x, t = Expr('x'), Expr('t')
        vec_x = bf(x)
        
        rho = Expr(r'\rho')
        phi = Expr(r'\phi')
        
        W = Expr(r'W')
        Y = Expr(r'Y')
        r = Expr(r'r')
        vec_r = bf(r)
        rdot = dot(r)
        vec_rdot = dot(bf(r))

        v = Expr('v')
        vec_v = bf(v)
        
        i, j = Expr(r'i'), Expr(r'j')
        a, b, c = Expr(r'a'), Expr(r'b'), Expr(r'c')
        
        f = Expr(r'f')
        var_e = Expr(r'\varepsilon ')

        p = Expr(r'p')
        vec_p = bf(p)
        P = Expr(r'P')

        L = Expr(r'L')
        H = Expr(r'H')
        U = Expr(r'U')
        phi = Expr(r'\phi ')
        
        delta = Expr(r'\delta ')
        
        lambda_ = Expr(r'\lambda ')
        
        h = Expr(r'h')
        dt = d(t)
        d_dt = lambda x: frac( d(x), dt )
        ddt = d_dt('')
        
        N = Expr(r'N')
        m = Expr(r'm')
        varphi = Expr(r'\varphi ')


        #doc.frame('test',
            #'double eq test',
            #Math( 
                #a == b == c
            #),
            #Math(
                #a == b and b == c
            #),
        #)
        #exit(0)
        
        doc.frame('kernel',
            'using a general kernel function',
            Math( 
                Int[vec_x]*W( vec_x - vec_r[a](t)) == 1 
            ),
            'in order to rewrite',
            Math( 
                Int[vec_x]*f( vec_x, t ) == Sum[a]*( f[a](t)/rho[a](t) ) 
            ),
            'reference density',
            Math( 
                rho( vec_x ) == Sum[a]*W( vec_x - vec_r[a]) 
            ),
            'which gives',
            Math( 
                Int[vec_x]*rho( vec_x ) == N 
            ),
        )

        doc.frame('derivatives',
            'useful computations',
            Math( rho( r[a] ) == Sum[b]*W( vec_r[a] - vec_r[b]) ),
            'the comoving derivative is',
            Math( 
                partial[vec_r[a] ]*rho[b] == (
                    ( Sum[c] * partial[vec_r[a]] * W[b*c] == Sum[c] * bf(Y)[b*c] * (delta[b*a] - delta[c*a]) ) ==
                    bf(Y)[a*b] + delta[b*a] * Sum[c] * bf(Y)[a*c]
                    )
            ),
            'where',
            Math(
                bf(Y)[b*c] == - bf(Y)[c*b]
            ),
            'a common form of term is',
            Math(
                Sum[b]* f[b] * partial[vec_r[a] ]*rho[b] == (
                    f[a] * Sum[c] * bf(Y)[a*c] + Sum[b]*f[b]* bf(Y)[a*b] == Sum[b] *( f[a] + f[b] ) * bf(Y)[a*b]
                    ),
            ),
            'which has the nice property that',
            Math(
                Sum[a] * Sum[b]* f[b] * partial[vec_r[a] ]*rho[b] == 0
            ),
        )
        
        doc.frame('lagrangean',
            'starting from the non-relativistic case',
            Math( 
                cal(L) == frac(1,2)* rho(vec_x)*vec_v(vec_x)**2 - var_e(vec_x)
            ),
            'the canonical momentum density is',
            Math( 
                cal(P) == ( 
                    frac( partial*cal(L), partial*vec_v ) == rho*vec_v
                    )
            ),            
            'the hamiltonian density is',
            Math( 
                cal(H) == ( 
                    vec_v * frac( partial*cal(L), partial*vec_v ) - cal(L) == frac(vec_p**2, 2*rho) + var_e
                    )
            ),
            'the discretization needs to be consistently made in both scenarios',
        )
            
        doc.frame('equations of motion',
            'in the lagrangean side, a straight forward paramaterization is',
            Math(
                L(vec_r, vec_rdot) == (
                    frac(1,2) * Int[vec_x]( rho(x)* v(x)**2 ) - Int[vec_x] * var_e(x),
                    frac(1,2) * Sum[a] * vec_v[a]**2 - Sum[a] * frac(var_e[a],rho[a]),
                    )
            ),
            Math(
                vec_p[a] == Sum[b]*v[b]*W[a*b]
            ),            
            'on the other side',
            Math(
                H(vec_r, vec_p) == (
                    frac(1,2) * Int[vec_x] * frac( p(x)**2, rho(x) ) + Int[vec_x] * var_e(x),
                    frac(1,2) * Sum[a] * frac( vec_p[a]**2, rho[a]**2) + Sum[a] * frac(var_e[a],rho[a]),
                    )
            ),
        )

        doc.frame('equations of motion',
            Math(
                vec_p[a] == Sum[b]*vec_rdot[b]*W[a*b]
            ),            
            Math(
                matrix([[partial[vec_rdot[a]]],[partial[vec_r[a]]]]) == 
                matrix([ [ partial[vec_rdot[a]]*p[b], 0 ], [ partial[vec_r[a]]*p[b], delta[a*b] ]]) 
                * matrix( [ [partial[vec_p[b]]], [partial[vec_r[b]]] ] )
            ),
        )

        doc.frame('equations of motions',
            'and the equation of motion is',
            Math( 
                frac( d(''), d(t) )*partial[vec_rdot[a] ]* L == partial[ vec_r[a] ]* L
            ),
            Math(
                (
                    frac( d( vec_r[a] ), d(t) ) == partial[vec_p[a] ]* H,
                    frac( d( vec_p[a] ), d(t) ) == - partial[vec_r[a] ]* H,
                )
            ),
            'the total derivative of the hmailtoninan is',
            Math(
                frac( d(''), d(t) )*H(vec_r, vec_p) == Sum[a] * frac( d( vec_p[a] ), d(t) ) * partial[vec_p[a] ]*H + Sum[a] * frac( d( vec_r[a] ), d(t) )*partial[vec_r[a] ]* H
            ),
            'the angular momentum is',
            Math(
                frac( d(''), d(t) )*Sum[a] * cross( vec_r, vec_p ) == Sum[a] * cross( frac( d(vec_r[a]), d(t) ), vec_p[a] ) + Sum[a] * cross( vec_r[a], frac( d(vec_p[a]), d(t) ) )
            ),            
        )
        
        F = Expr(r'F')
        doc.frame('equations of motions',
            'and the equation of motion is',
            Math(
                frac( d( vec_r[a] ), d(t) ) == frac(vec_p[a], rho[a]**2),
            ),
            Math((
                frac( d( vec_p[a] ), d(t) ) == - Sum[b]*( F[a] + F[b] )*bf(Y)[a*b],
                F[a] == frac( partial[rho](var_e)[a], rho[a] ) - var_e[a]/rho**2 - vec_p[a]**2/rho[a]**3
            )),
            'the angular momentum is',
            Math(
            Sum[a] * cross( vec_r, frac( d(vec_p), d(t) ) ) == (
                Sum[a] * Sum[b]*( F[a] + F[b] )*cross( vec_r[a], bf(Y)[a*b] ),
                #frac(1,2)*Sum[a*b]*( ( F[a] + F[b] )*cross( vec_r[a], bf(Y)[a*b] ) + ( F[b] + F[a] )*cross( vec_r[b], bf(Y)[b*a] )),
                frac(1,2)*Sum[a*b]*( F[a] + F[b] )*( cross( vec_r[a], bf(Y)[a*b] ) - cross( vec_r[b], bf(Y)[a*b] )),
                frac(1,2)*Sum[a*b]*( F[a] + F[b] )*( cross( vec_r[a] - vec_r[b], bf(Y)[a*b] ) ),
                )
            ),
        )
        
        M = Expr(r'M')
        doc.frame('continuity',
            Math(
                Int[vec_x] * rho(x) == Sum[a]* 1
            ),
            Math(
                Int[vec_x] * rho(x)*vec_v(x) == Sum[a] * v[a]
            ),
            'this gives the canonical momentum as',
            Math( 
                vec_p[a] == ( partial[vec_rdot[a] ]* L == vec_rdot[a] )
            ),
            'and the equation of motion is',
            Math( 
                frac( d(vec_p[a]), d(t) ) == ( 
                    partial[ vec_r[a] ]* L == - Sum[b] * (phi[a] + phi[b]) * bf(Y)[a*b]
                )
            ),
            Math(
                phi[a] == ( ( (partial[rho]*var_e[a])/rho[a] - var_e[a]/rho[a]**2 ) == frac( d(''), d(rho) )*frac(var_e,rho) )
            ),
        )
        
        doc.frame('reference density',
            'the hamiltonian density is',
            Math(
                cal(H) == frac(1,2) * frac( p(x)**2, rho(x) ) + var_e(x),
            ),
            'using a generic reference density',
            Math(
                phi[a] == Sum[b]*W[a*b]
            ),
            'the integral of the hamiltonian is',
            Math(
                H(vec_r, vec_p) == Sum[a]* frac(1,2) * frac( vec_p[a]**2, phi[a] * rho[a] ) +  Sum[a] * (var_e[a]/phi[a])
            ),
            'from the conservation of energy',
            Math(
                Sum[a]* frac(1,2) * ddt* frac( vec_p[a]**2, phi[a] * rho[a] ) == - Sum[a] * ddt*(var_e[a]/phi[a])
            ),            
        )

        doc.frame('reference density',
            'the velocity equation is',
            Math(
                (
                    d_dt(vec_r[a]) == frac( vec_p[a], rho[a]*phi[a] ),
                    d_dt( phi[a] ) == Sum[b]*( frac( vec_p[a], rho[a]*phi[a] ) - frac( vec_p[b], rho[b]*phi[b] ) )*cdot('',bf(Y)[a*b])
                )
            ),
            'from the conservation of energy',
            Math(
                cdot(vec_p[a], d_dt(vec_p[a]) ) == - phi[a] * rho[a] *frac(vec_p[a]**2,2) * ddt*( frac(1,(phi[a] * rho[a]) )) - phi[a] * rho[a] * ddt*(var_e[a]/phi[a])
            ),            
            Math(
                cdot(vec_p[a], d_dt(vec_p[a]) ) == frac(vec_p[a]**2,2) * ( frac(1,rho[a]) * d_dt(rho[a]) + frac(1, phi[a])*d_dt( phi[a] ) ) - rho[a]*d_dt(var_e[a]) + frac(rho[a]*var_e[a], phi[a])*d_dt(phi[a])
            ),            
            Math(
                cdot(vec_p[a], d_dt(vec_p[a]) ) == frac(vec_p[a]**2,2*rho[a]) * d_dt(rho[a]) - rho[a]*d_dt(var_e[a]) + ( frac(vec_p[a]**2,2) + rho[a]*var_e[a] ) * frac(1, phi[a]) * d_dt(phi[a])
            ),            
        )

        doc.frame('reference density',
            Math(
                ddt*Int[vec_x] * rho(vec_x) == ( ddt*Sum[a]*frac( rho[a], phi[a] ) == 0 )
            ),
            Math(
                Sum[a]*frac( 1, phi[a] )*d_dt(rho[a]) == Sum[a]*frac( rho[a], phi[a]**2 )*d_dt(phi[a])
            ),
            Math(
                d_dt(rho[a]) == frac( rho[a], phi[a] )*d_dt(phi[a])
            ),
            Math(
                cdot(vec_p[a], d_dt(vec_p[a]) ) == - rho[a]*d_dt(var_e[a]) + frac(1, phi[a]) *( vec_p[a]**2 + rho[a]*var_e[a] ) * d_dt(phi[a])
            ),            
            Math(
                d_dt(var_e[a]) == - frac(vec_p[a],rho[a])*cdot('', d_dt(vec_p[a]) ) + ( frac(vec_p[a]**2, phi[a]*rho[a]) + var_e[a]/phi[a] ) * d_dt(phi[a])
            ),
            'from the hamiltonian density',
            Math(
                d_dt(var_e(vec_x)) == - frac(vec_p(vec_x),rho(vec_x) )*cdot('', d_dt( vec_p(vec_x) ) ) + ( frac(vec_p(vec_x)**2, 2*rho(vec_x)**2) ) * d_dt( rho(vec_x) )
            ),            
        )

        doc.frame('reference density',
            'using',
            Math(
                phi[a] == var_e[a]
            ),
            'the equations reduce to',
            Math(
                d_dt(var_e[a]) == - frac(vec_p[a],rho[a])*cdot('', d_dt(vec_p[a]) ) + ( frac(vec_p[a]**2, var_e[a]*rho[a]) + 1 ) * d_dt(var_e[a])
            ),
            Math(
                 frac(vec_p[a]**2, var_e[a]*rho[a]) * d_dt(var_e[a]) == frac(vec_p[a],rho[a])*cdot('', d_dt(vec_p[a]) )
            ),
            Math(
                 d_dt(var_e[a]) == var_e[a] * frac(vec_p[a], vec_p[a]**2)*cdot('', d_dt(vec_p[a]) )
            ),
            'from the hamiltonian',
            Math(
                H(vec_r, vec_p) == 
                    frac(1,2) * Sum[a] * frac( vec_p[a]**2, var_e[a]*rho[a]) + Sum[a] * 1
            ),
        )

        doc.frame('reference density',
            Math(
                d_dt(vec_p[a]) == (
                    - frac(1,2) * Sum[b] *vec_p[b]**2* partial[vec_r[a]]*frac( 1, var_e[b]*rho[b]),
                    frac(1,2) * Sum[b] *vec_p[b]**2* frac( 1, var_e[b]*rho[b]**2) * partial[vec_r[a]]*rho[b] + frac(1,2) * Sum[b] *vec_p[b]**2*frac( 1, var_e[b]**2*rho[b])* partial[vec_r[a]]*var_e[b],
                    frac(1,2) * Sum[b] *vec_p[b]**2* frac( 1, var_e[b]**2*rho[b]**2)( var_e[b]* partial[var_e]*rho[b]  + rho[b] )* partial[vec_r[a]]*var_e[b],
                )
            ),
        )

        T = Expr(r'T')
        u = Expr(r'u')
        g = Expr(r'g')
        mu = Expr(r'\mu ')
        nu = Expr(r'\nu ')
        P = Expr(r'P')
        J = Expr(r'J')
        tau = Expr(r'\tau ')
        gamma = Expr(r'\gamma ')
        S = Expr(r'S')
        ast = ('', '*')
        dtau = d(tau)
        vec_u = bf(u)
        omega = Expr(r'\omega ')
        
        doc.frame('relativistic case',
            'the invariant length',
            Math(
                dtau**2 == g[mu*nu] * d(x['',mu]) * d(x['',nu])
            ),
            Math(
                dtau**2 == dt**2 - cdot( d(vec_x), d(vec_x) )
            ),
            Math(
                frac(dt**2, dtau**2) - cdot( frac(d( vec_x), dtau), frac( d(vec_x), dtau ) ) == ( u[mu]*u['',mu] == 1 )
            ),
            'where',
            Math(
                u[0] == ( gamma == frac(dt, dtau) )
            ),
            Math(
                vec_u == ( gamma*vec_v == ( frac(d(vec_x), dtau) == frac(dt, dtau) * frac(d( vec_x), dt) ))
            ),
        )
        
        doc.frame('relativistic case',
            'from the energy-momentum tensor',
            Math(
                T['',mu*nu] == ( omega*u['',mu] * u['',nu] - P*g['',mu*nu] == var_e * u['',mu] * u['',nu] + P*(u['',mu] * u['',nu] - g['',mu*nu] ) )
            ),
            'which comes from the lagrangean as',
            Math(
                T['',mu*nu] == u['',mu] * partial[ u[nu] ] * cal(L) - g['',mu*nu]*cal(L)
            ),
            'from that we identify that the canonical momentum is', 
            Math(
                partial[ u[nu] ] * cal(L) == omega*u['',nu]
            ),
            'and the lagrangean density',
            Math(
                cal(L) == P
            ),
            'in order to make these two conclusions compatible we make use of the 4-velocity constraint in the lagrangean',
            Math(
                cal(L) == ( P + omega*(u['',mu]*u[mu] - 1) == -var_e + omega*u['',mu]*u[mu] )
            ),
        )

        doc.frame('relativistic case',
            'the hamiltonian expressed in terms of the canonical momentum',
            Math((
                    p['',mu] == omega*u['',mu],
                    vec_p == omega*gamma*vec_v
                )
            ),
            'it is important to rewrite the $\gamma$ in terms of $p$',
            Math(
                (
                    vec_p**2/omega**2 == ( gamma**2*vec_v**2 == gamma**2 - 1 ),
                    gamma**2 == vec_p**2/omega**2 + 1
                )
            ),
            'from that we identify that the canonical momentum is', 
            Math(
                T['',mu*nu] == p['',mu] * p['',nu]/omega - P*g['',mu*nu]
            ),
            'the hamiltonian density is',
            Math(
                cal(H) == ( T[t*t] == ( p[t]**2/omega - P == vec_p**2/omega + var_e ) )
            ),
        )
        
        doc.frame('energy flux',
            'the energy flux is', 
            Math(
                T['',i*t] == p['',t] * p['',i]/omega == gamma*vec_p
            ),
            'the equations of motion are',
            Math(
                partial[mu]*T['',mu*nu] == 0
            ),
            'which integrates to',
            Math(
                Int[vec_x]*partial[mu]*T['',mu*nu] == ( frac(d(''),dt) * Int[vec_x]*T['',t*nu] == 0 )
            ),
            'the space part translates into',
            Math(
                frac(d(''),dt) * Int[vec_x]*omega*gamma**2*vec_v == ( frac(d(''),dt) * Int[vec_x] * sqrt( vec_p**2/omega**2 + 1 ) * vec_p == 0 )
            ),
        )
        
        doc.frame('energy flux',
            'the time part translates into',
            Math(
                frac(d(''),dt) * Int[vec_x]* (omega*gamma**2 - P) == ( frac(d(''),dt) * Int[vec_x] * (vec_p**2/omega + var_e) == 0 )
            ),
            'by imposing that any quantity can be integrated to',
            Math(
                frac(d(''),dt) * Int[vec_x]* f(vec_x) == frac(d(''),dt) * Sum[a]*(f[a]/h[a])
            ),
            'this allows us to define the reference density based on the hamiltonian density so that',            
            Math(
                h[a] == ( omega[a]*gamma[a]**2 - P[a] == vec_p[a]**2/omega[a] + var_e[a] )
            ),
            'therefore',
            Math(
                h[a] == Sum[b] * W[a*b]
            ),
        )
        
        doc.frame('integrals',
            'given a local density its integral',
            Math(
                cal(S) == Int* d(tau)*d(vec_x) * cal(L)
            ),
            'which can be converted to the common frame',
            Math(
                cal(S) == ( Int* d(t)*d(vec_x) * frac(d(tau),d(t)) * cal(L) == Int* d(t)*d(vec_x) * cal(L)[ast] )
            ),
            'where',
            Math(
                frac(dt, d(tau) ) == ( gamma == frac( 1, sqrt(1-vec_v**2) ) )
            ),
            'so the lagrangean density on the common frame and the canonical momentum',
            Math(
                (
                    cal(L)[ast] == -(var_e/gamma),
                    frac( partial*cal(L)['','*'], partial*vec_v ) == var_e/gamma**2*frac( partial*gamma, partial*vec_v )
                )
            ),
            'the gamma derivative is',
            Math(
                frac( partial*gamma, partial*vec_v ) ==  gamma**3 * vec_v
            ),
        )
        
        doc.frame('relativistic case',
            'canonical momentum',
            Math(
                vec_p == ( frac( partial*cal(L)[ast], partial*vec_v ) == var_e* gamma * vec_v )
            ),
            Math(
                ( 
                    frac(vec_p**2, var_e**2 ) == ( gamma**2 * vec_v**2 == gamma**2 - 1 ),
                    gamma == sqrt(frac(vec_p**2, var_e**2 ) + 1),
                )
            ),
            'the hamiltonian density in the common frame is',
            Math(
                cal(H)[ast] == ( cdot(vec_v,vec_p) - cal(L)[ast] == frac(vec_p**2, gamma*var_e) + var_e/gamma) 
            ),
            'the hamiltonian density in the local frame is',
            Math(
                cal(H) == frac(vec_p**2, var_e) + var_e
            ),
        )

        doc.frame('reference density',
            
        )
        
        doc.frame('sph',
            'the SPH integrals are defined in the common frame so',
            Math(
                Int[vec_x]* frac(f(vec_x),gamma(vec_x) ) == Sum[a]*( f[a]/(gamma[a]*rho[a]) )
            ),
            'so the common frame lagrangean is',
            Math(
                L[ast] == - Sum[a]*( var_e[a]/(gamma[a]*rho[a]) )
            ),
            'so the common frame hamiltonian is',
            Math(
                H[ast] == Sum[a]*( vec_p[a]**2/(gamma[a]*var_e[a]*rho[a]) ) + Sum[a]*( var_e[a]/(gamma[a]*rho[a]) )
            ),
        )

        doc.frame('eqs of motion',
            'the partial derivative in respect to the canonical momentum is',
            Math(
                partial[vec_p[a]]*(1/gamma[b]) == -frac(vec_p[a], gamma[a]**3 * var_e[a]**2)*delta[a*b]
            ),
            'the partial derivative in respect to the canonical momentum is',
            Math(
                frac( d(vec_r[a]), dt ) == (
                    partial[ vec_p[a] ]*H[ast],
                    2*vec_p[a]/(gamma[a]*var_e[a]*rho[a]) - vec_p[a]**2*vec_p[a]/(gamma[b]**3*var_e[b]**3*rho[b]) - vec_p[a]/(gamma[a]**3*var_e[a]*rho[a]),
                    ( 2*gamma[a]**2 - (gamma[a]**2-1) - 1 ) * frac( vec_p[a], gamma[a]**3*var_e[a] * rho[a] ),
                    frac( vec_p[a], gamma[a]*var_e[a] * rho[a] ),
                    )
            ),
        )
        
        doc.frame('kernel revisited',
            'given a normalized kernel',
            Math( Int[vec_x]*W( vec_x - vec_r[a](t)) == 1 ),
            'the relation is imposed',
            Math( Int[vec_x]*f( vec_x, t ) == Sum[a]*( F[a](t)/rho( vec_r[a](t), t ) ) ),
            'where $F$ are free weights to ensure the relationship.',
            Math(
                Int[vec_x]*W( vec_x - vec_r[a](t)) == (
                    Sum[b]*( W( vec_r[b](t) - vec_r[a](t) )/rho( vec_r[b](t), t ) ) == 1 )
            ),
            'the explicity dependence on the index $a$ has to be taken into account',
            Math(
                W(0)/rho( vec_r[a](t), t ) + Sum[b]*( W( vec_r[b](t) - vec_r[a](t) )/rho( vec_r[b](t), t ) ) == 1
            ),            
        )
        doc.frame('kernel revisited',
            Math(
                rho( vec_r[a](t), t ) == W(0)/(1 - Sum[b]*( W( vec_r[b](t) - vec_r[a](t) )/rho( vec_r[b](t), t ) ))
            ),            
            'the usual reference density can be derived by summing the kernel integral for all fluid elements',
            Math(
                Sum[a]*Sum[b]*( W( vec_r[b](t) - vec_r[a](t) )/rho(vec_r[b](t), t) ) == N
            ),
            'the commuatation of the summations yeilds',
            Math(
                Sum[b]*frac(1, rho(vec_r[b](t), t) )*( Sum[a]* W( vec_r[b](t) - vec_r[a](t)) ) == N
            ),
            'therefore',
            Math(
                rho(vec_r[b](t), t) == Sum[a]* W( vec_r[b](t) - vec_r[a](t))
            ),
        )
        
        doc.frame('kernel revisited',
            'the same proceedure can be done with constant weights',
            Math(
                Sum[a]* var_e[a] *Sum[b]*( W( vec_r[b](t) - vec_r[a](t))/rho(vec_r[b](t), t) ) == Sum[a]* var_e[a]
            ),
            Math(
                Sum[b]*frac(1, rho(vec_r[b](t),t) )*( Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) ) ) == Sum[a]* var_e[a]
            ),
            Math(
                rho(vec_r[b](t), t) == frac(1, var_e[b]) * Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) )
            ),
            Math(
                var_e[b]*rho(vec_r[b](t), t) == Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) )
            ),
        )

        doc.frame('conserved quantity',
            'the same proceedure can be done with constant weights',
            Math(
                ddt* Sum[a]* var_e[a] *Sum[b]*( W( vec_r[b](t) - vec_r[a](t))/rho(vec_r[b](t), t) ) == ( ddt*( Sum[a]* var_e[a] ) == 0 )
            ),
            Math(
                ddt* Sum[b]*frac(1, rho(vec_r[b](t),t) )*( Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) ) ) == Sum[a]* d_dt( var_e[a] )
            ),
            Math(
                rho(vec_r[b](t), t) == frac(1, var_e[b]) * Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) )
            ),
            Math(
                var_e[b]*rho(vec_r[b](t), t) == Sum[a]* var_e[a] * W( vec_r[b](t) - vec_r[a](t) )
            ),
        )


        doc.frame('attempt with fields',
            'a',
            Math(
                rho(vec_x) == Sum[a]*m[a](t) * W( vec_x - vec_r[a](t) )
            ),
            Math(
                varphi(vec_x) == Sum[a]*phi[a](t) * W( vec_x - vec_r[a](t) )
            ),
            'so',
            Math(
                vec_v(vec_x) == Sum[a]*phi[a](t) * partial[vec_x]*W( vec_x - vec_r[a](t) )
            ),
        )

        
        doc.frame('attempt with fields',
            'a',
            Math(
                delta*rho[a] == Sum[b]*delta*m[b] * W[a*b] + m[b] *(delta*r[a] - delta*r[b]) * bf(Y)[a*b]
            ),
            Math(
                delta*varphi[a] == Sum[b]*delta*phi[b] * W[a*b] + phi[b] *(delta*r[a] - delta*r[b]) * bf(Y)[a*b]
            ),
            'so',
            Math(
                cal(L)( rho, varphi ) == frac(1,2)*rho*(partial[i]*varphi)*(partial[i]*varphi) - var_e(rho)
            ),
        )

        doc.frame('attempt with fields',
            'lagrangean density using the fields',
            Math(
                cal(L)( rho, varphi ) == frac(1,2)*rho*(partial[i]*varphi)*(partial[i]*varphi) - var_e(rho)
            ),
            'variation',
            Math(
                delta*cal(L)( rho, varphi ) == partial[rho]*cal(L)*delta*rho + partial[partial[i]*varphi]*cal(L)*delta(partial[i]*varphi)
            ),
            'the canonical momentum is',
            Math(
                delta*cal(L)( rho, varphi ) == partial[rho]*cal(L)*delta*rho + hat(varphi)[i]*delta(partial[i]*varphi)
            ),            
            Math(
                delta*cal(L)( rho, varphi ) == partial[rho]*cal(L)*delta*rho + delta( hat(varphi)[i]*partial[i]*varphi ) - partial[i]*varphi * delta*hat(varphi)[i]
            ),            
            Math(
                delta(hat(varphi)[i]*partial[i]*varphi - cal(L)( rho, varphi )) == - partial[rho]*cal(L)*delta*rho + partial[i]*varphi * delta*hat(varphi)[i]
            ),            
        )
            
        doc.frame('attempt with fields',
            'hamiltonian density using the fields',
            Math(
                cal(H)( rho, varphi ) == frac(1,2)*rho*(partial[i]*varphi)*(partial[i]*varphi) + var_e(rho),
                cal(H)( rho, varphi ) == frac(1,2)*(hat(varphi)[i]*hat(varphi)[i]/rho) + var_e(rho),
            ),
            'the canonical momenta are',
            Math(
                delta*cal(H) == delta*rho * partial[rho]*cal(H) + delta*hat(varphi)[i] * partial[hat(varphi)[i]]*cal(H)
            ),
            r'since the only dependence are in $\rho$ and $\hat{\varphi}$',
            Math(
                delta*rho * partial[rho]*cal(H) + delta*hat(varphi)[i] * partial[hat(varphi)[i]]*cal(H) == - partial[rho]*cal(L)*delta*rho + partial[i]*varphi * delta*hat(varphi)[i]
            ),
            'the variations of the fields associate to',
            Math(
                ( partial[rho]*cal(H) == - partial[rho]*cal(L), partial[i]*varphi == partial[hat(varphi)[i]]*cal(H) )
            ),
        )
            
        doc.frame('attempt with fields',
            'from the Lagrangian equations of motion',
            Math(
                partial[i] * partial[partial[i]*rho]*cal(L) + partial[i] * partial[partial[i]*varphi]*cal(L) == - partial[rho]*cal(L) - partial[varphi]*cal(L)
            ),
            Math(
                ( 
                    partial[i] * partial[partial[i]*varphi]*cal(L) == - partial[rho]*cal(L), 
                    partial[i] * hat(varphi)[i] == - partial[rho]*cal(L), 
                )
            ),
            Math(
                ( partial[rho]*cal(H) == partial[i] * hat(varphi)[i], partial[i]*varphi == partial[hat(varphi)[i]]*cal(H) )
            ),            
        )

            
