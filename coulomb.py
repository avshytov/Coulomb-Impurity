from numpy import * 
from pylab import * 
from scipy import special


def coulombkernel(r): 
    N = len(r)
    M = zeros((N,N))
    for i in range (0,N):
        r_i = r[i]
        print "Coulomb kernel:", i+1, "/", N
        for j in range (0,N):
            r_j = r[j]
            if i == j:
                if (j < N - 1):
  	            Dr = 0.5 * (r[j + 1] - r_j)
	        else:
		    Dr = 0.5 * (r_j - r[j - 1])
                    Ldr    =  Dr * math.log(1.0 / Dr**2)
                    Lconst =  2.0 * Dr * math.log(4.0)
                    Lplus  =  2.0 * Dr

                    M[i,j] = 0.5 * (Ldr + Lplus + Lconst)
            else:
                if j == 0:
		     a = r_j
                     b = 0.5 * (r_j + r[j + 1])
                elif j == (N - 1):
                    a = 0.5 * (r_j + r[j - 1])
                    b = r_j
                else:
                    a = 0.5 * (r_j + r[j - 1])
                    b = 0.5 * (r_j + r[j + 1])
                d = b - a
                mu_top =  (4 * r_i * b) / ((r_i + b)**2)
                mu_bot =  (4 * r_i * a) / ((r_i + a)**2)
                ellip_top = special.ellipk(mu_top)
                ellip_bot = special.ellipk(mu_bot)
                alpha_top = b / (r_i + b)
                alpha_bot = a / (r_i + a)
                I_top = ellip_top * alpha_top
                I_bot = ellip_bot * alpha_bot
                # i have added the constant of 4 that was previously forgotten
                M[i,j] = 4.0 * 0.5 * (I_top + I_bot) * d
    return M


if __name__ == '__main__':
    N = 1000
    r_min = 0.001
    r_max = 25.0
    r_0 = 1.0

    r = zeros((N))  

    ivals = array(range(0, N))
    spaceri = math.log(r_max/r_min) * (1.0/(N - 1.0))
    r = r_min * exp(spaceri * ivals)

    M = coulombkernel(r)
    rho_j = 1.0 / (1.0 + r**2)**1.5				
    phi_i = dot(M, rho_j) 
    U = (2 * np.pi) / ( 1.0 + r**2)**0.5 

    plot(r, phi_i)
    plot(r, U)
    savefig("image1.pdf")

    figure()
    loglog(r, phi_i)
    loglog(r, U)
    savefig("image2.pdf")

    figure()
    plot(r, (phi_i / U))
    savefig("image3.pdf")

    show()




