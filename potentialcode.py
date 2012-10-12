from numpy import * 
from pylab import * 
from scipy import special

N = 1000
r_min = 0.001
r_max = 25.0
r_0 = 1.0

r = zeros((N))

ivals = array(range(0, N))
spaceri = math.log(r_max/r_min) * (1.0/(N - 1.0))
r = r_min * exp(spaceri * ivals)

M = zeros((N,N))
for i in range (0,N):
	r_i = r[i]
	print "Coulomb kernel:", i, "/", N
	for j in range (0,N):
		r_j = r[j]
		if i == j:
			if (j < N - 1):
  			    Dr = 0.5 * (r[j + 1] - r_j)
		        else:
			    Dr = 0.5 * (r_j - r[j - 1])
			# Definitely needs simplifications!
			# In particular, the sum of Lplus and Lminus
			# can be simplified by taking into account dr << r_j
			Ldr    = 2.0 * Dr * math.log(1.0 / Dr)
			Lplus  = (2.0 * r_j + Dr) * math.log(2.0 * r_j + Dr)
			Lminus = (2.0 * r_j - Dr) * math.log(2.0 * r_j - Dr)
			Lconst =  2.0 * Dr * math.log(4.0)
			M[i,j] = 0.5 * (Ldr + Lplus - Lminus + Lconst)
			#M[i,j] = 0.5 * ((2.0 * Dr * math.log(1.0 / Dr)) + ((2 * r_j + Dr) * math.log(2 * r_j + Dr)) - ((2 * r_j - Dr) * math.log(2 * r_j - Dr)) + (2 * Dr * math.log(4)))
		else:
			if j == 0:
			        b = 0.5 * (r_j + r[j + 1])
			        d = b - r_j
				mu_top =  (4 * r_i * b) / (r_i + b)**2
				mu_bot =  (4 * r_i * r_j) / (r_i + r_j)**2
				ellip_top = special.ellipk(mu_top)
				ellip_bot = special.ellipk(mu_bot)
				alpha_top = b / (r_i + b)
				alpha_bot = (r_j)/(r_i + r_j)
				I_top = ellip_top * alpha_top
				I_bot = ellip_bot * alpha_bot
				M[i,j] = 0.5 * (I_top + I_bot) * d
			elif j == (N - 1):
				a = 0.5 * (r_j + r[j - 1])
				d = r_j - a
				mu_top =  (4 * r_i * r_j) /  (r_i + r_j)**2
				mu_bot =  (4 * r_i * a)   /  (r_i + a)**2 
				ellip_top = special.ellipk(mu_top)
				ellip_bot = special.ellipk(mu_bot)
				alpha_top = r_j / (r_i + r_j)
				alpha_bot = a/(r_i + a)
				I_top = ellip_top * alpha_top
				I_bot = ellip_bot * alpha_bot
				M[i,j] = 0.5 * (I_top + I_bot) * d
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
				M[i,j] = 0.5 * (I_top + I_bot) * d


M0 = array(M)			      
import numpy
numpy.save("Mx.dat", M0)
rho_j = 1.0 / (1.0 + r**2)**1.5				
phi_i = dot(M, rho_j) 



U = 1.0 / ( 1.0 + r**2)**0.5 
#U = np.zeros((1000,1))

#for i in range (0,1000):
#	U[i] = 1 / (math.sqrt(1 + (r_i[i])**2))

plot(r, phi_i)
plot(r, U)
savefig("image1.pdf")

figure()
loglog(r, phi_i)
loglog(r, U)
savefig("image2.pdf")
show()
