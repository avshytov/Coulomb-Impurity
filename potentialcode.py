N = 1000
r_min = 0.001
r_max = 25.0
z = r_max / r_min
r_0 = 1.0
r_i = np.zeros((N, 1))
r_j = np.zeros((N, 1))
rho_j = np.zeros((N, 1))
M = np.zeros((N,N))

for i in range (0,N):
	spaceri = math.log(z) * (float(i)/(float(N - 1)))
	r_i[i] = r_min * math.exp(spaceri)
	for j in range (0,N):
		spacerj = math.log(z) * (float(j)/(float(N - 1)))
		r_j[j] = r_min * math.exp(spacerj)
		spacerjnext = math.log(z) * (float(j + 1)/(float(N - 1)))
		r_jnext = r_min * math.exp(spacerjnext)	
		spacerjprev = math.log(z) * (float(j-1)/(float(N-1)))
		r_jprev = r_min * math.exp(spacerjprev)
		rho_j[j] = (r_0)/(((r_0)**2 + (r_j[j])**2)**(1.5))
		
		if i == j:
			Dr = (0.5 * (r_jnext - r_j[j]))
			M[i,j] = 0.5 * ((Dr * math.log(1/((Dr)**2))) + ((2 * r_j[j] + Dr) * math.log(2 * r_j[j] + Dr)) - ((2 * r_j[j] - Dr) * math.log(2 * r_j[j] - Dr)) + (2 * Dr * math.log(4)))
		else:
			if j == 0:
				b = 0.5 * (r_j[j] + r_jnext)
				d = b - r_j[j]
				mu_top =  (4 * r_i[i] * b) / ((r_i[i] + b)**2)
				mu_bot =  (4 * r_i[i] * r_j[j]) / ((r_i[i] + r_j[j])**2)
				ellip_top = scipy.special.ellipk(mu_top)
				ellip_bot = scipy.special.ellipk(mu_bot)
				alpha_top = b / (r_i[i] + b)
				alpha_bot = (r_j[j])/(r_i[i] + r_j[j])
				I_top = ellip_top * alpha_top
				I_bot = ellip_bot * alpha_bot
				M[i,j] = 0.5 * (I_top + I_bot) * d
			elif j == (N - 1):
				a = 0.5 * (r_j[j] + r_jprev)
				d = r_j[j] - a
				mu_top =  (4 * r_i[i] * r_j[j]) / ((r_i[i] + r_j[j])**2)
				mu_bot =  (4 * r_i[i] * a) / ((r_i[i] + a)**2)
				ellip_top = scipy.special.ellipk(mu_top)
				ellip_bot = scipy.special.ellipk(mu_bot)
				alpha_top = r_j[j] / (r_i[i] + r_j[j])
				alpha_bot = (a)/(r_i[i] + a)
				I_top = ellip_top * alpha_top
				I_bot = ellip_bot * alpha_bot
				M[i,j] = 0.5 * (I_top + I_bot) * d
			else:
				a = 0.5 * (r_j[j] + r_jprev)
				b = 0.5 * (r_j[j] + r_jnext)
				d = b - a
				mu_top =  (4 * r_i[i] * b) / ((r_i[i] + b)**2)
				mu_bot =  (4 * r_i[i] * a) / ((r_i[i] + a)**2)
				ellip_top = scipy.special.ellipk(mu_top)
				ellip_bot = scipy.special.ellipk(mu_bot)
				alpha_top = b / (r_i[i] + b)
				alpha_bot = (a)/(r_i[i] + a)
				I_top = ellip_top * alpha_top
				I_bot = ellip_bot * alpha_bot
				M[i,j] = 0.5 * (I_top + I_bot) * d


phi_i = np.dot(M, rho_j)

U = np.zeros((1000,1))
for i in range (0,1000):
	U[i] = 1 / (math.sqrt(1 + (r_i[i])**2))

plot(r_i, phi_i)
plot(r_i, U)
savefig("image1.pdf")


loglog(r_i, phi_i)
loglog(r_i, U)
savefig("image2.pdf")

