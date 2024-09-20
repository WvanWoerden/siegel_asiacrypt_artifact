from sage.quadratic_forms.count_local_2 import CountAllLocalTypesNaive
from sage.arith.misc import valuation
from sage.rings.rational_field import QQ
from sage.quadratic_forms.genera.normal_form import p_adic_normal_form
from sage.misc.verbose import verbose

class QuadraticFormFixed(QuadraticForm):

	# This function introduces the new functionality to count local solutions more efficiently
	# Essentially this function will only be called for p=2 and k=3.
	# It runs in polynomial time in self.dim() and p^k in contrast to the function CountAllLocalTypesNaive which requires >= (p^k)^self.dim() operations.
	def CountAllLocalGoodTypesNormalForm(self, p, k, m, zvec, nzvec):
		r"""
		This is an internal routine, which is called by
		:meth:`sage.quadratic_forms.quadratic_form.QuadraticForm.count_congruence_solutions__good_type
		QuadraticForm.count_congruence_solutions__good_type`. See the documentation of
		that method for more details.

		INPUT:

		- ``Q`` -- quadratic form over `\ZZ` in local normal form at p with no zero blocks mod m_range
		- ``p`` -- prime number > 0
		- ``k`` -- integer > 0
		- ``m`` -- integer >= 0 (depending only on mod `p^k`)
		- ``zvec``, ``nzvec`` -- list of integers in ``range(Q.dim())``, or ``None``

		OUTPUT:

		an integer '\ge 0' representing the solutions of Good type.

		EXAMPLES::

			sage: from sage.quadratic_forms.count_local_2 import CountAllLocalGoodTypesNormalForm
			sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
			sage: Q_local_at2 = Q.local_normal_form(2)
			sage: Q_local_at3 = Q.local_normal_form(3)
			sage: CountAllLocalGoodTypesNormalForm(Q_local_at2, 2, 3, 3, None, None)
			64
			sage: CountAllLocalGoodTypesNormalForm(Q_local_at2, 2, 3, 3, [0], None)
			32
			sage: CountAllLocalGoodTypesNormalForm(Q_local_at3, 3, 2, 1, None, None)
			54
		"""
		n = self.dim()
		if n == 0:
			return 0
		m_range = p^k

		if zvec is None:
			zvec = []
		if nzvec is None:
			nzvec = []

		# determine local blocks
		blocks = []
		i=0
		while i < n-1:
			if self[i,i+1] != 0:
				blocks += [(i,i+1)]
				i+=2
			else:
				blocks += [(i,)]
				i+=1
		if i < n:
			blocks += [(i,)]

		solutions = [ [0,0] for _ in range(m_range) ] # [good, not good]
		solutions[0][1] = 1
		for b in blocks:
			Q_part = self.extract_variables(b)
			for i in range(Q_part.dim()):
				for j in range(Q_part.dim()):
					Q_part[i,j] = Q_part[i,j] % m_range
			
			zvec_local = range(len(b)) if (b[0] in zvec) else None 
			nzvec_local = range(len(b)) if (b[0] in nzvec) else None


			solutions_part = [ [0,0] for _ in range(m_range) ]
			for m_part in range(m_range):
				cnt = CountAllLocalTypesNaive(Q_part, p, k, m_part, zvec_local, nzvec_local)
				solutions_part[m_part][0] = cnt[1]
				solutions_part[m_part][1] = cnt[0] - cnt[1]

			# compute convolution of counts
			solutions_new = [ [0,0] for _ in range(8) ]
			for m1 in range(m_range):
				for m2 in range(m_range):
					total = (solutions[m1][0] + solutions[m1][1]) * (solutions_part[m2][0] + solutions_part[m2][1])
					good = total - solutions[m1][1] * solutions_part[m2][1]
					solutions_new[(m1+m2)%m_range][0] += good
					solutions_new[(m1+m2)%m_range][1] += total - good
			solutions = solutions_new

		return solutions[m%m_range][0]

	# count_congruence_solutions__good_type is a direct copy of count_congruence_solutions__good_type 
	# in sage.quadratic_forms.quadratic_form__count_local_2 with the following edits:
	# - Call new method CountAllLocalGoodTypesNormalForm instead of CountAllLocalTypesNaive
	def count_congruence_solutions__good_type(self, p, k, m, zvec, nzvec):
		r"""
		Count the good-type solutions of `Q(x) = m` (mod `p^k`) satisfying the
		additional congruence conditions described in
		assumes Q is in local normal form at p
		:meth:`QuadraticForm.count_congruence_solutions_as_vector`.

		INPUT:

		- ``p`` -- prime number > 0

		- ``k`` -- integer > 0

		- ``m`` -- integer (depending only on mod `p^k`)

		- ``zvec``, ``nzvec`` -- lists of integers up to dim(`Q`)

		EXAMPLES::

			sage: Q = DiagonalQuadraticForm(ZZ, [1,2,3])
			sage: Q.count_congruence_solutions__good_type(3, 1, 0, None, None)
			12
		"""
		return self.CountAllLocalGoodTypesNormalForm(p, k, m, zvec, nzvec)

	# local_density is a direct copy of local_density in sage.quadratic_forms.quadratic_form__local_density_interfaces with the following edits:
	# - Casting of local normal form from QuadraticForm to QuadraticFormFixed
	# - First use the more efficient p_adic_normal_form to get a (block)-diagonal normal form over Z_p, 
	#   then follow by QuadraticForm.local_normal_form to get the correct shape.
	def local_density(self, p, m):
		"""
		Return the local density.

		.. NOTE::

			This screens for imprimitive forms, and puts the quadratic
			form in local normal form, which is a *requirement* of the
			routines performing the computations!

		INPUT:

		- ``p`` -- a prime number > 0
		- ``m`` -- integer

		OUTPUT: a rational number

		EXAMPLES::

			sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])   # NOTE: This is already in local normal form for *all* primes p!
			sage: Q.local_density(p=2, m=1)
			1
			sage: Q.local_density(p=3, m=1)
			8/9
			sage: Q.local_density(p=5, m=1)
			24/25
			sage: Q.local_density(p=7, m=1)
			48/49
			sage: Q.local_density(p=11, m=1)
			120/121
		"""
		n = self.dim()
		if n == 0:
			raise TypeError("we do not currently handle 0-dim'l forms")

		# Find the local normal form and p-scale of Q 0    --  Note: This uses the valuation ordering of local_normal_form.
		#                                                     TO DO:  Write a separate p-scale and p-norm routines!
		if p == 2:
			G = self.matrix()
			verbose("Compute p_adic_normal_form \n")
			Q_local = QuadraticFormFixed(Matrix(ZZ, p_adic_normal_form(G, p)[0]))
			verbose("Compute local_normal_form \n")
			Q_local = QuadraticFormFixed(Q_local.local_normal_form(p).matrix())
			verbose("Done.\n")
		else:
			G = self.matrix()/2
			verbose("Compute p_adic_normal_form \n")
			Q_local = QuadraticFormFixed(2*Matrix(ZZ, p_adic_normal_form(G, p)[0]))
			verbose("Compute local_normal_form \n")
			Q_local = QuadraticFormFixed(Q_local.local_normal_form(p).matrix())
			verbose("Done.\n")
		if n == 1:
			p_valuation = valuation(Q_local[0, 0], p)
		else:
			p_valuation = min(valuation(Q_local[0, 0], p),
							  valuation(Q_local[0, 1], p))

		# If m is less p-divisible than the matrix, return zero
		if ((m != 0) and (valuation(m, p) < p_valuation)):   # Note: The (m != 0) condition protects taking the valuation of zero.
			return QQ(0)

		# If the form is imprimitive, rescale it and call the local density routine
		p_adjustment = QQ(1) / p**p_valuation
		m_prim = QQ(m) / p**p_valuation
		Q_prim = Q_local.scale_by_factor(p_adjustment)

		# Return the densities for the reduced problem
		return Q_prim.local_density_congruence(p, m_prim)

	# siegel_product is a direct copy of sage.quadratic_form.quadratic_form__siegel_product with the following edits:
	# - Removal of Q_normal = self.local_normal_form(p) (is already computed in local_density anyway)
	def siegel_product(self, u):
		"""
		Compute the infinite product of local densities of the quadratic
		form for the number `u`.

		EXAMPLES::

			sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
			sage: Q.theta_series(11)
			1 + 8*q + 24*q^2 + 32*q^3 + 24*q^4 + 48*q^5 + 96*q^6 + 64*q^7 + 24*q^8 + 104*q^9 + 144*q^10 + O(q^11)

			sage: Q.siegel_product(1)
			8
			sage: Q.siegel_product(2)      # This one is wrong -- expect 24, and the higher powers of 2 don't work... =(
			24
			sage: Q.siegel_product(3)
			32
			sage: Q.siegel_product(5)
			48
			sage: Q.siegel_product(6)
			96
			sage: Q.siegel_product(7)
			64
			sage: Q.siegel_product(9)
			104

			sage: Q.local_density(2,1)
			1
			sage: M = 4; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 1]) / M^3
			1
			sage: M = 16; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 1]) / M^3  # long time (2s on sage.math, 2014)
			1

			sage: Q.local_density(2,2)
			3/2
			sage: M = 4; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 2]) / M^3
			3/2
			sage: M = 16; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 2]) / M^3  # long time (2s on sage.math, 2014)
			3/2

		TESTS::

			sage: [1] + [Q.siegel_product(ZZ(a))  for a in range(1,11)] == Q.theta_series(11).list()  # long time (2s on sage.math, 2014)
			True
		"""
		# Protect u (since it fails often if it's an just an int!)
		u = ZZ(u)

		n = self.dim()
		d = self.det()       # ??? Warning: This is a factor of 2^n larger than it should be!

		# DIAGNOSTIC
		verbose("n = " + str(n))
		verbose("d = " + str(d))
		verbose("In siegel_product:  d = " + str(d) + "\n")

		# Product of "bad" places to omit
		S = 2 * d * u

		# DIAGNOSTIC
		verbose("siegel_product Break 1. \n")
		verbose(" u = " + str(u) + "\n")

		# Make the odd generic factors
		if n % 2:
			m = (n - 1) // 2
			d1 = fundamental_discriminant(((-1)**m) * 2*d * u)     # Replaced d by 2d here to compensate for the determinant
			f = abs(d1)

			# Make the ratio of factorials factor: [(2m)! / m!] * prod_{i=1}^m (2*i-1)
			factor1 = 1
			for i in range(1, m+1):
				factor1 *= 2*i - 1
			for i in range(m+1, 2*m + 1):
				factor1 *= i

			genericfactor = factor1 * ((u / f) ** m) \
				* QQ(sqrt((2 ** n) * f) / (u * d)) \
				* abs(QuadraticBernoulliNumber(m, d1) / bernoulli(2*m))

		# DIAGNOSTIC
		verbose("siegel_product Break 2. \n")

		# Make the even generic factor
		if ((n % 2) == 0):
			m = n // 2
			verbose("Compute fundamental_discriminant. \n")
			d1 = fundamental_discriminant(((-1)**m) * d)
			f = abs(d1)

			verbose("Compute generic factor." + str(m) + " " + str(d1) + "\n")
			genericfactor = m / QQ(sqrt(f*d)) \
				* ((u/2) ** (m-1)) * (f ** m) \
				/ abs(QuadraticBernoulliNumber(m, d1)) \
				* (2 ** m)  # This last factor compensates for using the matrix of 2*Q

		# Omit the generic factors in S and compute them separately
		omit = 1
		include = 1

		verbose("Start factoring.\n")
		S_divisors = prime_divisors(S)

		G = self.matrix()/2
		for p in S_divisors:
			# DIAGNOSTIC
			verbose(" p = " + str(p) + " and its Kronecker symbol (d1/p) = (" + str(d1) + "/" + str(p) + ") is " + str(kronecker_symbol(d1, p)) + "\n")

			omit *= 1 / (1 - (kronecker_symbol(d1, p) / (p**m)))

			# DIAGNOSTIC
			verbose(" omit = " + str(omit) + "\n")
			verbose(" p = " + str(p) + "\n")
			verbose(" u = " + str(u) + "\n")
			verbose(" include = " + str(include) + "\n")

			include *= self.local_density(p, u)

			# DIAGNOSTIC
			verbose("    ---  Exiting loop \n")

		# Return the final factor (and divide by 2 if n=2)
		if n == 2:
			return genericfactor * omit * include / 2
		return genericfactor * omit * include
