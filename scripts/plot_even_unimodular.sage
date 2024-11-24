load("utils.sage")	

import matplotlib.pyplot as plt
import numpy as np

plt.style.use('../plots/custom.mplstyle')

def bound_lemma(n):
	return 2*round(exp(-2*ball_log_vol(n)/n) * (3*zeta(n/4).n()/2.)^(2./n)/2.)

def bound_milnor(n):
	return 2*round(exp(-2*ball_log_vol(n)/n) * (5/3.)^(2./n)/2.)

def bound_smoothing(n, eps):
	return (17.1/eps)^(1./n)

def bound_gh(n):
	return exp(-2*ball_log_vol(n)/n)

def get_results():
	results = []
	for n in range(8,513,8):
		prec=n//2
		# average theta series for even unimodular lattice
		f=eisenstein_series_qexp(n//2, prec, normalization='constant').V(2)

		lambda1_sq = get_lambda1_sq(f)
		lambda1_sq = 2*ceil(lambda1_sq/2)

		lambda_gh_sq = 2*ceil(n / (4*pi*e))
		assert(lambda1_sq >= lambda_gh_sq)

		eps1 = 0.01
		smooth_eps1 = get_smoothing(f, eps1)

		eps2 = 1/sqrt(2^64 * 128)
		smooth_eps2 = get_smoothing(f, eps2)

		eps3 = 17.1 * exp(-n)
		smooth_eps3 = get_smoothing(f, eps3)

		result = (n, lambda1_sq, lambda_gh_sq, smooth_eps1, smooth_eps2, smooth_eps3)
		results += [result]
		print(result)
	return results

def print_results_packing(results):
	n_min = results[0][0]
	n_max = results[-1][0]

	fig = plt.figure(figsize=(10,4))

	bounds_lemma = [(n, bound_lemma(n)) for n in range(n_min, n_max+1)]
	bounds_lemma = list(zip(*bounds_lemma))
	plt.plot(bounds_lemma[0], np.sqrt(bounds_lemma[1]), linewidth=2, label="Lemma 5 (lower bound)")

	bounds_milnor = [(n, bound_milnor(n)) for n in range(n_min, n_max+1)]
	bounds_milnor = list(zip(*bounds_milnor))
	plt.plot(bounds_milnor[0], np.sqrt(bounds_milnor[1]), linewidth=2, color="red", linestyle='dashed', label="[MH+73] (lower bound)")

	bounds_gh = [(n, bound_gh(n)) for n in range(n_min, n_max+1)]
	bounds_gh = list(zip(*bounds_gh))
	plt.plot(bounds_gh[0], np.sqrt(bounds_gh[1]), linewidth=2, linestyle="dashdot", color="green", label=r'Gaussian Heuristic ($\omega_n^{-1/n} \approx \sqrt{n/2 \pi e}$)')

	res = list(zip(*results))
	plt.scatter(res[0], np.sqrt(res[1]), color="orange", marker="x", label="Concrete")

	plt.xlabel("Dimension (n)")
	plt.ylabel("First minimum $\\lambda_1(\\mathcal L )$")

	plt.legend()
	fig.tight_layout()
	plt.savefig("../plots/even_packing.pdf")	
	plt.close()

def print_results_smoothing(results):
	n_min = results[0][0]
	n_max = results[-1][0]

	fig = plt.figure(figsize=(10,4))

	eps1 = 0.01
	bounds_eps1 = [(n, bound_smoothing(n, eps1)) for n in range(n_min, n_max+1)]
	bounds_eps1 = list(zip(*bounds_eps1))
	plt.plot(bounds_eps1[0], bounds_eps1[1], color="red", label="Lemma 8 (upper bound)")

	eps2 = 1/sqrt(2^64 * 128)
	bounds_eps2 = [(n, bound_smoothing(n, eps2)) for n in range(n_min, n_max+1)]
	bounds_eps2 = list(zip(*bounds_eps2))
	plt.plot(bounds_eps2[0], bounds_eps2[1], color="blue")

	bounds_eps3 = [(n, bound_smoothing(n, 17.1*exp(-n))) for n in range(n_min, n_max+1)]
	bounds_eps3 = list(zip(*bounds_eps3))
	plt.plot(bounds_eps3[0], bounds_eps3[1], color="green")


	res = list(zip(*results))
	plt.scatter(res[0], res[3], color="red", marker=".", label=r'Concrete ($\varepsilon = 0.01$)')
	plt.scatter(res[0], res[4], color="blue", marker=".", label=r'Concrete ($\varepsilon = 2^{-71/2}$)')
	plt.scatter(res[0], res[5], color="green", marker=".", label=r'Concrete ($\varepsilon = 17.1 e^{-n}$)')

	plt.xlim((30, n_max))
	plt.ylim((0.9, np.e+0.1))

	plt.xlabel("Dimension (n)")
	plt.ylabel("Smoothing $\\eta_\\varepsilon(G)$")

	plt.legend()
	fig.tight_layout()
	plt.savefig("../plots/even_smoothing.pdf")	
	plt.close()

if __name__ == "__main__":
	results = get_results()
	print_results_packing(results)
	print_results_smoothing(results)