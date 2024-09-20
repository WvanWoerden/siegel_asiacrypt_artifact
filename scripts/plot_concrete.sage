load("utils.sage")	

import matplotlib.pyplot as plt
import numpy as np

plt.style.use('../plots/custom.mplstyle')

def bound_theorem(n, logdet):
	return round(exp(-2*ball_log_vol(n)/n + 2*logdet/n) * (7*zeta(3).n()/(9*zeta(2).n()))^(2./n))

def bound_gh(n, logdet):
	return exp(-2*ball_log_vol(n)/n + 2*logdet/n)


# values consists of pairs (index, N_genus[index])
def get_lambda1_sq_concrete(values):
	sm = 0.
	for j in range(0, values.shape[0]+1):
		sm += values[j,1]
		if sm >= 2:
			return values[j,0]
	return -1

def get_results():
	results = []

	q=521
	exeriment_cases_n = [16,32,48,64,80,96, 112, 128]

	for n in exeriment_cases_n:
		k = n//2
		logdet = RR((n-k) * log(q))
		values = np.loadtxt("../data/siegel_product_{}_{}_{}".format(n,k,q))
		lambda1_sq = get_lambda1_sq_concrete(values)
		
		lambda_gh_sq = bound_gh(n, logdet)
		assert(lambda1_sq >= lambda_gh_sq)

		result = (n, logdet, lambda1_sq/lambda_gh_sq)
		results += [result]
		print(result)
	return results

def print_results_packing(results):
	n_min = results[0][0]
	n_max = results[-1][0]

	q = 521

	fig = plt.figure(figsize=(10,4))

	bounds_theorem = [(n, bound_theorem(n, RR(n/2 * log(q)))/bound_gh(n, RR(n/2 * log(q)))) for n in range(n_min, n_max+1)]
	bounds_theorem = list(zip(*bounds_theorem))
	plt.plot(bounds_theorem[0], np.sqrt(bounds_theorem[1]), linewidth=2, color="red", label="Theorem 1 (lower bound)", zorder=2)

	bounds_gh = [(n, 1.) for n in range(n_min, n_max+1)]
	bounds_gh = list(zip(*bounds_gh))
	plt.plot(bounds_gh[0], np.sqrt(bounds_gh[1]), linewidth=2, linestyle="dashdot", color="green", label=r'Gaussian Heuristic', zorder=3)

	res = list(zip(*results))
	plt.scatter(res[0], np.sqrt(res[2]), color="orange", marker="x", label="Concrete", zorder=4)

	plt.xticks(list(range(16, 129, 16)))
	plt.xlabel("Dimension (n)")
	plt.ylabel("Ratio $\\lambda_1(\\mathcal L ) / gh(\\mathcal L )$")

	plt.legend()
	fig.tight_layout()
	plt.savefig("../plots/concrete_packing.pdf")	
	plt.close()

if __name__ == "__main__":
	results = get_results()
	print_results_packing(results)