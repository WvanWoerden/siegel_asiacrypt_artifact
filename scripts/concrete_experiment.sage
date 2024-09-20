import multiprocessing
import numpy as np
import sys

load("quadratic_form_fixed.sage")

n=16
k=8
q=521
values_start = 1
values_end = 800
multithreading_cores = 2

print(sys.argv)
if len(sys.argv) > 1:
	n = int(sys.argv[1])
if len(sys.argv) > 2:
	k = int(sys.argv[2])
if len(sys.argv) > 3:
	q = int(sys.argv[3])
if len(sys.argv) > 4:
    assert(len(sys.argv)>=6)
    values_start = int(sys.argv[4])
    values_end = int(sys.argv[5])
if len(sys.argv) > 6:
	multithreading_cores = int(sys.argv[6])
	assert(multithreading_cores >= 1)

def compute_siegel_product(n, k, q, val):
	B = block_matrix([[identity_matrix(k), 0],[0,q*identity_matrix(n-k)]])
	G = B*B.transpose()
	Q = QuadraticFormFixed(2*G)
	density = Q.siegel_product(val)
	return (val, density)

def map_compute(args):
        print(args)
        sol=compute_siegel_product(args[0], args[1], args[2], args[3])
        print(args, "done")
        return sol

pool = multiprocessing.Pool(multithreading_cores)
all_trials = {}
cases = [(n,k,q,val) for val in range(values_start, values_end)]
trials = list(pool.map(map_compute, cases))
all_trials = [(n, k, q, x[0], x[1]) for x in trials]
pool.close()

filename = "../data/siegel_product_"+str(n)+"_"+str(k)+"_"+str(q)
f = open(filename, "w")
for x in all_trials:
	f.write("{} {}\n".format(x[3], RR(x[4])))
f.close()
print(all_trials)

