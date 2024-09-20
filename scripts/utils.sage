from math import lgamma

def get_lambda1_sq(series):
	sm = 0.
	l = series.prec()
	sm = 0.
	for j in range(1, l):
		sm += series[j]
		if sm >= 2:
			return j

def get_smoothing(series, eps):
	l = series.prec()
	s_min = 0
	s_max = l

	while s_max - s_min > 1/100000:
		s = (s_min+s_max)/2
		sm = sum([series[j] * exp(-pi * s^2 * j) for j in range(1, l)])
		if sm < eps:
			s_max = s
		else:
			s_min = s
	return float(s)		

def ball_log_vol(n):
    return float((n/2.) * log(pi) - lgamma(n/2. + 1))