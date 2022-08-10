#!/usr/bin/env python
from random import seed
from random import random
import numpy as np
import matplotlib.pyplot as plt
import os
import sys


def main(r,p,num):
	a=[0,1]
	p1 = np.random.negative_binomial(r,p,num)
	np.savetxt('virus_prod.csv', p1, delimiter=',')
	return p1


if __name__ == "__main__":
	r=float(sys.argv[1])
	p=float(sys.argv[2])
	num=int(sys.argv[3])
	main(r,p,num)
