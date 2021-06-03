"""
This scripts solves equation B1 and generates the data in Figure 11 in the appendix. For each number of positive and negative
regulators of each node, it calculates how many negative edges we keep in the Boolean function of the same node in the simplified model.
"""


import math
import numpy as np




for p in range (1,26): # number of positive regulators
    for n in range (26): # number of negative regulators
        if p+n>25:
            break
        else:
            summation=0
            for j in range (1,p+1):
                j_of_p=math.factorial(p)//math.factorial((p-j))//math.factorial(j)
                sum=0
                if n >= 4 * j - 1:
                    limit = 4*j
                else:
                    limit = n + 1
                for i in range (limit):
                    sum=sum+(math.factorial(n)//math.factorial(n-i)//math.factorial(i))
                sum=sum*j_of_p
                summation=summation+sum

            sol=n-(np.log(summation/((2**p)-1))/np.log(2))
            x=n
            while True:
                if x<=sol:
                    break
                else:
                    x=x-1

            print('\nNumber of positive regulators (N_p): '+str(p))
            print('Total number of negative regulators (N_n): '+str(n))
            #print('the summation is: ' + str(summation))
            print('the minimum number of inactive negative regulators to be added in the Boolean function (x): '+str(x))



