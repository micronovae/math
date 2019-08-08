import numpy as np
import math

def SO3_mat(v,x):
    
    v = np.array(v)
    x = np.array(x)
    v = v/np.sqrt(sum(v*v))
    x = x/np.sqrt(sum(x*x))
#     print('ok')
    print(x*v)
#     theta = math.acos(x*v)
#     k = np.multiply(v,x)
#     k = k/np.sqrt(sum(k*k))
    
    return 0
    
    