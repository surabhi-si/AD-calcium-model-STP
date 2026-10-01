import numpy as np
import sys
from matplotlib.pyplot import *

# Distribution of VDCC with a standard deviation
"""
effDist = 1
while(abs(effDist) > 0.001):
    n = 80
    mu, sigma = 0, float(sys.argv[1])
    x = np.random.normal(mu, 0.020, n)
    y = np.random.normal(mu, sigma, n)
    
    effDist = y.sum()/n
    print(effDist)

fname = 'vdcc_disk_' + str(n) + '_' + str(int(sigma*1000)) + '.mdl'
print(fname)
f = open(fname, 'w')
f.write('''vdcc_disk RELEASE_SITE {
  SHAPE = LIST
    MOLECULE_POSITIONS {
''')

for i in range(n):
    f.write("    VDCC_PQ_C0' [%f, %f, 0]\n" % (x[i], y[i]))

f.write('''  }
  SITE_RADIUS = 0.042
}''')

plot(x,y,'.')
gca().set_aspect('equal', adjustable='box')
#show()
"""

# Random distribution of VDCC within a square
"""
n = 81
l = 0.05 # lenth in um
for l in np.arange(0.02, 0.13, 0.01):
    fname = 'vdcc_disk_' + str(n) + '_' + str(int(l*2000+0.1)) + 'l.mdl'
    print(fname)


    x = np.random.uniform(-l, l, n)
    y = np.random.uniform(-l, l, n)


    f = open(fname, 'w')
    f.write('''vdcc_disk RELEASE_SITE {
      SHAPE = LIST
        MOLECULE_POSITIONS {
    ''')

    for i in range(n):
        f.write("    VDCC_PQ_C0' [%f, %f, 0]\n" % (x[i], y[i]))

    f.write('''  }
      SITE_RADIUS = 0.042
    }''')
    
    
    plot(x,y,'.')
    gca().set_aspect('equal', adjustable='box')
    show()
""" 

# Random distribution of VDCC within a square

n = 81
l = 1.0 # lenth in um

for n in range(81,210,2):
    fname = 'vdcc_disk_' + str(n) + '.mdl'
    print(fname)


    x = np.random.uniform(-l, l, n)
    y = np.random.uniform(-l, l, n)

    
    f = open(fname, 'w')
    f.write('''vdcc_disk RELEASE_SITE {
    SHAPE = LIST
    MOLECULE_POSITIONS {\n''')

    for i in range(n):
        f.write("        VDCC_PQ_C0' [%f, %f, 0]\n" % (x[i], y[i]))

    f.write('''    }
    SITE_RADIUS = 0.042\n}''')
    
    """
    plot(x,y,'.')
    gca().set_aspect('equal', adjustable='box')
    show()
    """