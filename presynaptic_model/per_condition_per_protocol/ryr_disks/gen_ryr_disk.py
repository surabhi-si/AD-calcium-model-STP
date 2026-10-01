import numpy as np
from matplotlib.pyplot import *

n = 150 # number of RyRC1_L molecules
#for n in [2**i for i in range(11)]:
for n in [n]:
    print n
    pre = '''ryr_disk RELEASE_SITE {
      SHAPE = LIST
      MOLECULE_POSITIONS {\n'''

    post = '''  }
      SITE_RADIUS = 0.03
    }'''

    text = pre
    for i in range(n):
        s = np.random.uniform(-1,1,2)
        plot(s[0], s[1], 'o')
        text += '  RyRC1_L\' ['+str(s[0])+','+str(s[1])+',0]\n'
    text += post

    #print text
    f = open('ryr_disk_'+str(n)+'.mdl','w')
    f.write(text)
show()
