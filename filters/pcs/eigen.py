import sys, numpy
from numpy import linalg as LA

# t = [74.49476,  -6.79461, 160.36840, -79.94501, -63.83840]
t = [24.51397, 3.33062, -10.21228, -8.83543, -10.12881]
w, v = LA.eig(numpy.array([[t[0], t[1], t[2]], [t[1], t[3], t[4]], [t[2], t[4], -t[0] - t[3]]]))

print w
#print v

x = []
for i in range(3):
    x.append([abs(w[i]), w[i]])
x.sort()
for i in range(3):
    w[i] = x[i][1]
print 'Xax = ', w[2] - 0.5 * (w[0] + w[1])
print 'Xrh = ', w[0] - w[1]
