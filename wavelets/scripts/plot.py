
import matplotlib.pyplot as plt
import numpy as np
import operator as op

def binom(n, k):
    k = min(k, n-k)
    if k == 0: return 1
    numer = reduce(op.mul, xrange(n, n-k, -1))
    denom = reduce(op.mul, xrange(1, k+1))
    return numer//denom

def lagrange(p, x):
    L = []
    for i in range(1,2*p+1):
        P = 1.0
        for j in range(1,2*p+1):
            if i != j:
                P*= (x+p-j)/(i-j)
        L.append(P)
    
    return L

def deslaurierDubuc(m, x):
    H = []
    L = lagrange(m,x)

    #print "Lagrange: " + str(L)
    
    for k in range(-2*m+1,2*m):
        if k%2 == 0:
            H.append(1.0 if k==0 else 0.0)
        else:
            H.append(L[m + (k-1)//2])

    return H

def generateScalingFunction(m,levels):

    # !! normally support is 4*m-1 wide !! #
    length= 2**(levels)*(4*m-2)+1
    
    phi = np.zeros(shape=(2, length), dtype=float)
    phi[0][(length-1)//2] = 1.0

    H = deslaurierDubuc(m,0.5)
    print "Smoothing filter = " + str(H)

    for l in range(1, levels+1):
        for i in range(0, (length-1)//(2**(levels-l))):
            xk = i*2**(levels-l)

            phi[l%2][xk] = 0.0
            for j in range(-2*m+1,2*m):
                jj = j*2**(levels-l)
                if(xk+jj >= 0 and xk+jj < length):
                    phi[l%2][xk] += phi[(l+1)%2][xk+jj]*H[2*m-1+j] 

    return phi[levels%2]

def generateWavelet(j,k,pmax,maxLevels,interval,boundaryMode='truncate'):

  if not hasattr(generateWavelet, "scalingFuncs"):
     generateWavelet.scalingFuncs = []
    
  if(boundaryMode == 'truncate'):
      p = pmax
      colorId = j
  elif(boundaryMode == 'adapt'):
      p = pmax
      while k > 2**j - 2*p + 1 or k < 2*p - 1:
          p -= 1
          if(p==1): break
      colorId = p
  else:
      p = max(0,min(j,pmax))

  #generate scaling functions at maximum accuracy level (for j = 0)
  for i in range(len(generateWavelet.scalingFuncs)+1, p+2):
     print "generating scaling func " + str(i) + " !"
     generateWavelet.scalingFuncs.append(generateScalingFunction(i,maxLevels))

  #subsample the good scaling function to build the targetted wavelet
  funcLength= generateWavelet.scalingFuncs[p-1].size
  intLength= 2**(maxLevels)*(interval[1] - interval[0])+1

  offset = 2**(maxLevels-j)*k 

  phi = np.zeros(shape=(intLength),dtype=float)

  for i in range(0,max(intLength,funcLength)):

      ip = i
      if ip < 0 or ip >= intLength: 
          continue
      
      iw = 2**(j+1)*(i-offset) + 2**(maxLevels)*(2*p-1)
      if iw < 0 or iw >= funcLength:
          continue

      phi[ip] = generateWavelet.scalingFuncs[p-1][iw]

  return [phi,colorId]
      


plt.title("Wavelets level = 5")
jmax = 5; 
pmax = 5;

maxLevels = 10;
interval = [0,1]
colors = ['b','g','r','c','m','y','k']

nPoints = 2**(maxLevels)*(interval[1]-interval[0]) + 1
x=np.linspace(interval[0],interval[1],nPoints)

for j in range(0,jmax+1):
    for k in range(0,2**j+1):
        if j==0 or (j>0 and k%2==1):
            W = generateWavelet(j,k,pmax,maxLevels,interval,boundaryMode='adapt')
            if(j==jmax):
                plt.plot(x, W[0],color=colors[W[1]%len(colors)])


plt.show()

#jmax = 1; 
#levels = 0;

#pmax = jmax+1;
#interval = [-2*pmax+1,2*pmax-1]

#nPoints = 2**(levels)*(interval[1]-interval[0]) + 1

#x=np.linspace(interval[0],interval[1],nPoints)
#plt.plot(x, generateScalingFunction(pmax,levels),marker='x')

#plt.show()

