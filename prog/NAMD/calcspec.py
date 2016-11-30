#!/usr/bin/env python



import os, sys, getopt, glob
from math import *



fnin = "spec.log"
fnout = ""
fnrdf = "spec.rdf"
dordf = False
donc = False
nccut = 0
i0 = 0
every = 1
nn = 0



def showhelp():
  print "calcspec.py [Options] spec.log"
  print "Description:"
  print "  Compute quantities from the coordinats of special atoms"
  print "Options:"
  print "  --rdf            compute the radius distribution function"
  print "  --nc=            set the cutoff for the number of contacts"
  print "  --i0=            set the first frame"
  print "  -e, --every=     set the stride of frames"
  print "  --nn=            set the nearest neighbor"
  exit(1)



def doargs():
  global fnin, fnout, i0, every
  global dordf, donc, nccut, nn

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "ho:e:",
        ["help", "output=", "rdf",
          "nc=", "i0=", "every=", "nn="])
  except:
    print "Error parsing the command line"
    sys.exit(1)

  for o, a in opts:
    if o in ("-h", "--help",):
      showhelp()
    elif o in ("-o", "--output",):
      fnout = a
    elif o in ("--rdf",):
      dordf = True
    elif o in ("--nc",):
      donc = True
      nccut = float(a)
    elif o in ("--i0",):
      i0 = int(a)
    elif o in ("-e", "--every",):
      every = int(a)
    elif o in ("--nn",):
      nn = int(a)

  if len(args) == 0:
    fnins = glob.glob("spec*.log")
    if not fnins: showhelp()
    fnin = fnins[0]
  else:
    fnin = args[0]



def vdiff(a, b):
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

def vdot(a, b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vnorm(a):
  return sqrt(vdot(a, a))

def vdist(a, b):
  return vnorm(vdiff(a, b))



def ncontact(arr, cutoff):
  ''' compute the number of contact '''
  global nn
  n = len(arr)
  nc = 0
  for i in range(1, n):
    for j in range(i - nn):
      dis = vdist(arr[i], arr[j])
      if dis < cutoff:
        nc += 1
  return nc



class RDF:
  def __init__(self, dx, nn):
    self.dx = dx
    xmax1 = 50
    self.n = int( xmax1 / dx + 0.9999 )
    self.xmax = self.n * self.dx
    self.arr = [0] * self.n
    self.np = 0
    self.cnt = 0
    self.nn = nn

  def extend(self, xmax1):
    if xmax1 > self.xmax:
      m = int( (xmax1 - self.xmax) / self.dx ) + 1
      self.n += m
      self.arr += [0] * m
      self.xmax = self.n * self.dx

  def add(self, arr):
    if self.np <= 0:
      self.np = len(arr)
    for i in range(self.np):
      for j in range(i - self.nn):
        dis = vdist(arr[i], arr[j])
        if dis > self.xmax:
          self.extend(dis + 10)
        k = int( dis / self.dx )
        self.arr[k] += 1
    self.cnt += 1

  def save(self, fn):
    if self.cnt <= 0: return
    for imin in range(self.n):
      if self.arr[imin] > 0: break
    for imax in range(self.n, imin + 1, -1):
      if self.arr[imax - 1] > 0: break
    src = "# %s %s\n" % (self.np, self.cnt)
    norm = 1.0/self.cnt
    for i in range(imin, imax):
      x = self.dx * (i + 0.5)
      y = self.arr[i]
      src += "%s\t%s\t%s\n" % (x, y*norm, y)
    open(fn, "w").write(src)
    print "saving to the file %s" % fn



def calcspec(fn, fnout):
  global i0, every
  global dordf, donc, nccut, nn

  i = 0
  sout = ""
  rdf = RDF(0.1, nn)
  for ln in open(fn).readlines():
    s = ln.strip()
    i += 1
    if i <= i0 or i % every != 0: continue
    tok = s.split("\t")
    step = tok[0]
    arr = None
    for j in range(len(tok)):
      if tok[j][-1] == ";":
        arr = tok[j][:-1].split(";")
        break
    if not arr: continue
    # parse the array to an array of 3d vectors
    arr = [[float(x) for x in s.split(",")] for s in arr]
    xout = []
    # do some calculation here
    if nccut > 0:
      xout.append( str( ncontact(arr, nccut) ) )
    if dordf:
      rdf.add(arr)
    xout = tok[:j] + xout + tok[j+1:]
    ln = "\t".join(xout) + "\n"
    sout += ln

  if not fnout:
    fnout = os.path.splitext(fn)[0] + "_nc.log"
  open(fnout, "w").write(sout)

  rdf.save(fnrdf)
  print "saving results to %s" % fnout



if __name__ == "__main__":
  doargs()
  calcspec(fnin, fnout)
