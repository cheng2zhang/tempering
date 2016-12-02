#!/usr/bin/env python



import os, sys, getopt, glob
from math import *



fnin = "spec.log"
fnout = ""
fnrdf = "spec.rdf"
fndsep = "dsep.dat"
dordf = False
donc = False
dodsep = True
sepmin = 2
sepmax = 9
nccut = 0
nn = 0
i0 = 0
every = 1



def showhelp():
  print "calcspec.py [Options] spec.log"
  print "Description:"
  print "  Compute quantities from the coordinates of special atoms\n"
  print "Options:"
  print "  --nc=            set the cutoff for the number of contacts (turns on this calculation"
  print "  --rdf            compute the radial distribution function"
  print "  --fnrdf=         set the output file name for the radial distribution function"
  print "  --dsep           compute the distance vs. residue separation"
  print "  --sepmin=        set the minimal separation"
  print "  --sepmax=        set the maximal separation"
  print "  --fndsep=        set the output file name for distance vs. residue separation"
  print "  --nn=            set the exclusion limit for neighboring atoms, 1 excludes nearest neighbor, 2 excludes nearest and next-nearest neighbors"
  print "  --i0=            set the first frame"
  print "  -e, --every=     set the stride of frames"
  print "\n"
  print "Once `spec.log` is created with the coordinates of speical atoms,"
  print "one can postprocess the coordiates with the python script `calcspec.py`,"
  print "and compute the radial distribution function, the distribution of contacts, etc.\n"
  print "--- Radial distribution function ---\n"
  print "To compute the radial of distribution function"
  print "\n  calcspec.py spec.log -e 100 --rdf  \n"
  print "This command computes the radial distribution of the special atoms,"
  print "with a stride of 100 frames."
  print "The output is saved as `spec.rdf`.\n"
  print "--- Number of contacts ---\n"
  print "To compute the number of contacts with a cutoff of 7.5 angstrom, use"
  print "\n  calcspec.py spec.log --nc=7.5 --nn=3\n"
  print "The option `--nn=3` excludes the nearest neighbor (i -- i+1),"
  print "next-nearest neighbors (i -- i+2), and next-next-nearest neighbors (i -- i+3)."
  print "The column in `spec.log` for coordinates will be replaced by the number of contacts"
  print "in the output file `spec_nc.log`.\n"
  print "--- Number of separations ---\n"
  print "To compute the distance vs. the separation of residues, use"
  print "\n  calcspec.py spec.log --dsep\n"
  print "The output is `dsep.dat`."
  print "The first column is the residue separation,"
  print "the second column is the average distance,"
  print "the third column is the variance of distance,"
  print "the fourth column is the number of samples."
  print "From this file, we can compute the exponent nu."
  print "\n"

  exit(1)



def doargs():
  global fnin, fnout, fnrdf, fndsep, i0, every
  global dordf, donc, dodsep, nccut, nn, sepmin, sepmax

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "ho:e:",
        ["help", "output=", "rdf", "nc=",
          "dsep", "sepmin=", "sepmax=",
          "i0=", "every=", "nn="])
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
    elif o in ("--fnrdf",):
      dordf = True
      fnrdf = a
    elif o in ("--nc",):
      donc = True
      nccut = float(a)
    elif o in ("--dsep",):
      dodsep = True
    elif o in ("--fndsep",):
      dodsep = True
      fndsep = a
    elif o in ("--sepmin",):
      dosep = True
      sepmin = int(a)
    elif o in ("--sepmax",):
      dosep = True
      sepmax = int(a)
    elif o in ("--i0",):
      i0 = int(a)
    elif o in ("-e", "--every",):
      every = int(a)

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



def getdismat(arr):
  ''' compute the pairwise distance matrix '''
  n = len(arr)
  # initialize an nxn matrix of zeros
  dis = [[0]*n for i in range(n)]
  for i in range(n):
    for j in range(i):
      dij = vdist(arr[i], arr[j])
      dis[i][j] = dij
      dis[j][i] = dij
  return dis



def ncontact(dismat, cutoff):
  ''' compute the number of contact '''
  global nn
  n = len(dismat)
  nc = 0
  for i in range(1, n):
    for j in range(i - nn):
      if dismat[i][j] < cutoff:
        nc += 1
  return nc



class RDF:
  ''' radial distribution function '''
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

  def add(self, dis):
    if self.np <= 0:
      self.np = len(dis)
    for i in range(self.np):
      for j in range(i - self.nn):
        disij = dis[i][j]
        if disij > self.xmax:
          self.extend(disij + 10)
        k = int( disij / self.dx )
        self.arr[k] += 1
    self.cnt += 1

  def save(self, fn):
    if self.cnt <= 0: return
    for imin in range(self.n):
      if self.arr[imin] > 0: break
    for imax in range(self.n, imin + 1, -1):
      if self.arr[imax - 1] > 0: break
    src = "# %s %s\n" % (self.np, self.cnt)
    # normalization
    norm = 1.0/(self.cnt * self.dx)
    for i in range(imin, imax):
      x = self.dx * (i + 0.5)
      y = self.arr[i]
      src += "%s\t%s\t%s\n" % (x, y*norm, y)
    open(fn, "w").write(src)
    print "saving RDF to the file %s" % fn



class Ave:
  ''' accumulator for average and varaince '''
  def __init__(self):
    self.cnt = 0
    self.sx = 0
    self.sxx = 0

  def add(self, x, w = 1):
    self.cnt += w
    self.sx += w * x
    self.sxx += w * x * x

  def count(self):
    return self.cnt

  def ave(self):
    return self.sx / self.cnt if self.cnt > 0 else 0

  def var(self):
    if self.cnt <= 0: return 0
    av = self.ave()
    return self.sxx / self.cnt - av * av



class DisSep:
  ''' compute the average distance of atoms several residues apart
      currently assuming all atoms are CA atoms '''
  def __init__(self, dismat):
    self.n = len(dismat)
    self.dis = [Ave() for i in range(self.n)]
    self.frames = 0
    self.add(dismat)

  def add(self, dismat):
    if self.n != len(dismat):
      print "dimension mismatch %s vs %s" % (self.n, len(dismat))
      raise
    for i in range(1, self.n): # loop over the separations
      for j in range(self.n - i):
        self.dis[i].add( dismat[j][j+i] )
    self.frames += 1

  def save(self, fn):
    s = ""
    for i in range(1, self.n):
      cnt = self.dis[i].cnt
      ave = self.dis[i].ave()
      var = self.dis[i].var()
      s += "%d\t%s\t%s\t%s\n" % (i, ave, var, cnt)
    open(fn, "w").write(s)
    print "saving distance vs. separation results to", fn

  def computenu(self, imin, imax):
    ''' compute the exponent nu '''
    if imax > self.n - 1: imax = self.n - 1
    if imin >= imax: return
    s1 = 0
    sx = 0
    sy = 0
    sxx = 0
    sxy = 0
    for i in range(imin, imax + 1):
      x = log( 1.0 * i )
      y = log( self.dis[i].ave() )
      s1 += 1
      sx += x
      sy += y
      sxx += x * x
      sxy += x * y
    sx /= s1
    sy /= s1
    sxx = sxx / s1 - sx * sx
    sxy = sxy / s1 - sx * sy
    nu = sxy / sxx
    print "nu %s as in R ~ N^nu" % nu



def calcspec(fn, fnout):
  global i0, every
  global dordf, donc, dodsep, nccut, nn

  i = 0
  sout = ""
  rdf = RDF(0.1, nn)
  dsep = None
  xcnt = 0
  for ln in open(fn).readlines():
    s = ln.strip()
    if s == "" or s.startswith("#"): continue
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

    # do some calculation here
    dismat = getdismat(arr) # compute the distance matrix
    xout = []
    # compute the number of contacts
    if nccut > 0:
      xout.append( str( ncontact(dismat, nccut) ) )
    # compute RDF
    if dordf:
      rdf.add(dismat)
    # compute the distance vs separation
    if dodsep:
      if not dsep:
        dsep = DisSep(dismat)
      else:
        dsep.add(dismat)

    xcnt = len(xout)
    xout = tok[:j] + xout + tok[j+1:]
    ln = "\t".join(xout) + "\n"
    sout += ln

  if xcnt > 0:
    if not fnout: fnout = os.path.splitext(fn)[0] + "_nc.log"
    open(fnout, "w").write(sout)

  if dordf: rdf.save(fnrdf)
  if dodsep:
    dsep.computenu(sepmin, sepmax)
    dsep.save(fndsep)



if __name__ == "__main__":
  doargs()
  calcspec(fnin, fnout)

