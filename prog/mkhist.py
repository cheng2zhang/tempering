#!/usr/bin/env python



import sys, os, math, glob, getopt



fninp = "ene.log"
fnout = None
dE = 0.5
dT = 5
T = 300



def showhelp():
  print "mkhist.py [Options] ene.log"
  print "Options:"
  print "  -T, --tp=:     set the temperature"
  print "  --dE=:         set the energy bin size"
  print "  --dT=:         set the temperature tolerance"
  print "  -o, --output=: set the output file"



def doargs():
  global fninp, fnout, dE, dT, T

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hT:o:i:", ["dE=", "de=", "dT=", "dt=", "tp=", "input=", "output="])
  except:
    print "Error parsing the command line"
    sys.exit(1)

  for o, a in opts:
    if o == "-h":
      showhelp()
    elif o in ("-T", "--tp"):
      T = float(a)
    elif o in ("--dT", "--dt"):
      dT = float(a)
    elif o in ("--dE", "--de"):
      dE = float(a)
    elif o in ("-i", "--input"):
      fninp = a
    elif o in ("-o", "--output"):
      fnout = a

  if len( args ) > 0:
    fninp = args[0]
  elif not os.path.exists(fninp):
    fninp = glob.glob("e*.log")[0]
  #print args, opts, fninp, fnout, T, dT, dE



class Hist:
  def __init__(self, xmin1, xmax1, dx):
    if xmin1 >= xmax1:
      xmax1 = xmin1 + 10 * dx
      xmin1 -= 10 * dx
    self.dx = dx
    self.n = int( (xmax1 - xmin1) / dx + 0.9999 )
    self.xmin = round(xmin1 / dx) * dx
    self.xmax = self.xmin + self.n * self.dx
    self.arr = [0] * self.n

  def add(self, newx):
    if newx < self.xmin:
      # extend to the left
      lftx = newx - self.dx * 10 # leave some margin
      deln = int( (self.xmin - lftx) / self.dx )
      delx = deln * self.dx
      self.xmin -= delx
      self.n += deln
      self.arr = [0] * deln + self.arr
    
    if newx >= self.xmax - self.dx:
      # extend to the right
      rtx = newx + self.dx * 10
      deln = int( (rtx - self.xmax) / self.dx )
      delx = deln * self.dx
      self.xmax += delx
      self.n += deln
      self.arr += [0] * deln

    i = int( (newx - self.xmin) / self.dx + 0.5 )
    if i < 0 or i >= self.n:
      print i, self.n, newx, self.xmin, self.xmax
      raw_input()
    self.arr[i] += 1

  def save(self, fn):
    tot = math.fsum(self.arr)
    fac = 1.0 / (self.dx * tot)
    s = ""
    for i in range(self.n):
      dist = self.arr[i] * fac
      s += "%s\t%s\t%s\n" % (self.xmin + (i + 0.5) * self.dx,
          dist, self.arr[i])
    open(fn, "w").write(s)
    print "saving histogram file %s, total %s" % (fn, tot)



def mkhist2(s, fnout):
  n = len(s)
  hist = None
  for i in range(n):
    x = s[i].split()
    try:
      ene = float(x[1])
    except:
      break
    if i == 0:
      hist = Hist(ene, ene, dE)
    hist.add(ene)
  if hist:
    hist.save(fnout)



def mkhist3(s, fnout):
  n = len(s)
  hist = None
  for i in range(n):
    x = s[i].split()
    try:
      tp = float(x[2])
    except:
      break
    if tp < T - dT or tp > T + dT:
      continue
    ene = float(x[1])
    if not hist:
      hist = Hist(ene, ene, dE)
    hist.add(ene)
  if hist:
    hist.save(fnout)



def mkhist(fnin):
  global fnout
  s = open(fnin).readlines()
  if not fnout:
    fnout = os.path.splitext(fnin)[0] + ".his"

  tp = len( s[0].split() )
  if tp == 3:
    mkhist3(s, fnout)
  else:
    mkhist2(s, fnout)



if __name__ == "__main__":
  doargs()
  mkhist(fninp)
