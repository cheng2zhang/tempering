#!/usr/bin/env python



import os, sys, getopt, glob
from math import *



fnin = "spec.log"
fnout = ""


def showhelp():
  print "calcspec.py [Options] spec.log"
  print "Description:"
  print "  Compute quantities from the coordinats of special atoms"
  print "Options:"
  exit(1)



def doargs():
  global fnin, fnout

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "ho:",
        ["help","output="])
  except:
    print "Error parsing the command line"
    sys.exit(1)

  for o, a in opts:
    if o in ("-h", "--help",):
      showhelp()
    elif o in ("-o", "--output",):
      fnout = a

  if len(args) == 0:
    fnins = glob.glob("spec*.log")
    if not fnins: showhelp()
    fnin = fnins[0]
  else:
    fnin = args[0]


def vdiff(a, b):
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

def vsqr(a):
  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]

def vnorm(a):
  return sqrt(vsqr(a))

def vdist(a, b):
  return vnorm(vdiff(a, b))

def ncontact(arr, cutoff):
  # compute the number of contact
  n = len(arr)
  nc = 0
  for i in range(1, n):
    for j in range(i):
      dis = vdist(arr[i], arr[j])
      if dis < cutoff:
        nc += 1
  return nc


def calcspec(fn, fnout):
  i = 0
  sout = ""
  for s in open(fn).readlines():
    tok = s.strip().split("\t")
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
    nc = ncontact(arr, 10)
    i += 1
    sout += "%s\t%s\t%s\n" % (step, nc, "\t".join(tok[j+1:]))

  if not fnout:
    fnout = os.path.splitext(fn)[0] + "_nc.log"
  open(fnout, "w").write(sout)
  print "saving results to %s" % fnout
  
  

if __name__ == "__main__":
  doargs()
  calcspec(fnin, fnout)
