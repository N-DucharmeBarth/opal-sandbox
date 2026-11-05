#!/bin/sh

#  ---------
#   PHASE 9
#  ---------

mfclo64 alb.frq 08.par 09.par -switch 2 1 1 200 2 116 30    # 200 iterations, max F=0.3
mfclo64 alb.frq 09.par 09.par -switch 2 1 1 2000 2 116 300  # 2000 iterations, max F=3.0

# Parameter report
mfclo64 alb.frq 09.par report -switch 2 1 1 1 1 246 1  # produce indepvar.rpt file

# Hessian analysis
mfclo64 alb.frq 09.par hessian -switch 1 1 145 1         # compute Hessian
mfclo64 alb.frq 09.par hessian -switch 1 1 145 5         # produce .cor file
mfclo64 alb.frq 09.par hessian -switch 2 1 37 1 1 145 2  # compute derivatives for reference points
mfclo64 alb.frq 09.par hessian -switch 1 1 145 4         # produce .var file
