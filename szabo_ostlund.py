#!/usr/bin/env python
import math

pi = math.pi

def HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
  """
  DOES A HARTREE-FOCK CALCULATION FOR A TWO-ELECTRON DIATOMIC
  USING THE 1S MINIMAL STO-NG BASIS SET
  MINIMAL BASIS SET HAS BASIS FUNCTIONS 1 AND 2 ON NUCLEI A AND B

  IOP = 0 NO PRINTING WHATSOEVER (TO OPTIMIZE EXPONENTS, SAY)
  IOP = 1 PRINT ONLY CONVERGED RESULTS
  IOP = 2 PRINT EVERY INTERATION
  N       STO-NG CALCULATION (N=1,2 OR 3)
  R       BONDLENGTH (AU)
  ZETA1   SLATER ORBITAL EXPONENT (FUNCTION 1)
  ZETA2   SLATER ORBITAL EXPONENT (FUNCTION 2)
  ZA      ATOMIC NUMBER (ATOM A)
  ZB      ATOMIC NUMBER (ATOM B)
  """
  if (IOP != 0):
    print("STO-{} FOR ATOMIC NUMBERS {} AND {}".\
    format(N, ZA, ZB))

  # CALCULATE ALL THE ONE AND TWO-ELECTRON INTEGRALS
  S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111,\
      V2111, V2121, V2211, V2221, V2222 = INTGRL(IOP, N, R, ZETA1,\
      ZETA2, ZA, ZB)

  # BE INEFFICIENT AND PUT ALL INTEGRALS IN PRETTY ARRAYS
  S, H, X, XT, TT = COLECT(IOP, T11, T22, T12, V11A, V11B, V12A, V12B,\
      V22A, V22B, V1111, V2111, V2121, V2211, V2221, V2222, S12)

  # PERFORM THE SCF CALCULATION
  SCF(IOP, R, ZA, ZB, TT, H, X, S, XT)

def INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
  """
  CALCULATES ALL THE BASIC INTEGRALS NEEDED FOR SCF CALCULATION
  """

  # THESE ARE THE CONTRACTION COEFFICIENTS AND EXPONENTS FOR
  # A NORMALIZED 1S SLATER ORBITAL WITH EXPONENT 1.0 IN TERMS OF
  # NORMALIZED 1S PRIMITIVE GAUSSIANS
  COEF = [[1, 0.678914, 0.444635], [0, 0.430129, 0.535328],\
          [0, 0, 0.154329]]
  EXPON = [[0.270950, 0.151623, 0.109818], [0, 0.851819, 0.405771],\
          [0, 0, 2.22766]]
  D1 = [0]*3
  A1 = [0]*3
  D2 = [0]*3
  A2 = [0]*3
  R2 = R**2

  # SCALE THE EXPONENTS (A) OF PRIMITIVE GAUSSIANS
  # INCLUDE NORMALIZATION IN CONTRACTION COEFFICIENTS (D)
  for I in range(N):
    A1[I] = EXPON[I][N-1]*ZETA1**2
    D1[I] = COEF[I][N-1]*(2*A1[I]/pi)**0.75
    A2[I] = EXPON[I][N-1]*ZETA2**2
    D2[I] = COEF[I][N-1]*(2*A2[I]/pi)**0.75

  # D AND A ARE NOW CONTRACTION COEFFICIENTS AND EXPONENTS
  # IN TERMS OF UNNORMALIZED PRIMITIVE GAUSSIANS
  S12 = T11 = T12 = T22 = V11A = V12A = V22A = V11B = V12B = V22B = V1111 =\
      V2111 = V2121 = V2211 = V2221 = V2222 = 0

  # CALCULATE ONE-ELECTRON INTEGRALS
  # CENTER A IS FIRST ATOM, CENTER B IS SECOND ATOM
  # ORIGIN IS ON CENTER A
  # V12A = OFF-DIAGONAL NUCLEAR ATTRACTION TO CENTER A, ETC.
  for I in range(N):
    for J in range(N):

      # RAP2 = SQUARED DISTANCE BETWEEN CENTER A AND CENTER P, ETC.
      RAP = A2[J]*R/(A1[I]+A2[J])
      RAP2 = RAP**2
      RBP2 = (R-RAP)**2
      S12 = S12+S(A1[I], A2[J], R2)*D1[I]*D2[J]
      T11 = T11+T(A1[I], A1[J], 0)*D1[I]*D1[J]
      T12 = T12+T(A1[I], A2[J], R2)*D1[I]*D2[J]
      T22 = T22+T(A2[I], A2[J], 0)*D2[I]*D2[J]
      V11A = V11A+V(A1[I], A1[J], 0, 0, ZA)*D1[I]*D1[J]
      V12A = V12A+V(A1[I], A2[J], R2, RAP2, ZA)*D1[I]*D2[J]
      V22A = V22A+V(A2[I], A2[J], 0, R2, ZA)*D2[I]*D2[J]
      V11B = V11B+V(A1[I], A1[J], 0, R2, ZB)*D1[I]*D1[J]
      V12B = V12B+V(A1[I], A2[J], R2, RBP2, ZB)*D1[I]*D2[J]
      V22B = V22B+V(A2[I], A2[J], 0, 0, ZB)*D2[I]*D2[J]

  # CALCULATE TWO-ELECTRON INTEGRALS
  for I in range(N):
    for J in range(N):
      for K in range(N):
        for L in range(N):
          RAP = A2[I]*R/(A2[I]+A1[J])
          RBP = R-RAP
          RAQ = A2[K]*R/(A2[K]+A1[L])
          RBQ = R-RAQ
          RPQ = RAP-RAQ
          RAP2 = RAP**2
          RBP2 = RBP**2
          RAQ2 = RAQ**2
          RBQ2 = RBQ**2
          RPQ2 = RPQ**2
          V1111 = V1111+TWOE(A1[I], A1[J], A1[K], A1[L], 0, 0, 0)*\
              D1[I]*D1[J]*D1[K]*D1[L]
          V2111 = V2111+TWOE(A2[I], A1[J], A1[K], A1[L], R2, 0, RAP2)*\
              D2[I]*D1[J]*D1[K]*D1[L]
          V2121 = V2121+TWOE(A2[I], A1[J], A2[K], A1[L], R2, R2, RPQ2)*\
              D2[I]*D1[J]*D2[K]*D1[L]
          V2211 = V2211+TWOE(A2[I], A2[J], A1[K], A1[L], 0, 0, R2)*\
              D2[I]*D2[J]*D1[K]*D1[L]
          V2221 = V2221+TWOE(A2[I], A2[J], A2[K], A1[L], 0, R2, RBQ2)*\
              D2[I]*D2[J]*D2[K]*D1[L]
          V2222 = V2222+TWOE(A2[I], A2[J], A2[K], A2[L], 0, 0 ,0)*\
              D2[I]*D2[J]*D2[K]*D2[L]
  
  if (IOP != 0):
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format("R", "ZETA1", "ZETA2", "S12", "T11"))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format(R, ZETA1, ZETA2, S12, T11))
    print("\n{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format("T12", "T22", "V11A", "V12A", "V22A"))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format(T12, T22, V11A, V12A, V22A))
    print("\n{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format("V11B", "V12B", "V22B", "V1111", "V2111"))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format(V11B, V12B, V22B, V1111, V2111))
    print("\n{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format("V2121", "V2211", "V2221", "V2222"))
    print("{:<20}\t{:<20}\t{:<20}\t{:<20}".\
        format(V2121, V2211, V2221, V2222))

  return S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111,\
      V2111, V2121, V2211, V2221, V2222

def F0(ARG):
  """
  CALCULATES THE F FUNCTION
  F0 ONLY (S-TYPE ORBITALS)
  """
  if (ARG >= 1e-6):

    # F0 IN TERMS OF THE ERROR FUNCTION
    F0 = math.sqrt(pi/ARG)*DERF(math.sqrt(ARG))/2

  else:

    # ASYMPTOTIC VALUE FOR SMALL ARGUMENTS
    F0 = 1-ARG/3

  return F0

def DERF(ARG):
  """
  CALCULATES THE ERROR FUNCTION ACCORDING TO A RATIONAL
  APPROXIMATION FROM M. ABRAMOWITZ AND I.A. STEGUN,
  HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER.
  ABSOLUTE ERROR IS LESS THAN 1.6*10**(-7)
  CAN BE REPLACED BY A BUILT-IN FUNCTION ON SOME MACHINES
  """
  DERF = math.erf(ARG)

  return DERF

def S(A, B, RAB2):
  """
  CALCULATES OVERLAPS FOR UN-NORMALIZED PRIMITIVES
  """
  S = (pi/(A+B))**1.5*math.exp(-A*B*RAB2/(A+B))

  return S

def T(A, B, RAB2):
  """
  CALCULATES KINETIC ENERGY INTEGRALS FOR UN-NORMALIZED PRIMITVES
  """
  T = A*B/(A+B)*(3-2*A*B*RAB2/(A+B))*(pi/(A+B))**1.5*math.exp(-A*B*RAB2/(A+B))

  return T

def V(A, B, RAB2, RCP2, ZC):
  """
  CALCULATES UN-NORMALIZED NUCLEAR ATTRACTION INTEGRALS
  """
  V = 2*pi/(A+B)*F0((A+B)*RCP2)*math.exp(-A*B*RAB2/(A+B))
  V = -V*ZC

  return V

def TWOE(A, B, C, D, RAB2, RCD2, RPQ2):
  """
  CALCULATES TWO-ELECTRON INTEGRALS FOR UN-NORMALIZED PRIMITIVES
  A, B, C, D ARE THE EXPONENTS ALPHA,BETA, ETC.
  RAB2 EQUALS SQUARED DISTANCE BETWEEN CENTER A AND CENTER B, ETC.
  """
  TWOE = 2*(pi**2.5)/((A+B)*(C+D)*math.sqrt(A+B+C+D))*F0((A+B)*(C+D)*RPQ2/\
      (A+B+C+D))*math.exp(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D))

  return TWOE

def COLECT(IOP, T11, T22, T12, V11A, V11B, V12A, V12B, V22A, V22B, V1111,\
    V2111, V2121, V2211, V2221, V2222, S12):
  """
  THIS TAKES THE BASCI INTEGRALS FROM COMMON AND ASSEMBLES THE
  RELEVANT MATRICES, THAT IS S,H,X,XT, AND TWO-ELECTRON INTEGRALS
  """
  S = [[0, 0], [0, 0]] 
  X = [[0, 0], [0, 0]]
  XT = [[0, 0], [0, 0]]
  H = [[0, 0], [0, 0]]
  F = [[0, 0], [0, 0]]
  G = [[0, 0], [0, 0]]
  C = [[0, 0], [0, 0]]
  FPRIME = [[0, 0], [0, 0]]
  CPRIME = [[0, 0], [0, 0]]
  P = [[0, 0], [0, 0]]
  OLDP = [[0, 0], [0, 0]]
  E = [[0, 0], [0, 0]]
  TT = [[[[0 for i in range(2)] for j in range(2)] for k in range(2)]\
          for k in range(2)]

  # FORM CORE HAMILTONIAN
  H[0][0] = T11+V11A+V11B
  H[0][1] = H[1][0] = T12+V12A+V12B
  H[1][1] = T22+V22A+V22B

  # FORM OVERLAP MATRIX
  S[0][0] = S[1][1] = 1
  S[0][1] = S[1][0] = S12

  # USE CANONICAL ORTHOGONALIZATION
  X[0][0] = X[1][0] = 1/math.sqrt(2*(1+S12))
  X[0][1] = 1/math.sqrt(2*(1-S12))
  X[1][1] = -X[0][1]

  # TRANSPOSE OF TRANSFORMATION MATRIX
  XT[0][0] = X[0][0]
  XT[0][1] = X[1][0]
  XT[1][0] = X[0][1]
  XT[1][1] = X[1][1]

  # MATRIX OF TWO-ELECTRON INTEGRALS
  TT[0][0][0][0] = V1111
  TT[1][0][0][0] = TT[0][1][0][0] = TT[0][0][1][0] = TT[0][0][0][1] = V2111
  TT[1][0][1][0] = TT[0][1][1][0] = TT[1][0][0][1] = TT[0][1][0][1] = V2121
  TT[1][1][0][0] = TT[0][0][1][1] = V2211
  TT[1][1][1][0] = TT[1][1][0][1] = TT[1][0][1][1] = TT[0][1][1][1] = V2221
  TT[1][1][1][1] = V2222
  if (IOP != 0):
    MATOUT(S, "S")
    MATOUT(X, "X")
    MATOUT(H, "H")
    print("")
    for I in range(2):
      for J in range(2):
        for K in range(2):
          for L in range(2):
            print("({} {} {} {}) {}".format(I, J, K, L, TT[I][J][K][L]))

  return S, H, X, XT, TT

def SCF(IOP, R, ZA, ZB, TT, H, X, S, XT):
  """
  PERFORMS THE SCF ITERATIONS
  """
  P = [[0, 0], [0, 0]]
  F = [[0, 0], [0, 0]]
  OLDP = [[0, 0], [0, 0]]

  # CONVERGENCE CRITERION FOR DENSITY MATRIX
  CRIT = 1e-4

  # MAXIMUM NUMBER OF ITERATIONS
  MAXIT = 25

  # ITERATION NUMBER
  ITER = 0

  # USE CORE-HAMILTONIAN FOR INITIAL GUESS AT F, I.E. (P=0)
  for I in range(2):
    for J in range(2):
      P[I][J] = 0

  if (IOP >= 2):
    MATOUT(P, "P")

  # START OF ITERATION LOOP
  while True:
    ITER = ITER+1
    if (IOP >= 2):
      print("\nSTART OF ITERATION NUMBER = {}".format(ITER))

    # FORM TWO-ELECTRON PART OF FOCK MATRIX FROM P
    G = FORMG(P, TT)
    if (IOP >= 2):
      MATOUT(G, "G")

    # ADD CORE HAMILTONIAN TO GET FOCK MATRIX
    for I in range(2):
      for J in range(2):
        F[I][J] = H[I][J]+G[I][J]

    # CALCULATE ELECTRONIC ENERGY
    EN = 0
    for I in range(2):
      for J in range(2):
        EN = EN+0.5*P[I][J]*(H[I][J]+F[I][J])

    if (IOP >= 2):
      MATOUT(F, "F")
      print("\nELECTRONIC ENERGY = {}".format(EN))

    # TRANSFORM FOCK MATRIX USING G FOR TEMPORARY STORAGE
    G = MULT(F, X)
    FPRIME = MULT(XT, G)

    # DIAGONALIZE TRANSFORMED FOCK MATRIX
    CPRIME, E = DIAG(FPRIME)

    # TRANSFORM EIGENVECTORS TO GET MATRIX C
    C = MULT(X, CPRIME)

    # FORM NEW DENSITY MATRIX
    for I in range(2):
      for J in range(2):

        # SAVE PRESENT DENSITY MATRIX
        # BEFORE CREATING NEW ONE
        OLDP[I][J] = P[I][J]
        P[I][J] = 0
        for K in range(1):
          P[I][J] = P[I][J]+2*C[I][K]*C[J][K]

    if (IOP >= 2):
      MATOUT(FPRIME, "F'")
      MATOUT(CPRIME, "C'")
      MATOUT(E, "E")
      MATOUT(C, "C")
      MATOUT(P, "P")

    # CALCULATE DELTA
    DELTA = 0
    for I in range(2):
      for J in range(2):
        DELTA = DELTA+(P[I][J]-OLDP[I][J])**2
    DELTA = math.sqrt(DELTA/4)
    if (IOP != 0):
      print("DELTA(CONVERGENCE OF DENSITY MATRIX) = {}".format(DELTA))

    # CHECK FOR CONVERGENCE
    if (DELTA >= CRIT):

      # NOT YET CONVERGEND
      # TEST FOR MAXIMUM NUMBER OF ITERATIONS
      # IF MAXIMUM NUMBER NOT YET REACHED THEN
      # GO BACK FOR ANOTHER ITERATION
      if (ITER >= MAXIT):
        break

    else:
      break

  # CALCULATION CONVERGEND IF IT GOT HERE
  # ADD NUCLEAR REPULSION TO GET TOTAL ENERGY
  ENT = EN+ZA*ZB/R
  if (IOP != 0):
    print("\nCALCULATION CONVERGED")
    print("ELECTRONIC ENERGY = {}".format(EN))
    print("TOTAL ENERGY = {}".format(ENT))

  if (IOP == 1):

    # PRINT THE FINAL RESULTS IF
    # HAVE NOT DONE SO ALREADY
    MATOUT(G, "G")
    MATOUT(F, "F")
    MATOUT(E, "E")
    MATOUT(C, "C")
    MATOUT(P, "P")

  # PS MATRIX HAS MULLIKEN POPULATIONS
  OLDP = MULT(P, S)
  if (IOP != 0):
    MATOUT(OLDP, "PS")

def FORMG(P, TT):
  """
  CALCULATES THE G MATRIX FROM THE DENSITY MATRIX
  AND TWO-ELECTRON INTEGRALS
  """
  G = [[0, 0], [0, 0]]
  for I in range(2):
    for J in range(2):
      for K in range(2):
        for L in range(2):
          G[I][J] = G[I][J]+P[K][L]*(TT[I][J][K][L]-0.5*TT[I][L][K][J])

  return G

def DIAG(F):
  """
  DIAGONALIZES F TO GIVE EIGENVECTORS IN C AND EIGENVALUES IN E
  THETA IS THE ANGLE DESCRIBING SOLUTION
  """
  C = [[0, 0], [0, 0]]
  E = [[0, 0], [0, 0]]
  if (abs(F[0][0]-F[1][1]) <= 1e-20):

    # HERE IS SYMMETRY DETERMINED SOLUTION (HOMONUCLEAR DIATOMIC)
    THETA = pi/4
  else:

    # SOLUTION FOR HETERONUCLEAR DIATOMIC
    THETA = 0.5*math.atan(2*F[0][1]/(F[0][0]-F[1][1]))

  C[0][0] = math.cos(THETA)
  C[1][0] = math.sin(THETA)
  C[0][1] = math.sin(THETA)
  C[1][1] = math.cos(THETA)
  E[0][0] = F[0][0]*math.cos(THETA)**2+F[1][1]*math.sin(THETA)**2+F[0][1]*\
      math.sin(2*THETA)
  E[1][1] = F[1][1]*math.cos(THETA)**2+F[0][0]*math.sin(THETA)**2-F[0][1]*\
      math.sin(2*THETA)
  E[1][0] = E[0][1] = 0

  # ORDER EIGENVALUES AND EIGENVECTORS
  if (E[1][1] <= E[0][0]):
    TEMP = E[1][1]
    E[1][1] = E[0][0]
    E[0][0] = TEMP
    TEMP = C[0][1]
    C[0][1] = C[0][0]
    C[0][0] = TEMP
    TEMP = C[1][1]
    C[1][1] = C[1][0]
    C[1][0] = TEMP

  return C, E

def MULT(A, B):
  """
  MULTIPLIES TWO SQUARE MATRICES A AND B TO GET C
  """
  C = [[0, 0], [0, 0]]
  for I in range(2):
      for J in range(2):
          for K in range(2):
              C[I][J] = C[I][J]+A[I][K]*B[K][J]

  return C

def MATOUT(A, LABEL):
  """
  PRINT MATRICES OF SIZE M BY N
  """
  print("\nTHE {} ARRAY".format(LABEL))
  max_width = [max(len(str(row[i])) for row in A) for i in range(len(A[0]))]
  for i, row in enumerate(A):
    print("\t".join(str(row[i]).ljust(max_width[i]) for i in range(len(row))))

"""
MINIMAL BASIS STO-3G CALCULATION ON HEH+
THIS IS A LITTLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC
"""
if __name__ == "__main__":
  IOP = 2
  N = 3
  R = 1.4632
  ZETA1 = 2.0925
  ZETA2 = 1.24
  ZA = 2
  ZB = 1
  HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)