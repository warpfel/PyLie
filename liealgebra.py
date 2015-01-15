import sympy
import numpy as np

#Lie Algebra Class for Complex Lie-Algebras 

class LieStructureConstantTable(object):    
  def __init__(self, dimv):
    self.Dim = dimv #Dimension of the Algebra
    self.List = [[[0 for i in range(self.Dim)] for i in range(self.Dim)] for i in range(self.Dim)] #List with LieStructureConstants
    self.LieFlag = False #True if structure constants fullfilling Lie Algebra Axioms
  
  def __eq__(self, StrCons2):
    flag = True
    for i in range(self.Dim):
      for j in range(self.Dim):
	for k in range(self.Dim):
	  if self.List[i][j][k] != StrCons2.List[i][j][k]:
	    print self.List[i][j][k], StrCons2.List[i][j][k]
	    return False
    return True
  
  def SetEntry(self, i, j, ck): #For [X_i,X_j] = Sum_i ck[i] X_[i], [X_j, X_i] is defined via the negative, too
    if np.array(ck).shape == np.array(self.List[i][j]).shape :
      for k in range(self.Dim):
	self.List[i][j][k] = ck[k]
	self.List[j][i][k] = -ck[k]
    else :
      print("Shape-error")
      
  def AutoComplete(self): #Auto Completion, not specified structure constants are zero 
    zeroarray = [0 for i in range(self.Dim)]
    for i in range(self.Dim):
      self.SetEntry(i,i, zeroarray)
      
  def TestJacobiID(self): #Tests if the StructureConstantTable fulfills the Jacobi Identity:
    FlagJacobi = True
    for i in range(self.Dim):
      for k in range(self.Dim):
	for l in range(self.Dim):
	  for m in range(self.Dim):
	    jacobitestsum = 0
	    for n in range(self.Dim): #Summation
	      jacobitestsum = jacobitestsum + self.List[l][m][n]* self.List[k][n][i] + self.List[m][k][n]*self.List[l][n][i] + self.List[k][l][n]*self.List[m][n][i]
	    if sympy.simplify(jacobitestsum) == 0:
	      FlagJacobi = True
	    else :
	      FlagJacobi = False
	      break
	  if (FlagJacobi == False): break
	if (FlagJacobi == False): break
      if (FlagJacobi == False): break
    return FlagJacobi
  
  def TestAntiSymmetry(self): #The Antisymmetry for non diagonal Elements is implicitly given via SetEntry. This Tests if DiagonalElements are zero
    FlagDiag = True
    zeroarray = np.array([0 for i in range(self.Dim)])
    for i in range(self.Dim):
      if(np.array_equal(np.array(self.List[i][i]), zeroarray)):
	   FlagDiag = True
      else :
	FlagDiag = False
	break
      return FlagDiag
    
  def IsLie(self): #Tests if the StructureConstants.List defines a Lie-Algebra
    self.LieFlag = self.TestJacobiID() and self.TestAntiSymmetry()
    return self.LieFlag
  
class LieAlgebra(object):
  def __init__(self, dimv, Namev=""):
    self.Name = Namev
    self.Dim = dimv
    self.StructureConstants = LieStructureConstantTable(self.Dim)
    self.UnimodularFlag = False
    self.SemisimpleFlag = False
    self.SolvableFlag = False
    self.NilpotentFlag = False
    self.AdjReprIsFaithfulFlag = False
  
  """Constructing Lie Algebras: """
  
  def LieAlgebraByStructureConstants(self, SCTable): #Defining the Lie-Algebra via a LieStructureConstantTable SCTable
    self.StructureConstants = SCTable
    self.StructureConstants.AutoComplete()
    if not self.StructureConstants.IsLie():
      print("this is not a Lie Algebra, StructureConstants set to zero")
      self.StructureConstants = LieStructureConstantTable(self.Dim)
    else :
      print("%s is a Lie Algebra"% self.Name)

  def LinearLieAlgebraByMatrices(self, *Basev):
    Base = []
    for i in Basev:
      Base.append(sympy.Matrix(i))
    InterRepresentation = LieRepresentation(self.Dim, Base[0].shape[0], *Base)
    for i in range(self.Dim):
      for j in range(self.Dim):
	base1= [0 for i1 in range(self.Dim)]
	base2= [0 for i2 in range(self.Dim)]
	base1[i] = 1
	base2[j] = 1
	Comm = InterRepresentation.CoordinateCommutator(base1, base2)
	for k in range(self.Dim):
	  self.StructureConstants.List[i][j][k] = Comm[k]
    """
    if len(Basev) != self.Dim: return False
    #Testing if Basev is really a Base
    Base = []
    for i in Basev:
      Base.append(sympy.Matrix(i))
    M1 = sympy.zeros(Base[0].shape[0]*Base[0].shape[1], len(Base))
    for i in range(len(Base)):
      rownumber = 0
      for j in range(Base[i].shape[0]):
	for k in range(Base[i].shape[1]):
	  M1[rownumber,i] = Base[i][j,k]
	  rownumber = rownumber + 1
    if M1.nullspace() == []: #All are linear independent, now calculating the Structure Constants
      M2 = sympy.zeros(Base[0].shape[0]*Base[0].shape[1], len(Base)+1)
      for i in range(len(Base)):
	rownumber = 0
	for j in range(Base[i].shape[0]):
	  for k in range(Base[i].shape[1]):
	    M2[rownumber,i] = Base[i][j,k]
	    rownumber = rownumber + 1
      for base1 in range(len(Base)):
	for base2 in range(base1, len(Base)):
	  Commutator = Base[base1]*Base[base2]-Base[base2]*Base[base1]
	  rownumber = 0
	  for j in range(Base[i].shape[0]):
	    for k in range(Base[i].shape[1]):
	      #print base1, base2, j, k
	      M2[rownumber,len(Base)] = Commutator[j,k]
	      rownumber = rownumber + 1
	  l = sympy.symbols('x:%d' % self.Dim)
	  result = sympy.solve_linear_system(M2, *l)
	  vk = []
	  for k in range(self.Dim):
	     vk.append(result[l[k]])
	  self.StructureConstants.SetEntry(base1, base2, vk)
      self.StructureConstants.AutoComplete()
      print self.StructureConstants.List
      print "bla"
      self.PrintNonTrivialBrackets()
      return True
      """
  
  #def SimpleLieAlgebra():
  
  """The Adjoint Representation and the KillingForm:"""
  
  def AdjointMatrix(self, a): #Computes the representation matrix of the adjoint representation in the base {x_i}
    adbasematrix = [0 for i in range(self.Dim)] #list of base matrices ad_(x_i)
    AdjReprMatr = sympy.Matrix([[0 for i in range(self.Dim)] for i in range(self.Dim)])
    for i in range(self.Dim):
      adbasematrix[i]= sympy.Matrix(self.StructureConstants.List[i]).T
    for k in range(self.Dim):
      AdjReprMatr = AdjReprMatr + a[k]*adbasematrix[k]
    return AdjReprMatr  
  
  def KillingMatrix(self): #Computes the Coordinaterepresentation matrix of the killing form
    killing = [[0 for p in range(self.Dim)]for q in range(self.Dim)]
    for p in range(self.Dim):
      for q in range(self.Dim):
	for l in range(self.Dim):
	  for m in range(self.Dim):
	    killing[p][q] = killing[p][q] + self.StructureConstants.List[p][l][m]*self.StructureConstants.List[q][m][l]
    killing = sympy.Matrix(killing)
    return killing
  
  def IsAdjReprFaithful(self): #Tests if the adj. representation is faithful
    if self.LieCenter() == []:
      return True
    else:
      return False
 
  def AdjointRepresentation(self): #Returns the AdjointRepresentation as a LieRepresentation
    m = []
    for i in range(self.Dim):
      basev = [0 for j in range(self.Dim)]
      basev[i] = 1
      m.append(self.AdjointMatrix(basev))
    return LieRepresentation(self.Dim, self.Dim, *m)
    
  """Other Calculations:"""
  
  def LieBracket(self, a, b): # Computes LieBracket of [sum a[i]*x_i, sum b[i]*x[i]] for the base {x_i}
    brac = [0 for i in range(self.Dim)]
    for i in range(self.Dim):
      for j in range(self.Dim):
	for k in range(self.Dim):
	  brac[k] = brac[k] + a[i]*b[j]*self.StructureConstants.List[i][j][k]
    return brac
    
  def KillingForm(self, a, b): #Computes k(a, b), a, b should be lists
    M = self.KillingMatrix()
    result = 0
    for p in range(self.Dim):
      for q in range(self.Dim):
	result = result + a[p]*b[q]*M[p,q]
    return result
    
  def GLV_Action(self, a): #Computes the Action of an Element of GL(V) on the Variety of Structure Constants
    M = sympy.Matrix(a)
    if M.det == 0:
      print "no GL(V) Matrix"
    M_inv  = M**(-1)
    new = LieStructureConstants(4)
    for i in range(self.Dim):
      for j in range(self.Dim):
	for k in range(self.Dim):
	    for l in range(self.Dim):
	      for m in range(self.Dim):
		for n in range(self.Dim):
		  new.List[i][j][k] = new.List[i][j][k] + M_inv[k,n]*M[l,i]*M[m,j]*self.StructureConstants.List[l][m][n]
    self.LieAlgebraByStructureConstants(new)    
  
  def NonNilpotentElement(self): #Computes a NonNilpotentElement
    zerolist = [0 for i in range(self.Dim)]
    x = [sympy.zeros(self.Dim,self.Dim) for i in range(self.Dim)]
    for i in range(self.Dim):
      baseelement = zerolist
      baseelement[i] = 1
      x[i] = self.AdjointMatrix(baseelement)
      if not x[i].is_nilpotent():
	return x[i]
    for i in range(self.Dim):
      for j in range(self.Dim):
	if not (x[i]+x[j]).is_nilpotent():
	  return x[i]+x[j]
    return 0
  
  def NonNilpotentElementSubAlgebra(self, Base): #Computes a NonNilpotentElement
    zerolist = [0 for i in range(self.Dim)]
    x = [sympy.zeros(self.Dim,self.Dim) for i in range(self.Dim)]
    for i in range(len(Base)):
      baseelement = Base[i]
      x[i] = self.AdjointMatrix(baseelement)
      if not x[i].is_nilpotent():
	return x[i]
    for i in range(self.Dim):
      for j in range(self.Dim):
	if not (x[i]+x[j]).is_nilpotent():
	  return x[i]+x[j]
    return 0
  
  def EqualSubSpace(self, Base1, Base2): #Tests if the Subspaces spanned by Base1 and by Base2 are equal
    M1 =sympy.zeros(self.Dim, len(Base2)+1)
    l1 = sympy.symbols('l:%d' % self.Dim)
    M2 =sympy.zeros(self.Dim, len(Base1)+1)
    l2 = sympy.symbols('g:%d' % self.Dim)
    for i in range(len(Base1)):
      for k in range(self.Dim):
	for j in range(len(Base2)):
	  M1[k,j] = Base2[j][k]
	M1[k, len(Base2)] = Base1[i][k]
      result1 = sympy.solve_linear_system(M1, *l1)
      if result1 == None:
	return False
    for i in range(len(Base2)):
      for k in range(self.Dim):
	for j in range(len(Base1)):
	  M2[k,j] = Base1[j][k]
	M2[k, len(Base1)] = Base2[i][k]
      result2 = sympy.solve_linear_system(M2, *l2)
      if result2 == None:
	return False
    return True
  
  """Constructing Distinct Subalgebras, Subspaces and other shiny tings:"""
  
  #def LieNilradical(self):
    
  #def LieSolvableRadical(self):
  
  def FittingNullComponent(self, *nilpotentbase): #Computes FittingNullComponent in respect to a nilpotent sub-algebra spanned by nilpotentbase (coordinate-vector)
    base = []
    zerolist = [0 for i in range(self.Dim)]
    i = 1
    dim = 0
    dimfittingone = len(self.FittingOneComponent(*nilpotentbase))
    while dim < self.Dim - dimfittingone: #Abbruchbedingung: Wenn die Dimension von L0 + Dimension von L1 = dimension der lie algebra ist ist L0 fertig berechnet wg. direkte summe und so
      for x in nilpotentbase:
	nspace = ((self.AdjointMatrix(x))**i).nullspace() #Berechnung des Kern(adx^i)
	for j in nspace: #Anhaengen des Kern(adx^i) an das Erzeugendensystem
	  base.append(j) 
	Matrix = sympy.zeros(len(base), self.Dim)
	#Berechnung der rref um das Erzeugendensystem zu einer Basis zu machen
	for j1 in range(len(base)):
	  for j2 in range(self.Dim):
	    Matrix[j1, j2] = base[j1][j2]
	Matrix = Matrix.rref()[0]
	#Umwanldung der Matrix in reduced row echelon form in eine Liste aus sympy vektoren
	for j1 in range(len(base)):
	  newvector = []
	  for j2 in range(self.Dim):
	    newvector.append(Matrix[j1,j2])
	  if newvector != zerolist:
	    base.append(sympy.Matrix(newvector))
	dim = len(base) #dimension update 
      i = i+1 #exponenten updaten

    #nochmal Berechnung der reduced row echelon Form um das Erzeugendensystem schlussendlich zu einer Basis zu machen
    Matrix = sympy.zeros(len(base), self.Dim)
    for j1 in range(len(base)):
      for j2 in range(self.Dim):
	Matrix[j1, j2] = base[j1][j2]  
    Matrix = Matrix.rref()[0]
    #Umwandlung der rREF in eine Liste aus sympy vektoren
    newbase = []
    for j1 in range(len(base)):
      newvector = []
      for j2 in range(self.Dim):
	newvector.append(Matrix[j1,j2])
      if newvector != zerolist:
	newbase.append(sympy.Matrix(newvector))  
    #Umwandlung in eine normale Python Liste, da das der Standarduebergabewert in diesen Lie-Algebren klassen ist
    listbase = []
    for j1 in range(len(newbase)):
      newvec = []
      for j2 in range(self.Dim):
	newvec.append(newbase[j1][j2])
      listbase.append(newvec)
    return listbase
      
  def FittingOneComponent(self, *nilpotentbase): #Computes FittingOneComponent in respect to a nilpotent sub-algebra spanned by nilpotentbase (coordinate-vector)
    LieBase = []
    for i in range(self.Dim):
      zerolist = [0 for j in range(self.Dim)]
      zerolist[i] = 1
      LieBase.append(zerolist)
    OldBase = LieBase
    while True:
      NewBase = self.Productspace(nilpotentbase, LieBase)
      if self.EqualSubSpace(OldBase, NewBase):
	return OldBase
      OldBase = NewBase
     
  def CartanSubalgebra(self): #Noch unfertig
    CartanBase = []
    if self.IsLieNilpotent() == True:
      for i in range(self.Dim):
	basevec = [0 for j in range(self.Dim)]
	basevec[i] = 1
	CartanBase.append(basevec)
      return CartanBase
    x = self.NonNilpotentElement() 
    while True:
      print "a"
      P = self.FittingNullComponent(x)
      y = self.NonNilpotentElementSubAlgebra(P)
      if y == 0:
	return P
      for i in range(-self.Dim, self.Dim+1):
	print i
	print len(self.FittingNullComponent(x + i*(y-x)))
	if len(self.FittingNullComponent(x + i*(y-x))) < len(self.FittingNullComponent(x)):
	  x = x+i(y-x)
	  break
    
  #def QuotientAlgebra(self):
  
  def LieCentralizer(self, *Base): #Computes the Centralizer of a Subspace spanned by Base
    M = sympy.zeros(self.Dim*len(Base), self.Dim)
    for i in range(self.Dim):
      rowcounter = 0
      for l in range(len(Base)):
	for k in range(self.Dim):
	  lsum = 0
	  for j in range(self.Dim):
	    lsum = lsum + Base[l][j]*self.StructureConstants.List[i][j][k]
	  M[rowcounter, i] = lsum
	  rowcounter = rowcounter + 1
    result = M.nullspace()
    return result

  def LieCenter(self): #Computes the Center of the Lie Algebra
    Base = []
    for i in range(self.Dim):
      NewElement = []
      for j in range(self.Dim):
	if i == j:
	  NewElement.append(1)
	else:
	  NewElement.append(0)
      Base.append(NewElement)
    return self.LieCentralizer(*Base)
    
  def LieNormalizer(self, *Base): #Computes the Norlmalizer of a Subspace spanned by Base
    M = sympy.zeros(self.Dim*len(Base), self.Dim+len(Base)*len(Base))
    rowcounter = 0
    for l in range(len(Base)):
      for k in range(self.Dim):
	for i in range(self.Dim):
	  jsum = 0
	  for j in range(self.Dim):
	    jsum = jsum + Base[l][j]*self.StructureConstants.List[i][j][k]
	  M[rowcounter, i] = jsum
	for m in range(0, len(Base)):
	  M[rowcounter, l*len(Base)+m +self.Dim] = Base[m][k]  
	rowcounter = rowcounter + 1
    result = M.nullspace()
    resultlist = []
    for i in result:
      for k in range(len(Base)*len(Base)):
	i.row_del(self.Dim + k )
      newvec = []
      for k in range(self.Dim):
	newvec.append(i[k])
      resultlist.append(newvec)
    return resultlist
  
  def Productspace(self, base1, base2): #Computes the Productspace [<base1>,<base2>] of two subspaces spanned by base1 and base2
    newbase = []
    zerolist = [0 for i in range(self.Dim)]
    Matrix = sympy.zeros(len(base1)*len(base2), self.Dim)
    for i in base1:
      for j in base2:
	newbase.append(self.LieBracket(i, j))
    for i in range(len(base1)*len(base2)):
      for j in range(self.Dim):
	Matrix[i,j] = newbase[i][j]
    Matrix = Matrix.rref()[0]
    newbase = []
    for i in range(len(base1)*len(base2)):
      newvector = []
      for j in range(self.Dim):
	newvector.append(Matrix[i,j])
      if newvector != zerolist:
	newbase.append(newvector)
    return newbase
    
  def LowerCentralSeries(self, k): #Computes the Lower Central series L(k) of the liealg g: L(k=1) g = g, L(k= i+1)g = [g, L(i)g]
    fullbase = []
    for i in range(self.Dim):
      zerolist = [0 for j in range(self.Dim)]
      zerolist[j] = 1
      fullbase.append(zerolist)
    productbase = fullbase
    for i in range(k-1):
      productbase = self.Productspace(fullbase, productbase)
    return productbase
      
  def DerivedSeries(self, k): #Computes the derived series D(k) of the liealg g: D(k=1) g = g, D(k=i+1)g = [L(i)g, L(i)g]
    fullbase = []
    for i in range(self.Dim):
      zerolist = [0 for j in range(self.Dim)]
      zerolist[j] = 1
      fullbase.append(zerolist)
    productbase = fullbase
    for i in range(k-1):
      productbase = self.Productspace(productbase, productbase)
    return productbase
  
  #def DirectSumDecomposition(self):
  
  """Properties of the Lie-Algebra"""
  
  def IsLie(self): #Tests if the Lie Algebra is really a Lie Algebra
    return self.StructureConstants.IsLie()

  def IsLieNilpotent(self): #Tests if the Lie Algebra is nilpotent
    if self.NonNilpotentElement() == 0:
      return True
    else:
      return False
    """
    zerolist = [0 for i in range(self.Dim)]
    x = [sympy.zeros(self.Dim,self.Dim) for i in range(self.Dim)]
    for i in range(self.Dim):
      baseelement = zerolist
      baseelement[i] = 1
      x[i] = self.AdjointMatrix(baseelement)
      if not x[i].is_nilpotent():
	return False
    for i in range(self.Dim):
      for j in range(self.Dim):
	if not (x[i]+x[j]).is_nilpotent():
	  return False
    return True
    """
  
  def IsLieSolvable(self): #Tests if the Lie ALgebra is solvable
    for i in range(self.Dim):
      for j in range(self.Dim):
	for k in range(self.Dim):
	  x = [0 for p in range(self.Dim)]
	  x[i] = 1
	  y = [0 for q in range(self.Dim)]
	  y[j] = 1
	  z = [0 for r in range(self.Dim)]
	  z[k] = 1
	  if(self.KillingForm(self.LieBracket(x,y),z) != 0):
	    return False
    return True    
  
  def IsLieAbelian(self): #Tests if the Lie Algebra is abelian
    if self.StructureConstants == LieStructureConstantTable(self.Dim):
      return True
    return False
  
  def IsLieSemiSimple(self): #Tests if the Lie Algebra is semisimple
    if self.KillingMatrix().nullspace() == []:
      return True
    else:
      return False 
    
  def IsLieUnimodular(self): #Tests if the Lie Algebra is unimodular
    for i in range(self.Dim):
      x = [0 for p in range(self.Dim)]
      x[i] = 1
      #print self.ComputeAdjReprMatr(x)
      if((self.AdjointMatrix(x).trace())!=0):
	return False
    return True  
  
  """
  def IsLieCompact(self): #Tests if the Lie Algebra is Compact, definition Def12.1.1 in "Structure and Geometry of Lie Groups", Hilgert, Neeb
    A = sympy.zeros(self.Dim, self.Dim)
    M1 = sympy.zeros(self.Dim*self.Dim*self.Dim, self.Dim*self.Dim)
    M2 = sympy.zeros(self.Dim*self.Dim*self.Dim, self.Dim*self.Dim)
    rowcounter1 = 0
    colcounter1 = 0
    for a in range(self.Dim):
      for c in range(self.Dim):
	colcounter1 = 0
	for b in range(self.Dim):
	  for l in range(self.Dim):
	    M1[rowcounter1, colcounter1] = self.StructureConstants.List[c][a][l]
	    colcounter1 = colcounter1 + 1
	  rowcounter1 = rowcounter1 + 1
    rowcounter2 = 0
    colcounter2 = 0
    for c in range(self.Dim):
      for b in range(self.Dim):
	colcounter2 = 0
	for a in range(self.Dim):
	  for l in range(self.Dim):
	    M2[rowcounter2, colcounter2] = self.StructureConstants.List[c][b][l]
	    colcounter2 = colcounter2 + 1
	  rowcounter2 = rowcounter2 + 1
    M = M1 + M2
    result = M.nullspace()
    return len(result)
    """
  """Presentation:"""
  
  def PrintNonTrivialBrackets(self): #Prints all non trivial lie brackets
    for i in range(self.Dim):
      for j in range(i):
	zerolist = [0 for k in range(self.Dim)]
	if self.StructureConstants.List[i][j] != zerolist:
	  print i, j, self.StructureConstants.List[i][j] 
	  
  def PrintProperties(self): 
    print "#############"
    if not self.Name == "":
      print "Properties of Lie-Algebra %s" % self.Name
    self.PrintNonTrivialBrackets()
    if self.IsLieNilpotent():
      print "nilpotent"
      self.NilpotentFlag = True
    else:
      print "not nilpotent"
      self.NilpotentFlag = False
    if self.IsLieUnimodular():
      print "unimodular"
      self.UnimodularFlag = True
    else:
      print "not unimodular"
      self.UnimodularFlag = False
    if self.IsLieSolvable():
      print "solvable"
      self.SolvableFlag = True
    else:
      print "not solvable"
      self.SolvableFlag = False
    if self.IsLieSemiSimple():
      print "semisimple"
      self.SemisimpleFlag = True
    else:
      print "not semisimple"
      self.SemisimpleFlag = False
    if self.IsAdjReprFaithful():
      print "Adjoint Representation is faithful"
      AdjReprIsFaithfulFlag = True
    else:
      print "Adjoint Representation is not faithful"
      AdjReprIsFaithfulFlag = False
    print "#############"
 
class LieRepresentation(object):
  def __init__(self, LieDimv, VDimv, *Matricesv):
    self.LieDim = LieDimv #Dimension of the Lie Algebra
    self.VDim = VDimv #Dimension fo the Vectorspace
    self.lineardependenceflag = False
    if len(Matricesv) != self.LieDim:
      print "Number of Matrices should be the dimension of the Lie-Algebra"
      print "Representation set to trivial Representation"
      self.BaseMatrices = [sympy.zeros(self.VDim, self.VDim) for i in range(self.LieDim)]
    else:
      self.BaseMatrices = []
      for m in Matricesv:
	self.BaseMatrices.append(sympy.Matrix(m))
      for m in self.BaseMatrices:
	if m.shape != (self.VDim, self.VDim):
	  print "Matrix has wrong shape"
	  print "Representation set to trivial Representation"
	  self.BaseMatrices = [sympy.zeros(self.VDim, self.VDim) for i in range(self.LieDim)]
	  break
    if not self.IsBase():
      print "The Matrices are not linear independent"
      self.lineardependenceflag = True
      
  def IsBase(self): #Tests if the Matrices in BaseMatrices are linear independent, tests also if it is faithful
    TestBase = self.BaseMatrices
    M1 = sympy.zeros(TestBase[0].shape[0]*TestBase[0].shape[1], len(TestBase))
    for i in range(len(TestBase)):
      rownumber = 0
      for j in range(TestBase[i].shape[0]):
	for k in range(TestBase[i].shape[1]):
	  M1[rownumber,i] = TestBase[i][j,k]
	  rownumber = rownumber + 1    
    if M1.nullspace() == []: 
      return True
    else:
      return False
    
  def RepresentationMatrix(self, a):
    if len(a) == self.LieDim:
      A = sympy.zeros(self.VDim, self.VDim)
      for i in range(self.LieDim):
	A = A+a[i]*self.BaseMatrices[i]
      return A
    else:
	print "DimensionError"
    return False
    
  def MatrixCommutator(self, a, b):
    A = self.RepresentationMatrix(a)
    B = self.RepresentationMatrix(b)
    return A*B-B*A
  
  def CoordinateCommutator(self, a, b):
    if not self.lineardependenceflag:
      M2 = sympy.zeros(self.VDim*self.VDim, self.LieDim +1)
      for i in range(self.LieDim):
	rownumber = 0
	for j in range(self.VDim):
	  for k in range(self.VDim):
	    M2[rownumber,i] = self.BaseMatrices[i][j,k]
	    rownumber = rownumber + 1
      Commutator = self.MatrixCommutator(a,b)
      rownumber = 0
      for j in range(self.VDim):
	for k in range(self.VDim):
	  M2[rownumber,self.LieDim] = Commutator[j,k]
	  rownumber = rownumber + 1
      l = sympy.symbols('x:%d' % self.LieDim)
      result = sympy.solve_linear_system(M2, *l)
      vk = []
      for k in range(self.LieDim):
	vk.append(result[l[k]])
      return vk
    else:
      print "Commutator is not unique, base is not linear independent"
  
  def IsRepresentationLie(self, LieAlg):
    if not self.lineardependenceflag:
      TestLieAlg = LieAlgebra(self.LieDim, "TestLieAlgebra")
      TestLieAlg.LinearLieAlgebraByMatrices(*self.BaseMatrices)
      if TestLieAlg.StructureConstants == LieAlg.StructureConstants:
	return True
      else:
	return False
    else:
      leftside = sympy.zeros(self.VDim, self.VDim)
      rightside = sympy.zeros(self.VDim, self.VDim)
      for i in range(self.LieDim):
	for j in range(self.LieDim):
	  v1 = [0 for i1 in range(self.LieDim)]
	  v2 = [0 for i2 in range(self.LieDim)]
	  v1[i] = 1
	  v2[j] = 1
	  leftside = self.MatrixCommutator(v1, v2)
	  rightside = sympy.zeros(self.VDim, self.VDim)
	  for k in range(self.LieDim):
	    rightside = rightside + LieAlg.StructureConstants.List[i][j][k]*self.BaseMatrices[k]
	  if leftside != rightside:
	    return False
      return True
  
  def IsFaithful(self, LieAlg):
    return self.IsBase() and self.IsRepresentationLie(LieAlg)

  def PrintBase(self):
    for i in self.BaseMatrices:
      print i
  