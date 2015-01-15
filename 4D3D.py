from liealgebra import *


#L0  
L0Struct = LieStructureConstantTable(4)
L0Struct.AutoComplete()
L0 = LieAlgebra(4, "L0")
L0.LieAlgebraByStructureConstants(L0Struct)

#L1
L1Struct = LieStructureConstantTable(4)
L1Struct.SetEntry(0,1, [0,0,1,0])
L1Struct.AutoComplete()
L1 = LieAlgebra(4,"L1")
L1.LieAlgebraByStructureConstants(L1Struct)

#L2
L2Struct = LieStructureConstantTable(4)
L2Struct.SetEntry(0,1,[0,0,1,0])
L2Struct.SetEntry(0,2,[0,0,0,1])
L2Struct.AutoComplete()
L2 = LieAlgebra(4, "L2")
L2.LieAlgebraByStructureConstants(L2Struct)

#L3
L3Struct = LieStructureConstantTable(4)
L3Struct.SetEntry(0,1,[0,1,0,0])
L3Struct.SetEntry(0,2,[0,0,1,0])
L3Struct.SetEntry(0,3,[0,0,0,1])
L3 = LieAlgebra(4, "L3")
L3.LieAlgebraByStructureConstants(L3Struct)

#L4
a = sympy.Symbol('a')
L4Struct = LieStructureConstantTable(4)
L4Struct.SetEntry(0,1,[0,1,0,0])
L4Struct.SetEntry(0,2,[0,0,1,0])
L4Struct.SetEntry(0,3,[0,0,1,a])
L4 = LieAlgebra(4, "L4")
L4.LieAlgebraByStructureConstants(L4Struct)


#L4(infty)
L4inftyStruct = LieStructureConstantTable(4)
L4inftyStruct.SetEntry(0,1,[0,1,0,0])
L4infty = LieAlgebra(4, "L4infty")
L4infty.LieAlgebraByStructureConstants(L4inftyStruct)

#L5
L5Struct = LieStructureConstantTable(4)
L5Struct.SetEntry(0,1,[0,1,0,0])
L5Struct.SetEntry(0,2,[0,0,1,0])
L5Struct.SetEntry(0,3,[0,0,0,2])
L5Struct.SetEntry(1,2,[0,0,0,1])
L5 = LieAlgebra(4, "L5")
L5.LieAlgebraByStructureConstants(L5Struct)

#L6
L6Struct = LieStructureConstantTable(4)
L6Struct.SetEntry(0,1,[0,1,0,0])
L6Struct.SetEntry(0,2,[0,0,-1,0])
L6Struct.SetEntry(1,2,[1,0,0,0])
L6 = LieAlgebra(4, "L6")
L6.LieAlgebraByStructureConstants(L6Struct)

#L7
p = sympy.Symbol('p')
q = sympy.Symbol('q')
L7Struct = LieStructureConstantTable(4)
L7Struct.SetEntry(0,1,[0,1,0,0])
L7Struct.SetEntry(0,2,[0,1,p,0])
L7Struct.SetEntry(0,3,[0,0,1,q])
L7 = LieAlgebra(4, "L7")
L7.LieAlgebraByStructureConstants(L7Struct)

#L8
r = sympy.Symbol('r')
L8Struct = LieStructureConstantTable(4)
L8Struct.SetEntry(0,1,[0,1,0,0])
L8Struct.SetEntry(0,2,[0,1,r,0])
L8Struct.SetEntry(0,3,[0,0,0,r+1])
L8Struct.SetEntry(1,2,[0,0,0,1])
L8 = LieAlgebra(4, "L8")
L8.LieAlgebraByStructureConstants(L8Struct)

#L9
L9Struct = LieStructureConstantTable(4)
L9Struct.SetEntry(0,1,[0,1,0,0])
L9Struct.SetEntry(2,3,[0,0,0,1])
L9 = LieAlgebra(4, "L9")
L9.LieAlgebraByStructureConstants(L9Struct)
x1= L9.AdjointMatrix([1,0,0,0])
x2= L9.AdjointMatrix([0,1,0,0])


#Bianchi1
Bianchi1Struct = LieStructureConstantTable(3)
Bianchi1 = LieAlgebra(3, "B1")
Bianchi1.LieAlgebraByStructureConstants(Bianchi1Struct)

#Bianchi2
Bianchi2Struct = LieStructureConstantTable(3)
Bianchi2Struct.SetEntry(1,2,[1,0,0])
Bianchi2 = LieAlgebra(3, "B2")
Bianchi2.LieAlgebraByStructureConstants(Bianchi2Struct)

#Bianchi3
Bianchi3Struct = LieStructureConstantTable(3)
Bianchi3Struct.SetEntry(0,2,[1,0,0])
Bianchi3 = LieAlgebra(3, "B3")
Bianchi3.LieAlgebraByStructureConstants(Bianchi3Struct)

#Bianchi4
Bianchi4Struct = LieStructureConstantTable(3)
Bianchi4Struct.SetEntry(0,2,[1,0,0])
Bianchi4Struct.SetEntry(1,2,[1,1,0])
Bianchi4 = LieAlgebra(3, "B4")
Bianchi4.LieAlgebraByStructureConstants(Bianchi4Struct)

#Bianchi5
Bianchi5Struct = LieStructureConstantTable(3)
Bianchi5Struct.SetEntry(0,2,[1,0,0])
Bianchi5Struct.SetEntry(1,2,[0,1,0])
Bianchi5 = LieAlgebra(3, "B5")
Bianchi5.LieAlgebraByStructureConstants(Bianchi5Struct)

#Bianchi6
h = sympy.Symbol("h")
Bianchi6Struct = LieStructureConstantTable(3)
Bianchi6Struct.SetEntry(0,2,[1,0,0])
Bianchi6Struct.SetEntry(1,2,[0,h,0])
Bianchi6 = LieAlgebra(3, "B6")
Bianchi6.LieAlgebraByStructureConstants(Bianchi6Struct)

#Bianchi7_1
Bianchi7_1Struct = LieStructureConstantTable(3)
Bianchi7_1Struct.SetEntry(0,2,[0,1,0])
Bianchi7_1Struct.SetEntry(1,2,[-1,0,0])
Bianchi7_1 = LieAlgebra(3, "B7_1")
Bianchi7_1.LieAlgebraByStructureConstants(Bianchi7_1Struct)

#Bianchi7_2
g = sympy.Symbol("g")
Bianchi7_2Struct = LieStructureConstantTable(3)
Bianchi7_2Struct.SetEntry(0,2,[0,1,0])
Bianchi7_2Struct.SetEntry(1,2,[-1,g,0])
Bianchi7_2 = LieAlgebra(3, "B7_2")
Bianchi7_2.LieAlgebraByStructureConstants(Bianchi7_2Struct)

#Bianchi8
Bianchi8Struct = LieStructureConstantTable(3)
Bianchi8Struct.SetEntry(0,1,[1,0,0])
Bianchi8Struct.SetEntry(0,2,[0,2,0])
Bianchi8Struct.SetEntry(1,2,[0,0,1])
Bianchi8 = LieAlgebra(3, "B8")
Bianchi8.LieAlgebraByStructureConstants(Bianchi8Struct)
Bianchi8killing = Bianchi8.KillingMatrix()

#Bianchi9
Bianchi9Struct = LieStructureConstantTable(3)
Bianchi9Struct.SetEntry(0,1,[0,0,1])
Bianchi9Struct.SetEntry(1,2,[1,0,0])
Bianchi9Struct.SetEntry(2,0,[0,1,0])
Bianchi9 = LieAlgebra(3, "B9")
Bianchi9.LieAlgebraByStructureConstants(Bianchi9Struct)


print "Properties of 4D LieAlgebras"
L0.PrintProperties()
L1.PrintProperties()
L2.PrintProperties()
L3.PrintProperties()
L4.PrintProperties()
L4infty.PrintProperties()
L5.PrintProperties()
L6.PrintProperties()
L7.PrintProperties()
L8.PrintProperties()
L9.PrintProperties()

print "Properties of 3D LieAlgebras"
Bianchi1.PrintProperties()
Bianchi2.PrintProperties()
Bianchi3.PrintProperties()
Bianchi4.PrintProperties()
Bianchi5.PrintProperties()
Bianchi6.PrintProperties()
Bianchi7_1.PrintProperties()
Bianchi8.PrintProperties()
Bianchi9.PrintProperties()

print "Representations of 4D LieAlgebras"
L0Repr4D = LieRepresentation(4,4,[[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
#L1Repr4D = LieRepresentation(4,4, [[0,1,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]],[[0,0,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,0,0,0]])
L1Repr3D = LieRepresentation(4,3, [[0,1,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,1],[0,0,0]],[[0,0,1],[0,0,0],[0,0,0]],[[1,0,0],[0,1,0],[0,0,1]])
L2Repr = LieRepresentation(4, 4, [[0,1,0,0],[0,0,0,0],[0,0,0,0],[1,0,0,0]],[[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]],[[0,0,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,1,0]])
L3Repr = L3.AdjointRepresentation()
L4Repr = L4.AdjointRepresentation()
L4inftyRepr3D = LieRepresentation(4,3, [[1,0,0],[0,0,0],[0,0,0]],[[0,1,0],[0,0,0],[0,0,0]],[[1,0,0],[0,1,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,1]])
L5Repr = L5.AdjointRepresentation()
L6Repr3D = LieRepresentation(4, 3, [[0,0,0],[0,1,0],[0,0,-1]],[[0,0,1],[-1,0,0],[0,0,0]],[[0,-1,0],[0,0,0],[1,0,0]],[[1,0,0],[0,1,0],[0,0,1]])
L6Repr4D = LieRepresentation(4,4, [[0,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,0]],[[0,0,1,0],[-1,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,-1,0,0],[0,0,0,0],[1,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
L7Repr = L7.AdjointRepresentation()
L8Repr = L8.AdjointRepresentation()
L9Repr = L9.AdjointRepresentation()

print L0Repr4D.IsRepresentationLie(L0)
#print L1Repr4D.IsRepresentationLie(L1)
print L1Repr3D.IsRepresentationLie(L1)
print L2Repr.IsRepresentationLie(L2)
print L3Repr.IsRepresentationLie(L3)
print L4Repr.IsRepresentationLie(L4)
print L4inftyRepr3D.IsRepresentationLie(L4infty)
print L5Repr.IsRepresentationLie(L5)
print L6Repr3D.IsRepresentationLie(L6)
print L6Repr4D.IsRepresentationLie(L6)
print L7Repr.IsRepresentationLie(L7)
print L8Repr.IsRepresentationLie(L8)

print L9Repr.IsRepresentationLie(L9)

print "L0 4D"
L0Repr4D.PrintBase()
print "L1 3D"
L1Repr3D.PrintBase()
print "L2"
L2Repr.PrintBase()
print "L3"
L3Repr.PrintBase()
print "L4"
L4Repr.PrintBase()
print "L4infty"
L4inftyRepr3D.PrintBase()
print "L5"
L5Repr.PrintBase()
print "L6 3D"
L6Repr3D.PrintBase()
print "L6 4D"
L6Repr4D.PrintBase()
print "L7"
L7Repr.PrintBase()
print "L8"
L8Repr.PrintBase()
print "L9"
L9Repr.PrintBase()
"""
a1 = sympy.Symbol("a1")
a2 = sympy.Symbol("a2")
a3 = sympy.Symbol("a3")
a4 = sympy.Symbol("a4")

print "L0 4D"
print L0Repr4D.RepresentationMatrix([a1,a2,a3,a4]))
print "L1 3D"
print L1Repr3D.RepresentationMatrix([a1,a2,a3,a4]))
print "L2"
print L2Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L3"
print L3Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L4"
print L4Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L4infty"
print L4inftyRepr3D.Representation([a1,a2,a3,a4]))
print "L5"
print L5Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L6 3D"
print L6Repr3D.RepresentationMatrix([a1,a2,a3,a4]))
print "L6 4D"
print L6Repr4D.RepresentationMatrix([a1,a2,a3,a4]))
print "L7"
print L7Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L8"
print L8Repr.RepresentationMatrix([a1,a2,a3,a4]))
print "L9"
print L9Repr.RepresentationMatrix([a1,a2,a3,a4]))
"""
