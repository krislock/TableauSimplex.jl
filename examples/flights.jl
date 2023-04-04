using TableauSimplex
using LinearAlgebra

table = [
 -140  -100  -140  -100  -140  -100  0  0   1   1   1  0  0  0  -1220
    1     0     1     0     1     0  1  0   0   0   0  0  0  0      6
    0     1     0     1     0     1  0  1   0   0   0  0  0  0      8
  140   100     0     0     0     0  0  0  -1   0   0  1  0  0    300
    0     0   140   100     0     0  0  0   0  -1   0  0  1  0    700
    0     0     0     0   140   100  0  0   0   0  -1  0  0  1    220//1 ]

basis = [7, 8, 12, 13, 14]

vars        = ["A1", "B1", "A2", "B2", "A3", "B3"]
slacks      = ["sA", "sB", "r1", "r2", "r3"]
artificials = ["a1", "a2", "a3"]
varnames    = [vars; slacks; artificials]

T = Tableau(table, basis, varnames)
println(T)

T1 = simplex(T)
println(T1)

table = T1.table[:, [1:11;15]]
table[1,1:6] .= [3000, 2400, 2500, 2000, 2000, 1800]
basis = T1.basis
varnames = [ vars; slacks ]

T2 = Tableau(table, basis, varnames)
clean!(T2)
println(T2)

simplex!(T2)
println(T2)

Topt = solveip(T2, maxit=50)
println(Topt)
