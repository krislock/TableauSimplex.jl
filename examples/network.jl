using TableauSimplex

V = [1, 2, 3, 4, 5, 6]
E = [(1,2), (1,3), (1,5), (2,3), (2,4), (3,4), (3,5), (4,6), (5,6)]

b = [-1, -3, 0, 0, 0, 4]
c = [2, 3, 3, 2, 4, 1, 2, 3, 1]

m, n = length(V), length(E)

# Phase I
E0 = copy(E)
for i in V
    if b[i] < 0
        push!(E0, (i, m+1))
    else
        push!(E0, (m+1, i))
    end
end

A0 = zeros(eltype(b), m+1, length(E0))
for (k, e) in enumerate(E0)
    i, j = e
    A0[i,k] = -1
    A0[j,k] = 1
end

table = Rational{Int}[ zeros(1, n) ones(1, length(E0)-n) 0; 
                       A0[1:end-1,:] b ]
basis = collect(n+1:n+m)
varnames = ["x$i$j" for (i,j) in E]
artificials = ["a$i" for i in V]
T1 = Tableau(table, basis, [varnames; artificials])
println(T1)
println("\nCleaning tableau ...")
clean!(T1)
simplex!(T1, verbose=true)

# Phase II
inds = [i for i=1:m if T1.basis[i] <= n]
table = [
    c' 0
    T1.table[inds.+1, 1:n] T1.table[inds.+1, end] ]
basis = T1.basis[inds]
T2 = Tableau(table, basis, varnames)
println(T2)
println("\nCleaning tableau ...")
clean!(T2)
simplex!(T2, verbose=true)

println("Optimal solution:")
optsoln = [E[T2.basis] rhs(T2)]
