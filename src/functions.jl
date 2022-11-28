import Base.show
import Base.copy
using Printf, LinearAlgebra

__show_mode__ = :compact

struct Tableau{T <: Real}
    table::Matrix{T}
    basis::Vector{Int}
    varnames::Vector{String}
end

################
function Tableau(table::Matrix{T}, basis::Vector{Int}) where {T <: Real}
    n = size(table, 2) - 1
    varnames = ["x$i" for i = 1:n]
    return Tableau(table, basis, varnames)
end

################
function pivot!(A::Matrix{T}, p::Int, q::Int) where {T}

    m, n = size(A)

    A[p,:] /= A[p,q]

    for i = 1:m
        if i != p
            tmp = A[i,q]
            for j = 1:n
                A[i,j] -= tmp*A[p,j]
            end
        end
    end

    return A
end

################
function pivot!(Tab::Tableau{T}, entervar::Int, leavevar::Int) where {T}

    var = findfirst(isequal(leavevar), Tab.basis)
    isnothing(var) && error("x$leave is not a basic variable.")

    pivot!(Tab.table, var+1, entervar)
    Tab.basis[var] = entervar

    return Tab
end

################
function tableau(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}

    m, n = size(A)

    table = [-c' zeros(T, 1, m) zero(T); A I b]
    basis = collect(n+1:n+m)

    return Tableau(table, basis)
end

################
function bfs(Tab::Tableau{T}) where {T}

    x = zeros(T, nvars(Tab))
    x[Tab.basis] = rhs(Tab)

    return x
end

################
rhs(Tab::Tableau{T}) where {T} = Tab.table[2:end,end]
nvars(Tab::Tableau{T}) where {T} = size(Tab.table, 2) - 1
reducedcosts(Tab::Tableau{T}) where {T} = Tab.table[1,1:nvars(Tab)]
objval(Tab::Tableau{T}) where {T} = Tab.table[1,end]

################
function clean!(Tab::Tableau{T}) where {T}

    for var in Tab.basis
        pivot!(Tab, var, var)
    end

    return Tab
end

################
function simplex!(Tab::Tableau{T}; verbose=false) where {T}

    verbose && println(Tab)
    status = "Unsolved"

    while true
        # Optimality test
        cbar = reducedcosts(Tab)
        if all(cbar .>= 0)
            status = "Optimal"
            break
        end

        # Find entering variable using Dantzig's Rule
        _, entervar = findmin(cbar)

        # Unboundedness test
        entercol = Tab.table[2:end,entervar]
        if all(entercol .<= 0)
            status = "Unbounded"
            break
        end

        # Perform a simplex pivot
        simplexpivot!(Tab, entervar, verbose=verbose)
        verbose && println(Tab)
    end

    verbose && println("\n", status, "\n")

    return Tab
end

################
function simplexpivot!(
        Tab::Tableau{T}, entervar::String; verbose=false) where {T}

    ind = findfirst(isequal(entervar), Tab.varnames)

    return simplexpivot!(Tab, ind, verbose=verbose)
end

################
function simplexpivot!(Tab::Tableau{T}, entervar::Int; verbose=false) where {T}

    # Find leaving variable
    entercol = Tab.table[2:end,entervar]
    rows = entercol .> 0
    ratios = rhs(Tab)[rows]./entercol[rows]
    _, minratio = findmin(ratios)
    leavevar = Tab.basis[rows][minratio]

    # Perform pivot
    pivot!(Tab, entervar, leavevar)

    entername, leavename = Tab.varnames[entervar], Tab.varnames[leavevar]
    verbose && @printf("\n%s enters, %s leaves\n", entername, leavename)
    return Tab
end

################
function dualsimplex!(Tab::Tableau{T}; verbose=false) where {T}

    verbose && println(Tab)
    status = "Unsolved"

    while true
        # Optimality test
        bbar = rhs(Tab)
        if all(bbar .>= 0)
            status = "Optimal"
            break
        end

        # Find leaving variable
        _, leaveind = findmin(bbar)
        leavevar = Tab.basis[leaveind]

        # Infeasibility test
        leaverow = Tab.table[leaveind+1,1:end-1]
        if all(leaverow .>= 0)
            status = "Infeasible"
            break
        end

        # Perform a dual simplex pivot
        dualsimplexpivot!(Tab, leavevar, verbose=verbose)
        verbose && println(Tab)
    end

    verbose && println("\n", status, "\n")

    return Tab
end

################
function dualsimplexpivot!(
        Tab::Tableau{T}, leavevar::String; verbose=false) where {T}

    ind = findfirst(isequal(leavevar), Tab.varnames)

    return dualsimplexpivot!(Tab, ind, verbose=verbose)
end

################
function dualsimplexpivot!(
        Tab::Tableau{T}, leavevar::Int; verbose=false) where {T}

    # Find entering variable
    ind = findfirst(isequal(leavevar), Tab.basis)
    leaverow = Tab.table[ind+1,1:end-1]
    cols = findall(leaverow .< 0)
    ratios = reducedcosts(Tab)[cols]./leaverow[cols]
    _, maxratio = findmax(ratios)
    entervar = cols[maxratio]

    # Perform pivot
    pivot!(Tab, entervar, leavevar)

    verbose && @printf("\nx%d enters, x%d leaves\n", entervar, leavevar)
    return Tab
end

################
function addGomorycuts(
        Tab::Tableau{T}, vars::AbstractVector{Int}) where {T}

    inds = [findfirst(isequal(var), Tab.basis)
            for var in vars if var in Tab.basis]

    ncuts = length(inds)
    rows = Tab.table[inds.+1,:]

    # Generate the Gomory cuts
    newrows = floor.(rows) - rows

    # Allocate memory for an array with one more row and column
    table = zeros(T, size(Tab.table).+ncuts)

    # Fill in the entries of the new table
    table[1:end-ncuts,1:end-ncuts-1] = Tab.table[:,1:end-1]
    table[1:end-ncuts,end] = Tab.table[:,end]

    table[end-ncuts+1:end,1:end-ncuts-1] = newrows[:,1:end-1]
    table[end-ncuts+1:end,end-ncuts:end-1] = diagm(ones(T, ncuts))
    table[end-ncuts+1:end,end] = newrows[:,end]

    # Add the slack variable for the cut to the basis
    basis = [Tab.basis; nvars(Tab)+1:nvars(Tab)+ncuts]

    # Add varname for slack variable
    varnames = [Tab.varnames; ["s$i" for i=nvars(Tab)+1:nvars(Tab)+ncuts]]

    return Tableau(table, basis, varnames)
end

################
function addGomorycuts(Tab::Tableau{T}) where {T}

    inds = findall(.!isinteger.(rhs(Tab)))
    fracvars = Tab.basis[inds]

    return addGomorycuts(Tab, fracvars)
end

################
addGomorycut(Tab::Tableau{T}, var::Int) where {T} = addGomorycuts(Tab, [var])

################
copy(Tab::Tableau{T}) where {T} =
    Tableau(copy(Tab.table), copy(Tab.basis), copy(Tab.varnames))

################
simplex(Tab::Tableau{T}; kwargs...) where {T} = simplex!(copy(Tab); kwargs...)

################
simplexpivot(Tab::Tableau{T}, entervar::String; kwargs...) where {T} =
    simplexpivot!(copy(Tab), entervar; kwargs...)

################
simplexpivot(Tab::Tableau{T}, entervar::Int; kwargs...) where {T} =
    simplexpivot!(copy(Tab), entervar; kwargs...)

################
dualsimplex(Tab::Tableau{T}; kwargs...) where {T} =
    dualsimplex!(copy(Tab); kwargs...)

################
dualsimplexpivot(Tab::Tableau{T}, leavevar::String; kwargs...) where {T} =
    dualsimplexpivot!(copy(Tab), leavevar; kwargs...)

################
dualsimplexpivot(Tab::Tableau{T}, leavevar::Int; kwargs...) where {T} =
    dualsimplexpivot!(copy(Tab), leavevar; kwargs...)

################
function solveip(Tab::Tableau{T}; verbose=false) where {T}

    # Solve the linear relaxation
    myTab = simplex(Tab, verbose=verbose)

    while true
        # Find noninteger entries in the rhs
        inds = findall(.!isinteger.(rhs(myTab)))
        # Test for optimality of Integer Program (IP)
        isempty(inds) && break
        # Find all fractional basic variables
        fracvars = myTab.basis[inds]
        # Randomly select one of the fractional basic variables
        var = rand(fracvars, 1)[1]
        # Add the Gomory cut for the randomly selected variable
        myTab = addGomorycut(myTab, var)
        # Run the dual simplex method to return to optimality
        dualsimplex!(myTab, verbose=verbose)
    end

    return myTab
end

################
function solveip(
        A::Matrix{T}, b::Vector{T}, c::Vector{T}; verbose=false) where {T}

    all(isinteger.(A)) && all(isinteger.(b)) ||
        error("A and b must be integer")

    return solveip(tableau(A, b, c), verbose=verbose)
end


################
function num2string(num::T) where {T <: Real}

    # Pretty print integers and fractions
    if num == 0
        ns = ""
    elseif T == Rational{Int}
        ns = "$(numerator(num))"
        if denominator(num) != 1
            ns = "$ns/$(denominator(num))"
        end
    else
        ns = "$num"
    end

    return ns
end

################
function show(io::IO, Tab::Tableau{T}) where {T <: Real}

    m, n = size(Tab.table)

    # Label the columns with the non-basic variables
    vars = repeat(" ", 5)
    for j = 1:n-1
        if j ∉ Tab.basis || __show_mode__ == :full
            bvar = @sprintf(" %6s", Tab.varnames[j])
        else
            bvar = ""
        end
        vars = "$vars$bvar"
    end

    s = @sprintf("%-4s|", "Z")
    line = ""
    for i = 1:m

        # Print a horizontal line after the objective row
        if i == 2
            line = repeat("-", length(s) - 1)
            s = "$s$line\n"
        end

        # Print the basic variable
        if i != 1
            bvar = Tab.varnames[Tab.basis[i-1]]
            bvar = @sprintf("%-4s", bvar)
            s = "$s$bvar|"
        end

        # Print the array
        for j = 1:n
            tmp = Tab.table[i,j]

            # Before the last column, print a divider |
            if j == n
                s = "$s |"
            end

            # Only print the non-basic columns
            if j ∉ Tab.basis || __show_mode__ == :full
                ns = num2string(tmp)
                ns = @sprintf(" %6s", ns)
            else
                ns = ""
            end
            s = "$s$ns"
        end

        # Print a newline character
        s = "$s\n"
    end

    # Remove the last newline
    s = s[1:end-1]

    # Add vars to the first line
    s = "\n$vars\n$line\n$s"

    print(io, s)
end
