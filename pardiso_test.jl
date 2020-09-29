
#ENV["OMP_NUM_THREADS"] = 2

using SparseArrays, Pardiso

A = sparse(rand(10,10))
B = rand(10,2)
X = zeros(10,2)

ps = PardisoSolver()

key = 2 # real symmetric positive definite
#  set_matrixtype!(ps, key)

print(get_nprocs(ps))

pardiso(ps, X, A, B)

open("test.out", "w") do t
    write(X)
end
