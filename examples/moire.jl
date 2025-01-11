using MonteCarlo
using MonteCarloModels

tm = MonteCarlo.hopping_matrix(BTHubbardModel(3,3; tre=-4.348374610127436, tim=-3.1842927600873527*im, boundary=:open))