"""L2 norm squared"""
norm2(x) = (x' * x)[1]

"""The objective function value"""
flsa(x, y, D, λ) = 0.5 * norm2(y-x) + λ*norm(D*x, 1)
