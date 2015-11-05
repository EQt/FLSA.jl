macro inplace(ex)
    if ex.head == :+=
        # broadcast(+, $(ex.args[1]), $(ex.args[2]), $(ex.args[2]))
        :(BLAS.axpy!(1.0, $(ex.args[1]), $(ex.args[2])))
    else
        ex
    end
end
