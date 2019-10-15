using PyPlot

function plot_duplo()
    x = range(-3pi, stop = 3pi, length = 300)
    y = sin.(x)

    fig, ax = subplots(1, 2, figsize=(10,10))
    ax[1].plot(x, y, linestyle = "--", color="black")  # numbering respests julia indexing
    ax[1].set_xlim(-1.5,1.5)
    ax[2].plot(x, y, linestyle = "-.",label="y1")  # point sintaxe looks like python
    ax[2].set_ylim(-2.5,2.5)
    tight_layout()                      # there is no name space like "plt."
    return fig
end

plot_duplo()
savefig("teste.png")
show()

# k_vec = collect(1:3:10)
# M     = zeros(length(k_vec), length(x))
#
# for k = 1:length(k_vec)
#     M[k,:] = cos.(k_vec[k]*x)
# end
# println("size(M) = ", size(M))
# println("size(M') = ", size(M'))
#
# plot(x, M')
# xlim(-2,2)
# show()
