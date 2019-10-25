include("HamGaSb.jl")

import .HamGaSb

kx = range(-pi, stop=pi, length=201)./HamGaSb.aSystem * 0.05

HamGaSb.spectrum_inf_DS(kx, 0., 60., HamGaSb.params_97)
