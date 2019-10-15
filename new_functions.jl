using LinearAlgebra
using Plots


#*******************************************************************************
#                         Sistema com 97.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************
params97 = (
    EnC = 0.0320879,
    EnV = 0.0314029,
    gammaC = 36.917,
    gammaV = -22.4782,
    P = -0.108812,
    delta = -0.00019967,
    eta2 = -0.280174,
    eta3 = -1.16087
);

#*******************************************************************************
#                         Sistema com 103.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************

params103 = (
    EnC = 0.0313913,
    EnV = 0.0314028,
    gammaC = 39.9073,
    gammaV = -21.6738,
    P = 0.0972126,
    delta = -0.000201637,
    eta2 = -0.295556,
    eta3 = -1.17352
);

#*******************************************************************************
#                         Sistema com 110.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************

params110 = (
    EnC = 0.030775,
    EnV = 0.0314028,
    gammaC = 42.9456,
    gammaV = -15.4158,
    P = 0.0865811,
    delta = -0.000203489,
    eta2 = 0.309555,
    eta3 = 1.18696
);

#*******************************************************************************
#                         Hamiltonian Operators
#*******************************************************************************


function Hamil3x3(kx, ky, p)
    H_11 = p.EnC + (kx^2 + ky^2) * p.gammaC;
    H_12 = (kx^2 - ky^2) * p.eta2 + 1im * kx * ky * p.eta3;
    H_13 = - 1im * (kx + 1im * ky) * p.P;

    H_21 = (kx^2 - ky^2) * p.eta2 - 1im * kx * ky * p.eta3;
    H_22 = p.EnV + (kx^2 + ky^2) * p.gammaV;
    H_23 = -1im * p.delta ;

    H_31 = + 1im * (kx - 1im * ky) * p.P;;
    H_32 = 1im * p.delta ;
    H_33 = H_22;
    return [H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33]
end

function Hamil3Nyx3Ny(kx, p; dim = 3)
    k_y_first = (-1im/(2*dy)) * Tridiagonal(-ones(Ny-1), zeros(Ny), ones(Ny-1))
    k_y_second = (-1/(dy^2)) * SymTridiagonal(-2*ones(Ny), ones(Ny-1))

    H_11 = I*p.EnC + (I*kx^2 + k_y_second)*p.gammaC;
    H_12 = (I*kx^2 - k_y_second)*p.eta2 + 1im*kx*k_y_first*p.eta3;
    H_13 = -1im*(kx*I + 1im*k_y_first)*p.P;


    H_21 = transpose(conj(H_12));
    H_22 = I*p.EnV + (I*kx^2 + k_y_second)*p.gammaV;
    H_23 = -1im*p.delta*I;


    H_31 = transpose(conj(H_13));
    H_32 = 1im*p.delta*I;
    H_33 = H_22;

    if dim == 2
        return Matrix([H_11 H_13; H_31 H_33])
    else
        return Matrix([H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33])
    end
end

function spectrum_DS(k_vec, p; spinor_dim = 3)

    #********************************************************************#
    #           Essa função gera uma matriz em que                       #
    #   ▸ cada coluna corrresponde a um valor de kₓ;                     #
    #   ▸ cada linha corresponde uma banda de índice "n";                #
    #                                                                    #
    #                                                                    #
    #                                                                    #
    #       • "k_vec" =  Array com os valores de kₓ                      #
    #                                                                    #
    #       • "H_func" = função que gera a matriz Hamiltoniana           #
    #                                                                    #
    #       • "spinor_dim" = dimensão do spinor que descreve o           #
    #        autoestado da Hamiltoniana                                  #
    #                                                                    #
    #********************************************************************#

    E_vals = zeros(spinor_dim * Ny,length(k_vec));
    cont = 0;
    for k in k_vec
        cont += 1;
        E_vals[:, cont] = eigvals( Hamil3Nyx3Ny(k, p, dim = spinor_dim) );
    end
    return E_vals
end

function convertEnergy(energy)
    return  Ry * 1000 * energy
end

function convertMomentum(momentum)
    return momentum./(a₀*10^(-1))
end

function print_bandas(kx_span, M_vals)

    #******************************************************************#
    #                                                                  #
    #     Essa função imprime na tela a estrutura de bandas em pares   #
    #     • "M_vals" = matriz resultado da função "spectrum()"         #
    #                                                                  #
    #     • "kx_span" = Array de valores de kₓ                         #
    #                                                                  #
    #                                                                  #
    #******************************************************************#
    # print(size(M_vals))
    kx_span = convertMomentum(kx_span)
    M_vals  = convertEnergy(M_vals)
    ymin, ymax = [400 500]
    xmin, xmax = [kx_span[1] kx_span[end]]
    p = plot( kx_span , M_vals',
            xlims = (xmin, xmax),
            ylims = (ymin, ymax),
            color = "black",
            legend=false)
    return p
end


Ry = 13.6; # constante de Rydberg;
a₀ = 0.529167; # raio de Bohr

aSystem = 6.21; # parâmetro de rede em Å
aSystem = aSystem/a₀; # parâmetro de rede em unidades de raio de Bohr

Ny = 300; # número de pontos na dimensão "y"
Ly = 6000.0; # largura do sistema
y = range(-Ly/2, stop=Ly/2, length=Ny);
dy = y[2] - y[1];

porcent = 0.05 ; # fração da Zona de Brillouin (1 = tudo)
Nkx = 500 ; # número de pontos do espaço recíproco
vec_k_limited = range(-1,stop= 1, length=Nkx) * pi/aSystem * porcent;

chosen_dim = 3;

E_097 = spectrum_DS(vec_k_limited, params97, spinor_dim = chosen_dim);

# y_min,y_max = 410,470;
# x_min,x_max = -0.1,0.1;
#
# figure(1)
# plot(vec_k_limited./(a₀*10^(-1)), Ry * 1000 * E_097')
fig = print_bandas(vec_k_limited,  E_097)
# display(display(fig))
savefig(fig,"bands_97_old_params.png")
# ylim(y_min,y_max)
# xlim(x_min,x_max
# length(vec_k_limited)
# size(E_097')
# plot(vec_k_limited./(a₀*10^(-1)), E_097[1,:])
# show()
