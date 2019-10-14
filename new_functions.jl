using LinearAlgebra
using Plots


#*******************************************************************************
#                         Sistema com 97.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************

EnCB_97 = 0.0320879;
EnVB_97 = 0.0314029;
GammaCB_97 = 36.917;
GammaVB_97 = -22.4782;
P_97 = -0.108812;
Delta_97 = -0.00019967;
eta2_97 = -0.280174;
eta3_97 = -1.16087;


#*******************************************************************************
#                         Sistema com 103.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************


EnCB_103 = 0.0313913;
EnVB_103 = 0.0314028;
GammaCB_103 = 39.9073;
GammaVB_103 = -21.6738;
P_103 = 0.0972126;
Delta_103 = -0.000201637;
eta2_103 = -0.295556;
eta3_103 = -1.17352;


#*******************************************************************************
#                         Sistema com 110.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************

EnCB_110 = 0.030775;
EnVB_110 = 0.0314028;
GammaCB_110 = 42.9456;
GammaVB_110 = -15.4158;
P_110 = 0.0865811;
Delta_110 = -0.000203489;
eta2_110 = 0.309555;
eta3_110 = 1.18696;


#*******************************************************************************
#                         Hamiltonian Operators
#*******************************************************************************


function Hamil3x3(kx, ky, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ)
    H_11 = enCBR + (kx^2 + ky^2) * γCR;
    H_12 = (kx^2 - ky^2) * η2R + 1im * kx * ky * η3R;
    H_13 = - 1im * (kx + 1im * ky) * PXR;

    H_21 = (kx^2 - ky^2) * η2R - 1im * kx * ky * η3R;
    H_22 = enVBR + (kx^2 + ky^2) * γVR;
    H_23 = -1im * Δ ;

    H_31 = + 1im * (kx - 1im * ky) * PXR;;
    H_32 = 1im * Δ ;
    H_33 = H_22;
    return [H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33]
end

function Hamil3Nyx3Ny(kx, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ ; dim = 3)
    k_y_first = (-1im/(2*dy)) * Tridiagonal(-ones(Ny-1), zeros(Ny), ones(Ny-1))
    k_y_second = (-1/(dy^2)) * SymTridiagonal(-2*ones(Ny), ones(Ny-1))

    H_11 = I*enCBR + (I*kx^2 + k_y_second)*γCR;
    H_12 = (I*kx^2 - k_y_second)*η2R + 1im*kx*k_y_first*η3R;
    H_13 = -1im*(kx*I + 1im*k_y_first)*PXR;


    H_21 = transpose(conj(H_12));
    H_22 = I*enVBR + (I*kx^2 + k_y_second)*γVR;
    H_23 = -1im*Δ*I;


    H_31 = transpose(conj(H_13));
    H_32 = 1im*Δ*I;
    H_33 = H_22;

    if dim == 2
        return Matrix([H_11 H_13; H_31 H_33])
    else
        return Matrix([H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33])
    end
end

function spectrum_DS(k_vec, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ; spinor_dim = 3)

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
        E_vals[:, cont] = eigvals( Hamil3Nyx3Ny(k, enCBR, enVBR, γCR, γVR,
                            η2R, η3R, PXR, Δ, dim = spinor_dim) );
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
    display(p)
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

E_097 = spectrum_DS(vec_k_limited,
                EnCB_97, EnVB_97,
                GammaCB_97, GammaVB_97,
                eta2_97, eta3_97,
                P_97, Delta_97,
                spinor_dim = chosen_dim);

# y_min,y_max = 410,470;
# x_min,x_max = -0.1,0.1;
#
# figure(1)
# plot(vec_k_limited./(a₀*10^(-1)), Ry * 1000 * E_097')
print_bandas(vec_k_limited,  E_097)
# ylim(y_min,y_max)
# xlim(x_min,x_max
# length(vec_k_limited)
# size(E_097')
# plot(vec_k_limited./(a₀*10^(-1)), E_097[1,:])
# show()
