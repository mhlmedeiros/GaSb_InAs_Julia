
module HamGaSb

using PyPlot
using LinearAlgebra


## Constants :
const Ry = 13.6; # constante de Rydberg;
const a₀ = 0.529167; # raio de Bohr
const aSystem = 6.21/a₀; # parâmetro de rede em unidades de raio de Bohr

## Modules that are visible to "main"
export spectrumBulk
export hamilBulk
export params_97
export params_103
export params_110

function γC(gamma :: Float64)
    ## for 97 Å
    gamma
end

function γC(x::Float64, A::Float64, B::Float64)
    ## for 103 Å and 110 Å
    A + B*x^2
end

function γV(gamma::Float64)
    gamma ## for 97 Å
end

function γV(x::Float64, A::Float64, B::Float64, C::Float64)
    ## for 103 Å and 110 Å
    A + B*x + C*x^2
end

function Δγ(x::Float64, A::Float64, B::Float64, C::Float64)
    ## 97 Å
    A + B*exp(C*x)^2
end

function Δγ(x::Float64, A::Float64, B::Float64, C::Float64,
    D::Float64, F::Float64, G::Float64)
    ## 103 Å and 110 Å
    x < 6 ? A + B * exp(C*x) : D + F*x + G*x^2
end

function αC(x::Float64, A::Float64, B::Float64)
    # for 97 Å
    A + B * x
end

function αC(x::Float64, A::Float64, B::Float64, C::Float64)
    ## 103 and 110 Å
    A + B*x + C*x^2
end

function αV(x::Float64, A::Float64, B::Float64, C::Float64)
    ## for all systems
    A + B*exp(C*x)
end

function ΔR(x::Float64, A::Float64)
    ## for all systems
    A*x
end

function Px(P::Float64)
    P ## for 97 and 110 Å
end

function Px(x::Float64, A::Float64, B::Float64, C::Float64)
    A + B*x + C*x^2 ## for 103 Å
end


params_97 = Dict(
    :EC => 0.0320879,
    :EV => 0.0314029,
    :γCargs => [36.917],
    :γVargs => [-22.4782],
    :Pxargs => [-0.108812],
    :Δγargs => [-0.231075, 4.11407, -0.255277],
    :ΔRargs => [2.79538e-5],
    :αCargs => [3.24001e-5, 8.28089e-5],
    :αVargs => [0.0219519, -0.02344, -0.286293])

params_103 = Dict(
    :EC => 0.0313913,
    :EV => 0.0314028,
    :γCargs => [42.6859, -0.00694795],
    :γVargs => [-19.1031, -0.0203693, 0.00417352],
    :Pxargs => [0.0971718, 1.7474e-7, -5.92765e-7],
    :Δγargs => [1.31607,-3.7044,-0.140137,-0.158825,0.0213621,-0.0031922],
    :ΔRargs => [-2.72091e-5],
    :αCargs => [-5.41767e-5, 2.48155e-4, -3.66946e-7],
    :αVargs => [0.0164299,-0.0166324,-0.43601])

params_110 = Dict(
    :EC => 0.030775,
    :EV => 0.0314028,
    :γCargs => [39.3037, -8.92854e-3,],
    :γVargs => [-16.9182, 0.0434725, 2.84609e-3],
    :Pxargs => [0.125543],
    :Δγargs => [-0.0609416, -0.560249, -0.536596, 0.494559, -0.156445, -2.21618e-3],
    :ΔRargs => [-2.77126e-5],
    :αCargs => [6.47954e-5, 4.07089e-4, -4.18461e-6],
    :αVargs => [0.0206344, -0.0210344, -0.270751])



function preparargs(γCargs, γVargs, Pxargs, ΔRargs, Δγargs, αCargs, αVargs, εF)

    if length(γCargs) != 1
        pushfirst!(γCargs, εF)
        pushfirst!(γVargs, εF)
    end

    if length(Pxargs) != 1
        pushfirst!(Pxargs, εF)
    end

    pushfirst!(ΔRargs, εF)
    pushfirst!(Δγargs, εF)
    pushfirst!(αCargs, εF)
    pushfirst!(αVargs, εF)
end


function showparams(EC, EV, γCargs, γVargs,
                    Pxargs, Δγargs, ΔRargs, αCargs, αVargs)
    println("gammaC: ", γCargs)
    println("gammaV: ", γVargs)
    println("DeltaR: ", ΔRargs)
    println("alphaC: ", αCargs)
    println("alphaV: ", αVargs)
end

function hamilBulk(kx, ky, εF; EC, EV, γCargs, γVargs, Pxargs,
                    Δγargs, ΔRargs, αCargs, αVargs)

    if length(ΔRargs) == 1
        # show_params(EC, EV, γCargs, γVargs, Pxargs, Δγargs, ΔRargs, αCargs, αVargs)
        preparargs(γCargs, γVargs, Pxargs, ΔRargs, Δγargs, αCargs, αVargs, εF)
        # show_params(EC, EV, γCargs, γVargs, Pxargs, Δγargs, ΔRargs, αCargs, αVargs)
    end


    H_11 = EC + αC(αCargs...)*(kx + ky) + γC(γCargs...)*(kx^2 + ky^2)
    H_12 =  1im * (kx - 1im * ky) * Px(Pxargs...)
    H_21 = -1im * (kx + 1im * ky) * Px(Pxargs...)

    if length(γVargs) == 1
        H_22 = EV + ΔR(ΔRargs...) + (γV(γVargs...) + Δγ(Δγargs...)) * (kx^2 + ky^2) + αV(αVargs...) * (kx + ky)
        H_23 = 0
        H_32 = 0
        H_33 = EV - ΔR(ΔRargs...) + (γV(γVargs...) - Δγ(Δγargs...)) * (kx^2 + ky^2) + αV(αVargs...) * (kx + ky)
    else
        H_22 = EV + (γV(γVargs...) + Δγ(Δγargs...)) * (kx^2 + ky^2) + αV(αVargs...) * (kx + ky)
        H_23 =  1im*ΔR(ΔRargs...)
        H_32 = -1im*ΔR(ΔRargs...)
        H_33 = EV + (γV(γVargs...) - Δγ(Δγargs...)) * (kx^2 + ky^2) + αV(αVargs...) * (kx + ky)
    end

    [H_11 H_12  0; H_21 H_22 H_23;  0 H_32 H_33]
end


function spectrumBulk(k_vec_x, ky, εF, params)

    H_inf = [hamilBulk(kx, ky, εF; params...) for kx in k_vec_x];
    H_spectrum = [Ry * 1000 * eigvals( H_inf[i] ) for i = 1 : length(k_vec_x) ];

    H_1 = map(n -> n[1], H_spectrum);
    H_2 = map(n -> n[2], H_spectrum);
    H_3 = map(n -> n[3], H_spectrum);

    # plot(k_vec_x./(2*pi/aSystem), H_1,linestyle="",marker=".",markersize=0.5, color="red")
    # plot(k_vec_x./(2*pi/aSystem), H_2,linestyle="",marker=".",markersize=0.5, color="blue")
    # plot(k_vec_x./(2*pi/aSystem), H_3,linestyle="",marker=".",markersize=0.5, color="green")

    plot(k_vec_x./(2*pi/aSystem), H_1,linestyle="-",markersize=0.5, color="red")
    plot(k_vec_x./(2*pi/aSystem), H_2,linestyle="-",markersize=0.5, color="blue")
    plot(k_vec_x./(2*pi/aSystem), H_3,linestyle="-",markersize=0.5, color="green")
    show()
end


function print_parameters()
    println(params_97)

    println(params_103)

    println(params_110)
end

end # end module



# H = hamil3x3_up(0, 0, 0; params_103...)
# println("H :", H)


# function Hamil3Nyx3Ny(kx, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ ; dim = 3)
#     k_y_first = (-1im/(2*dy)) * full(Tridiagonal(-ones(Ny-1), zeros(Ny), ones(Ny-1)));
#     k_y_second = (-1/(dy^2)) * full(SymTridiagonal(-2*ones(Ny), ones(Ny-1)));
#
#     H_11 = eye(Ny)*enCBR + (eye(Ny)*kx^2 + k_y_second)*γCR;
#     H_12 = (eye(Ny)*kx^2 - k_y_second)*η2R + 1im*kx*k_y_first*η3R;
#     H_13 = -1im*(kx*eye(Ny) + 1im*k_y_first)*PXR;
#
#
#     H_21 = transpose(conj(H_12));
#     H_22 = eye(Ny)*enVBR + (eye(Ny)*kx^2 + k_y_second)*γVR;
#     H_23 = -1im*Δ*eye(Ny);
#
#
#     H_31 = transpose(conj(H_13));
#     H_32 = 1im*Δ*eye(Ny);
#     H_33 = H_22;
#
#     if dim == 2
#         return [H_11 H_13; H_31 H_33]
#     else
#         return [H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33]
#     end
# end
#
#
#
#
#
# function Hamil3Nyx3Ny_BHZ(kx, spinor_dim = 3)
#
#     #**********************************************
#     #                                             #
#     #   Teste com o hamiltoniano BHZ (spin up):   #
#     #                                             #
#     #**********************************************
#
#     A = 364.5; # meV.nm
#     B = -686.0; # meV.nm^2
#     C = 0; # meV
#     D = -512.0; # meV.nm^2
#     M = -10.0; # meV
#     hbar = 0.65821188926; # somente ħ em meV.ps
#
#
#     # matriz simetrica e tridiagonal p/ derivada de ordem 2:
#     k_y_second = -1/(dy^2) * full(SymTridiagonal(-2 * ones(Ny), ones(Ny - 1)));
#     # matriz tridiagonal p/ derivada de ordem 1
#     k_y_first = -1im/(2*dy) * full(Tridiagonal(-1 * ones(Ny - 1), zeros(Ny), ones(Ny - 1)));
#
#
#     H_11 = eye(Ny)*(C + M) - (D + B)*(eye(Ny)*kx^2 + k_y_second);
#     H_12 = A * (eye(Ny) * kx + 1im * k_y_first);
#     H_13 = A * (eye(Ny) * kx + 1im * k_y_first); # igual a H_12
#     H_21 = A * (eye(Ny) * kx - 1im * k_y_first);
#     H_22 = (C - M) * eye(Ny) - (D - B)*(eye(Ny)*kx^2 + k_y_second);
#     H_31 = H_21;
#     H_33 = H_22;
#
#     if spinor_dim == 2
#         return [H_11 H_12; H_21 H_22]
#     else
#         return [H_11 H_12 0*H_13; H_21 H_22 zeros(Ny,Ny); 0*H_31 zeros(Ny,Ny) 0*H_33]
#     end
# end
#
#
#
#
#
# function spectrum_BHZ(k_vec, H_func; spinor_dim = 3)
#     E_vals = Array{Float64,2}(spinor_dim * Ny,length(k_vec));
#     cont = 0;
#     for k in k_vec
#         cont += 1;
#         E_vals[:, cont] = eigvals(H_func(k, spinor_dim));
#     end
#     return E_vals
# end
#
#
#
#
#
# function spectrum_DS(k_vec, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ; spinor_dim = 3)
#
#     #********************************************************************#
#     #           Essa função gera uma matriz em que                       #
#     #   ▸ cada coluna corrresponde a um valor de kₓ;                     #
#     #   ▸ cada linha corresponde uma banda de índice "n";                #
#     #                                                                    #
#     #                                                                    #
#     #                                                                    #
#     #       • "k_vec" =  Array com os valores de kₓ                      #
#     #                                                                    #
#     #       • "H_func" = função que gera a matriz Hamiltoniana           #
#     #                                                                    #
#     #       • "spinor_dim" = dimensão do spinor que descreve o           #
#     #        autoestado da Hamiltoniana                                  #
#     #                                                                    #
#     #********************************************************************#
#
#     E_vals = Array{Float64,2}(spinor_dim * Ny,length(k_vec));
#     cont = 0;
#     for k in k_vec
#         cont += 1;
#         E_vals[:, cont] = eigvals( Hamil3Nyx3Ny(k, enCBR, enVBR, γCR, γVR,
#                             η2R, η3R, PXR, Δ, dim = spinor_dim) );
#     end
#     return E_vals
# end
#
#
#
#
#
# function print_bandas(M_vals , kx_span)
#
#     #******************************************************************#
#     #                                                                  #
#     #     Essa função imprime na tela a estrutura de bandas em pares   #
#     #     • "M_vals" = matriz resultado da função "spectrum()"         #
#     #                                                                  #
#     #     • "kx_span" = Array de valores de kₓ                         #
#     #                                                                  #
#     #                                                                  #
#     #******************************************************************#
#
#     for i = 1 : size(M_vals)[1]
#         plot( kx_span , M_vals[ i, : ], linestyle="",marker=".",markersize=0.8, color = "black");
#     end
#     return 0
# end
#
#
#
#
#
#
#
#
#
# function auto_vec(;larg="97", d = 2, kx_valor = 0.25e-2, n_banda = Ny)
#
#     # Retorna um vetor de dimensões (2Ny)×1 onde
#     #   ▸ a primeira metade é a auto função para o estado de Elétron, e
#     #   ▸ a segunda metade é a auto função para o estado de Buraco.
#
#     if larg == "97"
#       H_kx = Hamil3Nyx3Ny(kx_valor, enCBR097 , enVBR097 ,
#                   γCR097, γVR097,
#                   η2R097, η3R097, PXR097, dim = d);
#     elseif larg == "103"
#       H_kx = Hamil3Nyx3Ny(kx_valor, enCBR103 , enVBR103 ,
#                   γCR103, γVR103,
#                   η2R103, η3R103, PXR103, dim = d);
#     else
#       H_kx = Hamil3Nyx3Ny(kx_valor, enCBR110 , enVBR110 ,
#                   γCR110, γVR110,
#                   η2R110, η3R110, PXR110, dim = d);
#     end
#
#
#
#     m_vecs = eigvecs(H_kx);
#     m_vals = Ry * 1000 * eigvals(H_kx);
#     println(m_vals[n_banda])
#     return m_vecs[:,n_banda]
# end
