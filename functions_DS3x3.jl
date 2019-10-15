
using PyPlot
using LinearAlgebra


#*******************************************************************************
#                         Sistema com 97.0 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************

# enCBR097 = 0.032178964249849396;
# enVBR097 = 0.0318623742113142;
# γCR097 = -4.209592573429029;
# γVR097 = -1.491366115285701;
# PXR097 = 0.10231705910476735;
# η2R097 = -0.5780367713609732;
# η3R097 = -0.6573780493650576;

# Valor usado pelo Guilherme nos resultados obtidos
# e publicados na BWSP:
# Δ097 = - 0.000207954;

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

# enCBR103 = 0.03148921676629649;
# enVBR103 = 0.03186235744169688;
# γCR103 = -4.206355461827499;
# γVR103 = -1.4913678288955927;
# η2R103 = 0.6006715228212958;
# η3R103 = 0.6822720291713092;
# PXR103 = -0.09056275347820107;

# Valor usado pelo Guilherme nos resultados obtidos
# e publicados na BWSP:
# Δ103 = - 0.000209932;

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

# enCBR110 = 0.030881710177341544;
# enVBR110 = 0.03186237510679507;
# γCR110 = -4.201037543186996;
# γVR110 = -1.4913661225738042;
# η2R110 = -0.6233292515281043;
# η3R110 = - 0.7072054946588765;
# PXR110 = 0.07956522136229474;

# Valor usado pelo Guilherme nos resultados obtidos
# e publicados na BWSP:
# Δ110 = - 0.000211285;


EnCB_110 = 0.030775;
EnVB_110 = 0.0314028;
GammaCB_110 = 42.9456;
GammaVB_110 = -15.4158;
P_110 = 0.0865811;
Delta_110 = -0.000203489;
eta2_110 = 0.309555;
eta3_110 = 1.18696;






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
    k_y_first = (-1im/(2*dy)) * full(Tridiagonal(-ones(Ny-1), zeros(Ny), ones(Ny-1)));
    k_y_second = (-1/(dy^2)) * full(SymTridiagonal(-2*ones(Ny), ones(Ny-1)));

    H_11 = eye(Ny)*enCBR + (eye(Ny)*kx^2 + k_y_second)*γCR;
    H_12 = (eye(Ny)*kx^2 - k_y_second)*η2R + 1im*kx*k_y_first*η3R;
    H_13 = -1im*(kx*eye(Ny) + 1im*k_y_first)*PXR;


    H_21 = transpose(conj(H_12));
    H_22 = eye(Ny)*enVBR + (eye(Ny)*kx^2 + k_y_second)*γVR;
    H_23 = -1im*Δ*eye(Ny);


    H_31 = transpose(conj(H_13));
    H_32 = 1im*Δ*eye(Ny);
    H_33 = H_22;

    if dim == 2
        return [H_11 H_13; H_31 H_33]
    else
        return [H_11 H_12 H_13; H_21 H_22 H_23; H_31 H_32 H_33]
    end
end





function Hamil3Nyx3Ny_BHZ(kx, spinor_dim = 3)

    #**********************************************
    #                                             #
    #   Teste com o hamiltoniano BHZ (spin up):   #
    #                                             #
    #**********************************************

    A = 364.5; # meV.nm
    B = -686.0; # meV.nm^2
    C = 0; # meV
    D = -512.0; # meV.nm^2
    M = -10.0; # meV
    hbar = 0.65821188926; # somente ħ em meV.ps


    # matriz simetrica e tridiagonal p/ derivada de ordem 2:
    k_y_second = -1/(dy^2) * full(SymTridiagonal(-2 * ones(Ny), ones(Ny - 1)));
    # matriz tridiagonal p/ derivada de ordem 1
    k_y_first = -1im/(2*dy) * full(Tridiagonal(-1 * ones(Ny - 1), zeros(Ny), ones(Ny - 1)));


    H_11 = eye(Ny)*(C + M) - (D + B)*(eye(Ny)*kx^2 + k_y_second);
    H_12 = A * (eye(Ny) * kx + 1im * k_y_first);
    H_13 = A * (eye(Ny) * kx + 1im * k_y_first); # igual a H_12
    H_21 = A * (eye(Ny) * kx - 1im * k_y_first);
    H_22 = (C - M) * eye(Ny) - (D - B)*(eye(Ny)*kx^2 + k_y_second);
    H_31 = H_21;
    H_33 = H_22;

    if spinor_dim == 2
        return [H_11 H_12; H_21 H_22]
    else
        return [H_11 H_12 0*H_13; H_21 H_22 zeros(Ny,Ny); 0*H_31 zeros(Ny,Ny) 0*H_33]
    end
end





function spectrum_BHZ(k_vec, H_func; spinor_dim = 3)
    E_vals = Array{Float64,2}(spinor_dim * Ny,length(k_vec));
    cont = 0;
    for k in k_vec
        cont += 1;
        E_vals[:, cont] = eigvals(H_func(k, spinor_dim));
    end
    return E_vals
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

    E_vals = Array{Float64,2}(spinor_dim * Ny,length(k_vec));
    cont = 0;
    for k in k_vec
        cont += 1;
        E_vals[:, cont] = eigvals( Hamil3Nyx3Ny(k, enCBR, enVBR, γCR, γVR,
                            η2R, η3R, PXR, Δ, dim = spinor_dim) );
    end
    return E_vals
end





function print_bandas(M_vals , kx_span)

    #******************************************************************#
    #                                                                  #
    #     Essa função imprime na tela a estrutura de bandas em pares   #
    #     • "M_vals" = matriz resultado da função "spectrum()"         #
    #                                                                  #
    #     • "kx_span" = Array de valores de kₓ                         #
    #                                                                  #
    #                                                                  #
    #******************************************************************#

    for i = 1 : size(M_vals)[1]
        plot( kx_span , M_vals[ i, : ], linestyle="",marker=".",markersize=0.8, color = "black");
    end
    return 0
end





function spectrum_inf_DS(k_vec_x, ky, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ )

    H_inf = [Hamil3x3(kx, 0.0, enCBR, enVBR, γCR, γVR, η2R, η3R, PXR, Δ) for kx in k_vec_x];
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

    return 0
end




function auto_vec(;larg="97", d = 2, kx_valor = 0.25e-2, n_banda = Ny)

    # Retorna um vetor de dimensões (2Ny)×1 onde
    #   ▸ a primeira metade é a auto função para o estado de Elétron, e
    #   ▸ a segunda metade é a auto função para o estado de Buraco.

    if larg == "97"
      H_kx = Hamil3Nyx3Ny(kx_valor, enCBR097 , enVBR097 ,
                  γCR097, γVR097,
                  η2R097, η3R097, PXR097, dim = d);
    elseif larg == "103"
      H_kx = Hamil3Nyx3Ny(kx_valor, enCBR103 , enVBR103 ,
                  γCR103, γVR103,
                  η2R103, η3R103, PXR103, dim = d);
    else
      H_kx = Hamil3Nyx3Ny(kx_valor, enCBR110 , enVBR110 ,
                  γCR110, γVR110,
                  η2R110, η3R110, PXR110, dim = d);
    end



    m_vecs = eigvecs(H_kx);
    m_vals = Ry * 1000 * eigvals(H_kx);
    println(m_vals[n_banda])
    return m_vecs[:,n_banda]
end
