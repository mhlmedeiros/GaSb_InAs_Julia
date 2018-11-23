include("functions_DS3x3.jl")
include("save_the_data.jl")
include("read_the_data.jl")



Ry = 13.6; # constante de Rydberg;
a₀ = 0.529167; # raio de Bohr

aSystem = 6.21; # parâmetro de rede em Å
aSystem = aSystem/a₀; # parâmetro de rede em unidades de raio de Bohr

Ny = 300; # número de pontos na dimensão "y"
Ly = 6000.0; # largura do sistema
y = linspace(-Ly/2, Ly/2, Ny);
dy = y[2] - y[1];

porcent = 0.05 ; # fração da Zona de Brillouin (1 = tudo)
Nkx = 500 ; # número de pontos do espaço recíproco
vec_k_limited = linspace(-1, 1, Nkx) * pi/aSystem * porcent;

# Change directory to that where the code is
cd("/home/marcos/Dropbox/projetos/sipahi_dias/")









#*******************************************************************************
#                     Information to saving the data
#*******************************************************************************

nome_do_arquivo = "dados_parametros_novos.h5";

# Rank of the Hamiltonian matrix:
chosen_dim = 3;

nome_grupo = string(Ly,"_",porcent,"_",Nkx,"_",Ny,"_dim_",chosen_dim)

group_attr = Dict("Ly" => round(Ly,1),
                "Nkx" => Nkx,
                "Ny" => Ny,
                "BZ_percent" => porcent,
                "aSystem" => aSystem );

# path_dados = "/home/marcos/Documents/julia_resultados/resultados_Guilherme/";
path_dados = "/home/marcos/Documents/julia_resultados/resultados_Guilherme/novos_parametros/";




#*******************************************************************************
#                         Sistema com 97 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************


E_097 = spectrum_DS(vec_k_limited,
                EnCB_97, EnVB_97,
                GammaCB_97, GammaVB_97,
                eta2_97, eta3_97,
                P_97, 10*Delta_97,
                spinor_dim = chosen_dim);

nome_dados_97 = "E_097_Dim_3";

# save_coord(nome_do_arquivo, nome_grupo, E_097, nome_dados, group_attr)
# save_coord_path( path_dados, nome_do_arquivo, nome_grupo, E_097,
#                 nome_dados_97, vec_k_limited, group_attr);



#*******************************************************************************
#                         Sistema com 103 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************


E_103 = spectrum_DS(vec_k_limited,
                EnCB_103, EnVB_103,
                GammaCB_103, GammaVB_103,
                eta2_103, eta3_103,
                P_97, 10*Delta_97,
                spinor_dim = chosen_dim);

nome_dados_103 = "E_103_Dim_3";
# save_coord(nome_do_arquivo, nome_grupo, E_103, nome_dados, group_attr)
# save_coord_path( path_dados, nome_do_arquivo, nome_grupo, E_103,
#                 nome_dados_103, vec_k_limited, group_attr);



#*******************************************************************************
#                         Sistema com 110 Å (Angstrom):
#                       Retirado de Ham3x3-097-103-110.nb
#*******************************************************************************


E_110 = spectrum_DS(vec_k_limited,
                EnCB_110, EnVB_110,
                GammaCB_110, GammaVB_110,
                eta2_110, eta3_110,
                P_110, 10*Delta_110,
                spinor_dim = chosen_dim);

nome_dados_110 = "E_110_Dim_3";
# save_coord(nome_do_arquivo, nome_grupo, E_110, nome_dados, group_attr)
# save_coord_path( path_dados, nome_do_arquivo, nome_grupo, E_110,
#                 nome_dados_110, vec_k_limited, group_attr);


#*******************************************************************************

# h5open(nome_do_arquivo,"r+") do file
#     try
#         g = file[nome_grupo];
#         g["kx"] = collect(vec_k_limited);
#         println("\nvec_k salvo!\n")
#     catch
#         println("\nvec_k, não foi salvo ou não foi sobrescrito.\n")
#     end
# end


# lista_de_nomes = [nome_dados_97, nome_dados_103, nome_dados_110];
#
# read_and_plot_three_results_path(path_dados, nome_do_arquivo, nome_grupo, lista_de_nomes)


cd()
cd("Dropbox/projetos/sipahi_dias/")




y_min,y_max = 415,440;
x_min,x_max = -0.01,0.01;

figure(1)
print_bandas(Ry * 1000 * E_097 , vec_k_limited./(2*pi/aSystem))
ylim(y_min,y_max)
xlim(x_min,x_max)

figure(2)
print_bandas(Ry * 1000 * E_103 , vec_k_limited./(2*pi/aSystem))
ylim(y_min,y_max)
xlim(x_min,x_max)

figure(3)
print_bandas(Ry * 1000 * E_110 , vec_k_limited./(2*pi/aSystem))
ylim(y_min,y_max)
xlim(x_min,x_max)

#*******************************************************************************
#                               Sistema BHZ:
#
#*******************************************************************************

# Ny =100;
# Ly = 200;
# y = linspace(-Ly/2, Ly/2, Ny);
# dy = y[2] - y[1];
# aSystem = 1;
#
# porcent = 0.05 ; # fração da Zona de Brillouin (1 = tudo)
# Nkx = 501 ; # número de pontos do espaço recíproco
# vec_k_limited = linspace(-1, 1, Nkx) * pi/aSystem * porcent
#
# println(Hamil3Nyx3Ny_BHZ(vec_k_limited[1]))
#
# figure(4)
# E_BHZ = spectrum_BHZ(vec_k_limited, Hamil3Nyx3Ny_BHZ ,spinor_dim = 3)
# print_bandas(E_BHZ , vec_k_limited)
# ylim(-50,50)
#
#
# figure(5)
# E_BHZ = spectrum_BHZ(vec_k_limited, Hamil3Nyx3Ny_BHZ ,spinor_dim = 2);
# print_bandas(E_BHZ , vec_k_limited)
# ylim(-50,50)
