using HDF5;
using PyPlot;

PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
PyPlot.matplotlib[:rc]("font", family="serif")

include("functions_DS3x3.jl")
Ry = 13.6; # constante de Rydberg;

#*******************************************************************************


# todos_arquivos_da_pasta = readdir();
# print(todos_arquivos_da_pasta);
# hdf5_files = [];
#
# for i = 1:length(todos_arquivos_da_pasta)
#     if endswith(todos_arquivos_da_pasta[i], ".h5")
#         push!(hdf5_files, todos_arquivos_da_pasta[i])
#     end
# end
#
# println("hdf5_files = ")
# println(hdf5_files)


function read_and_plot_three_results(nome_arquivo, nome_grupo)

    cd("/home/marcos/Dropbox/projetos/sipahi_dias/")

    # nome_arquivo = "espectro_sist_confinado_marcos.h5";
    # nome_grupo = "1000.0_0.01_701";



    primeira_camada = h5open(nome_arquivo,"r") do file
        # "primeira_camada" = array de strings onde
        # cada entrada será o nome de um grupo pertencente
        # ao arquivo nomeado por "nome_arquivo".
        names(file);
    end


    segunda_camada = h5open(nome_arquivo,"r") do file
        # "segunda_camada" = array de strings onde
        # cada entrada será o nome de conj. de dados
        # pertencente ao grupo denotado por "nome_grupo".
        names(file[nome_grupo])
    end

    #=

        Para executar o "plot", precisamos saber a qual
        array de "k's" o conjunto de dados se refere. Essa
        informação é dada via "atributos" do "grupo" que
        contém o conjunto de dados.

        Porém, como os atributos são sempre dados por "strings"
        é mais fácil salvar o "vec_k_limited" e lê-lo inteiramente
        quando necessáio. Construir a  array de "k's" fica  complicada
        quando as informações para tal não estão num fomato
        numérico.

    =#


    vec_k_limited = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, "kx");
    end

    E_097 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, segunda_camada[1]);
    end

    E_103 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, segunda_camada[2]);
    end

    E_110 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, segunda_camada[3]);
    end


    y_min,y_max = 415,440;
    x_min,x_max = -0.01,0.01;

    figure(1)
    print_bandas(Ry * 1000 * E_097 , vec_k_limited)
    ylim(y_min,y_max)
    xlim(x_min,x_max)

    figure(2)
    print_bandas(Ry * 1000 * E_103 , vec_k_limited)
    ylim(y_min,y_max)
    xlim(x_min,x_max)

    figure(3)
    print_bandas(Ry * 1000 * E_110 , vec_k_limited)
    ylim(y_min,y_max)
    xlim(x_min,x_max)

    return 0
end


function read_and_plot_three_results_path(path, nome_arquivo, nome_grupo, nomes_dados)

    cd(path)

    # nome_arquivo = "espectro_sist_confinado_marcos.h5";
    # nome_grupo = "1000.0_0.01_701";



    primeira_camada = h5open(nome_arquivo,"r") do file
        # "primeira_camada" = array de strings onde
        # cada entrada será o nome de um grupo pertencente
        # ao arquivo nomeado por "nome_arquivo".
        names(file);
    end


    segunda_camada = h5open(nome_arquivo,"r") do file
        # "segunda_camada" = array de strings onde
        # cada entrada será o nome de conj. de dados
        # pertencente ao grupo denotado por "nome_grupo".
        names(file[nome_grupo])
    end

    #=

        Para executar o "plot", precisamos saber a qual
        array de "k's" o conjunto de dados se refere. Essa
        informação é dada via "atributos" do "grupo" que
        contém o conjunto de dados.

        Porém, como os atributos são sempre dados por "strings"
        é mais fácil salvar o "vec_k_limited" e lê-lo inteiramente
        quando necessáio. Construir a  array de "k's" fica  complicada
        quando as informações para tal não estão num fomato
        numérico.

    =#


    vec_k_limited = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, "kx");
    end

    E_097 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        # println(string("figure 1 :",segunda_camada[1]))
        read(group, nomes_dados[1]);
    end

    E_103 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        # println(string("figure 2 :",segunda_camada[2]))
        read(group, nomes_dados[2]);
    end

    E_110 = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        # println(string("figure 3 :",segunda_camada[3]))
        read(group, nomes_dados[3]);
    end



    y_min,y_max = 415,440;
    x_min,x_max = -0.01,0.01;

    figure(1)
    print_bandas(Ry * 1000 * E_097 , vec_k_limited./(2*pi/aSystem) )
    edit_fig( lim_x=[-0.01, 0.01], lim_y=[415+10,440+10],
                nome=L"$97$ \AA", label_x= L"$k_x~[2\pi/a_0]$",
                label_y=L"$\varepsilon$ $[meV]$")
    # ylim(y_min,y_max+20)
    # xlim(x_min,x_max)

    figure(2)
    print_bandas(Ry * 1000 * E_103 , vec_k_limited./(2*pi/aSystem) )
    edit_fig( lim_x=[-0.01, 0.01], lim_y=[415,440],
                nome=L"$103$ \AA", label_x=L"$k_x~[2\pi/a_0]$",
                label_y=L"$\varepsilon$ $[meV]$")
    # ylim(y_min,y_max)
    # xlim(x_min,x_max)

    figure(3)
    print_bandas(Ry * 1000 * E_110 , vec_k_limited./(2*pi/aSystem) )
    edit_fig( lim_x=[-0.01, 0.01], lim_y=[415,440],
                nome=L"$110$ \AA", label_x=L"$k_x~[2\pi/a_0]$",
                label_y=L"$\varepsilon$ $[meV]$")
    # ylim(y_min,y_max)
    # xlim(x_min,x_max)

    return 0
end


function read_by_name(path, nome_arquivo, nome_grupo, data_set)

    cd(path)
    println( string("\n Diretório atual: ", path,"\n"))

    if !(nome_arquivo in readdir())
        return println("\nArquivo hdf não encontrado!\n")
    end

    primeira_camada = h5open(nome_arquivo,"r") do file
        # "primeira_camada" = array de strings onde
        # cada entrada será o nome de um grupo pertencente
        # ao arquivo nomeado por "nome_arquivo".
        names(file);
    end

    if !(nome_grupo in primeira_camada)
        return "\nGrupo não encontrado.\n"
    end

    segunda_camada = h5open(nome_arquivo,"r") do file
        # "segunda_camada" = array de strings onde
        # cada entrada será o nome de conj. de dados
        # pertencente ao grupo denotado por "nome_grupo".
        names(file[nome_grupo])
    end

    if !(data_set in segunda_camada)
        return "\nData set não encontrado.\n"
    end


    #=

        Para executar o "plot", precisamos saber a qual
        array de "k's" o conjunto de dados se refere. Essa
        informação é dada via "atributos" do "grupo" que
        contém o conjunto de dados.

        Porém, como os atributos são sempre dados por "strings"
        é mais fácil salvar o "vec_k_limited" e lê-lo inteiramente
        quando necessáio. Construir a  array de "k's" fica  complicada
        quando as informações para tal não estão num fomato
        numérico.

    =#


    vec_k_limited = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, "kx");
    end

    Espectro_lido = h5open(nome_arquivo,"r") do file
        group = file[nome_grupo];
        read(group, data_set);
    end


    y_min,y_max = 415,440;
    x_min,x_max = -0.01,0.01;

    figure(1)
    print_bandas(Ry * 1000 * Espectro_lido , vec_k_limited)
    ylim(y_min,y_max)
    xlim(x_min,x_max)

    return 0
end


function edit_fig(; lim_x=[-0.2,0.2], lim_y=[-60,60],
        nome="", label_x="", label_y="")
    ax = gca();
    axis("tight");
    title(nome, fontsize=20);
    xlabel(label_x,fontsize=20);
    ylabel(label_y,fontsize=20);
    setp(ax[:get_xticklabels](), fontsize=15, color="black", family="serif"); # X Axis font formatting
    ax[:ticklabel_format](style="sci", axis="x", scilimits=(0,0))
    setp(ax[:get_yticklabels](), fontsize=15, color="black", family="serif"); # Y Axis font formatting
    xlim(lim_x);
    ylim(lim_y);
    tight_layout()
    return 0
end



# println("Grupos pertencentes ao arquivo especificado: ")
# print(primeira_camada)
#
#
# println("Dados pertencentes ao grupo especificado: ")
# print(primeira_camada)


#*******************************************************************************




#*******************************************************************************


#*******************************************************************************
