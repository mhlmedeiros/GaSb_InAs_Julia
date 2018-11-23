using HDF5;



function save_coord_try(file_name, group_name, data, data_name, group_attr)

    #**************************************************************************#
    #   Essa função é responsável por coordenar as demais funções responsáveis                                                                        #
    #  pela gravação dos dados calculados.
    #       • "file_name": é simplesmente o nome do arquivo a ser criado ou
    #         aberto;
    #       • "group_name": nome do grupo ao qual os dados pertencem; como uma
    #          pasta dentro do arquivo "file_name";
    #       • "ψ": dados;
    #       • "data_name": nome da coleção de dados;
    #       • "group_attr": dicionário com atributos do grupo de dados.
    #**************************************************************************#
    global file
    global g

    cd("/home/marcos/Dropbox/projetos/sipahi_dias/")


    try
        # Tentativa de abertura do arquivo:
        file = h5open(file_name,"r+");
        println("\nAbrindo arquivo existente...")
    catch
        # Caso não haja tal arquivo, essa linha o cria:
        file = h5open(file_name,"w");
        println("\nCriando arquivo novo...")
    end


    try
        # Tentativa de abrir grupo:
        g = file[group_name]; # handle
        println("\nAbrindo grupo existente...")
    catch
        # Caso tal grupo seja novo, esse bloco o cria:
        g = g_create(file, group_name);
        # atributos:
        for (key, value) in group_attr
            attrs(g)[key] = string(value);
        end
        println("\nCriando grupo novo...")
    end


    try
        # Tentativa de salvar dados:
        g[data_name] = data ;
        println("\nDados salvos!")

    catch
        # caso já exista tal data_set, simplesmente informa
        # usuário disso.
        println("\nJá existe um arquivo com esse nome.")
    end

    close(file)
    println("\nFim de execução (arquivo fechado)!\n
               **********************************\n\n")

end



function save_coord_path(path, file_name, group_name, data, data_name, kx, group_attr)

    #**************************************************************************#
    #   Essa função é responsável por coordenar as demais funções responsáveis                                                                        #
    #  pela gravação dos dados calculados.
    #       • "file_name": é simplesmente o nome do arquivo a ser criado ou
    #         aberto;
    #       • "group_name": nome do grupo ao qual os dados pertencem; como uma
    #          pasta dentro do arquivo "file_name";
    #       • "ψ": dados;
    #       • "data_name": nome da coleção de dados;
    #       • "group_attr": dicionário com atributos do grupo de dados.
    #**************************************************************************#


    # global file
    # global g


    cd(path)


    if file_name in readdir()
        # Abertura do arquivo:
        file = h5open(file_name,"r+");
        println("\nAbrindo arquivo existente...")
    else
        # Caso não haja tal arquivo, essa linha o cria:
        file = h5open(file_name,"w");
        println("\nCriando arquivo novo...")
    end


    if group_name in names(file)
        g = file[group_name]; # handle
        println("\nAbrindo grupo existente...")
    else
        # Caso tal grupo seja novo, esse bloco o cria:
        g = g_create(file, group_name);
        # atributos:
        for (key, value) in group_attr
            attrs(g)[key] = string(value);
        end
        println("\nCriando grupo novo...")
    end


    if data_name in names(g)
        # caso já exista o data_set, simplesmente informe
        # o usuário disso.
        println("\nJá existe um arquivo com esse nome.")
    else
        # Tentativa de salvar dados:
        g[data_name] = data ;
        println("\nDados salvos!")
    end


    if "kx" in names(g)
        println("\nJá existe a array com 'kx'.\n")
    else
        g["kx"] = collect(kx)
        println("\nArray 'kx' salva!\n")
    end


    close(file)
    println("\nFim de execução (arquivo fechado)!\n
               **********************************\n\n")

end
