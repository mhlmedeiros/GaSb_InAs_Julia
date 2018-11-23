include("functions_DS3x3.jl")

PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
PyPlot.matplotlib[:rc]("font", family="serif")


function edit_fig(; nome="", label_x="", label_y="")
    ax = gca();
    axis("tight");
    title(nome, fontsize=20);
    xlabel(label_x,fontsize=20);
    ylabel(label_y,fontsize=20);
    ylim(y_min,y_max)
    xlim(x_min,x_max)
    setp(ax[:get_xticklabels](), fontsize=15, color="black", family="serif"); # X Axis font formatting
    setp(ax[:get_yticklabels](), fontsize=15, color="black", family="serif"); # Y Axis font formatting
    tight_layout()
    show()
    return 0
end





Ry = 13.6; # constante de Rydberg
a₀ = 0.529167; # raio de Bohr

aSystem = 6.21;
aSystem = aSystem/a₀;

Ny = 201;
Ly = 200;
y = linspace(-Ly/2, Ly/2, Ny);
dy = y[2] - y[1];

porcent = 0.15 ; # fração da Zona de Brillouin (1 = tudo)
Nkx = 1001 ; # número de pontos do espaço recíproco
vec_k_limited = linspace(-1, 1, Nkx) * pi/aSystem * porcent;

x_min, x_max  = vec_k_limited[1], vec_k_limited[end];
y_min, y_max = 400, 450;


# figure(1)
# spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
#             EnCB_97, EnVB_97,
#             GammaCB_97, GammaVB_97,
#             eta2_97, eta3_97,
#             P_97, -2*Delta_97)
# edit_fig( nome=L"$97$ \AA", label_x=L"$k_x~[2\pi/a_0]$",
#             label_y=L"$\varepsilon$ $[meV]$")

figure(1)
spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
            EnCB_103, EnVB_103,
            GammaCB_103, GammaVB_103,
            eta2_103, eta3_103,
            P_103, +1.0*Delta_103)
edit_fig( nome= string(L"$\Delta = $",Delta_103," meV"), label_x=L"$k_x~[2\pi/a_0]$",
            label_y=L"$\varepsilon$ $[meV]$")

# figure(2)
# spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
#             EnCB_103, EnVB_103,
#             GammaCB_103, GammaVB_103,
#             eta2_103, eta3_103,
#             P_103, -2.0*Delta_103)
# edit_fig( nome= string(L"$\Delta = $",-2*Delta_103," meV"), label_x=L"$k_x~[2\pi/a_0]$",
#             label_y=L"$\varepsilon$ $[meV]$")
# #
# figure(3)
# spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
#             EnCB_103, EnVB_103,
#             GammaCB_103, GammaVB_103,
#             eta2_103, eta3_103,
#             P_103, -4.0*Delta_103)
# edit_fig( nome= string(L"$\Delta = $",-4*Delta_103," meV"), label_x=L"$k_x~[2\pi/a_0]$",
#             label_y=L"$\varepsilon$ $[meV]$")
# #
# figure(4)
# spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
#             EnCB_103, EnVB_103,
#             GammaCB_103, GammaVB_103,
#             eta2_103, eta3_103,
#             P_103, -8.0*Delta_103)
# edit_fig( nome= string(L"$\Delta = $",-8*Delta_103," meV"), label_x=L"$k_x~[2\pi/a_0]$",
#             label_y=L"$\varepsilon$ $[meV]$")

# figure(3)
# spectrum_inf_DS(vec_k_limited, vec_k_limited[2]-vec_k_limited[1],
#             EnCB_110, EnVB_110,
#             GammaCB_110, GammaVB_110,
#             eta2_110, eta3_110,
#             P_110, -2*Delta_110)
# edit_fig( nome=L"$110$ \AA", label_x=L"$k_x~[2\pi/a_0]$",
#             label_y=L"$\varepsilon$ $[meV]$")
