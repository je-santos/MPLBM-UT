import matplotlib.pyplot as plt


def plot_capillary_pressure_data(Sw, Pc,
                                 Pc_label=r'$P_c$', line_color='k', marker='o', line_style='-',
                                 label_font_size=18, title_font_size=18):
    # f1 = plt.figure()
    plt.plot(Sw, Pc, color=line_color, marker=marker, linestyle=line_style, linewidth=3, markersize=7, label=Pc_label)
    plt.xlim([0, 1])
    plt.grid()
    plt.xlabel(r'$S_{w}$', fontsize=label_font_size)
    plt.ylabel(r'$P_c$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Capillary Pressure', fontsize=title_font_size)
    # save_txt = save_dir + 'pc_graph.png'
    # f1.savefig(save_txt, dpi=resolution)

    return


def plot_rel_perm_data(Sw, krw, krnw,
                       label_krw=r'$k_{r\ w}$', line_color_krw='b', marker_krw='o', line_style_krw='-',
                       label_krnw=r'$k_{r\ nw}$', line_color_krnw='r', marker_krnw='o', line_style_krnw='-',
                       label_font_size=18, title_font_size=18):

    # f1 = plt.figure()
    plt.plot(Sw, krw, color=line_color_krw, marker=marker_krw, linestyle=line_style_krw,
             linewidth=3, markersize=7, label=label_krw)
    plt.plot(Sw, krnw, color=line_color_krnw, marker=marker_krnw, linestyle=line_style_krnw,
             linewidth=3, markersize=7, label=label_krnw)
    plt.xlim([0, 1])
    plt.grid()
    plt.xlabel(r'$S_{w}$', fontsize=label_font_size)
    plt.ylabel(r'$k_r$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Relative Permeability', fontsize=title_font_size)
    plt.legend(fontsize=label_font_size)
    # save_txt = save_dir + 'rel_perm.png'
    # f1.savefig(save_txt, dpi=resolution)

    return


def plot_pc_and_rel_perm(Sw, Pc, krw, krnw,
                         Pc_label=r'$P_c$', line_color_pc='k', marker_pc='o', line_style_pc='-',
                         label_krw=r'$k_{r\ w}$', line_color_krw='b', marker_krw='o', line_style_krw='-',
                         label_krnw=r'$k_{r\ nw}$', line_color_krnw='r', marker_krnw='o', line_style_krnw='-'):
    """Plot Rel perm and Pc"""

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    # f1 = plt.figure(dpi=120, figsize=[6, 10])
    # plt.suptitle(r'250$^3$ LBM Real Grains, $\phi = 10\%$, $\theta = 60^o$', fontsize=20)
    plt.subplot(2, 1, 1)
    plot_capillary_pressure_data(Sw, Pc, Pc_label=Pc_label,
                                 line_color=line_color_pc, marker=marker_pc, line_style=line_style_pc)

    plt.subplot(2, 1, 2)
    plot_rel_perm_data(Sw, krw, krnw,
                       label_krw=label_krw, line_color_krw=line_color_krw, marker_krw=marker_krw, line_style_krw=line_style_krw,
                       label_krnw=label_krnw, line_color_krnw=line_color_krnw, marker_krnw=marker_krnw, line_style_krnw=line_style_krnw)
    plt.tight_layout()

    return


def create_image_plate():
    print("Not implemented yet!")

    return
