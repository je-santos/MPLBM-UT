import matplotlib.pyplot as plt


def plot_capillary_pressure_data(Sw, Pc, label_font_size=18, title_font_size=18):
    # f1 = plt.figure()
    plt.plot(Sw, Pc, 'ko-', linewidth=3, markersize=7)
    plt.xlim([0, 1])
    plt.grid()
    plt.xlabel(r'$S_{w}$', fontsize=label_font_size)
    plt.ylabel(r'$P_c$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Capillary Pressure', fontsize=title_font_size)
    # save_txt = save_dir + 'pc_graph.png'
    # f1.savefig(save_txt, dpi=resolution)

    return


def plot_rel_perm_data(Sw, krw, krnw, label_font_size=18, title_font_size=18):

    # f1 = plt.figure()
    plt.plot(Sw, krw, 'bo-', linewidth=3, markersize=7)
    plt.plot(Sw, krnw, 'ro-', linewidth=3, markersize=7)
    plt.xlim([0, 1])
    plt.grid()
    plt.xlabel(r'$S_{w}$', fontsize=label_font_size)
    plt.ylabel(r'$k_r$ [LBM Units]', fontsize=label_font_size)
    plt.title(r'Relative Permeability', fontsize=title_font_size)
    plt.legend([r'$k_{r\ w}$', r'$k_{r\ nw}$'], fontsize=label_font_size)
    # save_txt = save_dir + 'rel_perm.png'
    # f1.savefig(save_txt, dpi=resolution)

    return


def plot_pc_and_rel_perm(Sw, Pc, krw, krnw):
    """Plot Rel perm and Pc"""

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    # f1 = plt.figure(dpi=120, figsize=[6, 10])
    # plt.suptitle(r'250$^3$ LBM Real Grains, $\phi = 10\%$, $\theta = 60^o$', fontsize=20)
    plt.subplot(2, 1, 1)
    plot_capillary_pressure_data(Sw, Pc)

    plt.subplot(2, 1, 2)
    plot_rel_perm_data(Sw, krw, krnw)
    plt.tight_layout()

    return


def create_image_plate():
    print("Not implemented yet!")

    return
