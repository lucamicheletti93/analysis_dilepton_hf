"""
Helper file with fit utils
"""

import sys
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter

def perform_fit(file_name, histo_name, fitter_name, sgn_func, bkg_func,
                mass_min, mass_max, init_mass, part_name,
                filen_name_corrbkg, histo_name_corrbkg, histo_name_signalmc,
                fix_frac_bkgcorr, fix_tailpars_from_mc, signal_pars_tofix={},
                rebin=1, frac_bkgcorr=1):
    """

    """
    bkg_funcs = [bkg_func]
    sgn_funcs = [sgn_func]
    labels_signal_pdf = []
    if part_name == "dplus":
        labels_signal_pdf.append(r"$\mathrm{D}^+$")
        if mass_max > 1.99:
            sgn_funcs.append("gaussian")
            labels_signal_pdf.append(
                r"$\mathrm{D}^{*+} + \mathrm{D}^+\rightarrow\mathrm{\pi\pi\pi}$")
    elif part_name == "dstar":
        labels_signal_pdf.append(r"$\mathrm{D}^{*+}$")
    elif part_name == "dzero":
        labels_signal_pdf.append(r"$\mathrm{D}^{0}$")
    labels_bkg_pdf = ["Comb. bkg."]

    data_hdl = DataHandler(file_name, histoname=histo_name, limits=[mass_min, mass_max], rebin=rebin)
    data_hdl_corrbkg, data_hdl_signalmc = None, None
    frac_from_mc = 0.
    if filen_name_corrbkg is not None:
        data_hdl_corrbkg = DataHandler(
            filen_name_corrbkg, histoname=histo_name_corrbkg,
            limits=[mass_min, mass_max])
        data_hdl_signalmc = DataHandler(
            filen_name_corrbkg, histoname=histo_name_signalmc,
            limits=[mass_min, mass_max])
        frac_from_mc = data_hdl_corrbkg.get_norm() / data_hdl_signalmc.get_norm()
        bkg_funcs.insert(0, "hist")
        if part_name == "dplus":
            labels_bkg_pdf.insert(0, r"$\mathrm{D}^+,\mathrm{D_s}^+\rightarrow\mathrm{KK\pi}$")
        elif part_name == "dzero":
            labels_bkg_pdf.insert(0, r"$\mathrm{K-\pi}$ reflected")
    else:
        if fix_tailpars_from_mc or fix_frac_bkgcorr and sgn_funcs[0] in ["doublecb", "doublecbsymm"]:
            print("ERROR, parameters cannot be fixed from MC if file from MC not provided")
            sys.exit()

    fitter_mc = None
    if fix_tailpars_from_mc and data_hdl_signalmc is not None and sgn_funcs[0] in ["doublecb", "doublecbsymm"]:
        fitter_mc = F2MassFitter(data_hdl_signalmc,
                                 name_signal_pdf=[sgn_funcs[0]],
                                 name_background_pdf=["nobkg"],
                                 name=f"{fitter_name}_mc", tol=0.1)
        fitter_mc.set_particle_mass(0, mass=init_mass, limits=[init_mass * 0.95, init_mass * 1.05])
        fitter_mc.set_signal_initpar(0, "sigma", 0.005, limits=[0.0001, 0.06])
        if sgn_funcs[0] == "doublecbsymm":
            fitter_mc.set_signal_initpar(0, "alpha", 1.5, limits=[1., 3.])
            fitter_mc.set_signal_initpar(0, "n", 50, limits=[0., 100.])
        else:
            fitter_mc.set_signal_initpar(0, "alphar", 1.5, limits=[1., 3.])
            fitter_mc.set_signal_initpar(0, "alphal", 1.5, limits=[1., 3.])
            fitter_mc.set_signal_initpar(0, "nr", 50, limits=[0., 100.])
            fitter_mc.set_signal_initpar(0, "nl", 50, limits=[0., 100.])
        fitter_mc.mass_zfit()
        if sgn_funcs[0] == "doublecbsymm":
            signal_pars_tofix[0] = {
                "alpha": fitter_mc.get_signal_parameter(0, "alpha")[0],
                "n": fitter_mc.get_signal_parameter(0, "n")[0]
            }
        else:
            signal_pars_tofix[0] = {
                "alphar": fitter_mc.get_signal_parameter(0, "alphar")[0],
                "alphal": fitter_mc.get_signal_parameter(0, "alphal")[0],
                "nr": fitter_mc.get_signal_parameter(0, "nr")[0],
                "nl": fitter_mc.get_signal_parameter(0, "nl")[0]
            }

    fitter = F2MassFitter(data_hdl,
                          name_signal_pdf=sgn_funcs,
                          name_background_pdf=bkg_funcs,
                          name=fitter_name, tol=0.1,
                          label_signal_pdf=labels_signal_pdf,
                          label_bkg_pdf=labels_bkg_pdf)
    ibkg = 0
    if data_hdl_corrbkg is not None:
        fitter.set_background_template(ibkg, data_hdl_corrbkg)
        fitter.set_background_initpar(ibkg, "frac", 0.01, limits=[0., 0.08])
        if fix_frac_bkgcorr:
            fitter.fix_bkg_frac_to_signal_pdf(ibkg, 0, frac_from_mc * frac_bkgcorr)
        ibkg += 1
    if part_name == "dstar":
        fitter.set_particle_mass(ibkg, mass=init_mass, limits=[init_mass * 0.95, init_mass * 1.05])
        fitter.set_signal_initpar(0, "sigma", 0.0005, limits=[0.0001, 0.002])
        fitter.set_signal_initpar(0, "frac", 0.1)
        fitter.set_signal_initpar(0, "alphal", 1.5, limits=[1., 3.])
        fitter.set_signal_initpar(0, "alphar", 1.5, limits=[1., 3.])
        fitter.set_signal_initpar(0, "nl", 50, limits=[30., 100.])
        fitter.set_signal_initpar(0, "nr", 50, limits=[30., 100.])
        fitter.set_background_initpar(ibkg, "power", 0.5)
        fitter.set_background_initpar(ibkg, "c1", -20)
        fitter.set_background_initpar(ibkg, "c2", 500)
        fitter.set_background_initpar(ibkg, "c3", -5000)
    elif part_name == "dplus":
        fitter.set_particle_mass(0, mass=init_mass, limits=[init_mass * 0.95, init_mass * 1.05])
        fitter.set_signal_initpar(0, "sigma", 0.005, limits=[0.0001, 0.04])
        fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 0.5])
        if mass_max > 1.99:
            fitter.set_particle_mass(1, mass=2.01, limits=[1.99, 2.03])
            fitter.set_signal_initpar(1, "sigma", 0.02, limits=[0.001, 0.08])
            fitter.set_signal_initpar(1, "frac", 0.0001, limits=[0., 0.005])
        fitter.set_background_initpar(ibkg, "c0", 0.4)
        fitter.set_background_initpar(ibkg, "c1", -0.2)
        fitter.set_background_initpar(ibkg, "c2", -0.01)
        fitter.set_background_initpar(ibkg, "c3", 0.01)
        fitter.set_background_initpar(ibkg, "lam", -1.)
    elif part_name == "dzero":
        fitter.set_particle_mass(0, mass=init_mass, limits=[init_mass * 0.95, init_mass * 1.05])
        fitter.set_signal_initpar(0, "sigma", 0.005, limits=[0.0001, 0.06])
        fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 0.5])
        if sgn_funcs[0] == "doublecbsymm":
            fitter.set_signal_initpar(0, "alpha", 1.5, limits=[1., 3.])
            fitter.set_signal_initpar(0, "n", 50, limits=[0., 100.])
        elif sgn_funcs[0] == "doublecb":
            fitter.set_signal_initpar(0, "alphal", 1.5, limits=[1., 3.])
            fitter.set_signal_initpar(0, "alphar", 1.5, limits=[1., 3.])
            fitter.set_signal_initpar(0, "nl", 50, limits=[0., 100.])
            fitter.set_signal_initpar(0, "nr", 50, limits=[0., 100.])
        fitter.set_background_initpar(ibkg, "c0", 0.4)
        fitter.set_background_initpar(ibkg, "c1", -0.2)
        fitter.set_background_initpar(ibkg, "c2", -0.01)
        fitter.set_background_initpar(ibkg, "c3", 0.01)
        fitter.set_background_initpar(ibkg, "lam", -1.)

    if len(signal_pars_tofix) > 0:
        for isig in signal_pars_tofix:
            for par in signal_pars_tofix[isig]:
                fitter.set_signal_initpar(isig, par, signal_pars_tofix[isig][par], fix=True)

    fitter.mass_zfit()

    return fitter, fitter_mc
