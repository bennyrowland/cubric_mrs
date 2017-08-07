import argparse
import numpy as np
import pyx
import suspect
import os

from cubric_mrs import voxel, table

# define some useful global parameters
voxel_outline_color = pyx.color.rgb(1, 1, 0)
voxel_outline_width = pyx.style.linewidth(1.0)
plot_linewidth = pyx.style.linewidth(0.01)
simple_plot_color = pyx.color.rgb(0, 0, 1)
sr_plot_color = pyx.color.rgb(1, 0, 0)
plot_linewidth_gaba = pyx.style.linewidth(0.03)
plot_color_gaba = pyx.color.rgb(0, 1, 0)
plot_color_residual = pyx.color.cmyk.Gray
plot_color_fit = pyx.color.rgb.red
plot_color_data = pyx.color.rgb.black
min_ppm = 0.2
max_ppm = 5.0


def analyse_mega(mega_path, t1_path=None, wref_path=None, out_path=None):
    mega = suspect.io.load_twix(mega_path)

    # start with channel combination, obviously
    channel_weights = suspect.processing.channel_combination.svd_weighting(np.mean(mega, axis=(0, 1)))
    cc_mega = suspect.processing.channel_combination.combine_channels(mega, channel_weights)

    # frequency correction next, as this affects some later stuff
    def correct_frequency_sr(target):
        def correct_fid(fid):
            frequency_shift, phase_shift = suspect.processing.frequency_correction.spectral_registration(fid, target)
            return fid.adjust_frequency(-frequency_shift).adjust_phase(-phase_shift)
        return correct_fid
    sr_aligned_on = np.apply_along_axis(correct_frequency_sr(cc_mega[108, 0]), 1, cc_mega[:, 0])
    sr_aligned_off = np.apply_along_axis(correct_frequency_sr(cc_mega[108, 1]), 1, cc_mega[:, 1])
    sr_on = np.mean(sr_aligned_on, axis=0)
    sr_off = np.mean(sr_aligned_off, axis=0)
    # do an additional SR to align the on with the off
    global_frequency_shift, global_phase_shift = suspect.processing.frequency_correction.spectral_registration(sr_off,
                                                                                                               sr_on)
    sr_off = sr_off.adjust_frequency(-global_frequency_shift).adjust_phase(-global_phase_shift)
    # align the biggest peak in the diff spectrum (NAA) to 2.01ppm
    diff = sr_on - sr_off
    sr_peak_index = np.argmax(np.abs(diff.spectrum()))
    sr_peak_frequency = mega.frequency_axis_ppm()[sr_peak_index]
    frequency_shift = mega.f0 * (sr_peak_frequency - 2.01)
    # we have to correct this frequency on a lot of pieces of data
    sr_on = sr_on.adjust_frequency(frequency_shift)
    sr_off = sr_off.adjust_frequency(frequency_shift)
    sr_aligned_on = sr_aligned_on.adjust_frequency(frequency_shift)
    sr_aligned_off = sr_aligned_off.adjust_frequency(frequency_shift)
    cc_mega = cc_mega.adjust_frequency(frequency_shift)

    # now we will do some water removal from sr_off, which makes the phase estimation better
    components = suspect.processing.water_suppression.hsvd(sr_off, 20)
    water_components = [component for component in components if component["frequency"] < 70]
    water_fid = suspect.processing.water_suppression.construct_fid(water_components, mega.time_axis())
    dry_off = sr_off - water_fid

    # now do the phase estimation
    zp, fp = suspect.processing.phase.mag_real(dry_off)

    sr_diff = (sr_on - sr_off).adjust_phase(zp, fp)

    # also get the un-frequency-corrected diff
    phased_mega = cc_mega.adjust_phase(zp, fp)
    simple_on = np.mean(phased_mega[:, 0], axis=0)
    simple_off = np.mean(phased_mega[:, 1], axis=0)
    simple_diff = simple_on - simple_off

    tarquin_results = suspect.io.tarquin.process(sr_diff, {"pul_seq": "mega_press",
                                                           "int_basis": "megapress_gaba",
                                                           "ref_signals": "1h_naa",
                                                           "start_pnt": "10",
                                                           })

    tarquin_off_results = suspect.io.tarquin.process(sr_off.adjust_phase(zp, fp))

    diff_comparison = plot_diff_comparison(simple_diff, sr_diff)
    drift_heatmap = plot_drift_heatmap(phased_mega[:, 1], sr_aligned_off.adjust_phase(zp, fp))

    off_spectrum_plot = plot_off_spectrum(tarquin_off_results["plots"]["data"],
                                          tarquin_off_results["plots"]["fit"] + tarquin_off_results["plots"]["baseline"])

    gaba_plot = plot_fitted_gaba_spectrum(tarquin_results["plots"]["data"],
                                          tarquin_results["plots"]["fit"] + tarquin_results["plots"]["baseline"],
                                          tarquin_results["plots"]["metabolites"]["GABA_A"] +
                                          tarquin_results["plots"]["metabolites"]["GABA_B"] +
                                          tarquin_results["plots"]["baseline"])

    if t1_path is not None:
        t1 = suspect.image.load_nifti(t1_path)
        voxel_canvases = voxel.get_voxel_slices(t1, mega)

    conc_table = table.metabolite_table(tarquin_results["metabolite_fits"],
                                  float(tarquin_off_results["metabolite_fits"]["TCr"]["concentration"]))

    qual_table = table.quality_table(tarquin_results["quality"])

    one_page_canvas = pyx.canvas.canvas()
    one_page_canvas.insert(table.file_table({
        "T1": os.path.basename(t1_path),
        "MEGAPRESS": os.path.basename(mega_path)
    }), [pyx.trafo.translate(0, 6.0)])
    one_page_canvas.insert(diff_comparison, [pyx.trafo.translate(0, -diff_comparison.height - 0.5)])
    one_page_canvas.insert(drift_heatmap, [pyx.trafo.translate(10, -diff_comparison.height - 0.5)])
    one_page_canvas.insert(gaba_plot, [pyx.trafo.translate(0,
                                                           - diff_comparison.height
                                                           - gaba_plot.height - 1.5)])
    one_page_canvas.insert(off_spectrum_plot, [pyx.trafo.translate(0,
                                                                   - off_spectrum_plot.height
                                                                   - diff_comparison.height
                                                                   - gaba_plot.height - 2.5)])

    if t1_path is not None:
        one_page_canvas.insert(voxel_canvases)

    one_page_canvas.insert(conc_table, [pyx.trafo.translate(10, -diff_comparison.height - 2)])

    one_page_canvas.insert(qual_table,
                           [pyx.trafo.translate(10, -diff_comparison.height - 2.5 - conc_table.bbox().height())])

    output_sheet = pyx.document.page(one_page_canvas,
                                     paperformat=pyx.document.paperformat.A4)
    output_document = pyx.document.document([output_sheet])

    output_document.writePDFfile(out_path)


def plot_diff_comparison(simple_diff, corrected_diff):
    raw_diff_plot = pyx.graph.graphxy(width=9, height=6,
                                      x=pyx.graph.axis.linear(min=max_ppm, max=min_ppm),
                                      y=pyx.graph.axis.linear(min=np.amin(corrected_diff.spectrum().real) * 0.2,
                                                              max=np.amax(corrected_diff.spectrum().real) * 1.2,
                                                              parter=None))
    raw_diff_plot.plot(pyx.graph.data.values(
        x=simple_diff.frequency_axis_ppm(),
        y=simple_diff.spectrum().real
    ),
        [pyx.graph.style.line([plot_linewidth, simple_plot_color])]
    )
    raw_diff_plot.plot(pyx.graph.data.values(
        x=corrected_diff.frequency_axis_ppm(),
        y=corrected_diff.spectrum().real
    ),
        [pyx.graph.style.line([plot_linewidth, sr_plot_color])]
    )
    return raw_diff_plot


def plot_drift_heatmap(simple_data, corrected_data):
    # plot heatmap of the spectra over the repetitions to evaluate
    # frequency drift
    ppm_slice = simple_data.spectrum().slice_ppm(3.12, 2.72)
    sr_drift_map = np.real(corrected_data.spectrum()[:, ppm_slice])
    simple_drift_map = np.real(simple_data.spectrum()[:, ppm_slice].reshape(-1))

    x_axis = np.zeros_like(sr_drift_map)
    x_axis[:] = np.arange(sr_drift_map.shape[0])[:, np.newaxis]
    x_axis = x_axis.reshape(-1)

    y_axis = np.zeros_like(sr_drift_map)
    y_axis[:] = simple_data.frequency_axis_ppm()[ppm_slice]
    y_axis = y_axis.reshape(-1)

    sr_drift_plot = pyx.graph.graphxy(height=3, width=9,
                                      x=pyx.graph.axis.linear(min=0, max=sr_drift_map.shape[0], title="Repetitions"),
                                      y=pyx.graph.axis.linear(min=simple_data.frequency_axis_ppm()[ppm_slice][-1],
                                                              max=simple_data.frequency_axis_ppm()[ppm_slice][0])
                                      )

    sr_drift_data = pyx.graph.data.points(
        np.stack((x_axis, y_axis, sr_drift_map.reshape(-1)), axis=-1),
        x=1, y=2, color=3)
    sr_drift_plot.plot(sr_drift_data,
                       [pyx.graph.style.density(gradient=pyx.color.gradient.Jet,
                                                keygraph=None)])

    simple_drift_plot = pyx.graph.graphxy(height=3, width=9,
                                          x=pyx.graph.axis.linkedaxis(sr_drift_plot.axes["x"]),
                                          y=pyx.graph.axis.linear(min=simple_data.frequency_axis_ppm()[ppm_slice][-1],
                                                                  max=simple_data.frequency_axis_ppm()[ppm_slice][0])
                                          )
    simple_drift_plot.plot(pyx.graph.data.points(
        np.stack((x_axis, y_axis, simple_drift_map), axis=-1),
        x=1, y=2, color=3),
        [pyx.graph.style.density(gradient=pyx.color.gradient.Jet,
                                 keygraph=None)]
    )

    combined_canvas = pyx.canvas.canvas()
    combined_canvas.insert(simple_drift_plot, [pyx.trafo.translate(0, 3)])
    combined_canvas.insert(sr_drift_plot)
    return combined_canvas


def plot_fitted_gaba_spectrum(data, fit, gaba):
    gaba_plot = pyx.graph.graphxy(width=9, height=6,
                                  x=pyx.graph.axis.linear(min=max_ppm, max=min_ppm),
                                  y=pyx.graph.axis.linear(
                                      min=np.amin(data.real) * 0.2,
                                      max=np.amax(data.real) * 1.2,
                                      parter=None
                                  ),
                                  key=pyx.graph.key.key(pos="tr")
                                  )

    gaba_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=data.spectrum().real,
        title="Data"

    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_data])]
    )

    gaba_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=fit.real,
        title="Fit"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_fit])]
    )

    gaba_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=(data.real - fit.real - np.amax(data.real) * 0.2),
        title="Residual"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_residual])]
    )

    gaba_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=gaba.real,
        title="GABA"
    ),
        [pyx.graph.style.line([plot_linewidth_gaba, plot_color_gaba])]
    )

    return gaba_plot


def plot_off_spectrum(data, fit):
    off_plot = pyx.graph.graphxy(width=9, height=6,
                                 x=pyx.graph.axis.linear(min=max_ppm, max=min_ppm, title="ppm"),
                                 y=pyx.graph.axis.linear(
                                     min=-np.amax(data.real) * 0.5,
                                     max=np.amax(data.real) * 1.2,
                                     parter=None),
                                 key=pyx.graph.key.key(pos="tr")
                                 )

    off_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=data.real,
        title="Data"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_data])]
    )

    off_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=fit.real,
        title="Fit"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_fit])]
    )

    off_plot.plot(pyx.graph.data.values(
        x=data.frequency_axis_ppm(),
        y=(data.real - fit.real - np.amax(data.real) * 0.2),
        title="Residual"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_residual])]
    )
    return off_plot


def megapress_script():
    parser = argparse.ArgumentParser()

    parser.add_argument("--t1",
                        help="path to the T1 structural image",
                        default=None)
    parser.add_argument("--mega",
                        help="path to the megapress twix file",
                        required=True)
    parser.add_argument("--out_pdf",
                        help="path to save the one-page pdf output")

    args = parser.parse_args()

    analyse_mega(args.mega, t1_path=args.t1, out_path=args.out_pdf)
