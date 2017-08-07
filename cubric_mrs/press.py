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


def analyse_press(data_path, t1_path=None, wref_path=None, out_path=None):
    data = suspect.io.load_twix(data_path)[:, 1]

    if wref_path is not None:
        wref = suspect.io.load_twix(wref_path)[:, 1]

        channel_weights = suspect.processing.channel_combination.svd_weighting(np.mean(wref, axis=0))
        cc_wref = suspect.processing.channel_combination.combine_channels(wref, channel_weights)
    else:
        channel_weights = suspect.processing.channel_combination.svd_weighting(np.mean(data, axis=0))

    cc_data = suspect.processing.channel_combination.combine_channels(data, channel_weights)

    fc_data = np.mean(cc_data, axis=0)
    print(fc_data.dtype)

    tarquin_results = suspect.io.tarquin.process(fc_data)

    one_page_canvas = pyx.canvas.canvas()
    files_dict = {
        "PRESS": os.path.basename(data_path)
    }
    if wref_path is not None:
        files_dict["WREF"] = os.path.basename(wref_path)
    if t1_path is not None:
        files_dict["T1"] = os.path.basename(t1_path)
    one_page_canvas.insert(table.file_table(files_dict))

    if t1_path is not None:
        t1 = suspect.image.load_dicom_volume(t1_path)
        voxel_canvases = voxel.get_voxel_slices(t1, data)
        one_page_canvas.insert(voxel_canvases,
                               [pyx.trafo.translate(0, -0.5 - voxel_canvases.bbox().height())])

    tarquin_plot = plot_fitted_spectrum(tarquin_results["plots"]["data"],
                                        tarquin_results["plots"]["fit"] + tarquin_results["plots"]["baseline"])

    one_page_canvas.insert(tarquin_plot,
                           [pyx.trafo.translate(0, -1 - voxel_canvases.bbox().height() - tarquin_plot.height)])

    conc_table = table.metabolite_table(tarquin_results["metabolite_fits"],
                                        float(tarquin_results["metabolite_fits"]["TCr"]["concentration"]))
    one_page_canvas.insert(conc_table,
                           [pyx.trafo.translate(10, -1 - voxel_canvases.bbox().height())])

    quality_table = table.quality_table(tarquin_results["quality"])
    print(conc_table.bbox())
    one_page_canvas.insert(quality_table,
                           [pyx.trafo.translate(10, conc_table.bbox().bottom() - 1.5 - voxel_canvases.bbox().height())])

    output_sheet = pyx.document.page(one_page_canvas,
                                     paperformat=pyx.document.paperformat.A4)
    output_document = pyx.document.document([output_sheet])

    output_document.writePDFfile(out_path)


def plot_fitted_spectrum(data, fit):

    gaba_plot = pyx.graph.graphxy(width=9, height=6,
                                  x=pyx.graph.axis.linear(min=max_ppm, max=min_ppm),
                                  y=pyx.graph.axis.linear(
                                      #min=np.amin(data.real) * 0.2,
                                      #max=np.amax(data.real) * 1.2,
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

    return gaba_plot


def press_script():
    parser = argparse.ArgumentParser()

    parser.add_argument("--t1",
                        help="path to the T1 structural image in Nifti or DICOM format",
                        default=None)
    parser.add_argument("--press",
                        help="path to the press twix file",
                        required=True)
    parser.add_argument("--wref",
                        help="path to the water reference twix file",
                        default=None)
    parser.add_argument("--out_pdf",
                        help="path to save the one-page pdf output")

    args = parser.parse_args()

    analyse_press(args.mega,
                  t1_path=args.t1,
                  wref_path=args.wref,
                  out_path=args.out_pdf)
