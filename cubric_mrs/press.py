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


def analyse_press(data_path, t1_path=None, wref_path=None, out_path=None, out_csv=None, subject_id=None):
    data = suspect.io.load_twix(data_path)[:, 1]

    if subject_id is None:
        subject_id = data.metadata["patient_id"]

    num_channels = data.shape[-2]
    noise = np.moveaxis(data[:, :, -250:], 2, 0).reshape(num_channels, -1)

    white_data = suspect.processing.channel_combination.whiten(data, noise)
    if wref_path is not None:
        wref = suspect.io.load_twix(wref_path)[:, 1]
        white_wref = suspect.processing.channel_combination.whiten(wref, noise)
        channel_weights = suspect.processing.channel_combination.svd_weighting(np.mean(white_wref, axis=0))
        cc_wref = suspect.processing.channel_combination.combine_channels(white_wref, channel_weights)
        # assume no frequency drift over water signal, and doesn't matter for singlet peak anyway
        wref_final = np.mean(cc_wref, axis=0)
        # calculate the eddy current from the water signal
        ec = np.unwrap(np.angle(wref_final))
        # smooth the eddy current signal
        ec_smooth = suspect.processing.denoising.sliding_gaussian(ec, 32)
        ecc = np.exp(-1j * ec_smooth)
        wref_final *= ecc
        white_data *= ecc
    else:
        channel_weights = suspect.processing.channel_combination.svd_weighting(np.mean(white_data, axis=0))

    cc_data = suspect.processing.channel_combination.combine_channels(white_data, channel_weights)

    # re-bin into phase cycles and average over each cycle
    pc_data = np.mean(cc_data.reshape(cc_data.shape[0] // 16, 16, -1), axis=1)

    def correct_frequency_sr(target):
        def correct_fid(fid):
            frequency_shift, phase_shift = suspect.processing.frequency_correction.spectral_registration(fid, target)
            return fid.adjust_frequency(-frequency_shift).adjust_phase(-phase_shift)

        return correct_fid

    sr_data = np.mean(np.apply_along_axis(correct_frequency_sr(pc_data[0]), 1, pc_data), axis=0)

    one_page_canvas = pyx.canvas.canvas()
    files_dict = {
        "PRESS": os.path.basename(data_path)
    }
    if wref_path is not None:
        files_dict["WREF"] = os.path.basename(wref_path)
    if t1_path is not None:
        files_dict["T1"] = os.path.basename(t1_path)
    file_table = table.file_table(files_dict)
    one_page_canvas.insert(file_table,
                           [pyx.trafo.translate(1, 28)])

    if t1_path is not None:
        t1_path = os.path.abspath(t1_path)
        if os.path.splitext(t1_path)[1].upper() in [".IMA", ".DCM"]:
            t1 = suspect.image.load_dicom_volume(t1_path)
        elif os.path.splitext(t1_path)[1].upper() in [".NII"]:
            t1 = suspect.image.load_nifti(t1_path)
        elif t1_path.upper().endswith(".NII.GZ"):
            t1 = suspect.image.load_nifti(t1_path)
        else:
            print("could not load t1 from {}".format(os.path.splitext(t1_path)[1].upper()))
            exit(-1)
        voxel_canvases = voxel.get_voxel_slices(t1, data)
        segmentation = voxel.segment_voxel(t1_path, suspect.image.create_mask(data, t1))
        one_page_canvas.insert(voxel_canvases,
                               [pyx.trafo.translate(1, 28 - 0.5 - voxel_canvases.bbox().height())])

    water = wref_final if wref_path is not None else None
    tarquin_results = suspect.io.tarquin.process(sr_data, water, options={
        "w_conc": segmentation["water_conc"],
        "w_att": segmentation["water_att"]
    })

    tarquin_plot = plot_fitted_spectrum(tarquin_results["plots"]["data"],
                                        tarquin_results["plots"]["fit"] + tarquin_results["plots"]["baseline"])

    one_page_canvas.insert(tarquin_plot,
                           [pyx.trafo.translate(1, 28 - 1 - voxel_canvases.bbox().height() - tarquin_plot.height)])

    conc_table = table.metabolite_table(tarquin_results["metabolite_fits"],
                                        float(tarquin_results["metabolite_fits"]["TCr"]["concentration"]),
                                        "mM" if wref_path is not None else "A.U.")
    one_page_canvas.insert(conc_table,
                           [pyx.trafo.translate(11, 28 - 1 - voxel_canvases.bbox().height())])

    quality_table = table.quality_table(tarquin_results["quality"])
    #print(conc_table.bbox())
    one_page_canvas.insert(quality_table,
                           [pyx.trafo.translate(11, 28 - conc_table.bbox().height() - 1.5 - voxel_canvases.bbox().height())])

    properties_table = table.voxel_properties_table(segmentation)
    one_page_canvas.insert(properties_table,
                           [pyx.trafo.translate(1, 28 - 2 - voxel_canvases.bbox().height() - tarquin_plot.height)])

    output_sheet = pyx.document.page(one_page_canvas,
                                     paperformat=pyx.document.paperformat.A4,
                                     centered=False)
    output_document = pyx.document.document([output_sheet])

    output_document.writePDFfile(out_path)

    if out_csv is not None:
        output_csv(out_csv, subject_id, tarquin_results["metabolite_fits"])


def plot_fitted_spectrum(data, fit):
    data_range = np.amax(data.real) - np.amin(data.real)

    gaba_plot = pyx.graph.graphxy(width=9, height=6,
                                  x=pyx.graph.axis.linear(min=max_ppm, max=min_ppm),
                                  y=pyx.graph.axis.linear(
                                      min=np.amin(data.real) - 0.3 * data_range,
                                      max=np.amax(data.real) + 0.1 * data_range,
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
        y=(data.real - fit.real - data_range * 0.15),
        title="Residual"
    ),
        [pyx.graph.style.line([plot_linewidth, plot_color_residual])]
    )

    return gaba_plot


def output_csv(filename, subject_id, concentrations):
    with open(filename, 'a') as fout:
        if fout.tell() == 0:
            # this is a new file, have to output the header line
            fout.write("SubjectID")
            for metabolite_name in sorted(concentrations.keys()):
                fout.write(", {0}, {0}_SD".format(metabolite_name))
        fout.write("\n")
        # now write the actual concentration values
        fout.write("{}".format(subject_id))
        for name, data in sorted(concentrations.items()):
            fout.write(", {}, {}".format(data["concentration"], data["sd"]))


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
                        help="path to save the one-page pdf output",
                        default="out.pdf")
    parser.add_argument("--out_csv",
                        help="path to save the concentration csv file",
                        default=None)
    parser.add_argument("--id",
                        help="override the patient id in output files",
                        default=None)

    args = parser.parse_args()

    analyse_press(args.press,
                  t1_path=args.t1,
                  wref_path=args.wref,
                  out_path=args.out_pdf,
                  out_csv=args.out_csv,
                  subject_id=args.id)
