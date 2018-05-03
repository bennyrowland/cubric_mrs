import numpy as np
import PIL
import pyx
import os
import nipype
from nipype.interfaces import fsl

import suspect

voxel_outline_color = pyx.color.rgb(1, 1, 0)
voxel_outline_width = pyx.style.linewidth(1.0)

beta_csf = 1
beta_gm = 0.78
beta_wm = 0.65

water_conc = 55556


def make_canvas(slice_array):
    scaling_factor = 255 / np.amax(slice_array)
    slice_image = PIL.Image.fromarray((slice_array * scaling_factor).astype(np.uint8))
    slice_canvas = pyx.canvas.canvas()
    slice_canvas.insert(pyx.bitmap.bitmap(0, 0, slice_image,
                                          height=slice_image.height),
                        [pyx.trafo.scale(1, -1, 0, slice_image.height / 2)])
    return slice_canvas


def get_voxel_path(volume, mrs):
    voxel_coords_mrs = [[-0.5, -0.5, -0.5],
                        [0.5, -0.5, -0.5],
                        [0.5, 0.5, -0.5],
                        [-0.5, 0.5, -0.5],
                        [-0.5, -0.5, -0.5],
                        [-0.5, -0.5, 0.5],
                        [0.5, -0.5, 0.5],
                        [0.5, 0.5, 0.5],
                        [-0.5, 0.5, 0.5],
                        [-0.5, -0.5, 0.5],
                        [0.5, -0.5, 0.5],
                        [0.5, -0.5, -0.5],
                        [0.5, 0.5, -0.5],
                        [0.5, 0.5, 0.5],
                        [-0.5, 0.5, 0.5],
                        [-0.5, 0.5, -0.5]]
    voxel_coords_volume = volume.from_scanner(mrs.to_scanner(voxel_coords_mrs))
    path_points = [pyx.path.moveto(*voxel_coords_volume[0, 0:2])] + \
                  [pyx.path.lineto(*point) for point in voxel_coords_volume[1:, 0:2]]
    path = pyx.path.path(*path_points)
    return path


def get_voxel_slices(volume, mrs):
    slices = [volume.resample(volume.coronal_vector,
                           volume.axial_vector,
                           shape=(1, 256, 256),
                           centre=(mrs.centre[0], volume.centre[1], volume.centre[2])),
              volume.resample(volume.sagittal_vector,
                              -volume.coronal_vector,
                              shape=(1, 256, 256),
                              centre=(volume.centre[0], volume.centre[1], mrs.centre[2])),
              volume.resample(volume.sagittal_vector,
                              volume.axial_vector,
                              shape=(1, 256, 256),
                              centre=(volume.centre[0], mrs.centre[1], volume.centre[2]))]

    scaling_factor = 255 / np.amax(volume)
    canvases = [make_canvas(slc * scaling_factor) for slc in slices]
    paths = [get_voxel_path(slc, mrs) for slc in slices]

    combined_canvas = pyx.canvas.canvas()
    for i in range(3):
        canvases[i].stroke(paths[i], [voxel_outline_color, voxel_outline_width])
        combined_canvas.insert(canvases[i], [pyx.trafo.scale(6 / 256).translated(6.5 * i, 0)])

    return combined_canvas


def _save_nifti(dicom_file, nifti_file=None):
    """
    This function is intended to be wrapped in a nipype Node and used to take a
    DICOM file (and other files from the same series in the same file), and
    save the volume in Nifti format instead. By default the output file will be
    in the folder created by nipype, which is usually temporary, but it can
    also be given an output path to save it somewhere specific.

    :param dicom_file:
    :param nifti_file:
    :return:
    """
    import suspect
    import os

    volume = suspect.image.load_dicom_volume(dicom_file)
    if nifti_file is None:
        nifti_file = os.path.join(os.getcwd(), "volume.nii")
    suspect.image.save_nifti(nifti_file, volume)
    return nifti_file


def classify_tissues(volume_path):
    ext = os.path.splitext(volume_path)[1]

    wf = nipype.Workflow(name="classify_wf")

    bet = nipype.Node(fsl.BET(frac=0.5,
                              robust=True),
                      name="bet")
    fast = nipype.Node(fsl.FAST(output_type="NIFTI",
                                number_classes=3),
                       name="fast")

    wf.connect([(bet, fast, [("out_file", "in_files")])])
    if ext.upper() in [".DCM", ".IMA"]:
        save_nifti = nipype.Node(nipype.Function(input_names=["dicom_file"],
                                                 output_names=["nifti_file"],
                                                 function=_save_nifti),
                                 name="save_nifti")
        save_nifti.inputs.dicom_file = volume_path
        wf.connect([(save_nifti, bet, [("nifti_file", "in_file")])])
    else:
        bet.inputs.in_file = volume_path

    result = wf.run()

    for node in result.nodes():
        if node.name == "fast":
            wm = suspect.image.load_nifti(node.result.outputs.partial_volume_files[2])
            gm = suspect.image.load_nifti(node.result.outputs.partial_volume_files[1])
            csf = suspect.image.load_nifti(node.result.outputs.partial_volume_files[0])

    return wm, gm, csf


def segment_voxel(volume_path, mask):
    wm, gm, csf = classify_tissues(volume_path)
    voxel_volume = np.sum(mask)

    f_csf = float(np.sum(csf * mask) / voxel_volume)
    f_gm = float(np.sum(gm * mask) / voxel_volume)
    f_wm = float(np.sum(wm * mask) / voxel_volume)

    tissue_beta = f_csf * beta_csf + f_gm * beta_gm + f_wm * beta_wm
    corrected_water_conc = int(water_conc * tissue_beta)

    return {
        "wm": f_wm,
        "gm": f_gm,
        "csf": f_csf,
        "volume": voxel_volume,
        "tissue_beta": tissue_beta,
        "water_conc": corrected_water_conc,
        "water_att": gasparovic_attenuation_factor(f_wm, f_gm, f_csf)
    }


def gasparovic_attenuation_factor(f_wm, f_gm, f_csf, te=30, tr=2000, t1_m=1400, t2_m=200):
    """
    Calculates the attenuation factor that needs to be applied to obtain
    milliMolal absolute concentration units according to Gasparovich 2006
    which accounts for changes in the size of the water reference peak due to
    differences in T1 and T2 for grey and white matter and CSF. It also
    corrects for the fact that metabolites only appear in GM/WM and not in CSF.

    T1 and T2 values for water and metabolites are taken from Gussew 2012
    :param f_wm:
    :param f_gm:
    :param f_csf:
    :param te:
    :param tr:
    :param t1_m:
    :param t2_m:
    :return:
    """
    t1_wm = 1080
    t2_wm = 70
    t1_gm = 1820
    t2_gm = 100
    t1_csf = 4160
    t2_csf = 500

    r_wm = np.exp(-te / t2_wm) * (1 - np.exp(-tr / t1_wm))
    r_gm = np.exp(-te / t2_gm) * (1 - np.exp(-tr / t1_gm))
    r_csf = np.exp(-te / t2_csf) * (1 - np.exp(-tr / t1_csf))
    r_m = np.exp(-te / t2_m) * (1 - np.exp(-tr / t1_m))

    water_norm = f_wm * beta_wm + f_gm * beta_gm + f_csf * beta_csf
    fw_wm = beta_wm * f_wm / water_norm
    fw_gm = beta_gm * f_gm / water_norm
    fw_csf = beta_csf * f_csf / water_norm

    return (fw_wm * r_wm + fw_gm * r_gm + fw_csf * r_csf) / (1 - fw_csf) / r_m
