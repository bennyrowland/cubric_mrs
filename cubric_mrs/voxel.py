import numpy as np
import PIL
import pyx

voxel_outline_color = pyx.color.rgb(1, 1, 0)
voxel_outline_width = pyx.style.linewidth(1.0)


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
