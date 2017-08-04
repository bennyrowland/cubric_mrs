import nibabel
import numpy as np
import PIL
import pyx
import scipy.interpolate
import suspect


def make_index_coords(shape, size):
    x_pixel_size = size[0] / shape[0]
    y_pixel_size = size[1] / shape[1]

    i_coords = np.linspace(-x_pixel_size * (shape[0] - 1) / 2,
                           x_pixel_size * (shape[0] - 1) / 2,
                           shape[0],
                           endpoint=True)
    j_coords = np.linspace(-y_pixel_size * (shape[1] - 1) / 2,
                           y_pixel_size * (shape[1] - 1) / 2,
                           shape[1],
                           endpoint=True)
    return np.meshgrid(i_coords, j_coords)


def make_canvas(slice_array):
    scaling_factor = 255 / np.amax(slice_array)
    slice_image = PIL.Image.fromarray((slice_array * scaling_factor).astype(np.uint8))
    slice_canvas = pyx.canvas.canvas()
    slice_canvas.insert(pyx.bitmap.bitmap(0, 0, slice_image,
                                          height=slice_image.height),
                        [pyx.trafo.scale(1, -1, 0, slice_image.height / 2)])
    return slice_canvas


def get_sagittal_slice(volume, x=0, shape=(256, 256), size=(256, 256)):

    volume_centre = volume.to_scanner((np.array(volume.shape[::-1]) - 1) / 2)

    II, JJ = make_index_coords(shape, size)

    slice_coords = II[..., np.newaxis] * volume.coronal_vector + JJ[..., np.newaxis] * volume.axial_vector + volume_centre
    slice_coords[:, :, 0] = x

    slice_array = scipy.interpolate.interpn([np.arange(dim) for dim in volume.shape],
                                            volume,
                                            volume.from_scanner(slice_coords)[..., ::-1],
                                            "linear",
                                            bounds_error=False,
                                            fill_value=0)

    return make_canvas(slice_array)


def get_axial_slice(volume, z=0, shape=(256, 256), size=(256, 256)):

    volume_centre = volume.to_scanner((np.array(volume.shape[::-1]) - 1) / 2)

    II, JJ = make_index_coords(shape, size)

    slice_coords = II[..., np.newaxis] * volume.sagittal_vector + JJ[..., np.newaxis] * volume.coronal_vector + volume_centre
    slice_coords[:, :, 2] = z

    slice_array = scipy.interpolate.interpn([np.arange(dim) for dim in volume.shape],
                                            volume,
                                            volume.from_scanner(slice_coords)[..., ::-1],
                                            "linear",
                                            bounds_error=False,
                                            fill_value=0)

    return make_canvas(slice_array)


def get_coronal_slice(volume, y=0, shape=(256, 256), size=(256, 256)):

    volume_centre = volume.to_scanner((np.array(volume.shape[::-1]) - 1) / 2)

    II, JJ = make_index_coords(shape, size)

    slice_coords = II[..., np.newaxis] * volume.sagittal_vector + JJ[..., np.newaxis] * volume.axial_vector + volume_centre
    slice_coords[:, :, 1] = y

    slice_array = scipy.interpolate.interpn([np.arange(dim) for dim in volume.shape],
                                            volume,
                                            volume.from_scanner(slice_coords)[..., ::-1],
                                            "linear",
                                            bounds_error=False,
                                            fill_value=0)

    return make_canvas(slice_array)
