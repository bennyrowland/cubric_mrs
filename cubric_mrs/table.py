import pyx
import numpy as np


def file_table(files):
    table_runner = pyx.text.LatexRunner()
    table_runner.preamble("\\usepackage{tabularx}")

    file_strings = ["{} & : & {}".format(label, filename.replace("_", "\\_")) for label, filename in files.items()]
    table_string = """
\\begin{{tabularx}}{{19cm}}{{l c X}}
{}
\\end{{tabularx}}
""".format("\\\\\n".join(file_strings))
    return table_runner.text(0, 0, table_string, [pyx.text.valign.bottom])


def metabolite_table(concentrations, ref_cr, units="AU"):
    table_runner = pyx.text.LatexRunner()
    table_runner.preamble("\\usepackage{tabularx}")

    if ref_cr <= 0:
        ref_cr = np.inf

    metabolite_strings = ["{} & {} & {:.3f} & {}\\\\\n".format(
        name.replace("_", "\\_"), value["concentration"], float(value["concentration"]) / ref_cr, value["sd"]) for
        name, value in sorted(concentrations.items())]

    table_string = "\\begin{{tabularx}}{{9cm}}{{X r r r}}\nMetabolite & Conc/{} & /Cr & SD\%\\\\\n\hline\n{}\\end{{tabularx}}".format(
        units,
        "".join(metabolite_strings))

    table_latex = table_runner.text(0, 0, table_string, [pyx.text.valign.top])
    return table_latex


def quality_table(quality_info):
    table_runner = pyx.text.LatexRunner()
    table_runner.preamble("\\usepackage{tabularx}")
    # now display the fit quality metrics
    quality_strings = ["{} & {}\\\\\n".format(
        param, value) for param, value in sorted(quality_info.items())]

    quality_table_string = "\\begin{{tabularx}}{{9cm}}{{X r}}\nFit Parameters&\\\\\n\hline\n{}\\end{{tabularx}}".format(
        "".join(quality_strings))

    quality_table_latex = table_runner.text(0, 0, quality_table_string, [pyx.text.valign.top])
    return quality_table_latex
