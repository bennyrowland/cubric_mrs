import pyx


def file_table(files):
    table_runner = pyx.text.LatexRunner()
    table_runner.preamble("\\usepackage{tabularx}")

    file_strings = ["{} & : & {}".format(label, filename.replace("_", "\\_")) for label, filename in files.items()]
    table_string = """
\\begin{{tabularx}}{{19cm}}{{l c X}}
{}
\\end{{tabularx}}
""".format("\\\\\n".join(file_strings))
    print(table_string)
    return table_runner.text(0, 0, table_string, [pyx.text.valign.bottom])