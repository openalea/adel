def Write_Table(Table, filename, first=True, sep=" "):
    """
    returns a string representation of the content of the table
    """

    header = [sep.join(list(Table.keys()))]
    rows = list(zip(*[Table[k] for k in list(Table.keys())]))
    lines = [sep.join(map(str, r)) for r in rows]
    mode = "a"
    if first:
        mode = "w"
    fout = open(filename, mode)
    if first:
        fout.write("\n".join(header + lines))
    else:
        fout.write("\n")
        fout.write("\n".join(lines))
    fout.close()

    return (filename,)
