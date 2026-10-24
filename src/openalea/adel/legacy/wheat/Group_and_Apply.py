def Group_and_Apply(table, properties, function, data):
    # write the node code here.
    keys = table.keys()

    prop = []

    if type(properties) is str:
        if properties.split(" ")[0] in keys:
            prop.append(properties.split(" ")[0])

    else:
        for p in range(0, len(properties)):
            if properties[p] in keys:
                prop.append(properties[p])

    dat = []

    if type(data) is str:
        if data.split(" ")[0] in keys:
            dat.append(data.split(" ")[0])

    else:
        for p in range(0, len(data)):
            if data[p] in keys:
                dat.append(data[p])

    by = {}
    n = len(table[keys[0]])
    for i in range(n):
        key = (prop[0], table[prop[0]][i])

    if len(prop) > 1:
        for p in range(1, len(prop)):
            key = (key, (prop[p], table[prop[p]][i]))
    by.setdefault(key, []).append(i)

    l = range(0, len(table[keys[0]]))
    for k, v in by.items():
        args = {}

    for i in v:
        for d in range(0, len(dat)):
            args.setdefault(dat[d], []).append(table[dat[d]][i])

        # args1.append(table[dat[0]][i])

        # args2.append(table[dat[1]][i])
    # res =  function(k,*args1,*args2)

    res = sum(args[dat[0]]) / sum(args[dat[1]])
    for index in v:
        l[index] = res

    # return outputs
    return l
