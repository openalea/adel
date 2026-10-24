def groupOutput_mais(caribu_dict, key="EiSup", sum_=False):
    # add a unique id per phytomer
    caribu_dict["id"], caribu_dict["av"] = [], []
    for i in range(len(caribu_dict[key])):
        opt = str(int(caribu_dict["Opt"][i]))
        phyto = str(int(caribu_dict["Elt"][i]))
        phyto = (4 - len(phyto)) * "0" + phyto
        # shoot = str(int(caribu_dict['Opak'][i]))
        # shoot = (3-len(shoot))*'0'+shoot
        pl = str(int(caribu_dict["Plt"][i]))
        pl = (3 - len(pl)) * "0" + pl
        lab = opt + pl + phyto  # label mais
        caribu_dict["id"].append(lab)

    # compute average for elements of key with the same id
    i = 0
    av = {}
    id = caribu_dict["id"][i]
    while i < len(caribu_dict["id"]):
        v = []
        while id == caribu_dict["id"][i]:
            v.append(caribu_dict[key][i])
            i = i + 1
            if i == len(caribu_dict["id"]):
                break

        if sum_ is False:
            av[id] = mean(v)
        else:
            av[id] = sum(v)

        if i == len(caribu_dict["id"]):
            break
        else:
            id = caribu_dict["id"][i]

    # add the average into the av column
    for i in range(len(caribu_dict[key])):
        caribu_dict["av"].append(av[caribu_dict["id"][i]])

    return caribu_dict["av"], av  # caribu_dict#['output, ']
