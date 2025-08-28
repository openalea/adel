def Sum_EnergyDensities(energy1, area1, energy2, area2):
    # write the node code here.
    energy = (energy1[0]) * (area1[0])
    area = area1[0]
    for i in range(1, len(energy1)):
        energy += energy1[i] * area1[i]
        area += area1[i]
    for j in range(1, len(energy2)):
        energy += energy2[j] * area2[j]
        area += area2[j]
    energy_density = energy / area
    # return outputs
    return energy_density
