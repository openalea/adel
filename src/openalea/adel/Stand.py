"""
Class interface for stand generation
"""

from math import sqrt
from random import sample
from scipy.interpolate import interp1d

from openalea.adel.stand.stand import agronomicplot, regular_plot


class AgronomicStand:
    def __init__(
        self,
        sowing_density=10,
        plant_density=10,
        inter_row=0.8,
        noise=0,
        density_curve_data=None,
    ):
        self.sowing_density = sowing_density
        self.inter_row = inter_row
        self.plant_density = plant_density
        self.inter_plant = 1.0 / inter_row / sowing_density
        self.noise = noise
        self.density_curve_data = density_curve_data
        df = density_curve_data
        if df is None:
            self.density_curve = None
        else:
            # hs_curve = interp1d(df['HS'], df['density'])
            TT_curve = interp1d(df["TT"], df["density"])
            # self.density_curve = {'hs_curve':hs_curve,'TT_curve':TT_curve}
            self.density_curve = TT_curve

    def plot_dimensions(self, nplants=1, aspect="square"):
        if aspect == "square":
            nsown = nplants * self.sowing_density / float(self.plant_density)
            side = sqrt(1.0 / self.sowing_density * nsown)
            nrow = max(1, round(side / self.inter_row))
            plant_per_row = max(1, round(side / self.inter_plant))
            plot_length = self.inter_plant * plant_per_row
            plot_width = self.inter_row * nrow
            return plot_length, plot_width
        elif aspect == "line":
            plot_width = self.inter_row
            plot_length = (
                nplants
                * self.inter_plant
                * self.sowing_density
                / float(self.plant_density)
                if self.plant_density > 0.0
                else 0.0
            )
            return plot_length, plot_width
        else:
            return 0.5, 0.5

    def smart_stand(self, nplants=1, at=None, convunit=100):
        """return an (almost) square stand that match inter-row, current density and nplants in the stand,
        but (dynamicaly) adjusting inter-plant to solve the problem
        """

        density = self.plant_density
        if at is not None:
            if self.density_curve is not None:
                density = self.density_curve(at)

        # find a square design for sowing
        nsown = nplants * 1.0 * self.sowing_density / density
        side = sqrt(1.0 / self.sowing_density * nsown)
        nrow = int(max(1, round(side / self.inter_row)))
        plant_per_row = int(max(1, round(float(nsown) / nrow)))
        while nplants > (nrow * plant_per_row):
            plant_per_row += 1
        domain_area = nrow * self.inter_row * plant_per_row * self.inter_plant
        # adjust inter_plant spacing so that n_emerged / domain_area match plant density
        n_emerged = int(round(domain_area * density))
        # assert(n_emerged >= nplants)
        n_emerged = nplants
        target_domain_area = 1.0 * n_emerged / density
        inter_plant = target_domain_area / (plant_per_row * nrow * self.inter_row)

        positions, domain, domain_area = regular_plot(
            inter_plant,
            self.inter_row,
            nrow,
            plant_per_row,
            noise=self.noise,
            convunit=convunit,
        )

        positions = sample(positions, nplants)
        return nplants, domain, positions, domain_area

    def stand(self, nplants=1, aspect="square", convunit=100):
        length, width = self.plot_dimensions(nplants, aspect)
        n_emerged, positions, domain, domain_area, _ = agronomicplot(
            length,
            width,
            self.sowing_density,
            self.plant_density,
            self.inter_row,
            noise=self.noise,
            convunit=convunit,
        )

        return n_emerged, domain, positions, length * width

    def plot(self, positions):
        import pandas

        df = pandas.DataFrame(positions)
        df.plot(0, 1, style="o")


def agronomicStand_node(
    sowing_density=10, plant_density=10, inter_row=0.8, noise=0, density_curve_data=None
):
    return AgronomicStand(**locals())
