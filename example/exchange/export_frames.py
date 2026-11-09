"""
Exporting growing wheat shoots as encoding3d files: a stand of plants, one frame every 10 days, the format the other FSPM
platforms (CPlantBox, GroIMP, VirtualPlantLab) and ParaView read.

    python export_frames.py [--days D] [--plants N]      frames in ./outputs/frames.fspm.zip: frame_<step>.glb, ...

The scene is Adel's (openalea.adel.development.shoot_scene: ScenarioShoot plants sown in rows, daily steps at 15 °C,
as openalea.adel.simulation.simulate runs them). Its recorder, given a frame_period, writes the whole scene every 10
days: every plant, axis, organ and element with its variables (declared units, kinds), and the leaves' and stems'
triangles as a mesh. At the end the frames folder is packed into frames.fspm.zip (the archive checked, then the
folder removed), and the sequence is checked with encoding3d's validator.

A second, filtered export shows logging: the first plant only (with all its axes, organs and elements) and its
lengths and areas, in ./outputs/log/frames.fspm.zip (outputs only: such frames are not reimported).

Open outputs/frames.fspm.zip in ParaView with encoding3d's reader plugin (encoding3d/paraview/plugin_encoding3d.py) to
play the days; import_frames.py reads the files back.
"""
import argparse
import os

from openalea.encoding3d import read_index
from openalea.encoding3d.validate import validate_sequence
from openalea.metafspm.scene.recorder import SceneRecorder

from openalea.adel.development import ScenarioShoot, random_delays, shoot_scene
from openalea.metafspm.scene.population import stand
from openalea.adel.simulation import DAY

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUTS = os.path.join(HERE, "outputs")
FRAME_PERIOD = "10 d"
# a log of one plant (by its name in the planting table, with everything it holds), its elements' lengths and areas
LOG = dict(entities=["p0"], variables=["length", "area", "green_area"])


def export(days: int = 80, plants: int = 3, temperature: float = 15., outputs: str = OUTPUTS,
           frame_options: dict = None) -> str:
    """Grow a stand of *plants* plants for *days* days at *temperature* °C, a frame every FRAME_PERIOD days (filtered
    by *frame_options*, e.g. LOG), packed at the end; the sequence's archive."""
    recorder = SceneRecorder(outputs, summaries=False, frame_period=FRAME_PERIOD,   # the scene's frames only
                             frame_options=frame_options, pack_at_finish=True)
    built = shoot_scene(stand(plants, rows=4, plant_spacing=0.04, jitter=0.005, seed=0), model=ScenarioShoot, time_step=DAY,
                        delays=random_delays(plants, seed=0), start=0., recorder=recorder,
                        forcings={"air_temperature": temperature})        # °C: its degree-days above 0 °C
    built.simulate(days + 1)                                             # day 0: the thermal time set to its start
    return os.path.join(outputs, "frames.fspm.zip")


def main():
    parser = argparse.ArgumentParser(description="Adel's wheat shoots as encoding3d frames, one every 10 days")
    parser.add_argument("--days", type=int, default=80)
    parser.add_argument("--plants", type=int, default=3)
    arguments = parser.parse_args()
    index = export(arguments.days, arguments.plants)
    report = validate_sequence(index)
    frames = read_index(index)["frames"]
    print(f"{len(frames)} frames, their parts and their index packed ({archive_size(index):.1f} MB): {index}")
    print(f"encoding3d validator: {len(report.errors)} error(s), {len(report.warnings)} warning(s)")
    log = export(arguments.days, arguments.plants, outputs=os.path.join(OUTPUTS, "log"), frame_options=LOG)
    report = validate_sequence(log)
    print(f"logged (plant {LOG['entities'][0]}, {len(LOG['variables'])} variables): {archive_size(log):.1f} "
          f"MB, {len(report.errors)} error(s)")


def archive_size(path: str) -> float:
    """The size of the archive *path* (MB)."""
    return os.path.getsize(path) / 1e6


if __name__ == "__main__":
    main()
