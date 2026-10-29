"""
Exporting growing wheat shoots as fspm-TF files: a stand of plants, one frame every 10 days, the format the other FSPM
platforms (CPlantBox, GroIMP, VirtualPlantLab) and ParaView read.

    python export_frames.py [--days D] [--plants N]      frames in ./outputs/frames: frame_<step>.glb, scene.fspm.json

The scene is Adel's (openalea.adel.development.shoot_scene: ScenarioShoot plants sown in rows, daily steps at 15 °C,
as openalea.adel.simulation.simulate runs them). Its recorder, given a frame_period, writes the whole scene every 10
steps: every plant, axis, organ and element with its variables (declared units, kinds), and the leaves' and stems'
triangles as a mesh. The sequence is then checked with fspm-TF's validator.

Open outputs/frames/scene.fspm.json in ParaView with fspm-TF's reader plugin (fspm-TF/paraview/fspm_tf_reader.py) to
play the days; import_frames.py reads the files back.
"""
import argparse
import os

import numpy as np

from openalea.fspm_tf.validate import validate_sequence
from openalea.metafspm.scene.scene import SceneRecorder

from openalea.adel.development import ScenarioShoot, random_delays, shoot_scene, sown_stand
from openalea.adel.simulation import DAY

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUTS = os.path.join(HERE, "outputs")
FRAME_PERIOD = 10                       # days


def export(days: int = 80, plants: int = 3, temperature: float = 15., outputs: str = OUTPUTS) -> str:
    """Grow a stand of *plants* plants for *days* days at *temperature* °C, a frame every FRAME_PERIOD days; the
    sequence's index path."""
    built = shoot_scene(sown_stand(plants, seed=0), model=ScenarioShoot, time_step=DAY,
                        delays=random_delays(plants, seed=0), start=0.)
    built.recorder = SceneRecorder(outputs, frame_period=FRAME_PERIOD)   # the scene's outputs, with its frames
    ds = built.populations[0].data_structure
    n_plants = len(ds.entity_ids("Plant"))
    built.run()                                                         # day 0: the thermal time set to its start
    for _ in range(days):
        ds.set("thermal_time_step", np.full(n_plants, temperature))     # °Cd per day, above a base of 0 °C
        built.run()
    return os.path.join(outputs, "frames", "scene.fspm.json")


def main():
    parser = argparse.ArgumentParser(description="Adel's wheat shoots as fspm-TF frames, one every 10 days")
    parser.add_argument("--days", type=int, default=80)
    parser.add_argument("--plants", type=int, default=3)
    arguments = parser.parse_args()
    index = export(arguments.days, arguments.plants)
    report = validate_sequence(index)
    frames = sorted(name for name in os.listdir(os.path.dirname(index)) if name.endswith(".glb"))
    size = sum(os.path.getsize(os.path.join(os.path.dirname(index), name)) for name in frames)
    print(f"{len(frames)} frames ({size / 1e6:.1f} MB) and their index: {index}")
    print(f"fspm-TF validator: {len(report.errors)} error(s), {len(report.warnings)} warning(s)")


if __name__ == "__main__":
    main()
