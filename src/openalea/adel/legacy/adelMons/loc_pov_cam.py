def loc_pov_cam(
    num_rang=10, num_plante=2, interrang=0.8, interplante=0.125, z_camera=2.0
):
    """define camera position from stand caracteristics"""
    # camera_position = None;
    # write the node code here.
    pos = [z_camera, num_rang * interrang + 0.0, num_plante * interplante]
    # return outputs
    return pos
