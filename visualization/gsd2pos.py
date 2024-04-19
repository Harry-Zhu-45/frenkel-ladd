"""
version = 2.0920.3 由王进编写
version = 2.1101 由张晔修正
可以将 GSD 文件转化为 POS 文件
目前支持 直方盒子 + 凸多面体/球体粒子 的转换
如果 GSD 文件不带形状，自行添加形状到 POS 文件中
"""
import sys

import gsd.hoomd
import numpy as np


def get_box_dimensions(traj_frame: gsd.hoomd.Frame) -> str:
    """Get box dimensions from trajectory frame.

    Args:
        traj_frame (gsd.hoomd.Frame): Trajectory frame.

    Returns:
        str: Box dimensions.
    """
    return ' '.join(map(str, traj_frame.configuration.box))


def get_particle_definition(shape: str, diameter: float = None, vertices: list = None) -> str:
    """Get particle definition for POS file.

    Args:
        shape (str): Particle shape.
        diameter (float, optional): particle diameter. Defaults to None.
        vertices (list, optional): particle vertices. Defaults to None.

    Returns:
        str: Particle definition for POS file.
    """
    if shape == 'Sphere':
        return f'sphere {diameter}'
    if shape == 'ConvexPolyhedron':
        vtks_flat = np.array(vertices).ravel()
        poly3d = f'poly3d {len(vertices)}'
        for vtks_i in range(3 * len(vertices)):
            poly3d += f' {vtks_flat[vtks_i]}'
        return poly3d

    sys.exit(f'Error: Unknown particle shape: {shape}')


def write_box_info(pos_file: str, box_dims: str) -> None:
    """Write box information to POS file.

    Args:
        pos_file (str): POS file.
        box_dims (str): Box dimensions.

    Returns:
        None
    """
    pos_file.write(f'box {box_dims}\n')


def write_particle_definition(pos_file: str, particle_type: str, definition: str) -> None:
    """Write particle definition to POS file.

    Args:
        pos_file (str): POS file.
        particle_type (str): Particle type.
        definition (str): Particle definition.

    Returns:
        None
    """
    pos_file.write(f'def {particle_type} "{definition}"\n')


def write_particle_placement(pos_file: str, traj_frame: gsd.hoomd.Frame, color_list: list) -> None:
    """Write particle placement to POS file.

    Args:
        pos_file (str): POS file.
        traj_frame (str): Trajectory frame.
        color_list (list): Color list.

    Returns:
        None
    """
    for N_i in range(traj_frame.particles.N):
        particle_type = traj_frame.particles.types[traj_frame.particles.typeid[N_i]]
        color = color_list[traj_frame.particles.typeid[N_i]]
        pos_arr = traj_frame.particles.position
        orientation = traj_frame.particles.orientation[N_i]
        pos_file.write(f'{particle_type} {color} {pos_arr[N_i][0]} {pos_arr[N_i][1]} {pos_arr[N_i][2]} ')
        if traj_frame.particles.type_shapes[traj_frame.particles.typeid[N_i]]:
            pos_file.write(f'{orientation[0]} {orientation[1]} {orientation[2]} {orientation[3]} ')
        pos_file.write('\n')


def gsd2pos(gsd_file_name: str, pos_file_name: str, frame_list: str | list = 'all'):
    """Convert GSD file to POS file.

    Args:
        gsd_file_name (str): GSD file name.
        pos_file_name (str): POS file name.
        frame_list (str | list, optional): Frame list to convert. Defaults to 'all'.

    Returns:
        None

    Example:
        gsd2pos('trajectory.gsd', 'trajectory.pos')
        gsd2pos('trajectory.gsd', 'trajectory.pos', [0:10])
    """
    color_list = ['ffff0000', 'ff0000ff', 'ff00ff00', 'ff000000', 'ff888888', 'ffff00ff']

    print('开始' + gsd_file_name + '的转换……')
    traj = gsd.hoomd.open(name=gsd_file_name, mode='r')

    with open(file=pos_file_name, mode='w', encoding='utf-8') as pos_file:
        if frame_list == "all":
            frame_list = range(len(traj))

        for frame in frame_list:
            traj_frame = traj[frame]
            box_dims = get_box_dimensions(traj_frame)
            write_box_info(pos_file, box_dims)

            # print(traj_frame.particles.types)
            # print(traj_frame.particles.type_shapes)

            for j, particle_type in enumerate(traj_frame.particles.types):
                if traj_frame.particles.type_shapes[j]:
                    definition = get_particle_definition(
                        traj_frame.particles.type_shapes[j]['type'],
                        diameter=traj_frame.particles.type_shapes[j].get('diameter'),
                        vertices=traj_frame.particles.type_shapes[j].get('vertices'))
                else:
                    print('注意！此文件可能含有未定义形状的粒子！默认设为体积为1的球形粒子！')
                    definition = 'sphere 1.128'
                write_particle_definition(pos_file, particle_type, definition)

            write_particle_placement(pos_file, traj_frame, color_list)
            pos_file.write('eof\n')

    print("转换完成！")


if __name__ == "__main__":
    # gsd2pos('fcc_initial.gsd', 'fcc_initial.pos')
    gsd2pos('trajectory.gsd', 'trajectory.pos')
