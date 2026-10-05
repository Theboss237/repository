#!/usr/bin/env python3
"""
Calcul des coordonnées absolues (dans le repère monde) d'un point défini
dans un repère parent, pour Stäubli Robotics Suite / VAL3.

Chaîne de repères :
    world[0]  ->  fCell[0] (frame)  ->  p30Nok[0] (pointRx)

Pose absolue :  T_world_point = T_world_fCell · T_fCell_point

Convention Stäubli (VAL3) pour (x, y, z, rx, ry, rz) :
    - translations en mm, angles en degrés
    - rotations successives autour des axes MOBILES : X puis Y' puis Z''
      =>  R = Rx(rx) · Ry(ry) · Rz(rz)
"""

import math


# --------------------------------------------------------------------------
# Petite algèbre matricielle (sans dépendance externe)
# --------------------------------------------------------------------------
def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def inverse_homogene(T):
    """Inverse d'une matrice homogène : [R t]^-1 = [R^T  -R^T·t]."""
    Rt = [[T[j][i] for j in range(3)] for i in range(3)]
    t = [-sum(Rt[i][k] * T[k][3] for k in range(3)) for i in range(3)]
    return [Rt[0] + [t[0]], Rt[1] + [t[1]], Rt[2] + [t[2]], [0.0, 0.0, 0.0, 1.0]]


# --------------------------------------------------------------------------
# Conversions pose <-> matrice homogène 4x4
# --------------------------------------------------------------------------
def rot_x(a):
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0], [0, c, -s], [0, s, c]]


def rot_y(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]


def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]


def pose_vers_matrice(x, y, z, rx, ry, rz):
    """Pose Stäubli (mm, degrés) -> matrice homogène 4x4."""
    rx, ry, rz = map(math.radians, (rx, ry, rz))
    R = mat_mul(mat_mul(rot_x(rx), rot_y(ry)), rot_z(rz))
    return [R[0] + [x], R[1] + [y], R[2] + [z], [0.0, 0.0, 0.0, 1.0]]


def matrice_vers_pose(T):
    """Matrice homogène 4x4 -> pose Stäubli (x, y, z, rx, ry, rz) en mm / degrés.

    R = Rx·Ry·Rz donne :
        R[0,2] =  sin(ry)
        R[1,2] = -sin(rx)·cos(ry)    R[2,2] = cos(rx)·cos(ry)
        R[0,1] = -cos(ry)·sin(rz)    R[0,0] = cos(ry)·cos(rz)
    """
    R = T
    sy = max(-1.0, min(1.0, R[0][2]))
    ry = math.asin(sy)
    if abs(sy) < 1.0 - 1e-12:
        rx = math.atan2(-R[1][2], R[2][2])
        rz = math.atan2(-R[0][1], R[0][0])
    else:
        # Blocage de cardan (ry = ±90°) : on fixe rz = 0
        rz = 0.0
        rx = math.atan2(R[2][1], R[1][1])
    x, y, z = T[0][3], T[1][3], T[2][3]
    return (x, y, z,
            math.degrees(rx), math.degrees(ry), math.degrees(rz))


# --------------------------------------------------------------------------
# Calcul principal
# --------------------------------------------------------------------------
def pose_absolue(pose_parent, pose_locale):
    """Pose du point dans le monde, à partir de la pose du repère parent
    (exprimée dans le monde) et de la pose du point (exprimée dans le parent)."""
    T = mat_mul(pose_vers_matrice(*pose_parent), pose_vers_matrice(*pose_locale))
    return matrice_vers_pose(T)


def pose_relative(pose_parent, pose_absolue_pt):
    """Opération inverse : pose dans le parent d'un point connu dans le monde."""
    T = mat_mul(inverse_homogene(pose_vers_matrice(*pose_parent)),
                pose_vers_matrice(*pose_absolue_pt))
    return matrice_vers_pose(T)


if __name__ == "__main__":
    # Screen 2 : <Data name="fCell" type="frame">  fatherId="world[0]"
    fCell = (362.0420374, -49.909946, -265.899982,
             -0.0000135538, -0.00000434751, 90.0000199962)

    # Screen 1 : <Data name="p30Nok" type="pointRx">  fatherId="fCell[0]"
    p30Nok = (169.837325, 700.0334453, 372.398311,
              -179.999244, -0.000734286, 88.845311)

    absolu = pose_absolue(fCell, p30Nok)

    noms = ("x", "y", "z", "rx", "ry", "rz")
    print("Coordonnées absolues (repère world) à appliquer au bloc :")
    for n, v in zip(noms, absolu):
        unite = "mm" if n in ("x", "y", "z") else "deg"
        print(f"  {n:>2} = {v:14.6f} {unite}")

    # Vérification : retour dans le repère fCell
    retour = pose_relative(fCell, absolu)
    err = max(abs(a - b) for a, b in zip(retour, p30Nok))
    print(f"\nVérification aller-retour (écart max) : {err:.2e}")
