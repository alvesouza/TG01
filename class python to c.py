from ctypes import *
import main

class a:
    def __idiv__(self, other):
        return 0

if __name__ == '__main__':
    v = main.Vector2D(2)
    h = main.Vector2D(2, 2)
    #h = h + v
    v = main.Vector3D(0, 10, 0)
    q = main.Quaternion(main.Vector3D(1, 0, 0), 30)
    h = main.Quaternion(main.Vector3D(1, 0, 0), 60)
    q = h.Multiply(q)
    q = q.Multiply(main.Quaternion(main.Vector3D(1, 0, 0), -90))
    v_linha = q.Inverse().Get_Rotated_Vector((q.Get_Rotated_Vector(v)))
    print(q.x, q.y, q.z, q.w)
    print(v_linha.x, v_linha.y, v_linha.z)