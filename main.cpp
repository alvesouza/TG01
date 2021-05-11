#include <iostream>
#include "boost/python.hpp"
#include "float.h"

using namespace boost::python;
/*
g++ -o2 -I/boost-1_73_0 -I/usr/include/python3.8/ main.cpp -fPIC -lboost_python38 -Wl,-soname,greet -shared -o hello.so
*/
typedef struct Vector2D Vector2D;
struct Vector2D{
    double x,y;

    //Constructor
    Vector2D(double x = 0.0, double y = 0.0):x(x), y(y){}

    //Operators overloading
    inline Vector2D operator=(const Vector2D& a)const{return {a.x,a.y};}
    inline Vector2D operator+(const Vector2D& a)const{return {a.x+x, a.y+y};}
    inline Vector2D operator+=(const Vector2D& a)const{ return *this + a;}
    inline Vector2D operator-=(const Vector2D& a)const{ return *this - a;}
    inline Vector2D operator-(const Vector2D& a) const{return {x-a.x, y-a.y};}
    inline bool operator==(const Vector2D& a)const{return (a.x==x && a.y==y);}
    inline bool operator!=(const Vector2D& a)const{return (a.x!=x || a.y!=y);}
    inline Vector2D operator*(const double a)const{return {a*x, y*a};}
    inline Vector2D operator*=(const double a)const{ return *this * a;}
    inline Vector2D operator/=(const double a)const{ return *this / a;}
    inline Vector2D operator/(const double a)const{return {x/a, y/a};}

    //Methods
    inline double Abs() const {return sqrt(x*x + y*y);}
    inline double Dot_Prod(const Vector2D& a) const{ return x*a.x + y*a.y;}
    inline double Angle_Between_In_Radian(const Vector2D& a) const {return acos(Dot_Prod(a)/(sqrt(Dot_Prod(*this)*(a.Dot_Prod(a)))));}


};

void Module_Add_Vector2D(){
    class_<Vector2D>("Vector2D", init<double, double>((
                        arg("x") = 0.0f,
                        arg("y") = 0.0f
                    )
                )
            )
            .def_readwrite("x", &Vector2D::x)
            .def_readwrite("y", &Vector2D::y)
            .def("Abs", &Vector2D::Abs)
            .def("Dot_Prod", &Vector2D::Dot_Prod)
            .def("Angle_Between_In_Radian", &Vector2D::Angle_Between_In_Radian)
            .def("__add__", &Vector2D::operator+)
            .def("__iadd__", &Vector2D::operator+=)
            .def("__sub__", &Vector2D::operator-)
            .def("__isub__", &Vector2D::operator-=)
            .def("__eq__", &Vector2D::operator==)
            .def("__ne__", &Vector2D::operator!=)
            .def("__mul__", &Vector2D::operator*)
            .def("__rmul__", &Vector2D::operator*)
            .def("__imul__", &Vector2D::operator*=)
            .def("__truediv__", &Vector2D::operator/)
            .def("__rtruediv__", &Vector2D::operator/)
            .def("__idiv__", &Vector2D::operator/=);
}

void Module_Add_2DLogic(){
    Module_Add_Vector2D();
}

typedef struct Vector3D Vector3D;
struct Vector3D{
    double x,y,z;

    //Constructor
    Vector3D(double x = 0.0, double y = 0.0, double z = 0.0):x(x), y(y),z(z){}

    //Operators overloading
    inline Vector3D operator=(const Vector3D& a)const{return {a.x,a.y,a.z};}
    inline Vector3D operator+(const Vector3D& a)const{return {a.x+x, a.y+y, a.z+z};}
    inline Vector3D operator+=(const Vector3D& a)const{ return *this + a;}
    inline Vector3D operator-=(const Vector3D& a)const{ return *this - a;}
    inline Vector3D operator-(const Vector3D& a) const{return {x-a.x, y-a.y, z-a.z};}
    inline bool operator==(const Vector3D& a)const{return (a.x==x && a.y==y && a.z==z);}
    inline bool operator!=(const Vector3D& a)const{return (a.x!=x || a.y!=y || a.z!=z);}
    inline Vector3D operator*(const double a)const{return {a*x, y*a, z*a};}
    inline Vector3D operator*=(const double a)const{ return *this * a;}
    inline Vector3D operator/=(const double a)const{ return *this / a;}
    inline Vector3D operator/(const double a)const{return {x/a, y/a, z/a};}

    //Methods
    inline double Abs() const {return sqrt(x*x + y*y + z*z);}
    inline double Dot_Prod(const Vector3D& a) const{ return x*a.x + y*a.y + z*a.z;}
    inline Vector3D Cross_Prod(const Vector3D& a) const{ return {y*a.z-z*a.y,
                                                                 z*a.x-x*a.z,
                                                                 x*a.y-y*a.x};}
    inline double Angle_Between_In_Radian(const Vector3D& a) const {return acos(Dot_Prod(a)/(sqrt(Dot_Prod(*this)*(a.Dot_Prod(a)))));}


};

void Module_Add_Vector3D(){
    class_<Vector3D>("Vector3D", init<double, double, double>((
                                                              arg("x") = 0.0f,
                                                              arg("y") = 0.0f,
                                                              arg("z") = 0.0f
                                                      )
                     )
    )
            .def_readwrite("x", &Vector3D::x)
            .def_readwrite("y", &Vector3D::y)
            .def_readwrite("z", &Vector3D::z)
            .def("Abs", &Vector3D::Abs)
            .def("Dot_Prod", &Vector3D::Dot_Prod)
            .def("Cross_Prod", &Vector3D::Cross_Prod)
            .def("Angle_Between_In_Radian", &Vector3D::Angle_Between_In_Radian)
            .def("__add__", &Vector3D::operator+)
            .def("__iadd__", &Vector3D::operator+=)
            .def("__sub__", &Vector3D::operator-)
            .def("__isub__", &Vector3D::operator-=)
            .def("__eq__", &Vector3D::operator==)
            .def("__ne__", &Vector3D::operator!=)
            .def("__mul__", &Vector3D::operator*)
            .def("__rmul__", &Vector3D::operator*)
            .def("__imul__", &Vector3D::operator*=)
            .def("__truediv__", &Vector3D::operator/)
            .def("__rtruediv__", &Vector3D::operator/)
            .def("__idiv__", &Vector3D::operator/=);
}

typedef struct Quaternion Quaternion;
struct Quaternion{
    double x,y,z,w;
    Quaternion(const Vector3D& axis, double angle){
        const double radians_2 = angle*(double)M_PI / 360;//Radians divided by two
        const Vector3D norm = axis/axis.Abs();
        w = cos(radians_2);
        x = norm.x * sin(radians_2);
        y = norm.y * sin(radians_2);
        z = norm.z * sin(radians_2);
        Normalize();
    }
    Quaternion(const double w = 1.0, const double x = 0.0, const double y = 0.0, const double z = 0.0):x(x),y(y),z(z),w(w){}
    Quaternion(const double w, const Vector3D& v):x(v.x),y(v.y),z(v.z),w(w){}

    //Operators overloading
    inline Quaternion operator=(const Quaternion& a)const{return {a.w,a.x,a.y,a.z};}
    inline bool operator==(const Quaternion& a)const{return (a.x==x && a.y==y && a.z==z && a.w==w);}
    inline bool operator!=(const Quaternion& a)const{return (a.x!=x || a.y!=y || a.z!=z || a.w != w);}

    //Methods

    inline Vector3D Get_Vector()const{return {x, y, z};}
    inline Quaternion Inverse()const{return {w, Vector3D(-x,-y,-z)};}
    inline Quaternion Multiply(const Quaternion& a)const{
        const Vector3D vector1 = Get_Vector();
        const Vector3D vector2 = a.Get_Vector();
        return Quaternion(
                w*a.w - vector1.Dot_Prod(vector2),
                vector2*w + vector1*a.w + vector1.Cross_Prod(vector2)
                );
    }

    inline Vector3D Get_Rotated_Vector(const Vector3D& v)const{
        //q * v * q-1
        return (
                (
                        Multiply(
                                Quaternion(0, v)
                                )
                                ).Multiply(
                                        Inverse()
                                        )
                                        ).Get_Vector();
    }
    inline double Abs()const{
        return sqrt(x*x+y*y+z*z+w*w);
    }
    void Normalize(){
        const double abs = Abs();
        x /= abs;
        y /= abs;
        z /= abs;
        w /= abs;
    }

};

void Module_Add_Quaternion(){
    class_<Quaternion, boost::shared_ptr<Quaternion>>("Quaternion", init<Vector3D&, double>())
            .def(init< double, Vector3D&>())
            .def(init<double, double, double, double>(
                    (
                            arg("w") = 1.0f,
                            arg("x") = 0.0f,
                            arg("y") = 0.0f,
                            arg("z") = 0.0f
                        )
                    )
            )
            .def_readwrite("x", &Quaternion::x)
            .def_readwrite("y", &Quaternion::y)
            .def_readwrite("z", &Quaternion::z)
            .def_readwrite("w", &Quaternion::w)
            .def("Abs", &Quaternion::Abs)
            .def("Normalize", &Quaternion::Normalize)
            .def("Multiply", &Quaternion::Multiply)
            .def("Inverse", &Quaternion::Inverse)
            .def("Get_Vector", &Quaternion::Get_Vector)
            .def("Get_Rotated_Vector", &Quaternion::Get_Rotated_Vector)
            .def("__eq__", &Quaternion::operator==)
            .def("__ne__", &Quaternion::operator!=);
}

void Module_Add_3DLogic(){
    Module_Add_Vector3D();
    Module_Add_Quaternion();
}

char const* greet(){
    return "hello, world";
}

BOOST_PYTHON_MODULE(main){
    def("greet", greet);
    Module_Add_2DLogic();
    Module_Add_3DLogic();
}
