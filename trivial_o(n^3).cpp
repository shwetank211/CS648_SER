#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Random.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <bits/stdc++.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Min_sphere_of_points_d_traits_2<K, double> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_circle;
typedef K::Point_2 Point;

typedef CGAL::Exact_circular_kernel_2 Circular_k;
typedef CGAL::Point_2<Circular_k> Point_2;
typedef CGAL::Circle_2<Circular_k> Circle_2;
typedef CGAL::Line_2<Circular_k> Line_2;
typedef CGAL::Vector_2<Circular_k> Vector_2;
typedef Circular_k::FT FT;

double calculate_angle(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    Vector_2 v1(p1, p2);
    Vector_2 v2(p3, p2);
    double s = (v1 * v2).to_double();
    double cosine = (s) / (std::sqrt(v1.squared_length().to_double()) * std::sqrt(v2.squared_length().to_double()));
    double angle = std::acos(cosine) * 180 * M_1_PI;
    return angle;
}

bool check(int n, const Point_2 *P, const Circle_2 &circle)
{
    for (int i = 0; i < n; i++)
    {
        if (circle.bounded_side(P[i]) == CGAL::ON_UNBOUNDED_SIDE)
        {
            return false;
        }
    }
    return true;
}

int main()
{
    const int n = 100;
    Point P[n];
    Point_2 P2[n];
    CGAL::Random r; // random number generator

    for (int i = 0; i < n; ++i)
    {
        P[i] = Point(r.get_double(), r.get_double());
        P2[i] = Point_2(P[i].x(), P[i].y());
    }

    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE;

    double radius = MAXFLOAT;
    Point_2 center;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            Circle_2 circle = Circle_2(P2[i], P2[j], ori);
            if (check(n, P2, circle))
            {
                if (circle.squared_radius().to_double() < radius * radius)
                {
                    radius = std::sqrt(circle.squared_radius().to_double());
                    center = circle.center();
                }
            }
            else
            {
                Point_2 p;
                double angle = 180.00;
                for (int k = 0; k < n; k++)
                {
                    if (k == i || k == j)
                        continue;
                    double angle1 = calculate_angle(P2[i], P2[k], P2[j]);
                    if (angle1 < angle)
                    {
                        angle = angle1;
                        p = P2[k];
                    }
                }
                if (angle < 90.00)
                {
                    circle = Circle_2(P2[i], P2[j], p);
                    if (circle.squared_radius().to_double() < radius * radius && check(n,P2,circle))
                    {
                        radius = std::sqrt(circle.squared_radius().to_double());
                        center = circle.center();
                    }
                }
            }
        }
    }
    std::cout << "Center: " << center.x().to_double() << " " << center.y().to_double() << std::endl;
    std::cout << "Radius: " << radius << std::endl;


    Min_circle mc(P, P + n);

    Min_circle::Cartesian_const_iterator ccib = mc.center_cartesian_begin(), ccie = mc.center_cartesian_end();
    std::cout << "center:";
    for (; ccib != ccie; ++ccib)
    {
        std::cout << " " << *ccib;
    }
    std::cout << std::endl
              << "radius: " << mc.radius() << std::endl;
    return 0;
}