// O(n^3) implementation
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

// A utility function to swap to integers 
void swap (Point_2 *a, Point_2 *b) 
{ 
    auto temp = *a; 
    *a = *b; 
    *b = temp; 
} 

double calculate_angle(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    Vector_2 v1(p1, p2);
    Vector_2 v2(p3, p2);
    double s = (v1 * v2).to_double();
    double cosine = (s) / (std::sqrt(v1.squared_length().to_double()) * std::sqrt(v2.squared_length().to_double()));
    double angle = std::acos(cosine) * 180 * M_1_PI;
    return angle;
}

void randomize (int n, Point_2 *P) 
{ 
    // Use a different seed value so that 
    // we don't get same result each time
    // we run this program 
    srand (time(NULL)); 
 
    // Start from the last element and swap 
    // one by one. We don't need to run for 
    // the first element that's why i > 0 
    for (int i = n - 1; i > 0; i--) 
    { 
        // Pick a random index from 0 to i 
        int j = rand() % (i + 1); 
 
        // Swap arr[i] with the element 
        // at random index 
        swap(&P[i], &P[j]); 
    } 
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
    randomize (n,P2); 

    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE;
    Circle_2 circle = Circle_2(P2[0], P2[1], ori);

    for (int i = 2; i < n; i++)
    {
        if (circle.bounded_side(P2[i]) == CGAL::ON_UNBOUNDED_SIDE)
        {
            Circle_2 circle1 = Circle_2(P2[i], P2[0], ori);

            for (int j = 1; j < i; j++)
            {
                if (circle1.bounded_side(P2[j]) == CGAL::ON_UNBOUNDED_SIDE)
                {
                    // Line_2 line = Line_2(P2[i],P2[j]);
                    Point_2 p;
                    double angle = 180.00;
                    for (int k = 0; k < j; k++)
                    {
                        // find first k such that angle subtended by P[K] on line is minimum and is acute
                        // if no point found, then return 0
                        double angle1 = calculate_angle(P2[i],P2[k],P2[j]);
                        if (angle1 < angle)
                        {
                            angle = angle1;
                            p = P2[k];
                        }
                    }
                    if (angle < 90.00)
                    {
                        circle1 = Circle_2(P2[i], P2[j], p);
                    }
                    else
                    {
                        circle1 = Circle_2(P2[i], P2[j], ori);
                    }
                }
            }
            circle = circle1;
        }
    }
    std::cout << "Center: " << circle.center().x().to_double() << " " << circle.center().y().to_double() << std::endl;
    std::cout << "Radius: " << std::sqrt(circle.squared_radius().to_double()) << std::endl;

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