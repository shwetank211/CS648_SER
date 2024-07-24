#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Random.h>
#include <bits/stdc++.h>


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <vector>

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

void swap(Point_2 *a, Point_2 *b)
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

void randomize(int n, Point_2 *P)
{
    // Use a different seed value so that
    // we don't get same result each time
    // we run this program
    srand(time(NULL));

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
    //     N = int(tmp[0])
    // int Z = 10000;
    std::cout << "iterations     10        20        30        50          100" << std::endl;
    for (int i = 1; i <= 5; i += 1)
    {
        int n = pow(10, i);
        Point P[n];
        Point_2 P2[n];
        CGAL::Random r; // random number generator

        double time = 0;
        double count = 0;
        // double total_time_2 = 0;
        // double total_time_3 = 0;

        int num_iterations = 200;
        std::vector<double> total_time_1(num_iterations, 0.0);

        for (int iter = 0; iter < num_iterations; ++iter)
        {
            for (int i = 0; i < n; ++i)
            {
                P[i] = Point(r.get_double(), r.get_double());
                P2[i] = Point_2(P[i].x(), P[i].y());
            }

            CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE;

            auto start_time_1 = std::chrono::high_resolution_clock::now();

            // Your algorithm code here... Expected O(n)
            Circle_2 circle = Circle_2(P2[0], P2[1], ori);

            for (int i = 2; i < n; i++)
            {
                if (circle.bounded_side(P2[i]) == CGAL::ON_UNBOUNDED_SIDE)
                {
                    count++;
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
                                double angle1 = calculate_angle(P2[i], P2[k], P2[j]);
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
            auto end_time_1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end_time_1 - start_time_1;
            time += duration.count();
            total_time_1[iter] = duration.count();


            // auto start_time_2 = std::chrono::high_resolution_clock::now();
            // // O(n^3)
            // double radius = MAXFLOAT;
            // Point_2 center;
            // for (int i = 0; i < n; i++)
            // {
            //     for (int j = i + 1; j < n; j++)
            //     {
            //         Circle_2 circle = Circle_2(P2[i], P2[j], ori);
            //         if (check(n, P2, circle))
            //         {
            //             if (circle.squared_radius().to_double() < radius * radius)
            //             {
            //                 radius = std::sqrt(circle.squared_radius().to_double());
            //                 center = circle.center();
            //             }
            //         }
            //         else
            //         {
            //             Point_2 p;
            //             double angle = 180.00;
            //             for (int k = 0; k < n; k++)
            //             {
            //                 if (k == i || k == j)
            //                     continue;
            //                 double angle1 = calculate_angle(P2[i], P2[k], P2[j]);
            //                 if (angle1 < angle)
            //                 {
            //                     angle = angle1;
            //                     p = P2[k];
            //                 }
            //             }
            //             if (angle < 90.00)
            //             {
            //                 circle = Circle_2(P2[i], P2[j], p);
            //                 if (circle.squared_radius().to_double() < radius * radius && check(n, P2, circle))
            //                 {
            //                     radius = std::sqrt(circle.squared_radius().to_double());
            //                     center = circle.center();
            //                 }
            //             }
            //         }
            //     }
            // }

            // auto end_time_2 = std::chrono::high_resolution_clock::now();
            // duration = end_time_2 - start_time_2;
            // // worst_time = std::max(worst_time, duration.count());
            // total_time_2 += duration.count();

            // // O(n^4)
            // auto start_time_3 = std::chrono::high_resolution_clock::now();
            // radius = MAXFLOAT;
            // for (int i = 0; i < n; i++)
            // {
            //     for (int j = i + 1; j < n; j++)
            //     {
            //         Circle_2 circle = Circle_2(P2[i], P2[j], ori);
            //         if (check(n, P2, circle))
            //         {
            //             if (circle.squared_radius().to_double() < radius * radius)
            //             {
            //                 radius = std::sqrt(circle.squared_radius().to_double());
            //                 center = circle.center();
            //             }
            //         }
            //     }
            // }
            // // checking all 3-point circles
            // for (int i = 0; i < n; i++)
            // {
            //     for (int j = i + 1; j < n; j++)
            //     {
            //         for (int k = j + 1; k < n; k++)
            //         {
            //             Circle_2 circle = Circle_2(P2[i], P2[j], P2[k]);
            //             if (check(n, P2, circle))
            //             {
            //                 if (circle.squared_radius().to_double() < radius * radius)
            //                 {
            //                     radius = std::sqrt(circle.squared_radius().to_double());
            //                     center = circle.center();
            //                 }
            //             }
            //         }
            //     }
            // }
            // auto end_time_3 = std::chrono::high_resolution_clock::now();
            // duration = end_time_3 - start_time_3;
            // // worst_time = std::max(worst_time, duration.count());
            // total_time_3 += duration.count();
        }

        double avg_time = time / num_iterations;
        double avg_count = count / num_iterations;
        int count_10=0,count_20=0,count_50=0,count_100=0,count_30=0;
        for(int i=0;i<num_iterations;i++)
        {
            if(total_time_1[i]>1.1*avg_time)
                count_10++;
            if(total_time_1[i]>1.2*avg_time)
                count_20++;
            if(total_time_1[i]>1.5*avg_time)
                count_50++;
            if(total_time_1[i]>2*avg_time)
                count_100++;
            if(total_time_1[i]>1.3*avg_time)
                count_30++;
        }
        // double out_avg = out_total / num_iterations;
        // std::cout << n << "           \t " << avg_time << " \t\t" << worst_time << "\t\t " << avg_count << std::endl;
        std::cout << n << " \t\t" << count_10 << "\t\t " << count_20 << "\t\t" << count_30 <<"\t\t " << count_50 << "\t\t " << count_100 << std::endl;
    }

    return 0;
}
