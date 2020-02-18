#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <map>

#include <SFML/Graphics.hpp>
#include "triangulation.h"


bool almost_equal(const double x, const double y, int ulp = 2)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * static_cast<double>(ulp)
        // unless the result is subnormal
        || std::abs(x - y) < std::numeric_limits<double>::min();
}


struct Vector2
{
    Vector2() = default;
    Vector2(const Vector2 &v) = default;
    Vector2(Vector2&&) = default;
    Vector2(const double vx, const double vy, const double vz, const int y, const double z_inter);

    double dist2(const Vector2 &v) const;
    double dist(const Vector2 &v) const;
    double norm2() const;

    Vector2 &operator=(const Vector2&) = default;
    Vector2 &operator=(Vector2&&) = default;
    bool operator ==(const Vector2 &v) const;
    bool operator <(const Vector2 &v) const;

    friend std::ostream &operator <<(std::ostream &str, const Vector2 &v);

    double x;
    double y;
public:
    double z;
    int year;
    double z_interpolation;
};


Vector2::Vector2(const double vx, const double vy, double vz = 0, int y = 2016, double z_inter = 0) :
    x(vx), y(vy), z(vz), year(y), z_interpolation(z_inter)
{}

double
Vector2::dist2(const Vector2 &v) const
{
    const double dx = x - v.x;
    const double dy = y - v.y;
    return dx * dx + dy * dy;
}

double
Vector2::dist(const Vector2 &v) const
{
    return hypot(x - v.x, y - v.y);
}

double
Vector2::norm2() const
{
    return x * x + y * y;
}

bool
Vector2::operator ==(const Vector2 &v) const
{
    return (this->x == v.x) && (this->y == v.y);
}

bool
Vector2::operator <(const Vector2 &v) const
{
    return std::make_pair(x, y) < std::make_pair(v.x, v.y);
}


std::ostream &
operator <<(std::ostream &str, const Vector2 &v)
{
    str.precision(20);
    return str << "Point x: " << v.x << " y: " << v.y << " z: " << v.z << " year: " << v.year;
}

bool almost_equal(const Vector2 &v1, const Vector2 &v2, int ulp = 2)
{
    return almost_equal(v1.x, v2.x, ulp) && almost_equal(v1.y, v2.y, ulp);
}

bool cmp(Vector2 a, Vector2 b) {
    return a.x < b.x || a.x == b.x && a.y < b.y;
}

bool cw(Vector2 a, Vector2 b, Vector2 c) {
    return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) < 0;
}

bool ccw(Vector2 a, Vector2 b, Vector2 c) {
    return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) > 0;
}

void convex_hull(std::vector<Vector2> &a) {
    if (a.size() == 1)  return;
    sort(a.begin(), a.end(), &cmp);
    Vector2 p1 = a[0], p2 = a.back();
    std::vector<Vector2> up, down;
    up.push_back(p1);
    down.push_back(p1);
    for (size_t i = 1; i < a.size(); ++i) {
        if (i == a.size() - 1 || cw(p1, a[i], p2)) {
            while (up.size() >= 2 && !cw(up[up.size() - 2], up[up.size() - 1], a[i]))
                up.pop_back();
            up.push_back(a[i]);
        }
        if (i == a.size() - 1 || ccw(p1, a[i], p2)) {
            while (down.size() >= 2 && !ccw(down[down.size() - 2], down[down.size() - 1], a[i]))
                down.pop_back();
            down.push_back(a[i]);
        }
    }
    a.clear();
    for (size_t i = 0; i < up.size(); ++i)
        a.push_back(up[i]);
    for (size_t i = down.size() - 2; i > 0; --i)
        a.push_back(down[i]);
}


struct ang {
    double a, b;
};

bool operator < (const ang & p, const ang & q) {
    if (p.b == 0 && q.b == 0)
        return p.a < q.a;
    return p.a * 1ll * q.b < p.b * 1ll * q.a;
}

long long sq(Vector2 & a, Vector2 & b, Vector2 & c) {
    return a.x * 1ll * (b.y - c.y) + b.x * 1ll * (c.y - a.y) + c.x * 1ll * (a.y - b.y);
}

bool point_in_polygon(int n, std::vector<Vector2> p, Vector2 q) {
    int zero_id = 0;
    for (int i = 0; i < n; ++i) {
        if (p[i].x < p[zero_id].x || p[i].x == p[zero_id].x && p[i].y < p[zero_id].y)
            zero_id = i;
    }
    Vector2 zero = p[zero_id];
    rotate(p.begin(), p.begin() + zero_id, p.end());
    p.erase(p.begin());
    --n;
    std::vector<ang> a(n);
    for (int i = 0; i < n; ++i) {
        a[i].a = p[i].y - zero.y;
        a[i].b = p[i].x - zero.x;
        if (a[i].a == 0)
            a[i].b = a[i].b < 0 ? -1 : 1;
    }

    bool in = false;
    if (q.x >= zero.x)
        if (q.x == zero.x && q.y == zero.y)
            in = true;
        else {
            ang my = { q.y - zero.y, q.x - zero.x };
            if (my.a == 0)
                my.b = my.b < 0 ? -1 : 1;
            std::vector<ang>::iterator it = upper_bound(a.begin(), a.end(), my);
            if (it == a.end() && my.a == a[n - 1].a && my.b == a[n - 1].b)
                it = a.end() - 1;
            if (it != a.end() && it != a.begin()) {
                int p1 = int(it - a.begin());
                if (sq(p[p1], p[p1 - 1], q) <= 0)
                    in = true;
            }
        }
        return in;
}


struct Edge
{
    using VertexType = Vector2;

    Edge() = default;
    Edge(const Edge&) = default;
    Edge(Edge&&) = default;
    Edge(const VertexType &v1, const VertexType &v2);

    Edge &operator=(const Edge&) = default;
    Edge &operator=(Edge&&) = default;
    bool operator ==(const Edge &e) const;

    friend std::ostream &operator <<(std::ostream &str, const Edge &e);

    const VertexType *v;
    const VertexType *w;

    double get_y(double x) const {
        if (v->x == w->x) return std::min(v->y, w->y);
        return v->y + (w->y - v->y) * (x - v->x) / (w->x - v->x);
    }
    bool isBad = false;
};


Edge::Edge(const VertexType &v1, const VertexType &v2) :
    v(&v1), w(&v2)
{}

bool
Edge::operator ==(const Edge &e) const
{
    return (*(this->v) == *e.v && *(this->w) == *e.w) ||
        (*(this->v) == *e.w && *(this->w) == *e.v);
}

std::ostream&
operator <<(std::ostream &str, const Edge &e)
{
    return str << "Edge " << *e.v << ", " << *e.w;
}

bool almost_equal(const Edge &e1, const Edge &e2)
{
    return	(almost_equal(*e1.v, *e2.v) && almost_equal(*e1.w, *e2.w)) ||
        (almost_equal(*e1.v, *e2.w) && almost_equal(*e1.w, *e2.v));
}

bool operator< (const Edge & a, const Edge & b) {
    double x;
    if (a.v->x == a.w->x) {
        x = a.v->x;
    }
    else if (b.v->x == b.w->x) {
        x = b.v->x;
    }
    else {
        x = (std::max(std::min(a.v->x, a.w->x), std::min(b.v->x, b.w->x)) +
            std::min(std::max(a.v->x, a.w->x), std::max(b.v->x, b.w->x))) / 2;
    }
    return a.get_y(x) < b.get_y(x);
}


struct Triangle
{
    using EdgeType = Edge;
    using VertexType = Vector2;

    Triangle() = default;
    Triangle(const Triangle&) = default;
    Triangle(Triangle&&) = default;
    Triangle(const VertexType &v1, const VertexType &v2, const VertexType &v3);

    bool containsVertex(const VertexType &v) const;
    bool circumCircleContains(const VertexType &v) const;

    Triangle &operator=(const Triangle&) = default;
    Triangle &operator=(Triangle&&) = default;
    bool operator ==(const Triangle &t) const;
    friend std::ostream &operator <<(std::ostream &str, const Triangle &t);

    const VertexType *a;
    const VertexType *b;
    const VertexType *c;
    bool isBad = false;
};


Triangle::Triangle(const VertexType &v1, const VertexType &v2, const VertexType &v3) :
    a(&v1), b(&v2), c(&v3), isBad(false)
{}

bool
Triangle::containsVertex(const VertexType &v) const
{
    // return p1 == v || p2 == v || p3 == v;
    return almost_equal(*a, v) || almost_equal(*b, v) || almost_equal(*c, v);
}

bool
Triangle::circumCircleContains(const VertexType &v) const
{
    const double ab = a->norm2();
    const double cd = b->norm2();
    const double ef = c->norm2();

    const double ax = a->x;
    const double ay = a->y;
    const double bx = b->x;
    const double by = b->y;
    const double cx = c->x;
    const double cy = c->y;

    const double circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
    const double circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));

    const VertexType circum((circum_x / 2), (circum_y / 2));
    const double circum_radius = a->dist2(circum);
    const double dist = v.dist2(circum);
    return dist <= circum_radius;
}

bool
Triangle::operator ==(const Triangle &t) const
{
    return	(*this->a == *t.a || *this->a == *t.b || *this->a == *t.c) &&
        (*this->b == *t.a || *this->b == *t.b || *this->b == *t.c) &&
        (*this->c == *t.a || *this->c == *t.b || *this->c == *t.c);
}

std::ostream&
operator <<(std::ostream &str, const Triangle &t)
{
    return str << "Triangle:" << "\n\t" <<
        *t.a << "\n\t" <<
        *t.b << "\n\t" <<
        *t.c << '\n';
}

bool almost_equal(const Triangle &t1, const Triangle &t2)
{
    return	(almost_equal(*t1.a, *t2.a) || almost_equal(*t1.a, *t2.b) || almost_equal(*t1.a, *t2.c)) &&
        (almost_equal(*t1.b, *t2.a) || almost_equal(*t1.b, *t2.b) || almost_equal(*t1.b, *t2.c)) &&
        (almost_equal(*t1.c, *t2.a) || almost_equal(*t1.c, *t2.b) || almost_equal(*t1.c, *t2.c));
}

class Delaunay
{
public:
    using TriangleType = Triangle;
    using EdgeType = Edge;
    using VertexType = Vector2;

    Delaunay() = default;

    const std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices);
    const std::vector<TriangleType>& getTriangles() const;
    const std::vector<EdgeType>& getEdges() const;
    const std::vector<VertexType>& getVertices() const;

    

private:
    std::vector<TriangleType> _triangles;
    std::vector<EdgeType> _edges;
    std::vector<VertexType> _vertices;
};

const std::vector<Delaunay::TriangleType>&
Delaunay::triangulate(std::vector<VertexType> &vertices)
{
    std::sort(vertices.begin(), vertices.end());
    std::vector<double> points;
    for (const auto &i : vertices) {
        points.push_back(i.x);
        //std::cout << points.back() << " ";
        points.push_back(i.y);
        //std::cout << points.back() << " ";
    }
    //std::cout << std::endl;

    std::cout << "Num vert " << vertices.size() << std::endl;

    delaunator::Delaunator d(points);

    std::cout << "Num trian " << d.triangles.size() << std::endl;

    _vertices = vertices;
    const auto &t = d.triangles;
    for (int i = 0; i < d.triangles.size(); i += 3) {
        const Vector2 &p1 = *std::lower_bound(vertices.begin(), vertices.end(), Vector2(d.coords[2 * d.triangles[i]], d.coords[2 * d.triangles[i] + 1]));
        const Vector2 &p2 = *std::lower_bound(vertices.begin(), vertices.end(), Vector2(d.coords[2 * d.triangles[i + 1]], d.coords[2 * d.triangles[i + 1] + 1]));
        const Vector2 &p3 = *std::lower_bound(vertices.begin(), vertices.end(), Vector2(d.coords[2 * d.triangles[i + 2]], d.coords[2 * d.triangles[i + 2] + 1]));


        _triangles.push_back({p1, p2, p3});
        _edges.push_back({ p1, p2 });
        _edges.push_back({ p2, p3 });
        _edges.push_back({ p3, p1 });
    }

    return _triangles;
}

const std::vector<Delaunay::TriangleType>&
Delaunay::getTriangles() const
{
    return _triangles;
}

const std::vector<Delaunay::EdgeType>&
Delaunay::getEdges() const
{
    return _edges;
}

const std::vector<Delaunay::VertexType>&
Delaunay::getVertices() const
{
    return _vertices;
}

struct event {
    double x;
    int tp, id;
    event() { }
    event(double x, int tp, int id)
        : x(x), tp(tp), id(id) { }

    bool operator< (const event & e) const {
        if (abs(x - e.x) > 1e-9)  return x < e.x;
        return tp < e.tp;
    }

};

std::ostream& operator<<(std::ostream &str, const event &e)
{
    return str << e.x << " " << e.tp << " " << e.id << std::endl;
}

std::set<Edge> s;
std::vector < std::set<Edge>::iterator > where;

inline std::set<Edge>::iterator prev(std::set<Edge>::iterator it) {
    return it == s.begin() ? s.end() : --it;
}

inline std::set<Edge>::iterator next(std::set<Edge>::iterator it) {
    return ++it;
}
bool point_in_edge(const Vector2 &p, const Edge &e)
{
    double x1 = p.x, y1 = p.y, x2 = e.v->x, y2 = e.v->y, x3 = e.w->x, y3 = e.w->y;

    if ((x1 == x2 && y1 == y2) || (x1 == x3 && y1 == y3)) {
        return true;
    }

    double xa, xb, ya, yb;
    xa = x1 - x2;
    ya = y1 - y2;
    xb = x3 - x2;
    yb = y3 - y2;

    double xa1, xb1, ya1, yb1;
    xa1 = x1 - x3;
    ya1 = y1 - y3;
    xb1 = x2 - x3;
    yb1 = y2 - y3;

    double pr = xa*yb - xb*ya;
    double pr1 = xa*xb + ya*yb;
    double pr2 = xa1*yb1 - xb1*ya1;
    double pr22 = xa1*xb1 + ya1*yb1;

    if ((abs(pr) <= 1e-9) && (pr1 >= 0) && (abs(pr2) <= 1e-9) && (pr22 >= 0))
        return true;
    return false;
}



std::vector < std::pair<Vector2, std::pair<Edge, Edge>>> scanline(Delaunay &triangulation, std::vector<Vector2> pnt) {


    std::vector<Edge> a1 = triangulation.getEdges(), a;
    using namespace std;
    cout << a1.size() << " ";


    auto it = a1.begin();
    while (it != a1.end()) {
        auto edge = *it;
        if (edge.v->x != edge.w->x) {
            a.push_back(edge);
        }
            it++;
    }

    cout << a.size() << " \n";
    std::vector < std::pair<Vector2, std::pair<Edge, Edge>>> ans;

    int n = a.size();

    std::vector<event> e;
    for (int i = 0; i < n; ++i) {
        e.push_back(event(std::min(a[i].v->x, a[i].w->x), +1, i));
        e.push_back(event(std::max(a[i].v->x, a[i].w->x), -1, i));
        //cout << a[i] << " id= " << i << endl;
    }
    for (int i = 0; i < pnt.size(); i++) {
        e.push_back(event(pnt[i].x, 2, n + i));
        a.push_back(Edge(pnt[i], pnt[i]));
    }
    sort(e.begin(), e.end());


    s.clear();
    where.resize(a.size());


    for (size_t i = 0; i < e.size(); ++i) {

        int id = e[i].id;
        if (e[i].tp == +1) {
            std::set<Edge>::iterator nxt = s.lower_bound(a[id]);
            if (nxt == s.end() || !(*nxt == a[id])) {
                if (prev(nxt) == s.end() || (prev(nxt) != s.end() && !(*prev(nxt) == a[id]))) {
                    where[id] = s.insert(nxt, a[id]);
                }
            }

        }
        else if (e[i].tp == -1) {

            auto tmp = s.find(a[id]);
            if (!(tmp == s.end() || !(*tmp == a[id]))) {
                s.erase(tmp);
            }

        }
        else {
            std::set<Edge>::iterator nxt = s.lower_bound(a[id]),
                prv = prev(nxt);


            if (nxt == s.end()) {
                auto nt = s.upper_bound(a[id]);
                if (nt != s.end()) {
                    ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*nt, *nt)));
                }
                continue;

            }

            if (point_in_edge(*a[id].v, *nxt)) {
                ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*nxt, *nxt)));
                continue;
            }

            if (prv == s.end()) {
                ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*nxt, *nxt)));
                continue;
            }
            if (point_in_edge(*a[id].v, *prv)) {
                ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*prv, *prv)));
                continue;
            }

            if (nxt != s.end() && prv != s.end()) {
                if (*(nxt->v) == *(prv->v) || *(nxt->v) == *(prv->w) || *(nxt->w) == *(prv->v) || *(nxt->w) == *(prv->w))

                    ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*nxt, *prv)));
                else {
                    //std::cout << "!!!\n";
                }
                continue;
            }

        }
    }
    return ans;
}

double interpolation(const Edge &a, const Edge &b, const Vector2 &p) {
    Vector2 v1 = *a.v, v2 = *a.w, v3 = *b.v, v4 = *b.w;
    if (v1 == v2) {
        return v1.z;
    }
    if (v3 == v4) {
        return v3.z;
    }
    Vector2 p1, p2, p0;
    if (v1 == v3 && v2 == v4 || v1 == v4 && v2 == v3) {
        p1 = v1;
        p2 = v2;
        if (v1.x != v2.x) {
            return p1.z + (p2.z - p1.z) * (p.x - p1.x) / (p2.x - p1.x);
        }
        else {
            return p1.z + (p2.z - p1.z) * (p.y - p1.y) / (p2.y - p1.y);
        }
    }

    if (v1 == v3) {
        p1 = v2;
        p2 = v4;
        p0 = v1;
    } else
    if (v1 == v4) {
        p1 = v2;
        p2 = v3;
        p0 = v1;
    } else
    if (v2 == v3) {
        p1 = v1;
        p2 = v4;
        p0 = v2;
    } else{
        p1 = v1;
        p2 = v3;
        p0 = v2;
    }
    double a1 = p1.x - p0.x, a2 = p2.x - p0.x;
    double b1 = p1.y - p0.y, b2 = p2.y - p0.y;
    double c1 = p1.z - p0.z, c2 = p2.z - p0.z;
    if (a1*b2 - a2 * b1 == 0) {
        double d1 = p1.x - p.x;
        std::cout << v1.x - v3.x << " " << v1.x - v4.x << " " << v2.x - v3.x << " " << v2.x - v4.x << " " << d1 <<  " !";
    }
    return p0.z + ((p.y - p0.y) * (a1 * c2 - a2 * c1) + (p.x - p0.x)*(c1*b2 - c2 * b1)) / (a1*b2 - a2 * b1);
}


double dist_point_to_plane(Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3) {
    double x0 = p0.x, x1 = p1.x, x2 = p2.x, x3 = p3.x;
    double y0 = p0.y, y1 = p1.y, y2 = p2.y, y3 = p3.y;
    double z0 = p0.z, z1 = p1.z, z2 = p2.z, z3 = p3.z;
    double a = y1*(z2 - z3) - y2*(z1 - z3) + y3*(z1 - z2);
    double b = -x1*(z2 - z3) + x2 * (z1 - z3) - x3*(z1 - z2);
    double c = y1*(-x2 + x3) - y2*(-x1 + x3) + y3*(-x1 + x2);
    double d = y1*(-z2*x3 + x2*z3) - y2*(-z1*x3 + x1*z3) + y3 * (-z1*x2 + x1*z2);
    if (a == 0) {
        return abs(z0 - z1);
    }
    return abs(a*x0 + b*y0 + c*z0 + d) / sqrt(a*a + b*b + c*c);
}


double tria_area(Vector2 a, Vector2 b, Vector2 c) {
    double x1 = a.x - b.x, y1 = a.y - b.y;
    double x2 = c.x - b.x, y2 = c.y - b.y;
    return (abs(x1*y2 - x2*y1) / 2) * (40075000 / 360 * sqrt(2) / 2) * (20004274 / 180);
}


Vector2 seg_intersect(Vector2 c_16, Vector2 b_16, Vector2 c_18, Vector2 b_18)
{
    //std::cout << c_16 << b_16;
    double x1 = c_18.x, x2 = b_18.x, y1 = c_18.y, y2 = b_18.y, z1 = c_18.z, z2 = b_18.z;
    double x3 = c_16.x, x4 = b_16.x, y3 = c_16.y, y4 = b_16.y, z3 = c_16.z, z4 = b_16.z;


    if (x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0) {
        double a1 = z2 - z1, b1 = x1 - x2, c1 = -x1*z2 + z1*x2;
        double a2 = z4 - z3, b2 = x3 - x4, c2 = -x3*z4 + z3*x4;


        double x = (b1*c2 - b2*c1) / (a1 * b2 - a2*b1);
        double z = (c1*a2 - c2*a1) / (a1 * b2 - a2*b1);

        double y;
        if (z2 != z1) {
            y = (z - z1) * (b_18.y - c_18.y) / (z2 - z1) + c_18.y;
        }
        else {
            y = (x - x1) * (b_18.y - c_18.y) / (x2 - x1) + c_18.y;
        }
        Vector2 ans = Vector2(x, y, z);
        return ans;
    }
    double a1 = z2 - z1, b1 = y1 - y2, c1 = -y1*z2 + z1*y2;
    double a2 = z4 - z3, b2 = y3 - y4, c2 = -y3*z4 + z3*y4;


    double y = (b1*c2 - b2*c1) / (a1 * b2 - a2*b1);
    double z = (c1*a2 - c2*a1) / (a1 * b2 - a2*b1);

    double x;
    if (z2 != z1) {
        x = (z - z1) * (b_18.x - c_18.x) / (z2 - z1) + c_18.x;
    }
    else {
        x = (y - y1) * (b_18.x - c_18.x) / (y2 - y1) + c_18.x;
    }
    Vector2 ans = Vector2(x, y, z);
    return ans;

}





double vol_diff(Triangle tria) {
    auto p1 = *tria.a, p2 = *tria.b, p3 = *tria.c;

    //std::cout << p1 << p2 << p3;
    Vector2 a_16, b_16, c_16, a_18, b_18, c_18;
    if (p1.year == 2016) {
        a_16 = Vector2(p1.x, p1.y, p1.z);
        a_18 = Vector2(p1.x, p1.y, p1.z_interpolation);
    }
    else {
        a_18 = Vector2(p1.x, p1.y, p1.z);
        a_16 = Vector2(p1.x, p1.y, p1.z_interpolation);
    }
    if (p2.year == 2016) {
        b_16 = Vector2(p2.x, p2.y, p2.z);
        b_18 = Vector2(p2.x, p2.y, p2.z_interpolation);
    }
    else {
        b_18 = Vector2(p2.x, p2.y, p2.z);
        b_16 = Vector2(p2.x, p2.y, p2.z_interpolation);
    }
    if (p3.year == 2016) {
        c_16 = Vector2(p3.x, p3.y, p3.z);
        c_18 = Vector2(p3.x, p3.y, p3.z_interpolation);
    }
    else {
        c_18 = Vector2(p3.x, p3.y, p3.z);
        c_16 = Vector2(p3.x, p3.y, p3.z_interpolation);
    }
    double a = a_18.z - a_16.z, b = b_18.z - b_16.z, c = c_18.z - c_16.z;

    //std::cout << " " << a << " " << b << " " << c << std::endl;


    if ((a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0)) {
        //std::cout << " fir ";
        double v = (abs(a) + abs(b) + abs(c)) / 3 * tria_area(p1, p2, p3);

        return v;
    }
    Vector2 m, n;
    bool f = false;
    if ((a > 0 && c > 0 && b < 0)) {
        f = true;
    }

    if ((a < 0 && c < 0 && b > 0)) {
        std::swap(a_18, a_16);
        std::swap(b_18, b_16);
        std::swap(c_18, c_16);
        f = true;
    }

    if ((b > 0 && c > 0 && a < 0)) {
        std::swap(b_18, a_18);
        std::swap(b_16, a_16);
        f = true;
    }

    if ((b < 0 && c < 0 && a > 0)) {
        std::swap(b_18, a_18);
        std::swap(b_16, a_16);
        std::swap(a_18, a_16);
        std::swap(b_18, b_16);
        std::swap(c_18, c_16);
        f = true;

    }

    if ((a > 0 && b > 0 && c < 0)) {
        std::swap(b_18, c_18);
        std::swap(b_16, c_16);
        f = true;

    }

    if ((a < 0 && b < 0 && c > 0)) {
        std::swap(b_18, c_18);
        std::swap(b_16, c_16);
        std::swap(a_18, a_16);
        std::swap(b_18, b_16);
        std::swap(c_18, c_16);
        f = true;
    }

    if (f) {
        //std::cout << " sec ";
        m = Vector2((a_16.x * abs(a) + b_16.x * abs(b)) / (abs(a) + abs(b)),
            (a_16.y * abs(a) + b_16.y * abs(b)) / (abs(a) + abs(b)),
            (a_16.z * abs(a) + b_16.z * abs(b)) / (abs(a) + abs(b)));
        n = Vector2((b_16.x * abs(b) + c_16.x * abs(c)) / (abs(b) + abs(c)),
            (b_16.y * abs(b) + c_16.y * abs(c)) / (abs(b) + abs(c)),
            (b_16.z * abs(b) + c_16.z * abs(c)) / (abs(b) + abs(c)));

        double v;
        double h_b_mnb1 = dist_point_to_plane(b_16, m, n, b_18);
        double s_mnb1 = tria_area(m, n, b_18);
        double h_m_aa1c = dist_point_to_plane(m, a_16, a_18, c_16);
        double s_aa1c = tria_area(a_16, a_18, c_16);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double h_m_ncc1 = dist_point_to_plane(m, n, c_16, c_18);
        double s_ncc1 = tria_area(n, c_16, c_18);

        v = (h_b_mnb1*s_mnb1 + h_m_aa1c*s_aa1c + h_m_ac1c*s_ac1c + h_m_ncc1*s_ncc1) / 3;
        return v;
    }


    a = -a, b = -b, c = -c;

    if (a == 0 && b < 0 && c > 0) {
        //std::cout << "a1\n";
        m = seg_intersect(c_16, b_16, c_18, b_18);
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;

    }


    if (a == 0 && b > 0 && c < 0) {
        //std::cout << "a2\n";
        std::swap(b_16, c_16);
        std::swap(b_18, c_18);

        std::swap(b_18, b_16);
        std::swap(c_18, c_16);


        m = seg_intersect(c_16, b_16, c_18, b_18);
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;

    }

    if (b == 0 && a < 0 && c > 0) {
        //std::cout << "b1\n";
        std::swap(a_16, b_16);
        std::swap(a_18, b_18);

        m = seg_intersect(c_16, b_16, c_18, b_18);
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;
    }

    if (b == 0 && c < 0 && a > 0) {
        //std::cout << "b1\n";
        std::swap(a_16, b_16);
        std::swap(a_18, b_18);

        std::swap(a_18, a_16);
        std::swap(c_18, c_16);

        m = seg_intersect(c_16, b_16, c_18, b_18);
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;
    }


    if (c == 0 && b < 0 && a >0) {
        //std::cout << "c\n";
        std::swap(a_16, c_16);
        std::swap(a_18, c_18);

        //std::cout << a_16 << a_18 << b_16 << b_18 << c_16 << c_18 << "\n";


        m = seg_intersect(c_16, b_16, c_18, b_18);
        //std::cout << m;
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;
    }

    if (c == 0 && b > 0 && a < 0) {
        //std::cout << "c2\n";

        std::swap(a_16, c_16);
        std::swap(a_18, c_18);

        std::swap(b_18, b_16);
        std::swap(a_18, a_16);

        m = seg_intersect(c_16, b_16, c_18, b_18);
        double h_b_amb1 = dist_point_to_plane(b_16, a_16, m, b_18);
        double s_amb1 = tria_area(a_16, m, b_18);
        double h_m_ac1c = dist_point_to_plane(m, a_16, c_18, c_16);
        double s_ac1c = tria_area(a_16, c_18, c_16);
        double v = (h_b_amb1*s_amb1 + h_m_ac1c*s_ac1c) / 3;
        return v;
    }

    return 0;
}



Delaunay triangulation(std::vector<Vector2> &points_16_in_conv, sf::RenderWindow &window, sf::Color cl)
{
    Delaunay triangulation_16;
    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<Triangle> triangles_16 = triangulation_16.triangulate(points_16_in_conv);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << triangles_16.size() << " triangles generated in " << diff.count() << "s\n";

    const std::vector<Edge> edges_16 = triangulation_16.getEdges();


    std::vector<std::array<sf::Vertex, 2> > lines;

    for (int i = 0; i < edges_16.size() - 1; i++) {
        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(edges_16[i].v->x + 2.),
            static_cast<float>(edges_16[i].v->y + 2.));

        line[0].color = cl;
        line[1].position = sf::Vector2f(
            static_cast<float>(edges_16[i].w->x + 2.),
            static_cast<float>(edges_16[i].w->y + 2.));
        line[1].color = cl;
        window.draw(line);

    }

    return triangulation_16;
}






int main(int argc, char * argv[])
{
    FILE *f;
    double x, y, z;
    std::vector<Vector2> v_16, v_18;
    double t7, t2, t3, t4, t5, t6;
    double min_x_16 = 10000, min_x_18 = 10000, min_y_16 = 10000, min_y_18 = 1000;
    double max_x_16 = 0, max_x_18 = 0, max_y_16 = 0, max_y_18 = 0;

    f = freopen("2018.txt", "r", stdin);
    while (fscanf(f, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &x, &y, &z, &t7, &t2, &t3, &t4, &t5, &t6) == 9) {
        v_18.push_back(Vector2(x, y, z, 2018));
        if (x < min_x_18)
            min_x_18 = x;
        if (x > max_x_18)
            max_x_18 = x;
        if (y < min_y_18)
            min_y_18 = y;
        if (y > max_y_18)
            max_y_18 = y;
    }
    std::cout << v_18.size() << std::endl;
    fclose(f);

    f = freopen("2016.txt", "r", stdin);
    while (fscanf(f, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &x, &y, &z, &t7, &t2, &t3, &t4, &t5, &t6) == 9) {
        v_16.push_back(Vector2(x, y, z, 2016));
        if (x < min_x_16)
            min_x_16 = x;
        if (x > max_x_16)
            max_x_16 = x;
        if (y < min_y_16)
            min_y_16 = y;
        if (y > max_y_16)
            max_y_16 = y;
    }

    std::vector<Vector2> v_16_delta, v_18_delta;

    double min_x = std::max(min_x_16, min_x_18), min_y = std::max(min_y_16, min_y_18);
    double max_x = std::min(max_x_16, max_x_18), max_y = std::min(max_y_16, max_y_18);
    double delta = sqrt(((max_x - min_x) *(max_y - min_y)) * 500 / 6206696);

    double a_x = 800 / (max_x - min_x), b_x = (-800 * min_x) / (max_x - min_x);
    double a_y = 600 / (max_y - min_y), b_y = (-600 * min_y) / (max_y - min_y);


    std::cout << "jjfjgng " << a_x << " " << b_x << "\n";
    x = min_x + (max_x - min_x) / 2, y = min_y + (max_y - min_y) / 2;

    for (auto i : v_16) {
        if (i.x >= x - delta && i.x <= x + delta  && i.y >= y - delta && i.y <= y + delta) {
            v_16_delta.push_back(i);
        }
    }

    for (auto i : v_18) {
        if (i.x >= x - delta && i.x <= x + delta  && i.y >= y - delta && i.y <= y + delta) {
            v_18_delta.push_back(i);
        }
    }

    std::cout << v_16.size() << "\n";
    std::cout << min_x_16 << " " << max_x_16 << " " << min_y_16 << " " << max_y_16 << std::endl;
    std::cout << min_x_18 << " " << max_x_18 << " " << min_y_18 << " " << max_y_18 << std::endl;

    std::cout << delta << std::endl;

    std::cout << v_16_delta.size() << " " << v_18_delta.size() << "\n";

    fclose(f);

    std::vector<Vector2> points_16 = v_16_delta;

    // SFML window
    sf::RenderWindow window(sf::VideoMode(850, 650), "Delaunay triangulation");
    window.setFramerateLimit(1);

    std::vector<Vector2> points_18 = v_18_delta;

    //int numberPoints = 200, numberPoints_18 = 200;

    //std::default_random_engine eng(std::random_device{}());
    //std::uniform_real_distribution<double> dist_w(0, 800);
    //std::uniform_real_distribution<double> dist_h(0, 600);
    //std::uniform_real_distribution<double> dist_z(0, 5);

    //std::cout << "Generating " << numberPoints << " " << numberPoints_18 << " random points" << std::endl;

    //std::vector<Vector2> points_16;
    //for (int i = 0; i < numberPoints; ++i) {
    //    points_16.push_back(Vector2{ dist_w(eng), dist_h(eng), dist_z(eng), 2016 });
    //}

    //// SFML window
    //sf::RenderWindow window(sf::VideoMode(850, 650), "Delaunay triangulation");
    //window.setFramerateLimit(1);
    //// Transform each points of each vector as a rectangle
    //for (const auto p : points_16) {
    //    sf::RectangleShape s{ sf::Vector2f(1, 1) };
    //    s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
    //    s.setFillColor(sf::Color::Cyan);
    //    window.draw(s);
    //}

    //std::vector<Vector2> points_18;
    //for (int i = 0; i < numberPoints_18; ++i) {
    //    points_18.push_back(Vector2{ dist_w(eng), dist_h(eng), dist_z(eng), 2018 });
    //}

    //// Transform each points of each vector as a rectangle
    //for (const auto p : points_18) {
    //    sf::RectangleShape s{ sf::Vector2f(1, 1) };
    //    s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
    //    s.setFillColor(sf::Color::Red);
    //    window.draw(s);
    //}





    //=============== Построение выпуклых оболочек ================


    std::vector<Vector2> conv_16(points_16.size());
    std::copy(points_16.begin(), points_16.end(), conv_16.begin());

    convex_hull(conv_16);

    std::reverse(conv_16.begin(), conv_16.end());




    std::vector<Vector2> conv_18(points_18.size());

    std::copy(points_18.begin(), points_18.end(), conv_18.begin());

    convex_hull(conv_18);


    std::reverse(conv_18.begin(), conv_18.end());


    // ============================== ОТБОР ТОЧЕК ====================================



    std::vector<Vector2> points_18_in_conv;


    for (auto p : points_18) {
        if (point_in_polygon(conv_16.size(), conv_16, p)) {
            points_18_in_conv.push_back(p);
        }
    }


    for (const auto p : points_18_in_conv) {
        sf::RectangleShape s{ sf::Vector2f(2, 2) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Red);
        window.draw(s);
    }

    std::vector<Vector2> points_16_in_conv;

    for (auto p : points_16) {
        if (point_in_polygon(conv_18.size(), conv_18, p)) {
            points_16_in_conv.push_back(p);
        }
    }
    for (const auto p : points_16_in_conv) {
        sf::RectangleShape s{ sf::Vector2f(2, 2) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Yellow);
        window.draw(s);
    }

    std::cout << " There are " << points_18_in_conv.size() << " points of " << points_18.size() << " in 2018\n";
    std::cout << " There are " << points_16_in_conv.size() << "  points of " << points_16.size() << "in 2016 \n";






    // =========================== ТРИАНГУЛЯЦИИ ==================================








    Delaunay triangulation_16((triangulation(points_16_in_conv, window, sf::Color::Cyan)));
    Delaunay triangulation_18((triangulation(points_18_in_conv, window, sf::Color::Yellow)));



    

    conv_18.push_back(conv_18[0]);
    conv_16.push_back(conv_16[0]);

    std::vector <Vector2> t;

    for (int i = 0; i < int(points_18_in_conv.size()); i++) {
        t.push_back(points_18_in_conv[i]);
    }
    auto point_18_in_trian_16_tmp = scanline(triangulation_16, t);
    std::vector <Vector2> point_18_in_trian_16;

    for (int i = 0; i < point_18_in_trian_16_tmp.size(); i++) {
        auto pt = point_18_in_trian_16_tmp[i].first;
        auto res = point_18_in_trian_16_tmp[i].second;

        double z = interpolation(res.first, res.second, pt);
        if (z < 0 || z > 100) {
            std::cout << abs(z) << std::endl;
        }

        point_18_in_trian_16_tmp[i].first.z_interpolation = z;

        point_18_in_trian_16.push_back(point_18_in_trian_16_tmp[i].first);
    }


    std::vector <Vector2> t1;


    for (int i = 0; i < int(points_16_in_conv.size()); i++) {
        t1.push_back(points_16_in_conv[i]);
    }


    auto point_16_in_trian_18_tmp = scanline(triangulation_18, t1);
    std::vector <Vector2> point_16_in_trian_18;

    for (int i = 0; i < point_16_in_trian_18_tmp.size(); i++) {
        auto pt = point_16_in_trian_18_tmp[i].first;
        auto res = point_16_in_trian_18_tmp[i].second;

        double z = interpolation(res.first, res.second, pt);
        point_16_in_trian_18_tmp[i].first.z_interpolation = z;

        if (z < 0 || z > 100) {
            std::cout << abs(z) << std::endl;
        }

        /* if ((abs(z) >= 25) || (abs(z) <= 23)) {
             std::cout << abs(z) << std::endl;
             std::cout << res.first << " " << res.second << " " << pt << "\n";

             std::cout << (*res.first.v == *res.second.v) << " " << (*res.first.v == *res.second.w) << " "
                 << (*res.first.w == *res.second.v) << " " << (*res.first.w == *res.second.w) << "\n";

            // std ::cout << point

         } */


        point_16_in_trian_18.push_back(point_16_in_trian_18_tmp[i].first);

        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(res.first.v->x + 2.),
            static_cast<float>(res.first.v->y + 2.));

        line[0].color = sf::Color::White;
        line[1].position = sf::Vector2f(
            static_cast<float>(res.first.w->x + 2.),
            static_cast<float>(res.first.w->y + 2.));
        line[1].color = sf::Color::White;
        window.draw(line);

        line[0].position = sf::Vector2f(
            static_cast<float>(res.second.v->x + 2.),
            static_cast<float>(res.second.v->y + 2.));

        line[0].color = sf::Color::White;
        line[1].position = sf::Vector2f(
            static_cast<float>(res.second.w->x + 2.),
            static_cast<float>(res.second.w->y + 2.));
        line[1].color = sf::Color::White;
        window.draw(line);

        sf::RectangleShape s{ sf::Vector2f(3, 3) };
        s.setPosition(static_cast<float>(pt.x), static_cast<float>(pt.y));
        s.setFillColor(sf::Color::White);
        window.draw(s);
    }


    // ================================ ОБЩАЯ ТРИАНГУЛЯЦИЯ ====================================
    std::cout << point_16_in_trian_18.size() << " " << point_18_in_trian_16.size() << " ";

    for (int i = 0; i < point_18_in_trian_16.size(); i++) {
        point_16_in_trian_18.push_back(point_18_in_trian_16[i]);
    }

    std::cout << " total_size: " << point_16_in_trian_18.size() << std::endl;

    auto points_total = point_16_in_trian_18;



    Delaunay triangulation_total(triangulation(points_total, window, sf::Color::Red));

    auto ar = triangulation_total.getTriangles();
    double vol_dif = 0, area = 0;
    for (auto v : ar) {


        double tmp = vol_diff(v);

        /*if (tria_area(*v.a, *v.b, *v.c) < tmp) {
            std::cout << tria_area(*v.a, *v.b, *v.c) << " " << tmp << std::endl;
            std::cout << v.a->z << " " << v.a->z_interpolation << "\n";
            std::cout << v.b->z << " " << v.b->z_interpolation << "\n";
            std::cout << v.c->z << " " << v.c->z_interpolation << "\n";

        }

        */

        vol_dif += tmp;
        area += tria_area(*v.a, *v.b, *v.c);
        //std::cout << "area "  << tria_area(*v.a, *v.b, *v.c) << " vol " << tmp << "\n";
    }
    std::cout << vol_dif << " " << area;
    window.display();

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
    }

    Vector2 v1(-3, 0, -2), v2(-4, -3, 1), v3(2, 3, 5), v4(4, 3, 10);
    std::cout << seg_intersect(v1, v2, v3, v4);


    return 0;

}
