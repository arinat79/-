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
#include <unordered_set>

#include <SFML/Graphics.hpp>

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
    Vector2(const double vx, const double vy, const double vz);

    double dist2(const Vector2 &v) const;
    double dist(const Vector2 &v) const;
    double norm2() const;

    Vector2 &operator=(const Vector2&) = default;
    Vector2 &operator=(Vector2&&) = default;
    bool operator ==(const Vector2 &v) const;
    friend std::ostream &operator <<(std::ostream &str, const Vector2 &v);

    double x;
    double y;
    double z;
};


Vector2::Vector2(const double vx, const double vy, double vz = 0) :
    x(vx), y(vy), z(vz)
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

std::ostream &
operator <<(std::ostream &str, const Vector2 &v)
{
    return str << "Point x: " << v.x << " y: " << v.y << " z: " << v.z;
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
    for (size_t i = 1; i<a.size(); ++i) {
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
    for (size_t i = 0; i<up.size(); ++i)
        a.push_back(up[i]);
    for (size_t i = down.size() - 2; i>0; --i)
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
    for (int i = 0; i<n; ++i) {
        if (p[i].x < p[zero_id].x || p[i].x == p[zero_id].x && p[i].y < p[zero_id].y)
            zero_id = i;
    }
    Vector2 zero = p[zero_id];
    rotate(p.begin(), p.begin() + zero_id, p.end());
    p.erase(p.begin());
    --n;
    std::vector<ang> a(n);
    for (int i = 0; i<n; ++i) {
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
        if (v->x == w->x) return v->y;
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
    Delaunay(const Delaunay&) = delete;
    Delaunay(Delaunay&&) = delete;

    const std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices);
    const std::vector<TriangleType>& getTriangles() const;
    const std::vector<EdgeType>& getEdges() const;
    const std::vector<VertexType>& getVertices() const;

    Delaunay& operator=(const Delaunay&) = delete;
    Delaunay& operator=(Delaunay&&) = delete;

private:
    std::vector<TriangleType> _triangles;
    std::vector<EdgeType> _edges;
    std::vector<VertexType> _vertices;
};

const std::vector<Delaunay::TriangleType>&
Delaunay::triangulate(std::vector<VertexType> &vertices)
{
    // Store the vertices locally
    _vertices = vertices;

    // Determinate the super triangle
    double minX = vertices[0].x;
    double minY = vertices[0].y;
    double maxX = minX;
    double maxY = minY;

    for (std::size_t i = 0; i < vertices.size(); ++i)
    {
        if (vertices[i].x < minX) minX = vertices[i].x;
        if (vertices[i].y < minY) minY = vertices[i].y;
        if (vertices[i].x > maxX) maxX = vertices[i].x;
        if (vertices[i].y > maxY) maxY = vertices[i].y;
    }

    const double dx = maxX - minX;
    const double dy = maxY - minY;
    const double deltaMax = std::max(dx, dy);
    const double midx = (minX + maxX) / 2;
    const double midy = (minY + maxY) / 2;

    const VertexType p1(midx - 20 * deltaMax, midy - deltaMax);
    const VertexType p2(midx, midy + 20 * deltaMax);
    const VertexType p3(midx + 20 * deltaMax, midy - deltaMax);

    // Create a list of triangles, and add the supertriangle in it
    _triangles.push_back(TriangleType(p1, p2, p3));

    for (auto p = begin(vertices); p != end(vertices); p++)
    {
        std::vector<EdgeType> polygon;

        for (auto & t : _triangles)
        {
            if (t.circumCircleContains(*p))
            {
                t.isBad = true;
                polygon.push_back(Edge{ *t.a, *t.b });
                polygon.push_back(Edge{ *t.b, *t.c });
                polygon.push_back(Edge{ *t.c, *t.a });
            }
        }

        _triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](TriangleType &t) {
            return t.isBad;
        }), end(_triangles));

        for (auto e1 = begin(polygon); e1 != end(polygon); ++e1)
        {
            for (auto e2 = e1 + 1; e2 != end(polygon); ++e2)
            {
                if (almost_equal(*e1, *e2))
                {
                    e1->isBad = true;
                    e2->isBad = true;
                }
            }
        }

        polygon.erase(std::remove_if(begin(polygon), end(polygon), [](EdgeType &e) {
            return e.isBad;
        }), end(polygon));

        for (const auto e : polygon)
            _triangles.push_back(TriangleType(*e.v, *e.w, *p));

    }

    _triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3](TriangleType &t) {
        return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
    }), end(_triangles));

    for (const auto t : _triangles)
    {
        _edges.push_back(Edge{ *t.a, *t.b });
        _edges.push_back(Edge{ *t.b, *t.c });
        _edges.push_back(Edge{ *t.c, *t.a });
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


    std::vector<Edge> a = triangulation.getEdges();
    using namespace std;

    std::vector < std::pair<Vector2, std::pair<Edge, Edge>>> ans;

    int n = a.size();
    std::vector<event> e;
    for (int i = 0; i<n; ++i) {
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


    for (size_t i = 0; i<e.size(); ++i) {

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
                    continue;
                }
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
                ans.push_back(std::make_pair(*(a[id].v), std::make_pair(*nxt, *prv)));
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
    }
    if (v1 == v4) {
        p1 = v2;
        p2 = v3;
        p0 = v1;
    }
    if (v2 == v3) {
        p1 = v1;
        p2 = v4;
        p0 = v2;
    }
    if (v2 == v4) {
        p1 = v1;
        p2 = v3;
        p0 = v2;
    }
    double a1 = p1.x - p0.x, a2 = p2.x - p0.x;
    double b1 = p1.y - p0.y, b2 = p2.y - p0.y;
    double c1 = p1.z - p0.z, c2 = p2.z - p0.z;
    return p0.z + ((p.y - p0.y) * (a1 * c2 - a2 * c1) + (p.x - p0.x)*(c1*b2 - c2 * b1)) / (a1*b2 - a2 * b1);
}


int main(int argc, char * argv[])
{
    int numberPoints = 20, numberPoints_18 = 20;

    std::default_random_engine eng(std::random_device{}());
    std::uniform_real_distribution<double> dist_w(0, 800);
    std::uniform_real_distribution<double> dist_h(0, 600);
    std::uniform_real_distribution<double> dist_z(0, 20);

    std::cout << "Generating " << numberPoints << " " << numberPoints_18 << " random points" << std::endl;

    std::vector<Vector2> points_16;
    for (int i = 0; i < numberPoints; ++i) {
        points_16.push_back(Vector2{ dist_w(eng), dist_h(eng), dist_z(eng) });
    }

    // SFML window
    sf::RenderWindow window(sf::VideoMode(850, 650), "Delaunay triangulation");
    window.setFramerateLimit(1);
    // Transform each points of each vector as a rectangle
    for (const auto p : points_16) {
        sf::RectangleShape s{ sf::Vector2f(1, 1) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Cyan);
        window.draw(s);
    }

    std::vector<Vector2> points_18;
    for (int i = 0; i < numberPoints_18; ++i) {
        points_18.push_back(Vector2{ dist_w(eng), dist_h(eng), dist_z(eng) });
    }

    // Transform each points of each vector as a rectangle
    for (const auto p : points_18) {
        sf::RectangleShape s{ sf::Vector2f(1, 1) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Red);
        window.draw(s);
    }




    //=============== Ïîñòðîåíèå âûïóêëûõ îáîëî÷åê ================


    std::vector<Vector2> conv_16(points_16.size());
    std::copy(points_16.begin(), points_16.end(), conv_16.begin());

    convex_hull(conv_16);

    std::reverse(conv_16.begin(), conv_16.end());




    std::vector<Vector2> conv_18(points_18.size());

    std::copy(points_18.begin(), points_18.end(), conv_18.begin());

    convex_hull(conv_18);


    std::reverse(conv_18.begin(), conv_18.end());


    // ============================== ÎÒÁÎÐ ÒÎ×ÅÊ ====================================



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





    // =========================== ÒÐÈÀÍÃÓËßÖÈÈ ==================================









    Delaunay triangulation_16;
    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<Triangle> triangles_16 = triangulation_16.triangulate(points_16);
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

        line[0].color = sf::Color::Cyan;
        line[1].position = sf::Vector2f(
            static_cast<float>(edges_16[i].w->x + 2.),
            static_cast<float>(edges_16[i].w->y + 2.));
        line[1].color = sf::Color::Cyan;
        window.draw(line);

    }



    Delaunay triangulation_18;
    start = std::chrono::high_resolution_clock::now();
    const std::vector<Triangle> triangles_18 = triangulation_18.triangulate(points_18);
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;

    std::cout << triangles_18.size() << " triangles generated in " << diff.count() << "s\n";

    const std::vector<Edge> edges_18 = triangulation_18.getEdges();

    for (int i = 0; i < edges_18.size() - 1; i++) {
        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(edges_18[i].v->x + 2.),
            static_cast<float>(edges_18[i].v->y + 2.));

        line[0].color = sf::Color::Yellow;
        line[1].position = sf::Vector2f(
            static_cast<float>(edges_18[i].w->x + 2.),
            static_cast<float>(edges_18[i].w->y + 2.));
        line[1].color = sf::Color::Yellow;
        window.draw(line);

    }
    

    conv_18.push_back(conv_18[0]);
    conv_16.push_back(conv_16[0]);
/*

    for (int i = 0; i < conv_16.size() - 1; i++) {
        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(conv_16[i].x + 2.),
            static_cast<float>(conv_16[i].y + 2.));

        line[0].color = sf::Color::Blue;
        line[1].position = sf::Vector2f(
            static_cast<float>(conv_16[i + 1].x + 2.),
            static_cast<float>(conv_16[i + 1].y + 2.));
        line[1].color = sf::Color::Blue;
        window.draw(line);

    }

    for (int i = 0; i < conv_18.size() - 1; i++) {
        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(conv_18[i].x + 2.),
            static_cast<float>(conv_18[i].y + 2.));

        line[0].color = sf::Color::Red;
        line[1].position = sf::Vector2f(
            static_cast<float>(conv_18[i + 1].x + 2.),
            static_cast<float>(conv_18[i + 1].y + 2.));
        line[1].color = sf::Color::Red;
        window.draw(line);

    }
    */





    std::vector <Vector2> t;

    for (int i = 0; i < int(points_18_in_conv.size()); i++) {
        t.push_back(points_18_in_conv[i]);
    }
    auto tmp = scanline(triangulation_16, t);

    for (int i = 0; i < tmp.size(); i++) {
        auto pt = tmp[i].first;
        auto res = tmp[i].second;

        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(res.first.v->x + 2.),
            static_cast<float>(res.first.v->y + 2.));

        line[0].color = sf::Color::Green;
        line[1].position = sf::Vector2f(
            static_cast<float>(res.first.w->x + 2.),
            static_cast<float>(res.first.w->y + 2.));
        line[1].color = sf::Color::Green;
        window.draw(line);

        line[0].position = sf::Vector2f(
            static_cast<float>(res.second.v->x + 2.),
            static_cast<float>(res.second.v->y + 2.));

        line[0].color = sf::Color::Green;
        line[1].position = sf::Vector2f(
            static_cast<float>(res.second.w->x + 2.),
            static_cast<float>(res.second.w->y + 2.));
        line[1].color = sf::Color::Green;
        window.draw(line);

        sf::RectangleShape s{ sf::Vector2f(3, 3) };
        s.setPosition(static_cast<float>(pt.x), static_cast<float>(pt.y));
        s.setFillColor(sf::Color::White);
        window.draw(s);
    }


    std::vector <Vector2> t1;


    for (int i = 0; i < int(points_16_in_conv.size()); i++) {
        t1.push_back(points_16_in_conv[i]);
    }
    tmp = scanline(triangulation_18, t1);

    for (int i = 0; i < tmp.size(); i++) {
        auto pt = tmp[i].first;
        auto res = tmp[i].second;

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

    return 0;

}