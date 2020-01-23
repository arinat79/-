#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>
#include <cmath>
#include <limits>

#include <SFML/Graphics.hpp>

bool almost_equal(const double x, const double y, int ulp = 2)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * static_cast<double>(ulp)
        // unless the result is subnormal
        || std::abs(x - y) < std::numeric_limits<double>::min();
}

double half(const double x)
{
    return 0.5 * x;
}


struct Vector2
{
    Vector2() = default;
    Vector2(const Vector2 &v) = default;
    Vector2(Vector2&&) = default;
    Vector2(const double vx, const double vy);

    double dist2(const Vector2 &v) const;
    double dist(const Vector2 &v) const;
    double norm2() const;

    Vector2 &operator=(const Vector2&) = default;
    Vector2 &operator=(Vector2&&) = default;
    bool operator ==(const Vector2 &v) const;
    friend std::ostream &operator <<(std::ostream &str, const Vector2 &v);

    double x;
    double y;
};


Vector2::Vector2(const double vx, const double vy) :
    x(vx), y(vy)
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
    return str << "Point x: " << v.x << " y: " << v.y;
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

    const VertexType circum(half(circum_x), half(circum_y));
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
    const double midx = half(minX + maxX);
    const double midy = half(minY + maxY);

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











int main(int argc, char * argv[])
{
    int numberPoints = 100;

    std::default_random_engine eng(std::random_device{}());
    std::uniform_real_distribution<double> dist_w(0, 800);
    std::uniform_real_distribution<double> dist_h(0, 600);

    std::cout << "Generating " << numberPoints << " random points" << std::endl;

    std::vector<Vector2> points_16;
    for (int i = 0; i < numberPoints; ++i) {
        points_16.push_back(Vector2{ dist_w(eng), dist_h(eng) });
    }

    // SFML window
    sf::RenderWindow window(sf::VideoMode(850, 650), "Delaunay triangulation");
    window.setFramerateLimit(1);
    // Transform each points of each vector as a rectangle
    for (const auto p : points_16) {
        sf::RectangleShape s{ sf::Vector2f(2, 2) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Cyan);
        window.draw(s);
    }


    std::vector<Vector2> conv_16(points_16.size());
    std::copy(points_16.begin(), points_16.end(), conv_16.begin());

    convex_hull(conv_16);
    std::cout << conv_16.size();

    // Transform each points of each vector as a rectangle
    for (const auto p : conv_16) {
        sf::RectangleShape s{ sf::Vector2f(4, 4) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::Blue);
        window.draw(s);
    }

    std::reverse(conv_16.begin(), conv_16.end());


    std::vector<Vector2> points_18;
    for (int i = 0; i < numberPoints; ++i) {
        points_18.push_back(Vector2{ dist_w(eng), dist_h(eng) });
    }

    // Transform each points of each vector as a rectangle
    for (const auto p : points_18) {
        sf::RectangleShape s{ sf::Vector2f(2, 2) };
        s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
        s.setFillColor(sf::Color::White);
        window.draw(s);
    }

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

    std::cout << points_18_in_conv.size() << std::endl;


    conv_16.push_back(conv_16[0]);

    std::vector<std::array<sf::Vertex, 2> > lines;
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



    Delaunay triangulation;
    const auto start = std::chrono::high_resolution_clock::now();
    const std::vector<Triangle> triangles = triangulation.triangulate(points_16);
    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff = end - start;

    std::cout << triangles.size() << " triangles generated in " << diff.count()
        << "s\n";
    //	return 0;
    const std::vector<Edge> edges = triangulation.getEdges();

    /*

    std::cout << " ========= ";

    std::cout << "\nPoints : " << points.size() << std::endl;
    for (const auto &p : points)
        std::cout << p << std::endl;

    std::cout << "\nTriangles : " << triangles.size() << std::endl;
    for (const auto &t : triangles)
        std::cout << t << std::endl;

    std::cout << "\nEdges : " << edges.size() << std::endl;
    for (const auto &e : edges)
        std::cout << e << std::endl;

   */

    for (int i = 0; i < edges.size() - 1; i++) {
        sf::VertexArray line(sf::Lines, 2);
        line[0].position = sf::Vector2f(
            static_cast<float>(edges[i].v->x + 2.),
            static_cast<float>(edges[i].v->y + 2.));

        line[0].color = sf::Color::Cyan;
        line[1].position = sf::Vector2f(
            static_cast<float>(edges[i + 1].w->x + 2.),
            static_cast<float>(edges[i + 1].w->y + 2.));
        line[1].color = sf::Color::Cyan;
        window.draw(line);

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
