#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <string_view>
#include <strstream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>
#include <memory>
#include <cstdio>
#include <string>
#include <format>
#include <ranges>
#include <cmath>

// minwindef.h defines these
// they are bollocks and its creator should be shot
#undef near
#undef far


using namespace std::chrono_literals;

const float epsilon = 0.000001;

struct Color
{
	uint8_t r;
	uint8_t g;
	uint8_t b;
};

namespace ConstColor 
{
	static const Color WHITE = Color{255,255,255};
	static const Color BLACK = Color{0,0,0};
}



class Window : public olc::PixelGameEngine
{
public:
	Window(const int width, const int height) :
	m_started{false}
	{
		sAppName = "Raytracer";
		this->Construct(width, height, 1, 1);
	}

public:
	bool OnUserCreate() override
	{
		m_started = true;
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		std::chrono::milliseconds elapsed = std::chrono::milliseconds{static_cast<long>(fElapsedTime * 1000)};
		std::this_thread::sleep_for(m_interval - elapsed);
		return true;
	}

	bool started()
	{
		return m_started;	
	}


private:
	const std::chrono::nanoseconds m_interval = 16666666ns;
	bool m_started;
};

class Screen {
public:
	Screen(const int width, const int height) :
		m_width{width}, m_height{height}
	{ 
		m_screen.resize(m_width * m_height);
	}

	Screen(const Screen&) = delete;
	Screen(Screen&&)      = delete;

	~Screen() 
	{
		m_window_thread.join();
	}

	void setPxl(const int x, const int y, const Color col)
	{
		getColor(x,y) = col;
		setWindowPxl(x,y);
	}

	void redrawWindow()
	{
		for(size_t x = 0; x < m_width; ++x) 
			for(size_t y = 0; y < m_height; ++y)
				setWindowPxl(x,y);
	}

	void showWindow() 
	{
		if(m_window)
			return;
		m_window = std::make_unique<Window>(m_width, m_height);
		Window* w = m_window.get();
		m_window_thread = std::thread([w](){
			w->Start();
		});
		while(!m_window->started()) {};
	}

	struct iterator
	{
		size_t x;
		const Screen& parent;

		bool operator==(const iterator& it)
		{
			return x == it.x;
		}
		bool operator!=(const iterator& it)
		{
			return !(*this == it);
		}

		iterator& operator++()
		{
			++x;
			return *this;
		}

		std::tuple<const size_t, const size_t> operator*() const 
		{
			return {x/parent.m_width, x%parent.m_width};
		}

	};
	
	iterator begin() const
	{
		return iterator{0, *this};
	}

	iterator end() const
	{
		return iterator{m_screen.size(), *this};
	}

private:
	Color& getColor(const int x, const int y)
	{
		return m_screen[x + m_width*y];
	}
	void setWindowPxl(const int x, const int y)
	{
		if(!m_window)
			return;
		const Color& pxl = getColor(x,y);
		m_window->Draw(x, y, olc::Pixel{pxl.r, pxl.g, pxl.b});
	}
private:
	const int m_width;
	const int m_height;
	std::vector<Color> m_screen;
	
	std::thread m_window_thread;
	std::unique_ptr<Window> m_window;
};

struct Point3
{
	float x,y,z;
};

struct Vect3
{
	float x,y,z;

	Vect3 operator-() const {
		return {-x, -y, -z};
	}
};

Vect3 operator-(const Point3& p, const Point3& q)
{
	return Vect3{p.x - q.x, p.y - q.y, p.z - q.z};
}

Vect3 operator+(const Vect3& v, const Vect3& w)
{
	return {v.x + w.x, v.y + w.y, v.z + w.z};
}

float operator*(const Vect3& v, const Vect3& w)
{
	return v.x*w.x + v.y*w.y + v.z*w.z;
}

Vect3 operator*(const float f, const Vect3 v)
{
	return {f*v.x, f*v.y, f*v.z};
}

Vect3 operator*(const Vect3 v, const float f)
{
	return f*v;
}

Vect3 operator^(const Vect3& v, const Vect3& w)
{
	return Vect3{
		v.y*w.z - v.z*w.y, 
		v.z*w.x - v.x*w.z, 
		v.x*w.y - v.y*w.x
	};
}

Vect3 operator/(const Vect3& v, float f)
{
	return Vect3{v.x / f, v.y / f, v.z / f};
}

Point3 operator+(const Point3& p, const Vect3& v)
{
	return {p.x + v.x, p.y + v.y, p.z + v.z};
}

Point3 operator+(const Vect3& v, const Point3& p)
{
	return p + v;
}

Point3 operator-(const Point3& p, const Vect3& v)
{
	return {p.x - v.x, p.y - v.y, p.z - v.z};
}

Vect3 operator-(const Vect3& v, const Vect3& w)
{
	return {v.x - w.x, v.y - w.y, v.z - w.z};
}

float length_sqr(const Vect3 v)
{
	return v*v;
}

Vect3 normalize(const Vect3 v)
{
	float l = std::sqrtf(length_sqr(v));
	return {v.x/l, v.y/l, v.z/l};
}

struct Triangle
{
	Point3 a,b,c;
};

struct Ray
{
	Point3 origin;
	Vect3 direction;
};

bool hit(const Triangle tr, const Ray r, float& t) 
{
	Point3 p0 = tr.a;
	Point3 p1 = tr.b;
	Point3 p2 = tr.c;

	Vect3 edge1, edge2, h, s, q;
	float a,f,u,v;

	edge1 = p1 - p0;
	edge2 = p2 - p0;
	h = r.direction ^ edge2;
	a = edge1 * h;

	if (a > -epsilon && a < epsilon)
		return false; // ray is parallel to this triangle
	
	f = 1.0 / a;
	s = r.origin - p0;
	u = f * s * h;

	if (u < 0.0 || u > 1.0)
		return false;

	q = s ^ edge1;
	v = f * r.direction * q;

	if(v < 0.0 || u + v > 1.0)
		return false;

	float thit = f * edge2 * q;

	if (thit < epsilon || thit > t)
		return false;
	t = thit;
	//Point3 intersect = r.origin + r.direction * t;
	return true;
}

// Axes Algigned Bounding Box
// Has 2 points: one near and one far from the origin {0,0,0}
struct BB
{
	Point3 near;
	Point3 far;
};

bool hit(const BB box, const Ray r, float& t) {
	float txmin = (box.near.x - r.origin.x) / r.direction.x; 
    float txmax = (box.far.x  - r.origin.x) / r.direction.x; 
	if (txmin > txmax) std::swap(txmin, txmax); 

	float tymin = (box.near.y - r.origin.y) / r.direction.y; 
    float tymax = (box.far.y  - r.origin.y) / r.direction.y; 
	if (tymin > tymax) std::swap(tymin, tymax); 

	if(txmin > tymax || tymin > txmax) return false;

	float tmin = std::max(txmin, tymin);
	float tmax = std::min(txmax, tymax);

	float tzmin = (box.near.z - r.origin.z) / r.direction.z; 
    float tzmax = (box.far.z  - r.origin.z) / r.direction.z;
	if (tzmin > tzmax) std::swap(tzmin, tzmax); 

	if(tzmin > tmax || tmin > tzmax) return false;

	float tmin = std::max(tzmin, tmin);
	float tmax = std::min(tzmax, tmax);

	float thit = tmin < 0 ? tmax : tmin;

	if (thit < epsilon || thit > t)
		return false;
	t = thit;
	return true;
}

class World {
public:
	std::vector<Triangle> m_triangles;	
	std::vector<BB> m_boxes;

	void insert(std::vector<Triangle>&& trs)
	{
		std::move(trs.begin(), trs.end(), std::back_inserter(m_triangles));
	}

	bool hit(const Ray ray) const 
	{
		float t = 100;
		return std::any_of(m_boxes.cbegin(), m_boxes.cend(), [&t, ray](const BB box) { return ::hit(box, ray, t); });
		// return std::any_of(m_triangles.cbegin(), m_triangles.cend(), [&t, ray](const Triangle tr) { return ::hit(tr, ray, t); });
	}
};

std::string to_string(float x, float y, float z) 
{
	return std::format("{{{:8.4f}, {:>8.4f}, {:>8.4f}}}", x, y ,z);
}

std::string to_string(const Vect3 v) 
{
	return to_string(v.x, v.y ,v.z);
}

std::string to_string(const Point3 v) 
{
	return to_string(v.x, v.y ,v.z);
}

void print(const Point3& p)
{
	std::printf("%s", to_string(p).c_str());
}

void print(const Vect3& v)
{
	std::printf("%s", to_string(v).c_str());
}

enum LineTag
{
	COMMENT,
	VERTEX,
	FACE,
	NORMAL,
	TEXTURE,
	OBJNAME,
	BLANK,
	LINE,
	UNKNOWN
};

LineTag parseLineTag(std::istream& in)
{
	const size_t bufsize = 3;
	char prefix[bufsize];
	in.getline(prefix, bufsize, ' ');
	
	if(prefix[0] == 'v')  
	{
		if(prefix[1] == '\0') return LineTag::VERTEX;
		if(prefix[1] == 't' && prefix[2] == '\0') return LineTag::TEXTURE;	
		if(prefix[1] == 'n' && prefix[2] == '\0') return LineTag::NORMAL;	
	}
	if(prefix[0] == 'f' && prefix[1] == '\0')  return LineTag::FACE;
	if(prefix[0] == 'l' && prefix[1] == '\0')  return LineTag::LINE;
	if(prefix[0] == 'o' && prefix[1] == '\0')  return LineTag::OBJNAME;
	if(prefix[0] == '#')  return LineTag::COMMENT;
	if(prefix[0] == '\0') 
	{
		in.clear(std::ios_base::failbit);
		in.clear(std::ios_base::badbit);
		in.clear(std::ios_base::eofbit);
		in.setstate(std::ios_base::goodbit); 
		return LineTag::BLANK;
	}

	in.setstate(std::ios_base::failbit);
	return LineTag::UNKNOWN;
}

size_t parseFaceElement(std::istream& in)
{
	size_t idx;
	in >> idx;
	return idx;
}

class Camera
{
public:
	Camera(float width, float height, float focus, Point3 position, Point3 lookat, size_t pixels_horizontal, size_t pixels_vertical) :
	m_pxls_horizontal{pixels_horizontal}, m_pxls_vertical{pixels_vertical},
	m_focus{focus}, m_width{width}, m_height{height}, m_pos{position}, m_lookat{lookat}
	{
		reset();
	}

	void reset() 
	{
		Vect3 view_dir = normalize(m_lookat - m_pos);
		Vect3 up = {0,1,0};
		Vect3 view_right = view_dir ^ up;
		Vect3 view_up    = view_right ^ view_dir;
		view_right = normalize(view_right);
		view_up    = normalize(view_up);
		m_cam_focus = m_pos - view_dir * m_focus;

		m_top_left  = m_pos + (view_up * m_height - view_right * m_width) / 2.f;
		m_bot_right = m_pos - (view_up * m_height - view_right * m_width) / 2.f;

		m_pxl_size_horizontal = (view_right * m_width ) / static_cast<float>(m_pxls_horizontal);
		m_pxl_size_vertical   = (view_up    * m_height) / static_cast<float>(m_pxls_vertical);

		std::printf("pos cam:   %s\n", to_string(m_pos).c_str());
		std::printf("cam focus: %s\n", to_string(m_cam_focus).c_str());
		std::printf("top left:  %s\n", to_string(m_top_left).c_str());
		std::printf("bot right: %s\n", to_string(m_bot_right).c_str());
		std::printf("pixel hor: %s\n", to_string(m_pxl_size_horizontal).c_str());
		std::printf("pixel ver: %s\n", to_string(m_pxl_size_vertical).c_str());
	}

	// struct iterator {
	// 	Camera* parent;
	// 	void 
	// };

	Ray toRay(size_t x, size_t y) const
	{
		Ray ray;
		ray.origin    = m_top_left + (m_pxl_size_horizontal * (x + 0.5) + m_pxl_size_vertical * (float(m_height)/2.f - y - 2.5f ));
		ray.direction = normalize(ray.origin - m_cam_focus); // Vect3{0,0,-1};
		return ray;
	}

public:
	size_t m_pxls_horizontal;
	size_t m_pxls_vertical;

	float m_focus;
	float m_width;
	float m_height;
	
	Point3 m_pos, m_lookat;
	Point3 m_top_left, m_bot_right;
	Vect3  m_pxl_size_horizontal, m_pxl_size_vertical; // pixel size in 3d world space
	Point3 m_cam_focus;
};

std::vector<Triangle> load(std::ifstream& file)
{
	std::string buff;
	long line_number = 0;
	std::vector<Point3>   points;
	std::vector<Triangle> triangles;
	while(std::getline(file, buff))
	{
		std::istringstream line(buff);
		++line_number;
		switch(parseLineTag(line))
		{
			case LineTag::COMMENT:
			case LineTag::TEXTURE:
			case LineTag::OBJNAME:
			case LineTag::BLANK:
			case LineTag::LINE:
				break;
			case LineTag::NORMAL:
				Vect3 v;
				line >> v.x;
				line >> v.y;
				line >> v.z;
				break;
			case LineTag::VERTEX:
				Point3 p;
				line >> p.x;
				line >> p.y;
				line >> p.z;
				points.emplace_back(p);
				break;
			case LineTag::FACE:
				size_t idxA,idxB,idxC;
				idxA = parseFaceElement(line) - 1;
				line.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
				idxB = parseFaceElement(line) - 1;
				line.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
				idxC = parseFaceElement(line) - 1;

				if(line.fail())
					break;
				
				if(!(idxA < points.size() && idxB < points.size() && idxC < points.size()))
					throw std::runtime_error(
						std::format(
R"(Parsing error at line {}, position {}: point out of range.
points registered: {}
points requested: {}, {}, {})",
						line_number, line.gcount(), points.size(), idxA+1, idxB+1, idxC+1)
						);
				Triangle tr;
				tr.a = points[idxA];
				tr.b = points[idxB];
				tr.c = points[idxC];
				triangles.emplace_back(tr);
				break;
			case LineTag::UNKNOWN:
				throw std::runtime_error(std::format("Parsing error at line {}, position {}: unknown tag.", line_number, line.gcount()));
				break;
		}
		if(line.fail())
			throw std::runtime_error(std::format("Parsing error at line {}, position {}.", line_number, line.gcount()));
	}
	return triangles;
}

int main()
{
	auto start = std::chrono::steady_clock::now();

	const int width  = 1000;
	const int height = 1000;

	// float width, float height, float focus, Point3 position, Point3 lookat, size_t pixels_horizontal, size_t pixels_vertical
	Camera camera{8,8, 7, Point3{1,1,-2}, Point3{0,0,0}, width, height};
	World world;

	std::ifstream file;
	file.open("cube.obj");
	if(!file.is_open()) 
	{
		std::printf("File couldn't be opened.");
		return 1;
	}
	try
	{
		world.insert(load(file));
		std::printf("Number of triangles: %llu\n", world.m_triangles.size());
	} 
	catch(std::runtime_error& e)
	{
		std::printf("%s", e.what());
		return 1;
	}


	// std::for_each(world.m_triangles.begin(), world.m_triangles.end(), 
	// 	[](const Triangle& tr) {std::printf("%s %s %s\n", to_string(tr.a).c_str(), to_string(tr.b).c_str(), to_string(tr.c).c_str());}
	// );
	

	Screen screen(width, height);
	screen.showWindow();

	Ray r;
	for(auto[x,y] : screen) 
	{
		Ray r = camera.toRay(x,y);
		// std::printf("[%d %d]   %s %s\n", x, y, to_string(r.origin).c_str(), to_string(r.direction).c_str());
		if(world.hit(r))
			screen.setPxl(x,y, Color{static_cast<uint8_t>(255*x/width), static_cast<uint8_t>(255*y/height), 0});
		else
			screen.setPxl(x,y, ConstColor::WHITE);
	}

	screen.redrawWindow();
	auto duration = std::chrono::steady_clock::now() - start;
	std::cout << "Total duration: " << floor<std::chrono::milliseconds>(duration) << "\n";
	return 0;
}
