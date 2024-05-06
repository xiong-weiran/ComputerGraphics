#include "myglwidget.h"
#include <GL/glew.h>
#include <algorithm>

MyGLWidget::MyGLWidget(QWidget *parent)
	:QOpenGLWidget(parent)
{
}

MyGLWidget::~MyGLWidget()
{
	delete[] render_buffer;
	delete[] temp_render_buffer;
	delete[] normal_buffer;
	delete[] pos_mv_buffer;
	delete[] temp_z_buffer;
	delete[] z_buffer;
}

void MyGLWidget::resizeBuffer(int newW, int newH) {
	delete[] render_buffer;
	delete[] temp_render_buffer;
	delete[] temp_z_buffer;
	delete[] z_buffer;
	delete[] normal_buffer;
	delete[] pos_mv_buffer;
	WindowSizeW = newW;
	WindowSizeH = newH;
	render_buffer = new vec3[WindowSizeH*WindowSizeW];
	temp_render_buffer = new vec3[WindowSizeH*WindowSizeW];
	normal_buffer = new vec3[WindowSizeH * WindowSizeW];
	pos_mv_buffer = new vec3[WindowSizeH * WindowSizeW];
	temp_z_buffer = new float[WindowSizeH*WindowSizeW];
	z_buffer = new float[WindowSizeH*WindowSizeW];
}

void MyGLWidget::initializeGL()
{
	WindowSizeW = width();
	WindowSizeH = height();
	glViewport(0, 0, WindowSizeW, WindowSizeH);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glDisable(GL_DEPTH_TEST);
	offset = vec2(WindowSizeH / 2, WindowSizeW / 2);
	// 对定义的数组初始化
	render_buffer = new vec3[WindowSizeH*WindowSizeW];
	temp_render_buffer = new vec3[WindowSizeH*WindowSizeW];
	normal_buffer = new vec3[WindowSizeH * WindowSizeW];
	pos_mv_buffer = new vec3[WindowSizeH * WindowSizeW];
	temp_z_buffer = new float[WindowSizeH*WindowSizeW];
	z_buffer = new float[WindowSizeH*WindowSizeW];
	for (int i = 0; i < WindowSizeH*WindowSizeW; i++) {
		render_buffer[i] = vec3(0, 0, 0);
		temp_render_buffer[i] = vec3(0, 0, 0);
		normal_buffer[i] = vec3(0, 0, 0);
		pos_mv_buffer[i] = vec3(0, 0, 0);
		temp_z_buffer[i] = MAX_Z_BUFFER;			
		z_buffer[i] = MAX_Z_BUFFER;
	}
}

void MyGLWidget::keyPressEvent(QKeyEvent *e) {
	
	switch (e->key()) {
		case Qt::Key_0: scene_id = 0;update(); break;
		case Qt::Key_1: scene_id = 1;update(); break;
		case Qt::Key_9: degree += 35;update(); break;
	}
}

void MyGLWidget::paintGL()
{
	switch (scene_id) {
		case 0:scene_0(); break;
		case 1:scene_1(); break;
	}
}
void MyGLWidget::clearBuffer(vec3* now_buffer) {
	for (int i = 0; i < WindowSizeH*WindowSizeW; i++) {
		now_buffer[i] = vec3(0,0,0);
	}
}

void MyGLWidget::clearBuffer(int* now_buffer) {
	memset(now_buffer, 0, WindowSizeW * WindowSizeH * sizeof(int));
}


void MyGLWidget::clearZBuffer(float* now_buffer) {
	std::fill(now_buffer,now_buffer+WindowSizeW * WindowSizeH, MAX_Z_BUFFER);
}


// 窗口大小变动后，需要重新生成render_buffer等数组
void MyGLWidget::resizeGL(int w, int h)
{
	resizeBuffer(w, h);
	offset = vec2(WindowSizeH / 2, WindowSizeW / 2);
	clearBuffer(render_buffer);
	/*clearBuffer(normal_buffer);
	clearBuffer(pos_mv_buffer);*/
}

void MyGLWidget::scene_0()
{
	// 选择要加载的model
	objModel.loadModel("./objs/singleTriangle.obj");

	// 自主设置变换矩阵
	camPosition = vec3(100 * sin(degree * 3.14 / 180.0) + objModel.centralPoint.y, 100 * cos(degree * 3.14 / 180.0) + objModel.centralPoint.x, 10+ objModel.centralPoint.z);
	camLookAt = objModel.centralPoint;     // 例如，看向物体中心
	camUp = vec3(0, 1, 0);         // 上方向向量
	projMatrix = glm::perspective(radians(20.0f), 1.0f, 0.1f, 2000.0f);

	// 单一点光源，可以改为数组实现多光源
	lightPosition = objModel.centralPoint + vec3(0,100,100);
	clearBuffer(render_buffer);
	clearZBuffer(z_buffer);

	for (int i = 0; i < objModel.triangleCount; i++) {
		Triangle nowTriangle = objModel.getTriangleByID(i);
		drawTriangle(nowTriangle);
	}
	glClear(GL_COLOR_BUFFER_BIT);
	renderWithTexture(render_buffer,WindowSizeH,WindowSizeW);
}


void MyGLWidget::scene_1()
{
	auto start_time = std::chrono::high_resolution_clock::now();
	// 选择要加载的model
	objModel.loadModel("./objs/teapot_600.obj");
	//objModel.loadModel("./objs/teapot_8000.obj");
	//objModel.loadModel("./objs/rock.obj");
	//objModel.loadModel("./objs/cube.obj");
	//objModel.loadModel("./objs/singleTriangle.obj");
	
	// 自主设置变换矩阵
	camPosition = vec3(100 * sin(degree * 3.14 / 180.0) + objModel.centralPoint.y, 100 * cos(degree * 3.14 / 180.0) + objModel.centralPoint.x, 10+ objModel.centralPoint.z);
	/*camPosition = vec3(50 * sin(degree * 3.14 / 180.0) + objModel.centralPoint.y,
		50 * cos(degree * 3.14 / 180.0) + objModel.centralPoint.x,
		10 + objModel.centralPoint.z);*/
	camLookAt = objModel.centralPoint;     // 例如，看向物体中心
	camUp = vec3(0, 1, 0);         // 上方向向量
	projMatrix = glm::perspective(radians(20.0f), 1.0f, 0.1f, 2000.0f);

	// 单一点光源，可以改为数组实现多光源
	lightPosition = objModel.centralPoint + vec3(0,100,100);
	//lightPosition = objModel.centralPoint + vec3(50, 150, 50);
	clearBuffer(render_buffer);
	clearZBuffer(z_buffer);
	/*clearBuffer(normal_buffer);
	clearBuffer(pos_mv_buffer);*/
	for (int i = 0; i < objModel.triangleCount; i++) {
		Triangle nowTriangle = objModel.getTriangleByID(i);
		drawTriangle(nowTriangle);
	}
	glClear(GL_COLOR_BUFFER_BIT);
	renderWithTexture(render_buffer, WindowSizeH, WindowSizeW);
	auto end_time = std::chrono::high_resolution_clock::now();

	// 计算时间差并输出
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << "Time taken by scene_1: " << duration.count() << " milliseconds" << std::endl;
}

void MyGLWidget::original_set_pixel(int x, int y, const vec3& color, float depth) {
	/*if (x < 0 || x >= WindowSizeW || y < 0 || y >= WindowSizeH) {

	}*/
	if (x >= 0 && x < WindowSizeW && y >= 0 && y < WindowSizeH) {
		int index = y * WindowSizeW + x;
		temp_render_buffer[index] = color;
		temp_z_buffer[index] = depth;
	}
}

void MyGLWidget::set_pixel(int x, int y, const vec3& color, float depth, vec3 normal, vec3 pos_mv) {
	// 检查像素是否在窗口范围内
	if (x < 0 || x >= WindowSizeW || y < 0 || y >= WindowSizeH) {
		// 像素位于窗口之外，因此不执行任何操作
		
	}
	if (x >= 0 && x < WindowSizeW && y >= 0 && y < WindowSizeH) {
		// 计算像素在帧缓冲数组中的线性索引
		int index = y * WindowSizeW + x;
		// 在缓冲中设置像素的颜色、深度、法线、model_view中位置
		temp_render_buffer[index] = color;
		temp_z_buffer[index] = depth;
		normal_buffer[index] = normal;
		pos_mv_buffer[index] = pos_mv;
	}
}


void MyGLWidget::drawTriangle(Triangle triangle) {
	// 三维顶点映射到二维平面
	vec3* vertices = triangle.triangleVertices;
	vec3* normals = triangle.triangleNormals;
	FragmentAttr transformedVertices[3];
	float ka, kd, ks;
	ka = 0.2;
	kd = 0.6;
	ks = 0.3;
	clearBuffer(this->temp_render_buffer);
	clearZBuffer(this->temp_z_buffer);
	/*clearBuffer(this->normal_buffer);
	clearBuffer(this->pos_mv_buffer);*/
	mat4 viewMatrix = glm::lookAt(camPosition, camLookAt, camUp);

	for (int i = 0; i < 3; ++i) {
		vec4 ver_mv = viewMatrix * vec4(vertices[i], 1.0f);
		float nowz = glm::length(camPosition - vec3(ver_mv));
		vec4 ver_proj = projMatrix * ver_mv;
		transformedVertices[i].x = ver_proj.x + offset.x;
		transformedVertices[i].y = ver_proj.y + offset.y;
		transformedVertices[i].z = nowz;
		transformedVertices[i].pos_mv = ver_mv;
		mat3 normalMatrix = mat3(viewMatrix);
		vec3 normal_mv = normalMatrix * normals[i];
		transformedVertices[i].normal = normal_mv;

		transformedVertices[i].color = phong_model(transformedVertices[i].normal, transformedVertices[i].pos_mv,transformedVertices[i].color ,ka,kd, ks);
	}

	// 将当前三角形渲染在temp_buffer中

	/*DDA(transformedVertices[0], transformedVertices[1], 1);
	DDA(transformedVertices[1], transformedVertices[2], 2);
	DDA(transformedVertices[2], transformedVertices[0], 3);*/

	// HomeWork: 1、绘制三角形三边
	bresenham(transformedVertices[0], transformedVertices[1], 1);
	bresenham(transformedVertices[1], transformedVertices[2], 2);
	bresenham(transformedVertices[2], transformedVertices[0], 3);

	// HomeWork: 2: 用edge-walking填充三角形内部到temp_buffer中
	int firstChangeLine = edge_walking();

	// 合并temp_buffer 到 render_buffer, 深度测试
	// 从firstChangeLine开始遍历，可以稍快
	for (int h = firstChangeLine; h < WindowSizeH; h++) {
		auto render_row = &render_buffer[h * WindowSizeW];
		auto temp_render_row = &temp_render_buffer[h * WindowSizeW];
		auto z_buffer_row = &z_buffer[h * WindowSizeW];
		auto temp_z_buffer_row = &temp_z_buffer[h * WindowSizeW];
		for (int i = 0; i < WindowSizeW; i++) {
			if (z_buffer_row[i] < temp_z_buffer_row[i])
				continue;
			else
			{
				z_buffer_row[i] = temp_z_buffer_row[i];
				render_row[i] = temp_render_row[i];
			}
		}

	}
}

vec3 MyGLWidget::phong_model(const vec3& normal, const vec3& pos_mv, const vec3& original_color, float ka, float kd, float ks)
{
	// 定义光源颜色
	vec3 light_source(1.0f, 1.0f, 1.0f);

	// 对法向量进行标准化
	vec3 normal_normalized = normalize(normal);

	// 计算视线方向
	vec3 view_direction = normalize(camPosition - pos_mv);

	// 计算光线方向
	vec3 light_direction = normalize(lightPosition - pos_mv);

	// 计算环境光照
	vec3 ambient = ka * light_source * original_color;

	// 计算漫反射光照
	vec3 diffuse = kd * max(dot(normal_normalized, light_direction), 0.0f) * light_source * original_color;

	// 计算镜面反射光照
	vec3 reflect_direction = reflect(-light_direction, normal_normalized);
	float s = pow(max(dot(view_direction, reflect_direction), 0.0f), 32);
	vec3 specular = ks * s * light_source;

	// 合成最终颜色并限制在合理范围
	vec3 phong_model_color = ambient + diffuse + specular;
	return clamp(phong_model_color, 0.0f, 1.0f);
}



void MyGLWidget::DDA(FragmentAttr& start, FragmentAttr& end, int id) {
	int x1 = start.x;
	int y1 = start.y;
	int x2 = end.x;
	int y2 = end.y;

	int dx = x2 - x1, dy = y2 - y1;  // 计算 x 和 y 方向的增量
	int steps;
	float delta_x, delta_y, x = x1, y = y1;
	vec3 color(1.0f, 0.0f, 0.0f);  // 设置像素颜色为红色

	steps = (abs(dx) > abs(dy)) ? (abs(dx)) : (abs(dy));  // 计算步数，取增量较大的方向
	delta_x = dx / (float)steps;  // 计算 x 方向每步的增量
	delta_y = dy / (float)steps;  // 计算 y 方向每步的增量

	original_set_pixel(round(x), round(y), color,0);  // 设置起点像素颜色

	for (int i = 0; i < steps; ++i) {
		x += delta_x;  // 更新 x 坐标
		y += delta_y;  // 更新 y 坐标
		original_set_pixel(round(x), round(y), color,0);  // 设置当前像素颜色
	}
}

// 使用Bresenham算法在屏幕上绘制直线
void MyGLWidget::bresenham(FragmentAttr& start, FragmentAttr& end, int id) {
	// 提取起点和终点坐标
	int x1 = start.x;
	int y1 = start.y;
	int x2 = end.x;
	int y2 = end.y;
	vec3 black(0.0f, 0.0f, 0.0f);
	// 初始化当前像素的坐标
	int x = x1, y = y1;
	FragmentAttr interpolate;
	FragmentAttr temp;
	// 计算 x 和 y 方向上的增量
	int dx = abs(x2 - x1);
	int dy = abs(y2 - y1);

	// 判断斜率的变化是否大于等于1
	int k = (dx >= dy) ? 1 : 0;

	// 计算 x 和 y 方向的增量方向
	int sx = (x2 > x1) ? 1 : -1;
	int sy = (y2 > y1) ? 1 : -1;

	// 初始化 Bresenham 算法中的判别式参数 p
	int p = (k) ? (2 * dy - dx) : (2 * dx - dy);

	// 根据斜率情况循环绘制直线
	if (k) {
		while (x != x2 || y != y2) {
			// 根据判别式的值的正负选择不同的更新规则
			
			if (p <= 0) {
				x += sx;
				//t = fabs(float(x-x1)/(x1-x2));

				//线性插值计算边上点的信息
				interpolate = getLinearInterpolation(start, end, x,y);
				/*if (interpolate.color == black)
				{
					std::cout << "black" << endl;
				}*/
				//normal_buffer[interpolate.y * WindowSizeW + x] = interpolate.normal;
				//pos_mv_buffer[interpolate.y * WindowSizeW + x] = interpolate.pos_mv;
				set_pixel(x, y, interpolate.color, interpolate.z, interpolate.normal,interpolate.pos_mv);
				//if (std::isnan(interpolate.color.r)) {
				//	//std::cout << "left:" << temp_z_buffer[i * WindowSizeW + left_boarder-1] << std::endl;
				//	for (int j = 0; j < WindowSizeW; j++)
				//	{
				//		std::cout << j << temp_z_buffer[y * WindowSizeW + j] << std::endl;
				//	}
				//	std::cout << "right:" << temp_z_buffer[y * WindowSizeW + x] << std::endl;
				//}
				p = p + 2 * dy;
			}
			else {
				x += sx;
				y += sy;
				//t = fabs(float(x - x1) / (x1 - x2));
				interpolate = getLinearInterpolation(start, end, x,y);
				/*if (interpolate.color == black)
				{
					std::cout << "black" << endl;
				}*/
				//normal_buffer[interpolate.y * WindowSizeW + x] = interpolate.normal;
				//pos_mv_buffer[interpolate.y * WindowSizeW + x] = interpolate.pos_mv;
				set_pixel(x, y, interpolate.color, interpolate.z, interpolate.normal, interpolate.pos_mv);
				//if (std::isnan(interpolate.color.r)) {
				//	//std::cout << "left:" << temp_z_buffer[i * WindowSizeW + left_boarder-1] << std::endl;
				//	for (int j = 0; j < WindowSizeW; j++)
				//	{
				//		std::cout << j << temp_z_buffer[y * WindowSizeW + j] << std::endl;
				//	}
				//	std::cout << "right:" << temp_z_buffer[y * WindowSizeW + x] << std::endl;
				//}
				p = p + 2 * (dy - dx);
			}
		}
	}
	else {
		while (x != x2 || y != y2) {
			// 根据判别式的值的正负选择不同的更新规则
			if (p <= 0) {
				y += sy;
				//t = fabs(float(x - x1) / (x1 - x2));
				interpolate = getLinearInterpolation(start, end, x,y);
				/*if (interpolate.color == black)
				{
					std::cout << "black" << endl;
				}*/
				//normal_buffer[interpolate.y * WindowSizeW + x] = interpolate.normal;
				//pos_mv_buffer[interpolate.y * WindowSizeW + x] = interpolate.pos_mv;
				set_pixel(x, y, interpolate.color, interpolate.z, interpolate.normal, interpolate.pos_mv);
				/*if (std::isnan(interpolate.color.r)) {
					std::cout << "y:" << y << std::endl;
					for (int j = 0; j < WindowSizeW; j++)
					{
						std::cout << j << temp_z_buffer[y * WindowSizeW + j] << std::endl;
					}
					std::cout << "right:" << temp_z_buffer[y * WindowSizeW + x] << std::endl;
				}*/
				p = p + 2 * dx;
			}
			else {
				x += sx;
				y += sy;
				//t = fabs(float(x - x1) / (x1 - x2));
				interpolate = getLinearInterpolation(start, end, x,y);
				/*if (interpolate.color == black)
				{
					std::cout << "black" << endl;
				}*/
				//normal_buffer[interpolate.y * WindowSizeW + x] = interpolate.normal;
				//pos_mv_buffer[interpolate.y * WindowSizeW + x] = interpolate.pos_mv;
				set_pixel(x, y, interpolate.color, interpolate.z, interpolate.normal, interpolate.pos_mv);
				//if (std::isnan(interpolate.color.r)) {
				//	//std::cout << "left:" << temp_z_buffer[i * WindowSizeW + left_boarder-1] << std::endl;
				//	for (int j = 0; j < WindowSizeW; j++)
				//	{
				//		std::cout << j << temp_z_buffer[y * WindowSizeW + j] << std::endl;
				//	}
				//	std::cout << "right:" << temp_z_buffer[y * WindowSizeW + x] << std::endl;
				//}
				p = p + 2 * (dx - dy);
			}
		}
	}
}


int MyGLWidget::edge_walking() {
	// 初始化第一个颜色变化的水平线位置为窗口高度
	int firstChangeLine = WindowSizeH;
	FragmentAttr interpolate,temp1,temp2;
	float t;
	vec3 black(0.0f, 0.0f, 0.0f);

	// 从底部向上遍历每一行
	for (int i = WindowSizeH - 1; i >= 0; i--) {
		int left_boarder = -1;
		int right_boarder = -1;

		// 遍历当前行的每个像素
		for (int x = 0; x <= WindowSizeW; x++) {
			// 如果当前像素不是黑色
			if (temp_render_buffer[i * WindowSizeW + x] != black) {
				firstChangeLine = i;

				// 如果左边界已经记录，说明当前这个非黑色点对应的是右边界，记录右边界并结束循环
				if (left_boarder != -1) {
					right_boarder = x ;
					break;
				}

				// 如果左边界未记录，记录左边界
				if (temp_render_buffer[i * WindowSizeW + x + 1] == black) {
					left_boarder = x ;
				}
			}
		}

		// 如果左边界已记录而右边界未记录，则表示这一行只有一个三角形顶点，不需要用edge_walking上色，所以将left_boarder设为-1
		if (left_boarder != -1 && right_boarder == -1) {
			//right_boarder = left_boarder;
			left_boarder = -1;
		}

		// 如果左边界已记录，则将该行在三角形范围内的像素设置为灰色
		if (left_boarder != -1) {
			//获取该横行上在三角形上的两个点的信息
			temp1.x = left_boarder;
			temp1.y = i;
			temp1.color = temp_render_buffer[i * WindowSizeW + left_boarder];
			temp1.z = temp_z_buffer[i * WindowSizeW + left_boarder];
			temp1.normal = normal_buffer[i * WindowSizeW + left_boarder];
			temp1.pos_mv = pos_mv_buffer[i * WindowSizeW + left_boarder];
			temp2.x = right_boarder;
			temp2.y = i;
			temp2.color = temp_render_buffer[i * WindowSizeW + right_boarder];
			temp2.z = temp_z_buffer[i * WindowSizeW + right_boarder];
			temp2.normal = normal_buffer[i * WindowSizeW + right_boarder];
			temp2.pos_mv = pos_mv_buffer[i * WindowSizeW + right_boarder];
			for (int x = left_boarder; x <= right_boarder; x++) {
			//Gouraud shading		
				//interpolate = getLinearInterpolation(temp1, temp2, x,i);
				////printf("-%.2lf %.2lf %.2lf-\n", interpolate.color.r, interpolate.color.g, interpolate.color.b);
				////if (std::isnan(interpolate.color.r)) {
				////	//std::cout << "left:" << temp_z_buffer[i * WindowSizeW + left_boarder-1] << std::endl;
				////	for (int j = 0; j < WindowSizeW; j++)
				////	{
				////		std::cout << j<<temp_z_buffer[i * WindowSizeW + j] << std::endl;
				////	}
				////	std::cout <<"right:"<<temp_z_buffer[i * WindowSizeW + right_boarder] << std::endl;
				////}
				//set_pixel(x, i, interpolate.color, interpolate.z, interpolate.normal, interpolate.pos_mv);
				// 
				// 
				// 
				//phong模型
				interpolate = getLinearInterpolation(temp1, temp2, x,i);
				//printf("-%.2lf %.2lf %.2lf-\n", temp.color.r, temp.color.b, temp.color.g);
				interpolate.color= vec3(1.0f, 1.0f, 0.0f);
				//vec3 phong_color = PhoneShading(interpolate);
				//BlinnPhong模型
				vec3 phong_color = BlinnPhongShading(interpolate);
				interpolate.color = phong_color;
				set_pixel(x, i, interpolate.color, interpolate.z, interpolate.normal, interpolate.pos_mv); // setPixel函数用于在特定位置绘制像素
			
			}

		}
	}

	// 返回第一个发生颜色变化的水平线位置
	return firstChangeLine;
}


vec3 MyGLWidget::PhoneShading(FragmentAttr& nowPixelResult) {
	// 定义光源颜色
	vec3 light(1.0f, 1.0f, 1.0f);

	// 定义环境光、镜面反射的系数和幂指数
	float ka, s, alpha;

	// 计算观察方向
	vec3 view_direction = normalize(camPosition - nowPixelResult.pos_mv);

	// 计算光线方向
	vec3 light_direction = normalize(lightPosition - nowPixelResult.pos_mv);

	// 计算法线的规范化向量
	vec3 normal_normalized = normalize(nowPixelResult.normal);

	// 计算反射方向
	vec3 reflect_direction = reflect(-light_direction, normal_normalized);

	// 设置环境光系数和镜面反射的幂指数
	ka = 0.2;
	alpha = 100;

	// 计算环境光分量
	vec3 ambient = ka * light * nowPixelResult.color;

	// 计算漫反射分量
	vec3 diffuse = max(dot(normal_normalized, light_direction), 0.0f) * light * nowPixelResult.color;

	// 计算镜面反射分量
	s = pow(max(dot(view_direction, reflect_direction), 0.0f), alpha);
	vec3 specular = s * light; // 计算镜面反射的颜色分量

	// 最终颜色是环境光、漫反射和镜面反射三者相加
	vec3 finalColor = ambient + diffuse + specular;

	// 限制颜色范围在 [0, 1] 内
	return clamp(finalColor, 0.0f, 1.0f);
}


vec3 MyGLWidget::BlinnPhongShading(FragmentAttr& nowPixelResult) {
	// 定义光源颜色
	vec3 light(1.0f, 1.0f, 1.0f);

	// 定义环境光、镜面反射的系数和幂指数
	float ka, s, alpha;

	// 计算观察方向
	vec3 view_direction = normalize(camPosition - nowPixelResult.pos_mv);

	// 计算光线方向
	vec3 light_direction = normalize(lightPosition - nowPixelResult.pos_mv);

	// 计算法线的规范化向量
	vec3 normal_normalized = normalize(nowPixelResult.normal);

	// 计算半向量
	vec3 h = normalize(light_direction + view_direction);

	// 设置环境光系数和镜面反射的幂指数
	ka = 0.2;
	alpha = 100;

	// 计算环境光分量
	vec3 ambient = ka * light * nowPixelResult.color;

	// 计算漫反射分量
	vec3 diffuse = max(dot(normal_normalized, light_direction), 0.0f) * light * nowPixelResult.color;

	// 计算镜面反射分量
	s = pow(max(dot(normal_normalized, h), 0.0f), alpha);
	vec3 specular = s * light; // 计算镜面反射的颜色分量

	// 最终颜色是环境光、漫反射和镜面反射三者相加
	vec3 finalColor = ambient + diffuse + specular;

	// 限制颜色范围在 [0, 1] 内
	return clamp(finalColor, 0.0f, 1.0f);
}






