
#include "feature.h"
#include "raster.h"
#include "templates.hpp"
#include <type_traits>

// Matrix
GD::Matrix<float> AverageKernel(size_t size) {
	if ((size & 1) == 0) {
		std::cout << "size must be odd" << std::endl;
	}

	GD::Matrix<float> res(size, size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; j++) {
			res.Data[i][j] = 1.0 / (size * size);
		}
	}
	return res;
}

GD::Matrix<float> SobelKernelx() {
	return GD::Matrix<float>(
		{ {-1., 0.,  1.},
		  {-2., 0.,  2.},
		  {-1., 0.,  1.} });
}

GD::Matrix<float> SobelKernely() {
	return GD::Matrix<float>(
		{ {-1., -2.,  -1.},
		  { 0., 0.,  0.},
		  { 1., 2.,  1.} });
}

GD::Matrix<float> LaplaceKernel_8() {
	return GD::Matrix<float>(
		{ {1,  1, 1},
		  {1, -8, 1},
		  {1,  1, 1} }
	);
}

GD::Matrix<float> LaplaceKernel_4() {
	return GD::Matrix<float>(
		{ {0,  1, 0},
		  {1, -4, 1},
		  {0,  1, 0} }
	);
}

GD::Matrix<float> LaplaceSharpening_4() {
	return GD::Matrix<float>(
		{ {0,  1, 0},
		  {1, -5, 1},
		  {0,  1, 0} }
	);
}

GD::Matrix<float> LaplaceSharpening_8() {
	return GD::Matrix<float>(
		{ {1,  1, 1},
		  {1, -9, 1},
		  {1,  1, 1} }
	);
}

GD::Matrix<float> SobelKernel() {
	return GD::Matrix<float> ({
		{0, -1, 0},
		{-1, 0, 1},
		{0, 1, 0}
	});
}

GD::Matrix<float> GaussianKernel(size_t size, float stdX, float stdY) {
	GD::Matrix<float> ret(size, size);
	float sum = 0;
	for (int i = 0; i < size; ++i) {
		float x = i - size / 2;
		for (int j = 0; j < size; ++j) {
			float y = j - size / 2;
			ret.Data[i][j] = exp(-(x * x / (2 * stdX * stdX) + y * y / (2 * stdY * stdY)));
			sum += ret.Data[i][j];
		}
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; j++) {
			ret.Data[i][j] /= sum;
		}
	}
	return ret;
}

template<typename _Ty>
GD::Matrix<_Ty>::Matrix(size_t width, size_t height) :m_width(width), m_height(height) {
	Data = std::vector<_Ty*>(height, nullptr);
	for (int i = 0; i < height; ++i) {
		Data[i] = new _Ty[width];
		memset((void*)Data[i], 0, width * sizeof(_Ty));
	}
	int centerx = width / 2, centery = height / 2;
	Data[centerx][centery] = 1;
}
template<typename _Ty>
GD::Matrix<_Ty>::Matrix(std::initializer_list<std::initializer_list<_Ty>> const& init_list) {
	size_t listW = init_list.begin()->size(), listH = init_list.size();
	m_width = listW;
	m_height = listH;
	if (listW != m_width || listH != m_height) {
		std::cerr << "Matrix init error" << std::endl;
		exit(1);
	}
	int i = 0;
	Data = std::vector<_Ty*>(m_height);
	for (auto const& list : init_list) {
		if (list.size() != m_width) {
			std::cerr << "Matrix init error" << std::endl;
			exit(1);
		}
		int j = 0;
		Data[i] = new _Ty[m_width];
		memset((void*)Data[i], 0, m_width * sizeof(_Ty));
		for (auto e : list) {
			Data[i][j++] = e;
		}
		++i;
	}
}
// RasterBase
GD::RasterBase::RasterBase(const std::string& filepath) {}
GD::RasterBase::RasterBase(const char* filepath) {}
GD::Array_features GD::RasterBase::get(size_t i, size_t j) const { return Array_features(); }
float GD::RasterBase::get(size_t i, size_t j, size_t channel) const { return 0; }
size_t GD::RasterBase::width() const { return 0; }
size_t GD::RasterBase::height() const { return 0; }
size_t GD::RasterBase::bands() const { return 0; }
void* GD::RasterBase::data() const { return nullptr; }
GDALDataType GD::RasterBase::type() const { return GDT_Unknown; }
void GD::RasterBase::init(
	GDALDataset* pDatasetRead,
	size_t lWidth, size_t lHeight, size_t nBands,
	GDALDataType datatype) {};
void GD::RasterBase::save(std::string const&, std::string const& DriverName)const {};
std::vector<std::array<size_t, 256>> GD::RasterBase::histgram() const { return  std::vector< std::array<size_t, 256>>(); }
GD::RasterPtr GD::RasterBase::equalizeHist(std::vector< std::array<size_t, 256>> histgram) { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::equalizeHist() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::channel_extraction(size_t channel) { return GD::RasterPtr(nullptr); };
GD::RasterPtr GD::RasterBase::filter(Matrix<float> const& kernel) { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::MedianFilter() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::addSalt(size_t num) { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Sobel() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Laplace_4() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Laplace_8() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::GaussBlur(size_t size, float stdx, float stdy) { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::addGaussNoise(size_t size, float stdy) { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Laplace_4_Sharpening() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Laplace_8_Sharpening() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::Sobel_Sharpening() { return GD::RasterPtr(nullptr); }
GD::RasterPtr GD::RasterBase::toGray() { return GD::RasterPtr(nullptr); }
GD::RasterBase::~RasterBase() {};
void GD::RasterBase::resize(size_t width, size_t height) {};
void GD::RasterBase::set(size_t x, size_t y, GD::Array_features const& arr) {};
void GD::RasterBase::set(size_t x, size_t y, Eigen::Vector3f const& arr) {};
void GD::RasterBase::set(size_t x, size_t y, size_t channel, double const value) {};
void GD::RasterBase::add(size_t i, size_t j, size_t channel, double const value) {};
void GD::RasterBase::add(size_t i, size_t j, double const value) {};
void GD::RasterBase::add(size_t i, size_t j, size_t channel, GD::Array_features const& arr) {};
size_t GD::RasterBase::type_size() const { return 0; }
size_t GD::RasterBase::Area()const { return 0; }
GD::RasterPtr GD::RasterBase::copy() const { return nullptr; }
void GD::RasterBase::setnullptr() {};
// !RasterBase

// GD::Raster<T>
template<typename T>
std::vector<float> GD::Raster<T>::conv(long long pos, GD::Matrix<float> const& mat) {
	std::vector<float> res;
	long long newpos = pos % area;
	size_t x = newpos % lWidth, y = newpos / lWidth;
	long long n = pos / area;
	size_t matw = mat.m_width, math = mat.m_height;
	long long matx = matw / 2, maty = math / 2;
	for (int j = -maty; j <= maty; ++j) {
		int yj = y + j;
		for (int i = -matx; i <= matx; ++i) {
			int xi = x + i;
			if (xi < 0 || xi >= lHeight || yj < 0 || yj >= lWidth) {
				res.push_back(0);
				continue;
			}
			long long index = yj * lWidth + xi + n * area;
			res.push_back(((float)pData[index]) * ((float)mat.Data[maty + j][matx + i]));
		}
	}
	return res;
}
template<typename T>
GD::RasterPtr GD::Raster<T>::Laplace(Matrix<float> const& kernel) {
	GD::RasterPtr pLaplace{ filter(kernel) };
	GD::RasterPtr ret{ std::make_shared<GD::Raster<T>>(GD::Raster<T>(lWidth, lHeight, 1, datatype)) };
	ret->init(nullptr, lWidth, lHeight, 1, datatype);
	long long total = lWidth * lHeight * nBands;
	for (size_t i = 0; i < total; ++i) {
		/*double sum{ 0.0 };
		for (long long channel = 0; channel < nBands; ++channel) {
			long long index = i + channel * area;
			sum += ((T*)pLaplace->data())[index];
		}
		sum /= 3;
		((T*)ret->data())[i] = (T)sum;*/
		((T*)ret->data())[i%area] += ((T*)pLaplace->data())[i]/3;
	}
	return ret;
}
//
template<typename T>
GD::Array_features GD::Raster<T>::get(size_t i, size_t j) const {
	Array_features ret;
	for (int channel = 0; channel < this->nBands; ++channel) {
		auto pos = POS(j, i, channel, lWidth, area);
		ret.push_back(float(this->pData[pos]));
	}
	return ret;
}
//
template<typename T>
float GD::Raster<T>::get(size_t i, size_t j, size_t channel) const {
	auto pos = POS(j, i, channel, lWidth, area);
	return float(this->pData[pos]);
}
//
template<typename T>
size_t GD::Raster<T>::width() const { return lWidth; }
template<typename T>
size_t GD::Raster<T>::height() const { return lHeight; }
template<typename T>
size_t GD::Raster<T>::bands() const { return nBands; }
template<typename T>
void* GD::Raster<T>::data() const { return (void*)pData; }
template<typename T>
GDALDataType GD::Raster<T>::type() const { return datatype; }
template<typename T>
void GD::Raster<T>::init(
	GDALDataset* pDatasetRead,
	size_t lWidth, size_t lHeight, size_t nBands,
	GDALDataType datatype)
{
	this->pDatasetRead = pDatasetRead;
	this->lWidth = lWidth; this->lHeight = lHeight; this->nBands = nBands;
	this->datatype = datatype;
	this->area = lWidth * lHeight;
	T* p = GD::Malloc<T>(lWidth, lHeight, nBands);
	if (this->pDatasetRead)
		pDatasetRead->RasterIO(GF_Read, 0, 0, lWidth, lHeight, (void*)p, lWidth, lHeight, datatype, nBands, NULL, 0, 0, 0);//读取图像数据
	this->pData = p;
}
template<typename T>
void GD::Raster<T>::save(std::string const& filepath, std::string const& DriverName) const {
	GDALDataset* pDatasetSave;
	auto pDriver = GetGDALDriverManager()->GetDriverByName(DriverName.c_str());
	pDatasetSave = pDriver->Create(filepath.c_str(), lWidth, lHeight, nBands, datatype,
		nullptr);
	assert(pDatasetSave != nullptr);
	pDatasetSave->RasterIO(GF_Write, 0, 0, lWidth, lHeight, (void*)pData, lWidth, lHeight, datatype, nBands, NULL, 0, 0, 0);//创建图像数据
	GDALClose(pDatasetSave);
}

template<typename T>
std::vector<std::array<size_t, 256>> GD::Raster<T>::histgram() const {
	std::vector<std::array<size_t, 256>> ret(nBands);
	for (int i = 0; i < nBands; ++i) {
		ret[i].fill(0);
	}
	long long total = this->lWidth * this->lHeight;
	for (int channel = 0; channel < nBands; ++channel) {
		auto p = pData + channel * total;
		for (int i = 0; i < total; ++i) {
			++ret[channel][p[i]];
		}
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::equalizeHist(std::vector<std::array<size_t, 256>> histgram) {
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	size_t total = lHeight * lWidth;
	memcpy((T*)ret->data(), (T*)this->data(), nBands * total * sizeof(T));
	for (int channel = 0; channel < nBands; ++channel) {
		T* srcbegin = (T*)this->data() + lHeight * lWidth * channel;
		T* tarbegin = (T*)ret->data() + lHeight * lWidth * channel;
		for (int i = 1; i < 256; ++i) {
			histgram[channel][i] += histgram[channel][i - 1];
		}
		for (int i = 0; i < total; ++i) {
			tarbegin[i] = static_cast<T>(255.0 * histgram[channel][srcbegin[i]] / total);
		}
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::equalizeHist() {
	return equalizeHist(histgram());
}

template<typename T>
GD::RasterPtr GD::Raster<T>::channel_extraction(size_t channel) {
	assert(channel < nBands);
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	size_t total = lHeight * lWidth;
	T* srcbegin = pData + lHeight * lWidth * channel;
	T* tarbegin = (T*)ret->data() + lHeight * lWidth * channel;
	memset(ret->data(), 0, total * sizeof(T));
	for (int i = 0; i < total; ++i) {
		tarbegin[i] = srcbegin[i];
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::filter(Matrix<float> const& kernel) {
	if (kernel.m_width % 2 == 0 || kernel.m_height % 2 == 0) {
		std::cerr << "kernel size must be odd" << std::endl;
		exit(1);
	}
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	size_t total = lHeight * lWidth;
	auto conv_ = [&](size_t begin) {
		for (long long i = begin; i < total + begin; ++i) {
			auto res = conv(i, kernel);
			double s = sum(res);
			if (s < 0) s = fabs(s);
			if (s > GD::Type_Max<T>()) s = GD::Type_Max<T>();
			((T*)ret->data())[i] = (T)(s);
		}
	};
	std::vector<std::thread> threads;
	for(int i = 0; i < bands(); ++ i){
		threads.emplace_back(conv_, i * total);
	}
	for (auto& t : threads) {
		t.join();
	}

	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::MedianFilter() {
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	size_t total = lHeight * lWidth * nBands;
	Matrix<float> kernel{ {1.f, 1.f, 1.f} ,{1.f, 1.f, 1.f}, {1.f, 1.f, 1.f} };
#pragma omp parallel for num_threads(10)
	for (long long i = 0; i < total; ++i) {
		auto res = conv(i, kernel);
		((T*)ret->data())[i] = qSelect(res, res.size() / 2);
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::addSalt(size_t num) {
	size_t limit = lWidth * lHeight;
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	memcpy(ret->data(), this->data(), limit * nBands);
	for (long long i = 0; i < num; ++i) {
		int randomIndex = (rand() * rand()) % limit;
		bool is_white = rand() % 2;
		for (int channel = 0; channel < nBands; ++channel) {
			((T*)ret->data())[randomIndex + channel * limit] = is_white ? 255 : 0;
		}
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Sobel() {
	size_t total = lWidth * lHeight;
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, 1, datatype);
	GD::RasterPtr Sobelx;
	auto t1 = [&]() {
		Sobelx = filter(SobelKernelx());
	};
	std::thread t(t1);
	GD::RasterPtr Sobely(filter(SobelKernely()));
	t.join();
	total *= nBands;
#pragma omp parallel for num_threads(4)
	for (long long i = 0; i < total; ++i) {
		long long index = i % area;
		double x = ((T*)Sobelx->data())[i];
		double y = ((T*)Sobely->data())[i];
		double s = (x + y)/2 + ((T*)ret->data())[index];
		if (s > GD::Type_Max<T>()) s = GD::Type_Max<T>();
		((T*)ret->data())[index] = (T)(s);
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::GaussBlur(size_t size, float stdx, float stdy) {
	return filter(GaussianKernel(size, stdx, stdy));
}
template<typename T>
GD::RasterPtr GD::Raster<T>::addGaussNoise(size_t size, float std) {
	size_t total = lHeight * lWidth * nBands;
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	memcpy(ret->data(), this->data(), total * sizeof(T));
	std::default_random_engine e;
	std::normal_distribution<float> x(0, std);
	for (long long i = 0; i < size; ++i) {
		int randomIndex = (rand() * rand()) % total;
		double s = ((T*)ret->data())[randomIndex] + x(e);
		if (s < 0) s = 0;
		if (s > GD::Type_Max<T>()) s = GD::Type_Max<T>();
		((T*)ret->data())[randomIndex] = (T)(s);
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Laplace_4() {
	return Laplace(LaplaceKernel_4());
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Laplace_8() {
	return Laplace(LaplaceKernel_8());
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Laplace_4_Sharpening() {
	return filter(LaplaceSharpening_4());
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Laplace_8_Sharpening() {
	return filter(LaplaceSharpening_8());
}

template<typename T>
GD::RasterPtr GD::Raster<T>::Sobel_Sharpening() {
	size_t total = lWidth * lHeight * nBands;
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	GD::RasterPtr Sobelx(filter(SobelKernelx()));
	GD::RasterPtr Sobely(filter(SobelKernely()));
	for (long long i = 0; i < total; ++i) {
		T x = ((T*)Sobelx->data())[i];
		T y = ((T*)Sobely->data())[i];
		double s = sqrt(x * x + y * y) + pData[i];
		if (s < 0) s = fabs(s);
		if (s > GD::Type_Max<T>()) s = GD::Type_Max<T>();
		((T*)ret->data())[i] = (T)(s);
	}
	return ret;
}

template<typename T>
GD::RasterPtr GD::Raster<T>::toGray() {
	GD::RasterPtr ret = std::make_shared<GD::Raster<T>>();
	ret->init(nullptr, lWidth, lHeight, 1, datatype);
	std::vector<float> alphas(nBands);
	for (int i = 0; i < nBands; ++i) {
		alphas[i] = 1.f / nBands;
	}
	for (size_t i = 0; i < area; ++i) {
		double gray{ 0.0 };
		for (long long channel = 0; channel < nBands; ++channel) {
			long long index = i + channel * area;
			gray += ((T*)pData)[index] * alphas[channel];
		}
		gray += 0.5;
		if (gray > 255) gray = 255;
		else if (gray < 0) gray = 0;
		((T*)ret->data())[i] = (T)gray;
	}
	return ret;
}


GD::RasterPtr GD::getByType(GDALDataType datatype) {
	GD::RasterPtr ret = nullptr;
	switch (datatype) {
	case GDT_Byte:
		ret = std::make_shared<GD::Raster<unsigned char>>();
		break;
	case GDT_UInt16:
		ret = std::make_shared<GD::Raster<unsigned short>>();
		break;
	case GDT_Int16:
		ret = std::make_shared<GD::Raster<short>>();
		break;
	case GDT_UInt32:
		ret = std::make_shared<GD::Raster<unsigned int>>();
		break;
	case GDT_Int32:
		ret = std::make_shared<GD::Raster<int>>();
		break;
	case GDT_Float32:
		ret = std::make_shared<GD::Raster<float>>();
		break;
	case GDT_Float64:
		ret = std::make_shared<GD::Raster<double>>();
		break;
	}
	assert(ret != nullptr);
	return ret;
}

std::shared_ptr<GD::RasterBase> GD::open(const std::string& filepath) {
	GDALDataset* pDatasetRead{ nullptr };
	size_t lWidth{ 0 }, lHeight{ 0 }, nBands{ 0 };
	GDALDataType datatype{ GDT_Unknown };
	try {
		pDatasetRead = (GDALDataset*)GDALOpen(filepath.c_str(), GA_ReadOnly);
		lWidth = pDatasetRead->GetRasterXSize();//获取图像宽度
		lHeight = pDatasetRead->GetRasterYSize();//获取图像高度
		nBands = pDatasetRead->GetRasterCount();//获取图像波段数
		datatype = pDatasetRead->GetRasterBand(1)->GetRasterDataType();//获取图像格式类型
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		exit(1);
	}
	GD::RasterPtr ptr(GD::getByType(datatype));
	ptr->init(pDatasetRead, lWidth, lHeight, nBands, datatype);
	return ptr;
}

GD::RasterPtr GD::operator-(GD::RasterBase const& lhs, GD::RasterBase const& rhs) {
	assert(lhs.width() == rhs.width() && lhs.height() == rhs.height() && lhs.bands() == rhs.bands());
	GD::RasterPtr ret{ new Raster<unsigned char>(lhs.width(), lhs.height(), lhs.bands(), lhs.type()) };
	size_t total = lhs.width() * lhs.height() * lhs.bands();
	for (size_t i = 0; i < total; ++i) {
		double diff = ((unsigned char*)lhs.data())[i] - ((unsigned char*)rhs.data())[i];
		if (diff < 0) diff = 0;
		((unsigned char*)ret->data())[i] = (unsigned char)diff;
	}
	return ret;
}

template<typename T>
void GD::Raster<T>::resize(size_t newWidth, size_t newHeight) {
	size_t width = this->lWidth, height = this->lHeight, nBands = this->nBands;
	if (newWidth == width && newHeight == height) return;
	T* newdata = Malloc<T>(width, height, nBands);
	for (size_t i = 0; i < newHeight; ++i) {
		for (size_t j = 0; j < newWidth; ++j) {
			size_t oldi = i * height / newHeight;
			size_t oldj = j * width / newWidth;
			for (size_t channel = 0; channel < nBands; ++channel) {
				newdata[i * newWidth + j + channel * newWidth * newHeight] =
					pData[oldi * width + oldj + channel * width * height];
			}
		}
	}
	delete[] this->pData;
	this->pData = newdata;
}

template<typename T>
void GD::Raster<T>::set(size_t x, size_t y, GD::Array_features const& arr) {
	if (arr.size() < bands()) {
		std::cerr << "Array length ERROR!" << std::endl;
		return;
	}
	for (size_t channel = 0; channel < bands(); ++channel) {
		set(x, y, channel, arr[channel]);
	}
}

template<typename T>
void GD::Raster<T>::set(size_t x, size_t y, Eigen::Vector3f const& arr) {
	if (arr.size() < nBands) {
		std::cerr << "Array length ERROR!" << std::endl;
		return;
	}
	for (size_t channel = 0; channel < bands(); ++channel) {
		set(x, y, channel, arr[channel]);
	}
}

template<typename T>
void GD::Raster<T>::set(size_t i, size_t j, size_t channel, double value) {
	// if (value < 0.f || value > 255.f) return;
	pData[POS(j, i, channel, lWidth, area)] = static_cast<T>(value);
}
template<typename T>
void GD::Raster<T>::add(size_t i, size_t j, size_t channel, double const value) {
	float res = pData[POS(j, i, channel, lWidth, area)] + static_cast<T>(value);
	set(i, j, channel, res);
}

template<typename T>
void GD::Raster<T>::add(size_t i, size_t j, double const value) {
	for (size_t channel = 0; channel < bands(); ++channel) {
		float res = pData[POS(j, i, channel, lWidth, area)] + static_cast<T>(value);
		set(i, j, channel, res);
	}
}

template<typename T>
void GD::Raster<T>::add(size_t i, size_t j, size_t channel, GD::Array_features const& arr) {
	if (arr.size() <= bands()) {
		std::cerr << "Array length ERROR!" << std::endl;
		return;
	}
	for (size_t channel = 0; channel < bands(); ++channel) {
		pData[POS(j, i, channel, lWidth, area)] += static_cast<T>(arr[channel]);
	}
}

template<typename T>
GD::RasterPtr GD::Raster<T>::copy()const {
	GD::RasterPtr ret{ new GD::Raster<T> };
	ret->init(nullptr, lWidth, lHeight, nBands, datatype);
	size_t total = lHeight * lWidth * nBands;
	memcpy(ret->data(), this->data(), total * sizeof(T));
	return ret;
}

template<typename T>
size_t GD::Raster<T>::Area()const {
	return this->area;
}

template<typename T>
void GD::Raster<T>::setnullptr()
{
	pData = nullptr;
}

template<typename T>
size_t GD::Raster<T>::type_size() const {
	return sizeof(T);
}


void GD::HistgramMatch(GD::RasterPtr src, GD::RasterPtr dst) {
	// src -> dst
	assert(src->bands() == dst->bands());
	auto srchist = src->histgram();
	auto dsthist = dst->histgram();
	for (int channel = 0; channel < src->bands(); ++channel) {
		for (int i = 1; i < 256; ++i) {
			srchist[channel][i] += srchist[channel][i - 1];
			dsthist[channel][i] += dsthist[channel][i - 1];
		}
	}
	size_t area = src->Area();
	for (size_t i = 0; i < area; ++i) {
		auto p = src->get(i % src->width(), i / src->width());
		for (size_t channel = 0; channel < src->bands(); ++channel) {
			p[channel] = 255.0 * dsthist[channel][p[channel]] / srchist[channel][255];
		}
		src->set(i % src->width(), i / src->width(), p);
	}
}

template<typename T>
GD::Raster<T>::Raster() {
	pData = nullptr;
	this->lHeight = -1;
	this->lWidth = -1;
	pDatasetRead = nullptr;
	pDriver = nullptr;
	this->area = 0;
	this->datatype = GDT_Byte;
	this->nBands = 0;
}

template<typename T>
GD::Raster<T>::Raster(size_t width, size_t height, size_t nBands, GDALDataType datatype) {
	pData = nullptr;
	this->init(nullptr, width, height, nBands, datatype);
	this->area = width * height;
}

template<typename T>
GD::Raster<T>::~Raster() {
	if (pDatasetRead)
		GDALClose(pDatasetRead);
	if (pData)
		delete[] pData;
}

const Eigen::Matrix3f T = (
	Eigen::Matrix3f() <<
	1 / 3.f, 1 / 3.f, 1 / 3.f,
	-sqrt(2) / 6, -sqrt(2) / 6, sqrt(2) / 3,
	1 / sqrt(2), -1 / sqrt(2), 0
	).finished();

const Eigen::Matrix3f T_inverse = T.inverse();

void GD::RGB2IHS(Eigen::Vector3f const& rgb, Eigen::Vector3f& ihs) {
	ihs = T * rgb;
}

void GD::IHS2RGB(Eigen::Vector3f const& ihs, Eigen::Vector3f& rgb) {
	rgb = T_inverse * ihs;
}