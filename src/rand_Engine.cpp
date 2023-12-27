#include "rand_Engine.h"

float GD::UniformGenerator::Get(int dim, int threadnum) {
    return ((float)rand()) / RAND_MAX;
}

GD::LdsGenerator::LdsGenerator(int pixelNum, int dim) {
    Build(pixelNum, dim);
}

void GD::LdsGenerator::Build(int pixelNum, int dim) {
    uint32_t bit = 1;
    mBits = std::vector<uint32_t>(32);
    for (int i = 0; i < 32; i++) {
        mBits[i] = bit;
        bit = bit << 1;
    }

    int idxOffs = 1;
    mDirNumV_.resize(dim);

    // �����0��ά������m��Ϊ1
    std::vector<uint32_t> arrV(32);
    for (int k = 1; k <= 32; k++) {
        arrV[k - idxOffs] = 1u << (32 - k);
    }
    mDirNumV_[0] = arrV;

    // ��������ά��
    for (int d = 1; d < dim; d++) {
        int sj = mMSet[d].size();
        // ��ȫm
        std::vector<uint32_t> arrM = mMSet[d];
        arrM.resize(32);

        for (int k = sj + 1; k <= 32; k++) {        // k~[sj + 1,32]
            uint32_t mk = arrM[k - sj - idxOffs];            // ���һ��m[k-sj]
            mk ^= ((1u << sj) * arrM[k - sj - idxOffs]);    // �����ڶ���(2^sj)*m[k-sj]
            for (int i = 1; i < sj; i++) {        // i~[1, sj-1]i - idxOffs
                mk ^= ((1u << i) * arrM[k - i - idxOffs] * (mASet[d] & mBits[sj - i - idxOffs]));    // ������
            }
            arrM[k - idxOffs] = mk;
        }
        // ��v
        std::vector<uint32_t> arrV(32);
        for (int k = 1; k <= 32; k++) {
            arrV[k - idxOffs] = arrM[k - idxOffs] << (32 - k);
        }
        mDirNumV_[d] = arrV;
    }

    // ��ʼ�����������ͼ�¼
    mIndex.resize(dim);
    mDigit32.resize(dim);
    for (int i = 0; i < dim; i++) {
        mIndex[i] = std::vector<uint32_t>(pixelNum, 0);
        mDigit32[i] = std::vector<uint32_t>(pixelNum, 0);
    }
}

float GD::LdsGenerator::Get(int dim, int pixelNum) {
    int idxOffs = 1;
    uint32_t i = mIndex[dim][pixelNum];
    uint32_t ci = 1u;   // �ҵ������һ����0��λ
    while ((i & 1u) != 0u) {
        i = i >> 1u;
        ci++;
    }
    mIndex[dim][pixelNum]++;
    mDigit32[dim][pixelNum] = (mDigit32[dim][pixelNum] ^ mDirNumV_[dim][ci - idxOffs]);
    return float(mDigit32[dim][pixelNum]) * mFactor;
}

void GD::LdsGenerator::Reset() {

}

void GD::LdsGenerator::Reset(int dim, int thread) {
    mIndex[dim][thread] = 0;
    mDigit32[dim][thread] = 0;
}