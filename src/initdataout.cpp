#include "kadath_polar.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Kadath;

int main(int argc, char** argv) {

    // 1. 只需要一个参数：Kadath 的二进制结果文件
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.bin\n";
        return 1;
    }

    // 2. 打开输入文件
    const char* input_name = argv[1];
    std::cerr << "Attempting to read data from file: " << input_name << std::endl;

    FILE* fin = fopen(input_name, "r");
    if (!fin) {
        std::cerr << "Cannot open input file " << input_name << std::endl;
        return 1;
    }

    // 3. 读入空间与物理场（按 Kadath 的标准顺序）
    Space_polar space(fin);
    Scalar rsint(space, fin);
    Scalar phi(space,fin);

    fclose(fin);

    std::cerr << "Transformations complete.\n";

    // 5. 打开二维输出文本文件
    const char* output_name = "initfields_2d.txt";
    std::ofstream fout(output_name);
    if (!fout) {
        std::cerr << "Cannot open output file " << output_name << std::endl;
        return 1;
    }

    // 6. 在 (r, theta) 网格上采样
    // 假设 Point(2): 第 1 个分量是 r，第 2 个是 theta
    int Nr  = 500;        // 径向采样点数
    int Nth = 500;        // theta 采样点数
    double rmax = 10.0;   // 径向最大半径，可根据解调整

    Point P(2);

    for (int it = 0; it <= Nth; ++it) {
        double theta = M_PI * it / Nth;  // theta ∈ [0, pi]
        P.set(2) = theta;

        for (int ir = 0; ir <= Nr; ++ir) {
            double r = rmax * ir / Nr;   // r ∈ [0, rmax]
            P.set(1) = r;

            double phival   = phi.val_point(P);
            double rsintval = rsint.val_point(P);
 

            // 一行：r theta phi lapse A B bt
            fout << r      << " "
                 << theta  << " "
                 << rsintval << " "
                 << phival << "\n ";
                 
        }

        // 用空行分隔不同的 theta 扫描（可选，方便某些画图工具）
        fout << "\n";
    }

    fout.close();

    std::cerr << "2D data written to " << output_name << std::endl;
    return 0;
}
