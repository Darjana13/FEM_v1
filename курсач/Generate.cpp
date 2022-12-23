#include "Generate.h"

void Make_grid(std::string path)
{
    std::ofstream out;
    out.precision(15);  
    std::vector<double> all_R, all_Z;
    std::ifstream in(path + "grid.txt");
    double R, Z, kr, kz;
    int Nr, Nz;
    int count_r, count_z;
    in >> count_r >> count_z;
    all_R.resize(count_r);
    all_Z.resize(count_z);

    in >> all_R[0] >> all_Z[0];
    for (int curr_count_r = 0; curr_count_r < count_r - 1; )
    {
        in >> R >> Nr >> kr;
        double hx;
        if (kr == 1)
        {
            hx = (R - all_R[curr_count_r]) / Nr;
            for (int p = 1; p < Nr; p++)
            {
                all_R[curr_count_r + p] = all_R[curr_count_r] + hx * p;
            }
            curr_count_r += Nr;
        }
        else
        {
            hx = (R - all_R[curr_count_r]) * (kr - 1) / (pow(kr, Nr) - 1);
            for (int p = 0; p < Nr - 1; curr_count_r++, p++)
            {
                all_R[curr_count_r + 1] = all_R[curr_count_r] + hx * pow(kr, p);
            }
            curr_count_r++;
        }
        all_R[curr_count_r] = R;
    }
    for (int curr_count_z = 0; curr_count_z < count_z - 1; )
    {
        in >> Z >> Nz >> kz;
        double hy;
        if (kz == 1)
        {
            hy = (Z - all_Z[curr_count_z]) / Nz;
            for (int p = 1; p < Nz; p++)
            {
                all_Z[curr_count_z + p] = all_Z[curr_count_z] + hy * p;
            }
            curr_count_z += Nz;
        }
        else
        {
            hy = (Z - all_Z[curr_count_z]) * (kz - 1) / (pow(kz, Nz) - 1);
            for (int p = 0; p < Nz - 1; curr_count_z++, p++)
            {
                all_Z[curr_count_z + 1] = all_Z[curr_count_z] + hy * pow(kz, p);
            }
            curr_count_z++;
        }
        all_Z[curr_count_z] = Z;
    }
    in.close();
    out.open("rz.txt");
    for (int i = 0; i < count_z; i++)
    {
        for (int j = 0; j < count_r; j++)
        {
            out << all_R[j] << "\t" << all_Z[i] << "\n";
        }
    }
    out.close();

    // input area
    // lambda, sigma same in all area
    out.open("elem.txt");
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_r - 1; j++)
        {
            out << i * count_r + j << " " << i * count_r + j + 1 << " "
                << (i+1) * count_r + j << " " << (i + 1) * count_r + j + 1
                << " 0 0 \n";
        }
    }
    out.close();

    // bounders
    std::vector<int> bottom(count_r), right(count_z), left(count_z), top(count_r);

    // bottom
    for (int j = 0; j < count_r; j++)
    {
        bottom[j] = j;
    }
    // right
    for (int i = 1; i < count_z + 1; i++)
    {
        right[i-1] = i * count_r - 1;       
    }
    // left
    for (int i = 0; i < count_z; i++)
    {
        left[i] = i * count_r;
    }
    // top
    for (int j = 0; j < count_r; j++)
    {
        top[j] = (count_z - 1) * count_r + j;
    }

    // справа вторые, по остальным - первые
    out.open("S1.txt");
    out << top.size() + left.size() + bottom.size() - 2 << " 0\n";
    for (int j = 0; j < top.size(); j++)
    {
        out << top[j]<< " ";
    }
    for (int j = 0; j < left.size() - 1; j++)
    {
        out << left[j] << " ";
    }
    for (int j = 1; j < bottom.size(); j++)
    {
        out << bottom[j] << " ";
    }
    out.close();
    out.open("S2_r.txt");
    out << 0;
    out.close();
    out.open("S2_z.txt");
    out << 1 << "\n";
    out << right.size() << " 1\n";
    for (int j = 0; j < right.size(); j++)
    {
        out << right[j] << " ";
    }
    out.close();
    out.open("S3_r.txt");
    out << 0;
    out.close();
    out.open("S3_z.txt");
    out << 0;
    out.close();

    // везде первые
    /*out.open("S1.txt");
    out << top.size() + left.size() + bottom.size() + right.size() - 4 << " 0\n";
    for (int j = 0; j < top.size(); j++)
    {
        out << top[j] << " ";
    }
    for (int j = 0; j < left.size() - 1; j++)
    {
        out << left[j] << " ";
    }
    for (int j = 1; j < bottom.size(); j++)
    {
        out << bottom[j] << " ";
    }
    for (int j = 1; j < right.size() - 1; j++)
    {
        out << right[j] << " ";
    }
    out.close();
    out.open("S2_r.txt");
    out << 0;
    out.close();
    out.open("S2_z.txt");
    out << 0;   
    out.close();
    out.open("S3_r.txt");
    out << 0;
    out.close();
    out.open("S3_z.txt");
    out << 0;
    out.close();*/
}

void Make_grid2(std::string path)
{
    std::ofstream out;
    out.precision(15);
    std::vector<double> all_R, all_Z;
    std::ifstream in(path + "grid.txt");
    double R, Z, kr, kz;
    int Nr, Nz;
    int count_r, count_z;
    in >> count_r >> count_z;
    all_R.resize(count_r);
    all_Z.resize(count_z);

    in >> all_R[0] >> all_Z[0];
    for (int curr_count_r = 0; curr_count_r < count_r - 1; )
    {
        in >> R >> Nr >> kr;
        double hx;
        if (kr == 1)
        {
            hx = (R - all_R[curr_count_r]) / Nr;
            for (int p = 1; p < Nr; p++)
            {
                all_R[curr_count_r + p] = all_R[curr_count_r] + hx * p;
            }
            curr_count_r += Nr;
        }
        else
        {
            hx = (R - all_R[curr_count_r]) * (kr - 1) / (pow(kr, Nr) - 1);
            for (int p = 0; p < Nr - 1; curr_count_r++, p++)
            {
                all_R[curr_count_r + 1] = all_R[curr_count_r] + hx * pow(kr, p);
            }
            curr_count_r++;
        }
        all_R[curr_count_r] = R;
    }
    for (int curr_count_z = 0; curr_count_z < count_z - 1; )
    {
        in >> Z >> Nz >> kz;
        double hy;
        if (kz == 1)
        {
            hy = (Z - all_Z[curr_count_z]) / Nz;
            for (int p = 1; p < Nz; p++)
            {
                all_Z[curr_count_z + p] = all_Z[curr_count_z] + hy * p;
            }
            curr_count_z += Nz;
        }
        else
        {
            hy = (Z - all_Z[curr_count_z]) * (kz - 1) / (pow(kz, Nz) - 1);
            for (int p = 0; p < Nz - 1; curr_count_z++, p++)
            {
                all_Z[curr_count_z + 1] = all_Z[curr_count_z] + hy * pow(kz, p);
            }
            curr_count_z++;
        }
        all_Z[curr_count_z] = Z;
    }
    in.close();
    out.open("rz.txt");
    out << count_z * count_r << "\n";
    for (int i = 0; i < count_z; i++)
    {
        for (int j = 0; j < count_r; j++)
        {
            out << all_R[j] << "\t" << all_Z[i] << "\n";
        }
    }
    out.close();

    // input area
    // lambda, sigma same in all area
    in.open("area.txt");
    double water_bottom, water_right;
    double water_lambda, water_gamma;
    double pot_lambda, pot_gamma;
    double r_center, z_center;

    in >>  water_right >> water_bottom;
    in >> water_lambda >> water_gamma;
    in >> pot_lambda >> pot_gamma;
    in.close();

    out.open("elem.txt");
    out << (count_z - 1) * (count_r - 1) << "\n";

    for (int i = 0; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_r - 1; j++)
        {
            out << i * count_r + j << " " << i * count_r + j + 1 << " "
                << (i + 1) * count_r + j << " " << (i + 1) * count_r + j + 1;

            r_center = (all_R[j] + all_R[j + 1]) / 2;
            z_center = (all_Z[i] + all_Z[i + 1]) / 2;

            if(r_center >= water_right || z_center <= water_bottom)
                out << " 1 0\n";
            else
                out << " 0 0\n";
        }
    }
    out.close();

    // bounders
    std::vector<int> bottom(count_r), right(count_z), left(count_z), top(count_r);

    // bottom
    for (int j = 0; j < count_r; j++)
    {
        bottom[j] = j;
    }
    // right
    for (int i = 1; i < count_z + 1; i++)
    {
        right[i - 1] = i * count_r - 1;
    }
    // left
    for (int i = 0; i < count_z; i++)
    {
        left[i] = i * count_r;
    }
    // top
    for (int j = 0; j < count_r; j++)
    {
        top[j] = (count_z - 1) * count_r + j;
    }

    // справа вторые, по остальным - первые
    out.open("S1.txt");
    out << 0;
    out.close();
    out.open("S2_z.txt");
    out << 0;
    out.close();
    out.open("S2_r.txt");
    out << 1 << "\n";
    out << bottom.size() << " 0\n";
    for (int j = 0; j < bottom.size(); j++)
    {
        out << bottom[j] << " ";
    }
    out.close();
    out.open("S3_r.txt");
    out << 0 << "\n";

    //out << 1 << "\n";
    out << right.size() << " 1\n";
    //for (int j = 0; j < left.size(); j++) // null 2nd
    //{
    //    out << left[j] << " ";
    //}
    for (int j = 0; j < right.size(); j++)
    {
        out << right[j] << " ";
    }
    out.close();
    out.open("S3_z.txt");
    out << 0 << "\n";

    //out << 1 << '\n';
    out << top.size() << " 2\n";
    for (int j = 0; j < top.size(); j++)
    {
        out << top[j] << " ";
    }
    out.close();

    // везде первые
    /*out.open("S1.txt");
    out << 1 << "\n";
    out << top.size() + left.size() + bottom.size() + right.size() - 4 << " 1\n";
    for (int j = 0; j < top.size(); j++)
    {
        out << top[j] << " ";
    }
    for (int j = 0; j < left.size() - 1; j++)
    {
        out << left[j] << " ";
    }
    for (int j = 1; j < bottom.size(); j++)
    {
        out << bottom[j] << " ";
    }
    for (int j = 1; j < right.size() - 1; j++)
    {
        out << right[j] << " ";
    }
    out.close();
    out.open("S2_r.txt");
    out << 0;
    out.close();
    out.open("S2_z.txt");
    out << 0;
    out.close();
    out.open("S3_r.txt");
    out << 0;
    out.close();
    out.open("S3_z.txt");
    out << 0;
    out.close();*/
}


void Create_time_grid()
{
    std::ofstream out;
    out.precision(15);
    std::ifstream in;
    // time grid
    in.open("time_grid.txt");
    std::vector<double> time_grid;
    double T, kt;
    int Nt;
    int count_t;
    in >> count_t;
    time_grid.resize(count_t);

    in >> time_grid[0];
    for (int curr_count_t = 0; curr_count_t < count_t - 1; )
    {
        in >> T >> Nt >> kt;
        double ht;
        if (kt == 1)
        {
            ht = (T - time_grid[curr_count_t]) / Nt;
            for (int p = 1; p < Nt; p++)
            {
                time_grid[curr_count_t + p] = time_grid[curr_count_t] + ht * p;
            }
            curr_count_t += Nt;
        }
        else
        {
            ht = (T - time_grid[curr_count_t]) * (kt - 1) / (pow(kt, Nt) - 1);
            double pow_kt = 1;
            for (int p = 0; p < Nt - 1; curr_count_t++, p++)
            {
                time_grid[curr_count_t + 1] = time_grid[curr_count_t] + ht * pow_kt;
                pow_kt *= kt;
            }
            curr_count_t++;
        }
        time_grid[curr_count_t] = T;
    }
    in.close();
    out.open("time.txt");
    out << time_grid.size() << "\n";
    for (int i = 0; i < count_t; i++)
    {
        out << time_grid[i] << " ";
    }
    out.close();
}

void Make_triangle_grid(std::string path)
{
    std::ofstream out;
    out.precision(15);
    std::vector<double> all_R, all_Z;
    std::ifstream in(path + "grid.txt");
    double R, Z, kr, kz;
    int Nr, Nz;
    int count_r, count_z;
    in >> count_r >> count_z;
    all_R.resize(count_r);
    all_Z.resize(count_z);

    in >> all_R[0] >> all_Z[0];
    for (int curr_count_r = 0; curr_count_r < count_r - 1; )
    {
        in >> R >> Nr >> kr;
        double hx;
        if (kr == 1)
        {
            hx = (R - all_R[curr_count_r]) / Nr;
            for (int p = 1; p < Nr; p++)
            {
                all_R[curr_count_r + p] = all_R[curr_count_r] + hx * p;
            }
            curr_count_r += Nr;
        }
        else
        {
            hx = (R - all_R[curr_count_r]) * (kr - 1) / (pow(kr, Nr) - 1);
            for (int p = 0; p < Nr - 1; curr_count_r++, p++)
            {
                all_R[curr_count_r + 1] = all_R[curr_count_r] + hx * pow(kr, p);
            }
            curr_count_r++;
        }
        all_R[curr_count_r] = R;
    }
    for (int curr_count_z = 0; curr_count_z < count_z - 1; )
    {
        in >> Z >> Nz >> kz;
        double hy;
        if (kz == 1)
        {
            hy = (Z - all_Z[curr_count_z]) / Nz;
            for (int p = 1; p < Nz; p++)
            {
                all_Z[curr_count_z + p] = all_Z[curr_count_z] + hy * p;
            }
            curr_count_z += Nz;
        }
        else
        {
            hy = (Z - all_Z[curr_count_z]) * (kz - 1) / (pow(kz, Nz) - 1);
            for (int p = 0; p < Nz - 1; curr_count_z++, p++)
            {
                all_Z[curr_count_z + 1] = all_Z[curr_count_z] + hy * pow(kz, p);
            }
            curr_count_z++;
        }
        all_Z[curr_count_z] = Z;
    }
    in.close();
    out.open("rz.txt");
    out << count_z * count_r;
    for (int i = 0; i < count_z; i++)
    {
        for (int j = 0; j < count_r; j++)
        {
            out << all_R[j] << "\t" << all_Z[i] << "\n";
        }
    }
    out.close();

    // input area
    // lambda, sigma same in all area
    out.open("elem.txt");
    out << 2 * (count_z - 1) * (count_r - 1);
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int j = 0; j < count_r - 1; j++)
        {
            out << i * count_r + j << " " << i * count_r + j + 1 << " " << (i + 1) * count_r + j + 1 << " 0 \n";
            out << i * count_r + j << " " << (i + 1) * count_r + j << " " << (i + 1) * count_r + j + 1 << " 0 \n";
        }
    }
    out.close();

    // bounders
    std::vector<int> bottom(count_r), right(count_z), left(count_z), top(count_r);

    // bottom
    for (int j = 0; j < count_r; j++)
    {
        bottom[j] = j;
    }
    // right
    for (int i = 1; i < count_z + 1; i++)
    {
        right[i - 1] = i * count_r - 1;
    }
    // left
    for (int i = 0; i < count_z; i++)
    {
        left[i] = i * count_r;
    }
    // top
    for (int j = 0; j < count_r; j++)
    {
        top[j] = (count_z - 1) * count_r + j;
    }

    // справа вторые, по остальным - первые
    /*out.open("S1.txt");
    out << top.size() + left.size() + bottom.size() - 2 << " 0\n";
    for (int j = 0; j < top.size(); j++)
    {
        out << top[j] << " ";
    }
    for (int j = 0; j < left.size() - 1; j++)
    {
        out << left[j] << " ";
    }
    for (int j = 1; j < bottom.size(); j++)
    {
        out << bottom[j] << " ";
    }
    out.close();
    out.open("S2_r.txt");
    out << 0;
    out.close();
    out.open("S2_z.txt");
    out << 1 << "\n";
    out << right.size() << " 1\n";
    for (int j = 0; j < right.size(); j++)
    {
        out << right[j] << " ";
    }
    out.close();
    out.open("S3_r.txt");
    out << 0;
    out.close();
    out.open("S3_z.txt");
    out << 0;
    out.close();*/

    // везде первые
    out.open("S1.txt");
    out << top.size() + left.size() + bottom.size() + right.size() - 4 << " 0\n";
    for (int j = 0; j < top.size() - 1; j++)
    {
        out << top[j] << " " << top[j + 1] << "\n";
    }
    for (int j = 0; j < left.size() - 1; j++)
    {
        out << left[j] << " " << left[j + 1] << "\n";
    }
    for (int j = 0; j < bottom.size(); j++)
    {
        out << bottom[j] << " " << bottom[j + 1] << "\n";
    }
    for (int j = 0; j < right.size() - 1; j++)
    {
        out << right[j] << " " << right[j + 1] << "\n";
    }
    out.close();
}