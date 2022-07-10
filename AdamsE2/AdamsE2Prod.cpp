#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "groebner_steenrod_const.h"

#include <sstream>
#include <cstring>

// std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
int main(int argc, char** argv)
{
    int t_max = 392;
    std::string db_in = "AdamsE2.db";
    std::string db_out = "AdamsE2Prod.db";

    if (argc >= 2 && strcmp(argv[1], "-h") == 0) {
        std::cout << "Usage:\n  AdamsE2Export <t_max> <db_in> <db_out>\n\n";
        std::cout << "Default values:\n  t_max = " << t_max << "\n  db_in = " << db_in << "\n  db_out = " << db_out << "\n\n";
        std::cout << "Version:\n  1.0 (2022-7-7)" << std::endl;
        return 0;
    }

    if (argc >= 2) {
        std::istringstream ss(argv[1]);
        if (!(ss >> t_max) || t_max < 0) {
            std::cerr << "Invalid t_max: " << argv[1] << '\n';
            return 1;
        }
    }

    if (argc >= 3) {
        std::istringstream ss(argv[2]);
        if (!(ss >> db_in)) {
            std::cerr << "Invalid db_in: " << argv[2] << '\n';
            return 2;
        }
    }

    if (argc >= 4) {
        std::istringstream ss(argv[3]);
        if (!(ss >> db_out)) {
            std::cerr << "Invalid db_out: " << argv[3] << '\n';
            return 3;
        }
    }

    using namespace steenrod;
    bench::Timer timer;


    compute_products_ind(t_max, db_in, db_out);
    compute_products(t_max, db_in, db_out);

    // bench::Counter::print();
    return 0;
}