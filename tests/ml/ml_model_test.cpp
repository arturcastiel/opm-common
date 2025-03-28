/*
  Copyright (c) 2016 Robert W. Rose
  Copyright (c) 2018 Paul Maevskikh
  Copyright (c) 2024 NORCE
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/common/ErrorMacros.hpp>
#include <opm/ml/ml_model.hpp>

#include <tests/ml/ml_tools/include/test_dense_10x1.hpp>
#include <tests/ml/ml_tools/include/test_dense_10x10.hpp>
#include <tests/ml/ml_tools/include/test_dense_10x10x10.hpp>
#include <tests/ml/ml_tools/include/test_dense_1x1.hpp>
#include <tests/ml/ml_tools/include/test_dense_2x2.hpp>
#include <tests/ml/ml_tools/include/test_dense_relu_10.hpp>
#include <tests/ml/ml_tools/include/test_dense_tanh_10.hpp>
#include <tests/ml/ml_tools/include/test_relu_10.hpp>
#include <tests/ml/ml_tools/include/test_scalingdense_10x1.hpp>

#include <cstdio>

#include <fmt/format.h>

namespace Opm {

using namespace ML;

template <class Evaluation>
bool tensor_test()
{
    std::printf("TEST tensor_test\n");

    {
        const int i = 3;
        const int j = 5;
        const int k = 10;
        Tensor<Evaluation> t(i, j, k);

        Evaluation c = 1.f;
        for (int ii = 0; ii < i; ii++) {
            for (int jj = 0; jj < j; jj++) {
                for (int kk = 0; kk < k; kk++) {
                    t(ii, jj, kk) = c;
                    c += 1.f;
                }
            }
        }

        c = 1.f;
        int cc = 0;
        for (int ii = 0; ii < i; ii++) {
            for (int jj = 0; jj < j; jj++) {
                for (int kk = 0; kk < k; kk++) {
                    OPM_ERROR_IF(fabs(t(ii, jj, kk).value() - c.value()) > 1e-9,
                                 fmt::format("\n Expected "
                                             "{}"
                                             "got "
                                             "{}",
                                             c.value(),
                                             t(ii, jj, kk).value()));
                    OPM_ERROR_IF(fabs(t.data_[cc].value() - c.value()) > 1e-9,
                                 fmt::format("\n Expected "
                                             "{}"
                                             "got "
                                             "{}",
                                             c.value(),
                                             t.data_[cc].value()));

                    c += 1.f;
                    cc++;
                }
            }
        }
    }

    {
        const int i = 2;
        const int j = 3;
        const int k = 4;
        const int l = 5;
        Tensor<Evaluation> t(i, j, k, l);

        Evaluation c = 1.f;
        for (int ii = 0; ii < i; ii++) {
            for (int jj = 0; jj < j; jj++) {
                for (int kk = 0; kk < k; kk++) {
                    for (int ll = 0; ll < l; ll++) {
                        t(ii, jj, kk, ll) = c;
                        c += 1.f;
                    }
                }
            }
        }

        c = 1.f;
        int cc = 0;
        for (int ii = 0; ii < i; ii++) {
            for (int jj = 0; jj < j; jj++) {
                for (int kk = 0; kk < k; kk++) {
                    for (int ll = 0; ll < l; ll++) {
                        OPM_ERROR_IF(fabs(t(ii, jj, kk, ll).value() - c.value()) > 1e-9,
                                     fmt::format("\n Expected "
                                                 "{}"
                                                 "got "
                                                 "{}",
                                                 c.value(),
                                                 t(ii, jj, kk, ll).value()));
                        OPM_ERROR_IF(fabs(t.data_[cc].value() - c.value()) > 1e-9,
                                     fmt::format("\n Expected "
                                                 "{}"
                                                 "got "
                                                 "{}",
                                                 c.value(),
                                                 t.data_[cc].value()));
                        c += 1.f;
                        cc++;
                    }
                }
            }
        }
    }

    {
        Tensor<Evaluation> a(2, 2);
        Tensor<Evaluation> b(2, 2);

        a.data_ = {1.0, 2.0, 3.0, 5.0};
        b.data_ = {2.0, 5.0, 4.0, 1.0};

        Tensor<Evaluation> result = a + b;
        OPM_ERROR_IF(result.data_ != std::vector<Evaluation>({3.0, 7.0, 7.0, 6.0}),
                     "Vector add failed");
    }

    {
        Tensor<Evaluation> a(2, 2);
        Tensor<Evaluation> b(2, 2);

        a.data_ = {1.0, 2.0, 3.0, 5.0};
        b.data_ = {2.0, 5.0, 4.0, 1.0};

        Tensor<Evaluation> result = a.multiply(b);
        OPM_ERROR_IF(result.data_ != std::vector<Evaluation>({2.0, 10.0, 12.0, 5.0}),
                     "Vector multiply failed");
    }

    {
        Tensor<Evaluation> a(2, 1);
        Tensor<Evaluation> b(1, 2);

        a.data_ = {1.0, 2.0};
        b.data_ = {2.0, 5.0};

        Tensor<Evaluation> result = a.dot(b);
        OPM_ERROR_IF(result.data_ != std::vector<Evaluation>({2.0, 5.0, 4.0, 10.0}),
                     "Vector dot failed");
    }

    return true;
}

} // namespace Opm

int main()
{
    using Evaluation = Opm::DenseAd::Evaluation<double, 1>;

    Evaluation load_time = 0.0;
    Evaluation apply_time = 0.0;

    try {
        tensor_test<Evaluation>();
        test_dense_1x1<Evaluation>(&load_time, &apply_time);
        test_dense_10x1<Evaluation>(&load_time, &apply_time);
        test_dense_2x2<Evaluation>(&load_time, &apply_time);
        test_dense_10x10<Evaluation>(&load_time, &apply_time);
        test_dense_10x10x10<Evaluation>(&load_time, &apply_time);
        test_relu_10<Evaluation>(&load_time, &apply_time);
        test_dense_relu_10<Evaluation>(&load_time, &apply_time);
        test_dense_tanh_10<Evaluation>(&load_time, &apply_time);
        test_scalingdense_10x1<Evaluation>(&load_time, &apply_time);
    }
    catch(...) {
        return 1;
    }

    return 0;
}
