/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __CTP_GAUSS_HERMITE__H
#define	__CTP_GAUSS_HERMITE__H


#include <boost/math/constants/constants.hpp>


namespace votca {
    namespace ctp {
    namespace ub = boost::numeric::ublas;
    class GaussHermite {
        public:

            GaussHermite() {};
            ~GaussHermite(){};

            std::vector< ub::matrix<double> > integrate_shellpair(Shell* shell1,Shell* shell2);
          

        private:
            
            
            struct quadgrid{
                std::vector<double> _points;
                std::vector<double> _weights;
            };
            

            inline void fillgrid(int N){
                if (N == 0) {
                    quadgrid._points[0] = 0.00000000000000e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                } else if (N == 1) {
                    quadgrid._points[0] = -7.07106781186548e-01;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = 7.07106781186547e-01;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                } else if (N == 2) {
                    quadgrid._points[0] = -1.22474487139159e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = 0.00000000000000e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = 1.22474487139159e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                } else if (N == 3) {
                    quadgrid._points[0] = -1.65068012388578e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -5.24647623275290e-01;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = 5.24647623275290e-01;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = 1.65068012388578e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                } else if (N == 4) {
                    quadgrid._points[0] = -2.02018287045608e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -9.58572464613819e-01;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = 0.00000000000000e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = 9.58572464613817e-01;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = 2.02018287045608e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                } else if (N == 5) {
                    quadgrid._points[0] = -2.35060497367449e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -1.33584907401370e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -4.36077411927616e-01;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = 4.36077411927616e-01;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = 1.33584907401369e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 2.35060497367449e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                } else if (N == 6) {
                    quadgrid._points[0] = -2.65196135683523e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -1.67355162876747e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -8.16287882858965e-01;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = 0.00000000000000e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = 8.16287882858964e-01;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 1.67355162876747e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 2.65196135683523e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                } else if (N == 7) {
                    quadgrid._points[0] = -2.93063742025724e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -1.98165675669584e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -1.15719371244678e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -3.81186990207322e-01;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = 3.81186990207321e-01;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 1.15719371244678e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 1.98165675669584e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 2.93063742025724e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                } else if (N == 8) {
                    quadgrid._points[0] = -3.19099320178152e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -2.26658058453184e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -1.46855328921666e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -7.23551018752838e-01;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = 0.00000000000000e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 7.23551018752836e-01;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 1.46855328921666e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 2.26658058453184e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 3.19099320178152e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                } else if (N == 9) {
                    quadgrid._points[0] = -3.43615911883773e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -2.53273167423278e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -1.75668364929988e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -1.03661082978951e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -3.42901327223704e-01;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 3.42901327223704e-01;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 1.03661082978951e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 1.75668364929988e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 2.53273167423278e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 3.43615911883773e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                } else if (N == 10) {
                    quadgrid._points[0] = -3.66847084655957e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -2.78329009978165e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -2.02594801582575e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -1.32655708449493e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -6.56809566882100e-01;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = 0.00000000000000e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 6.56809566882098e-01;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 1.32655708449493e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 2.02594801582575e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 2.78329009978165e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 3.66847084655958e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                } else if (N == 11) {
                    quadgrid._points[0] = -3.88972489786977e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -3.02063702512088e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -2.27950708050106e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -1.59768263515260e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -9.47788391240164e-01;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -3.14240376254359e-01;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 3.14240376254358e-01;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 9.47788391240162e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 1.59768263515260e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 2.27950708050105e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 3.02063702512088e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 3.88972489786978e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                } else if (N == 12) {
                    quadgrid._points[0] = -4.10133759617864e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -3.24660897837240e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -2.51973568567823e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -1.85310765160151e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -1.22005503659075e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -6.05763879171060e-01;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = 0.00000000000000e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 6.05763879171059e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 1.22005503659074e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 1.85310765160151e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 2.51973568567823e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 3.24660897837240e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 4.10133759617863e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                } else if (N == 13) {
                    quadgrid._points[0] = -4.30444857047362e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -3.46265693360226e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -2.74847072498540e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -2.09518325850771e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -1.47668273114114e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -8.78713787329399e-01;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -2.91745510672562e-01;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 2.91745510672561e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 8.78713787329397e-01;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 1.47668273114114e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 2.09518325850771e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 2.74847072498539e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 3.46265693360226e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 4.30444857047362e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                } else if (N == 14) {
                    quadgrid._points[0] = -4.49999070730939e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -3.66995037340445e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -2.96716692790559e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -2.32573248617385e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -1.71999257518649e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -1.13611558521092e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -5.65069583255576e-01;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = 0.00000000000000e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 5.65069583255574e-01;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 1.13611558521092e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 1.71999257518648e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 2.32573248617385e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 2.96716692790560e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 3.66995037340444e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                    quadgrid._points[14] = 4.49999070730938e+00;
                    quadgrid._weights[14] = 1.99532420590461e-02;
                } else if (N == 15) {
                    quadgrid._points[0] = -4.68873893930580e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -3.86944790486012e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -3.17699916197995e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -2.54620215784747e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -1.95178799091625e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -1.38025853919888e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -8.22951449144656e-01;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = -2.73481046138152e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 2.73481046138152e-01;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 8.22951449144654e-01;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 1.38025853919888e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 1.95178799091625e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 2.54620215784747e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 3.17699916197995e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                    quadgrid._points[14] = 3.86944790486011e+00;
                    quadgrid._weights[14] = 1.99532420590461e-02;
                    quadgrid._points[15] = 4.68873893930581e+00;
                    quadgrid._weights[15] = 4.53000990550896e-03;
                } else if (N == 16) {
                    quadgrid._points[0] = -4.87134519367440e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -4.06194667587546e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -3.37893209114149e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -2.75776291570388e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -2.17350282666661e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -1.61292431422123e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -1.06764872574345e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = -5.31633001342655e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = 0.00000000000000e+00;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 5.31633001342653e-01;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 1.06764872574345e+00;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 1.61292431422123e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 2.17350282666661e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 2.75776291570388e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                    quadgrid._points[14] = 3.37893209114148e+00;
                    quadgrid._weights[14] = 1.99532420590461e-02;
                    quadgrid._points[15] = 4.06194667587546e+00;
                    quadgrid._weights[15] = 4.53000990550896e-03;
                    quadgrid._points[16] = 4.87134519367439e+00;
                    quadgrid._weights[16] = 1.57067320322857e-01;
                } else if (N == 17) {
                    quadgrid._points[0] = -5.04836400887446e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -4.24811787356812e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -3.57376906848625e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -2.96137750553160e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -2.38629908916668e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -1.83553160426162e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -1.30092085838962e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = -7.76682919267412e-01;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = -2.58267750519097e-01;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 2.58267750519096e-01;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 7.76682919267409e-01;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 1.30092085838961e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 1.83553160426162e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 2.38629908916668e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                    quadgrid._points[14] = 2.96137750553159e+00;
                    quadgrid._weights[14] = 1.99532420590461e-02;
                    quadgrid._points[15] = 3.57376906848626e+00;
                    quadgrid._weights[15] = 4.53000990550896e-03;
                    quadgrid._points[16] = 4.24811787356811e+00;
                    quadgrid._weights[16] = 1.57067320322857e-01;
                    quadgrid._points[17] = 5.04836400887446e+00;
                    quadgrid._weights[17] = 7.24629595224393e-01;
                } else if (N == 18) {
                    quadgrid._points[0] = -5.22027169053746e+00;
                    quadgrid._weights[0] = 1.77245385090552e+00;
                    quadgrid._points[1] = -4.42853280660377e+00;
                    quadgrid._weights[1] = 8.86226925452758e-01;
                    quadgrid._points[2] = -3.76218735196401e+00;
                    quadgrid._weights[2] = 8.86226925452758e-01;
                    quadgrid._points[3] = -3.15784881834759e+00;
                    quadgrid._weights[3] = 2.95408975150920e-01;
                    quadgrid._points[4] = -2.59113378979453e+00;
                    quadgrid._weights[4] = 1.18163590060368e+00;
                    quadgrid._points[5] = -2.04923170985061e+00;
                    quadgrid._weights[5] = 2.95408975150920e-01;
                    quadgrid._points[6] = -1.52417061939353e+00;
                    quadgrid._weights[6] = 8.13128354472457e-02;
                    quadgrid._points[7] = -1.01036838713431e+00;
                    quadgrid._weights[7] = 8.04914090005513e-01;
                    quadgrid._points[8] = -5.03520163423888e-01;
                    quadgrid._weights[8] = 8.04914090005513e-01;
                    quadgrid._points[9] = 0.00000000000000e+00;
                    quadgrid._weights[9] = 8.13128354472455e-02;
                    quadgrid._points[10] = 5.03520163423886e-01;
                    quadgrid._weights[10] = 1.99532420590462e-02;
                    quadgrid._points[11] = 1.01036838713431e+00;
                    quadgrid._weights[11] = 3.93619323152241e-01;
                    quadgrid._points[12] = 1.52417061939353e+00;
                    quadgrid._weights[12] = 9.45308720482942e-01;
                    quadgrid._points[13] = 2.04923170985061e+00;
                    quadgrid._weights[13] = 3.93619323152243e-01;
                    quadgrid._points[14] = 2.59113378979453e+00;
                    quadgrid._weights[14] = 1.99532420590461e-02;
                    quadgrid._points[15] = 3.15784881834760e+00;
                    quadgrid._weights[15] = 4.53000990550896e-03;
                    quadgrid._points[16] = 3.76218735196401e+00;
                    quadgrid._weights[16] = 1.57067320322857e-01;
                    quadgrid._points[17] = 4.42853280660376e+00;
                    quadgrid._weights[17] = 7.24629595224393e-01;
                    quadgrid._points[18] = 5.22027169053747e+00;
                    quadgrid._weights[18] = 7.24629595224393e-01;
                }

            }

        };

    }
}
#endif	/* LEBEDEV_H */