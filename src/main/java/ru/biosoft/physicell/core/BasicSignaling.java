package ru.biosoft.physicell.core;

import java.util.List;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2023, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

public class BasicSignaling
{
    //    #include "./PhysiCell_basic_signaling.h"
    //
    //    using namespace BioFVM; 
    //
    //    namespace PhysiCell{

    public static double Hill_response_function(double s, double half_max, double hill_power)
    {
        // newer. only one expensive a^b operation. 45% less computationl expense. 

        // give an early exit possibility to cut cost on "empty" rules
        if( s < 1e-16 ) // maybe also try a dynamic threshold: 0.0001 * half_max 
        {
            return 0.0;
        }

        // operations to reduce a^b operations and minimize hidden memory allocation / deallocation / copy operations. 
        // Hill = (s/half_max)^hill_power / ( 1 + (s/half_max)^hill_power  )
        double temp = s; // s 
        temp /= half_max; // s/half_max 
        double temp1 = Math.pow( temp, hill_power ); // (s/half_max)^h 
        temp = temp1; // (s/half_max)^h 
        temp += 1; // (1+(s/half_max)^h ); 
        temp1 /= temp; // (s/half_max)^h / ( 1 + s/half_max)^h) 
        return temp1;
    }


    public static double linear_response_function(double s, double s_min, double s_max)
    {
        if( s <= s_min )
        {
            return 0.0;
        }
        if( s >= s_max )
        {
            return 1.0;
        }
        s -= s_min; // overwrite s with s - s_min 
        s_max -= s_min; // overwrite s_max with s_max - s_min 
        s /= s_max; // now we have (s-s_min)/(s_max-s_min
        return s;
    }

    public static double decreasing_linear_response_function(double s, double s_min, double s_max)
    {
        if( s <= s_min )
        {
            return 1.0;
        }
        if( s >= s_max )
        {
            return 0.0;
        }
        // (smax-s)/(smax-smin); 
        // = -(s-smax)/(smax-smin)
        s -= s_max; // replace s by s-s_max 
        s_max -= s_min; // replace s_max = s_max - s_min 
        s /= s_max; // this is (s-s_max)/(s_max-s_min)
        s *= -1; // this is (s_max-s)/(s_max-s_min)
        return s;
    }

    public static double interpolate_behavior(double base_value, double max_changed_value, double response)
    {
        double output = max_changed_value; // bM
        output -= base_value; // (bM-b0); 
        output *= response; // R*(bM-b0); 
        output += base_value; // b0 + (bM-b0)*R; 
        return output;
    }

    public static double multivariate_Hill_response_function(List<Double> signals, List<Double> half_maxes, List<Double> hill_powers)
    {
        double temp1 = 0.0;
        double temp2 = 0.0;
        double temp3 = 0.0;
        // create the generalized (s^h), stored in temp1; 
        for( int j = 0; j < signals.size(); j++ )
        {
            temp2 = signals.get( j ); // s
            temp2 /= half_maxes.get( j ); // s/s_half 
            temp3 = Math.pow( temp2, hill_powers.get( j ) ); // (s/s_half)^h 
            temp1 += temp3;
        }
        temp2 = temp1; // numerator (S^h)
        temp1 += 1.0; // denominator (1+S^h)
        temp2 /= temp1; // numerator/denominator = S^h / (1+S^h)
        return temp2;
    }

    public static double multivariate_linear_response_function(double[] signals, double[] min_thresholds, double[] max_thresholds)
    {
        double output = 0.0;

        for( int j = 0; j < signals.length; j++ )
        {
            output += linear_response_function( signals[j], min_thresholds[j], max_thresholds[j] );
        }

        if( output > 1.0 )
        {
            return 1.0;
        }

        return output;
    }

    public static double[] linear_response_to_Hill_parameters(double s0, double s1)
    {
        double tol = 0.1;
        double param1 = ( 1 - tol ) / tol;
        double param2 = Math.log( param1 );

        // half max, then hill power 
        double hm = 0.5 * ( s0 + s1 );

        // hp so that H(s1) ~ (1-tol)
        double hp = Math.round( param2 / Math.log( s1 / hm ) );

        return new double[] {hm, hp};
        //    std::vector<double>output=
        //    { hm , hp };
        //
        //    return output;
    }

    public static double[] Hill_response_to_linear_parameters(double half_max, double Hill_power)
    {
        double tol = 0.1;
        double param1 = ( 1 - tol ) / tol;
        double param2 = Math.pow( param1, 1.0 / Hill_power );

        // s1 such that H(s1) ~ (1-tol)
        double s1 = half_max * param2;

        // s0 for symmetry
        double s0 = 2 * half_max - s1;
        if( s0 < 0 )
        {
            s0 = 0.0;
        }
        return new double[] {s0, s1};
        //    std::vector<double>output=
        //    {s0,s1};
        //
        //    return output;
    }
}