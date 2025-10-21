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

    /**
     * Hill response function
     * @param s
     * @param halfMax
     * @param hillPower
     * @return S / (S+1), where S = (s / hm)^hp
     */
    public static double hillResponse(double s, double halfMax, double hillPower)
    {
        // give an early exit possibility to cut cost on "empty" rules
        if( s < 1e-16 ) // maybe also try a dynamic threshold: 0.0001 * half_max 
            return 0.0;

        double result = Math.pow( s / halfMax, hillPower );
        return result / ( result + 1 );
    }


    /**
     * Linear response function
     */
    public static double linearResponse(double s, double sMin, double sMax)
    {
        if( s <= sMin )
            return 0.0;
        if( s >= sMax )
            return 1.0;
        return ( s - sMin ) / ( sMax - sMin );
    }

    public static double decreasingLinearResponse(double s, double sMin, double sMax)
    {
        if( s <= sMin )
            return 1.0;
        if( s >= sMax )
            return 0.0;
        return ( sMax - s ) / ( sMax - sMin );
    }

    public static double interpolateBehavior(double baseValue, double maxChangedValue, double response)
    {
        return baseValue + ( maxChangedValue - baseValue ) * response + baseValue;
    }

    /**
     * Multivariate Hill response function
     * @return S / (1 + S), where S = is sum of Hill responses 
     */
    public static double multivariateHillResponse(List<Double> signals, List<Double> halfMaxes, List<Double> hillPowers)
    {
        double result = 0.0;
        for( int j = 0; j < signals.size(); j++ )
            result += Math.pow( signals.get( j ) / halfMaxes.get( j ), hillPowers.get( j ) );
        return result / ( result + 1 );
    }

    /**
     * Multivariate Linear response function
     * @return sum of responses (signal - min)/ (max - min)
     */
    public static double multivariateLinearResponse(double[] signals, double[] min, double[] max)
    {
        double output = 0.0;
        for( int j = 0; j < signals.length; j++ )
            output += linearResponse( signals[j], min[j], max[j] );
        if( output > 1.0 )
            return 1.0;
        return output;
    }

    /**
     * Converts parameters of linear response to parameters of Hilll response
     * @param s0
     * @param s1
     * @return half max, and hill power
     */
    public static double[] linearToHill(double s0, double s1)
    {
        double tol = 0.1;

        // half max, then hill power 
        double hm = 0.5 * ( s0 + s1 );

        // hp so that H(s1) ~ (1-tol)
        double param2 = Math.log( ( 1 - tol ) / tol );
        double hp = Math.round( param2 / Math.log( s1 / hm ) );

        return new double[] {hm, hp};
    }

    /**
     * Converts parameters of Hill response to parameters of Linear response
     * @param halfMax
     * @param hillPower
     * @return s0, and s1
     */
    public static double[] hillToLinear(double halfMax, double hillPower)
    {
        double tol = 0.1;

        // s1 such that H(s1) ~ (1-tol)
        double s1 = halfMax * Math.pow( ( 1 - tol ) / tol, 1.0 / hillPower );

        // s0 for symmetry
        double s0 = 2 * halfMax - s1;
        if( s0 < 0 )
            s0 = 0.0;

        return new double[] {s0, s1};
    }
}