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

public class IntegratedSignal
{
    double base_activity;
    double max_activity;

    List<Double> promoters;
    List<Double> promoter_weights;
    double promoters_Hill;
    double promoters_half_max;

    List<Double> inhibitors;
    List<Double> inhibitor_weights;
    double inhibitors_Hill;
    double inhibitors_half_max;

    public IntegratedSignal()
    {
        base_activity = 0.0;
        max_activity = 1.0;

        promoters.clear();
        promoter_weights.clear();

        promoters_half_max = 0.1;
        promoters_Hill = 4;

        inhibitors.clear();
        inhibitor_weights.clear();

        inhibitors_half_max = 0.1;
        inhibitors_Hill = 4;
    }

    void reset()
    {
        promoters.clear();
        promoter_weights.clear();

        inhibitors.clear();
        inhibitor_weights.clear();
    }

    double compute_signal()
    {
        double pr = 0.0;
        double w = 0.0;
        for( int k = 0; k < promoters.size(); k++ )
        {
            pr += promoters.get( k );
            w += promoter_weights.get( k );
        }
        w += 1e-16;
        pr /= w;

        double inhib = 0.0;
        w = 0.0;
        for( int k = 0; k < inhibitors.size(); k++ )
        {
            inhib += inhibitors.get( k );
            w += inhibitor_weights.get( k );
        }
        w += 1e-16;
        inhib /= w;

        double Pn = Math.pow( pr, promoters_Hill );
        double Phalf = Math.pow( promoters_half_max, promoters_Hill );

        double In = Math.pow( inhib, inhibitors_Hill );
        double Ihalf = Math.pow( inhibitors_half_max, inhibitors_Hill );

        double P = Pn / ( Pn + Phalf );
        double I = 1.0 / ( In + Ihalf );

        double output = max_activity;
        output -= base_activity; //(max-base)
        output *= P; // (max-base)*P 
        output += base_activity; // base + (max-base)*P 
        output *= I; // (base + (max-base)*P)*I; 

        return output;
    };

    void add_signal(char signal_type, double signal, double weight)
    {
        if( signal_type == 'P' || signal_type == 'p' )
        {
            promoters.add( signal );
            promoter_weights.add( weight );
            return;
        }
        if( signal_type == 'I' || signal_type == 'i' )
        {
            inhibitors.add( signal );
            inhibitor_weights.add( weight );
            return;
        }
    }

    void add_signal(char signal_type, double signal)
    {
        add_signal( signal_type, signal, 1.0 );
    }

}
