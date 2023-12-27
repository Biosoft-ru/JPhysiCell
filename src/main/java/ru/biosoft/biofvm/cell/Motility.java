package ru.biosoft.biofvm.cell;

import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.VectorUtil;

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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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
public class Motility implements Cloneable
{
    boolean is_motile;

    double persistence_time; // mean time to keep going in one direction 
    // before resampling for a new direction. 
    double migration_speed; // migration speed along chosen direction, 
    // in absence of all other adhesive / repulsive forces 

    double[] migration_bias_direction = new double[0];; // a unit vector
    // random motility is biased in this direction (e.g., chemotaxis)
    double migration_bias; // how biased is motility
    // if 0, completely random. if 1, deterministic along the bias vector 

    boolean restrict_to_2D;
    // if true, set random motility to 2D only. 

    double[] motility_vector = new double[0];;

    int chemotaxis_index;
    int chemotaxis_direction;

    // advanced chemotaxis 
    double[] chemotactic_sensitivities = new double[0];

    public Motility()
    {
        is_motile = false;

        persistence_time = 1.0;
        migration_speed = 1.0;

        migration_bias_direction = new double[3];
        migration_bias = 0.0;

        restrict_to_2D = false;

        // update_migration_bias_direction = NULL; 

        motility_vector = new double[3];

        chemotaxis_index = 0;
        chemotaxis_direction = 1;
    }


    void sync(Microenvironment m)
    {
        chemotactic_sensitivities = VectorUtil.resize( chemotactic_sensitivities, m.number_of_densities(), 0 );
        //        chemotactic_sensitivities.resize( pNew_Microenvironment.number_of_densities() , 0.0 ); 
        return;
    }

    //double is ref
    double chemotactic_sensitivity(String name)
    {
        //        int n = microenvironment.find_density_index(name); 
        //        return chemotactic_sensitivities[n]; 
        return 0;
    }

    @Override
    public Motility clone()
    {
        try
        {
            Motility result = (Motility)super.clone();
            result.chemotactic_sensitivities = this.chemotactic_sensitivities.clone();
            result.motility_vector = this.motility_vector.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
