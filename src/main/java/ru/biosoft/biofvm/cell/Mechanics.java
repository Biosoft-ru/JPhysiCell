package ru.biosoft.biofvm.cell;

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
public class Mechanics implements Cloneable
{
    double cell_cell_adhesion_strength;
    double cell_BM_adhesion_strength;

    double cell_cell_repulsion_strength;
    double cell_BM_repulsion_strength;

    double[] cell_adhesion_affinities;

    // this is a multiple of the cell (equivalent) radius
    double relative_maximum_adhesion_distance;
    // double maximum_adhesion_distance; // needed? 

    /* for spring attachments */
    int maximum_number_of_attachments;
    double attachment_elastic_constant;

    double attachment_rate;
    double detachment_rate;

    /* to be deprecated */
    double relative_maximum_attachment_distance;
    double relative_detachment_distance;
    double maximum_attachment_rate;

    public Mechanics()
    {
        cell_cell_adhesion_strength = 0.4;
        cell_BM_adhesion_strength = 4.0;

        cell_cell_repulsion_strength = 10.0;
        cell_BM_repulsion_strength = 100.0;

        cell_adhesion_affinities = new double[] {1};

        // this is a multiple of the cell (equivalent) radius
        relative_maximum_adhesion_distance = 1.25;
        // maximum_adhesion_distance = 0.0; 

        /* for spring attachments */
        maximum_number_of_attachments = 12;
        attachment_elastic_constant = 0.01;

        attachment_rate = 0; // 10.0 prior ot March 2023
        detachment_rate = 0;

        /* to be deprecated */
        relative_maximum_attachment_distance = relative_maximum_adhesion_distance;
        relative_detachment_distance = relative_maximum_adhesion_distance;

        maximum_attachment_rate = 1.0;
    }

    public void sync_to_cell_definitions()
    {
        int number_of_cell_defs = CellDefinition.getDefinitionsCount();
        if( cell_adhesion_affinities.length != number_of_cell_defs )
            VectorUtil.resize( cell_adhesion_affinities, number_of_cell_defs, 1.0 );
    }

    double cell_adhesion_affinity(String type_name)
    {
        int n = CellDefinition.getCellDefinition( type_name ).type;
        return cell_adhesion_affinities[n];
    }

    void set_fully_heterotypic()
    {
        //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
        int number_of_cell_defs = CellDefinition.getDefinitionsCount();
        //        cell_adhesion_affinities.assign( number_of_cell_defs, 1.0);
        cell_adhesion_affinities = VectorUtil.assign( number_of_cell_defs, 1.0 );
    }

    void set_fully_homotypic(Cell pC)
    {
        //        extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
        int number_of_cell_defs = CellDefinition.getDefinitionsCount();
        //        cell_adhesion_affinities.assign( number_of_cell_defs, 0.0);
        cell_adhesion_affinities = new double[number_of_cell_defs];
        // now find my type and set to 1 
        //  cell_adhesion_affinity( pC->type_name ) = 1.0; 
    }


    // new on July 29, 2018
    // change the ratio without changing the repulsion strength or equilibrium spacing 
    void set_relative_maximum_adhesion_distance(double new_value)
    {
        // get old equilibrium spacing, based on equilibriation of pairwise adhesive/repulsive forces at that distance. 

        // relative equilibrium spacing (relative to mean cell radius)
        double s_relative = 2.0;

        double temp1 = cell_cell_adhesion_strength;
        temp1 /= cell_cell_repulsion_strength;
        temp1 = Math.sqrt( temp1 );

        double temp2 = 1.0;
        temp2 -= temp1; //  1 - sqrt( alpha_CCA / alpha_CCR );


        s_relative *= temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ); 

        temp1 /= relative_maximum_adhesion_distance; // sqrt( alpha_CCA / alpha_CCR)/f;
        temp2 = 1.0;
        temp2 -= temp1; // 1 - sqrt( alpha_CCA / alpha_CCR )/f;

        s_relative /= temp2; // 2*( 1 - sqrt( alpha_CCA / alpha_CCR ) ) / ( 1-1/f) ; 

        // now, adjust the relative max adhesion distance 

        relative_maximum_adhesion_distance = new_value;

        // adjust the adhesive coefficient to preserve the old equilibrium distance

        temp1 = s_relative;
        temp1 /= 2.0;

        temp2 = 1.0;
        temp2 -= temp1; // 1 - s_relative/2.0 

        temp1 /= relative_maximum_adhesion_distance; // s_relative/(2*relative_maximum_adhesion_distance); 
        temp1 *= -1.0; // -s_relative/(2*relative_maximum_adhesion_distance); 
        temp1 += 1.0; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance); 

        temp2 /= temp1;
        temp2 *= temp2;

        cell_cell_adhesion_strength = cell_cell_repulsion_strength;
        cell_cell_adhesion_strength *= temp2;

        return;
    }

    // new on July 29, 2018
    // set the cell-cell equilibrium spacing, accomplished by changing the 
    // cell-cell adhesion strength, while leaving the cell-cell repulsion 
    // strength and the maximum adhesion distance unchanged 
    void set_relative_equilibrium_distance(double new_value)
    {
        if( new_value > 2.0 )
        {
            //            std::cout << "**** Warning in function " << __FUNCTION__ << " in " << __FILE__ << " : " << std::endl 
            //                << "\tAttempted to set equilibrium distance exceeding two cell radii." << std::endl
            //                << "\tWe will cap the equilibrium distance at 2.0 cell radii." << std::endl 
            //                << "****" << std::endl << std::endl; 

            new_value = 2.0;
        }

        // adjust the adhesive coefficient to achieve the new (relative) equilibrium distance

        double temp1 = new_value;
        temp1 /= 2.0;

        double temp2 = 1.0;
        temp2 -= temp1; // 1 - s_relative/2.0 

        temp1 /= relative_maximum_adhesion_distance; // s_relative/(2*relative_maximum_adhesion_distance); 
        temp1 *= -1.0; // -s_relative/(2*relative_maximum_adhesion_distance); 
        temp1 += 1.0; // 1.0 -s_relative/(2*relative_maximum_adhesion_distance); 

        temp2 /= temp1;
        temp2 *= temp2;

        cell_cell_adhesion_strength = cell_cell_repulsion_strength;
        cell_cell_adhesion_strength *= temp2;

        return;
    }

    void set_absolute_equilibrium_distance(Phenotype phenotype, double new_value)
    {
        set_relative_equilibrium_distance( new_value / phenotype.geometry.radius );
    }

    // void Mechanics::set_absolute_maximum_adhesion_distance( double new_value );
    // void 
    @Override
    public Mechanics clone()
    {
        try
        {
            Mechanics result = (Mechanics)super.clone();
            result.cell_adhesion_affinities = cell_adhesion_affinities.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}