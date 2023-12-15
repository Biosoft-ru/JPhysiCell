package ru.biosoft.biofvm;

/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/
public class MicroenvironmentOptions
{
    Microenvironment pMicroenvironment;
    String name;

    String time_units;
    String spatial_units;
    double dx;
    double dy;
    double dz;

    boolean outer_Dirichlet_conditions;
    double[] Dirichlet_condition_vector;
    boolean[] Dirichlet_activation_vector;

    /* new in PhysiCell 1.7.0 to enable setting Dirichlet conditions 
       on a boundary-by-boundary basis */
    boolean[] Dirichlet_all = new boolean[0];

    //  boolean[] Dirichlet_interior; 
    boolean[] Dirichlet_xmin = new boolean[0];
    boolean[] Dirichlet_xmax = new boolean[0];
    boolean[] Dirichlet_ymin = new boolean[0];
    boolean[] Dirichlet_ymax = new boolean[0];
    boolean[] Dirichlet_zmin = new boolean[0];
    boolean[] Dirichlet_zmax = new boolean[0];

    double[] Dirichlet_xmin_values = new double[0];
    double[] Dirichlet_xmax_values = new double[0];
    double[] Dirichlet_ymin_values = new double[0];
    double[] Dirichlet_ymax_values = new double[0];
    double[] Dirichlet_zmin_values = new double[0];
    double[] Dirichlet_zmax_values = new double[0];

    double[] initial_condition_vector = new double[0];

    public boolean simulate_2D;
    double[] X_range;
    double[] Y_range;
    double[] Z_range;

    public boolean calculate_gradients;

    boolean use_oxygen_as_first_field;

    public boolean track_internalized_substrates_in_each_agent;

    MicroenvironmentOptions(Microenvironment microenvironment)
    {
        use_oxygen_as_first_field = true;

        //        if( Microenvironment.get_default_microenvironment() != null )
        //        {
        //            pMicroenvironment = Microenvironment.get_default_microenvironment();
        //        }
        //        else
        //        {
        //            pMicroenvironment = new Microenvironment();//TODO: check 
        //            Microenvironment.set_default_microenvironment( pMicroenvironment );
        //        }
        name = "microenvironment";
        pMicroenvironment = microenvironment;
        time_units = "min";
        spatial_units = "micron";
        dx = 20;
        dy = 20;
        dz = 20;

        outer_Dirichlet_conditions = false;
        Dirichlet_condition_vector = new double[pMicroenvironment.number_of_densities()];
        for( int i = 0; i < Dirichlet_condition_vector.length; i++ )
            Dirichlet_condition_vector[i] = 1;
        //        Dirichlet_condition_vector.assign( pMicroenvironment->number_of_densities() , 1.0 ); 
        //        Dirichlet_activation_vector.assign( pMicroenvironment->number_of_densities() , false ); 
        Dirichlet_activation_vector = new boolean[pMicroenvironment.number_of_densities()];

        initial_condition_vector = new double[0];
        //        initial_condition_vector.resize(0); //  = Dirichlet_condition_vector; 
        //        
        // set a far-field value for oxygen (assumed to be in the first field)
        Dirichlet_condition_vector[0] = 38.0;

        simulate_2D = false;

        X_range = new double[] { -500, 500};
        //        X_range.resize(2,500.0); 
        //        X_range[0] *= -1.0;

        Y_range = new double[] { -500, 500};
        //        Y_range.resize(2,500.0); 
        //        Y_range[0] *= -1.0;
        Z_range = new double[] { -500, 500};
        //        Z_range.resize(2,500.0); 
        //        Z_range[0] *= -1.0;

        calculate_gradients = false;

        track_internalized_substrates_in_each_agent = false;

        VectorUtil.push_back( Dirichlet_all, true );
        VectorUtil.push_back( Dirichlet_xmin, false );
        VectorUtil.push_back( Dirichlet_xmax, false );
        VectorUtil.push_back( Dirichlet_ymin, false );
        VectorUtil.push_back( Dirichlet_ymax, false );
        VectorUtil.push_back( Dirichlet_zmin, false );
        VectorUtil.push_back( Dirichlet_zmax, false );

        VectorUtil.push_back( Dirichlet_activation_vector, calculate_gradients );
        VectorUtil.push_back( Dirichlet_activation_vector, calculate_gradients );
        //        Dirichlet_all.push_back( true ); 
        //    //  Dirichlet_interior.push_back( true ); 
        //        Dirichlet_xmin.push_back( false ); 
        //        Dirichlet_xmax.push_back( false ); 
        //        Dirichlet_ymin.push_back( false ); 
        //        Dirichlet_ymax.push_back( false ); 
        //        Dirichlet_zmin.push_back( false ); 
        //        Dirichlet_zmax.push_back( false ); 

        //        default_microenvironment_options.Dirichlet_xmin_values.push_back( 1.0 ); 
        //        default_microenvironment_options.Dirichlet_xmax_values.push_back( 1.0 ); 
        //        default_microenvironment_options.Dirichlet_ymin_values.push_back( 1.0 ); 
        //        default_microenvironment_options.Dirichlet_ymax_values.push_back( 1.0 ); 
        //        default_microenvironment_options.Dirichlet_zmin_values.push_back( 1.0 ); 
        //        default_microenvironment_options.Dirichlet_zmax_values.push_back( 1.0 );


        //TODO: check default options
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_xmin_values, 1 );
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_xmax_values, 1 );
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_ymin_values, 1 );
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_ymax_values, 1 );
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_zmin_values, 1 );
        //        VectorUtil.push_back( default_microenvironment_options.Dirichlet_zmax_values, 1 );
    }


}
