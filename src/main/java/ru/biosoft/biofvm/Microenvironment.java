package ru.biosoft.biofvm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

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
public class Microenvironment
{
    private static Microenvironment defaultMicroenvironment = null;

    public static void set_default_microenvironment(Microenvironment m)
    {
        defaultMicroenvironment = m;
    }

    public static Microenvironment get_default_microenvironment()
    {
        return defaultMicroenvironment;
    }

    public MicroenvironmentOptions options;
    DiffusionDecaySolver solver;

    public void setSolver(DiffusionDecaySolver solver)
    {
        this.solver = solver;
    }

    /*! For internal use and accelerations in solvers */
    //        std::vector< std::vector<double> > temporary_density_vectors1;
    double[][] temporary_density_vectors1;
    /*! For internal use and accelerations in solvers */
    //        std::vector< std::vector<double> > temporary_density_vectors2; 
    double[][] temporary_density_vectors2;
    /*! for internal use in bulk source/sink solvers */
    //        std::vector< std::vector<double> > bulk_source_sink_solver_temp1;
    double[][] bulk_source_sink_solver_temp1;
    //        std::vector< std::vector<double> > bulk_source_sink_solver_temp2;
    double[][] bulk_source_sink_solver_temp2;
    //        std::vector< std::vector<double> > bulk_source_sink_solver_temp3;
    double[][] bulk_source_sink_solver_temp3;
    //        bool bulk_source_sink_solver_setup_done; 
    boolean bulk_source_sink_solver_setup_done;

    /*! stores pointer to current density solutions. Access via operator() functions. */
    //        std::vector< std::vector<double> >* p_density_vectors;
    double[][] p_density_vectors;

    //        std::vector< std::vector<gradient> > gradient_vectors;
    double[][][] gradient_vectors;
    //        std::vector<bool> gradient_vector_computed; 
    boolean[] gradient_vector_computed;

    /*! helpful for solvers -- resize these whenever adding/removing substrates */
    //        std::vector<double> one; 
    double[] one;
    //        std::vector<double> zero;
    double[] zero;
    //        std::vector<double> one_half; 
    double[] one_half;
    //        std::vector<double> one_third; 
    double[] one_third;

    /*! for internal use in diffusion solvers : these make the solvers safe across microenvironments ""*/
    //        std::vector< std::vector<double> > thomas_temp1; 
    double[][] thomas_temp1;
    //        std::vector< std::vector<double> > thomas_temp2; 
    double[][] thomas_temp2;
    //        std::vector<double> thomas_constant1x;
    double[] thomas_constant1x;
    //        std::vector<double> thomas_constant1y; 
    double[] thomas_constant1y;
    //        std::vector<double> thomas_constant1z;
    double[] thomas_constant1z;
    //        std::vector<double> thomas_neg_constant1x; 
    double[] thomas_neg_constant1x;
    //        std::vector<double> thomas_neg_constant1y; 
    double[] thomas_neg_constant1y;
    //        std::vector<double> thomas_neg_constant1z;
    double[] thomas_neg_constant1z;

    //        bool thomas_setup_done; 
    boolean thomas_setup_done;
    int thomas_i_jump;
    int thomas_j_jump;
    int thomas_k_jump;
    //        std::vector<double> thomas_constant1; 
    double[] thomas_constant1;
    //        std::vector<double> thomas_constant1a; 
    double[] thomas_constant1a;
    //        std::vector<double> thomas_constant2;
    double[] thomas_constant2;
    //        std::vector<double> thomas_constant3;
    double[] thomas_constant3;
    //        std::vector<double> thomas_constant3a;
    double[] thomas_constant3a;
    //        std::vector< std::vector<double> > thomas_denomx;
    double[][] thomas_denomx;
    //        std::vector< std::vector<double> > thomas_cx;
    double[][] thomas_cx;
    //        std::vector< std::vector<double> > thomas_denomy;
    double[][] thomas_denomy;
    //        std::vector< std::vector<double> > thomas_cy;
    double[][] thomas_cy;
    //        std::vector< std::vector<double> > thomas_denomz;
    double[][] thomas_denomz;
    //        std::vector< std::vector<double> > thomas_cz;
    double[][] thomas_cz;
    boolean diffusion_solver_setup_done;

    // on "resize density" type operations, need to extend all of these 

    /*
    std::vector<int> dirichlet_indices; 
    std::vector< std::vector<double> > dirichlet_value_vectors; 
    std::vector<bool> dirichlet_node_map; 
    */
    //        std::vector< std::vector<double> > dirichlet_value_vectors; 
    double[][] dirichlet_value_vectors;
    //        std::vector<bool> dirichlet_activation_vector; 
    boolean[] dirichlet_activation_vector;
    /* new in Version 1.7.0 -- activation vectors can be specified 
       on a voxel-by-voxel basis */

    //        std::vector< std::vector<bool> > dirichlet_activation_vectors; 
    boolean[][] dirichlet_activation_vectors;
    //     public:

    /*! The mesh for the diffusing quantities */
    public CartesianMesh mesh;
    AgentContainer agent_container = new AgentContainer();
    public String spatial_units;
    public String time_units;
    public String name;
    //
    //        // diffusing entities 
    //        std::vector< std::string > density_names;
    String[] density_names;
    //        std::vector< std::string > density_units; 
    String[] density_units;
    //     
    //        // coefficients 
    //        std::vector< double > diffusion_coefficients; 
    public double[] diffusion_coefficients;
    //        std::vector< double > decay_rates; 
    public double[] decay_rates;
    //        
    //        std::vector< std::vector<double> > supply_target_densities_times_supply_rates; 
    double[][] supply_target_densities_times_supply_rates;
    //        std::vector< std::vector<double> > supply_rates; 
    double[][] supply_rates;
    //        std::vector< std::vector<double> > uptake_rates; 
    double[][] uptake_rates;


    public Microenvironment()
    {
        if( defaultMicroenvironment == null )
        {
            defaultMicroenvironment = this;
        }


        name = "unnamed";
        spatial_units = "none";
        time_units = "none";

        bulk_source_sink_solver_setup_done = false;
        thomas_setup_done = false;
        diffusion_solver_setup_done = false;

        solver = new ConstantCoefficientsLOD3D();

        mesh = new CartesianMesh();
        mesh.resize( 1, 1, 1 );

        one = new double[] {1};
        zero = new double[] {0};

        temporary_density_vectors1 = new double[mesh.voxels.length][1];
        temporary_density_vectors2 = new double[mesh.voxels.length][1];
        p_density_vectors = temporary_density_vectors1;

        gradient_vectors = new double[mesh.voxels.length][1][3];
        //        gradient_vectors.resize( mesh.voxels.size() ); 
        //        for( unsigned int k=0 ; k < mesh.voxels.size() ; k++ )
        //        {
        //            gradient_vectors[k].resize( 1 ); 
        //            (gradient_vectors[k])[0].resize( 3, 0.0 );
        //        }
        //        gradient_vector_computed.resize( mesh.voxels.size() , false ); 
        gradient_vector_computed = new boolean[mesh.voxels.length];

        //        bulk_supply_rate_function = zero_function;
        //        bulk_supply_target_densities_function = zero_function;
        //        bulk_uptake_rate_function = zero_function;

        density_names = new String[] {"unnamed"};
        density_units = new String[] {"none"};
        //        density_names.assign( 1, "unnamed" );
        //        density_units.assign( 1, "none" );

        //        diffusion_coefficients.assign( number_of_densities(), 0.0 );
        //        decay_rates.assign( number_of_densities(), 0.0 );
        diffusion_coefficients = new double[number_of_densities()];
        decay_rates = new double[number_of_densities()];

        //        one_half = one;
        //        one_half *= 0.5;
        //        one_third = one;
        //        one_third /= 3.0;

        one_half = new double[] {0.5};
        one_third = new double[] {1.0 / 3.0};

        //        dirichlet_activation_vector.assign( 1, false );
        //        dirichlet_value_vectors.assign( mesh.voxels.size(), one );
        //        dirichlet_activation_vectors.assign( 1, dirichlet_activation_vector );
        dirichlet_value_vectors = VectorUtil.assign( mesh.voxels.length, one );
        dirichlet_activation_vector = new boolean[] {false};
        dirichlet_activation_vectors = VectorUtil.assign( 1, dirichlet_activation_vector );

        options = new MicroenvironmentOptions( this );
        options.Dirichlet_all = new boolean[] {true};
        options.Dirichlet_xmin = new boolean[] {false};
        options.Dirichlet_xmax = new boolean[] {false};
        options.Dirichlet_ymin = new boolean[] {false};
        options.Dirichlet_ymax = new boolean[] {false};
        options.Dirichlet_zmin = new boolean[] {false};
        options.Dirichlet_zmax = new boolean[] {false};

        options.Dirichlet_xmin_values = new double[] {1.0};
        options.Dirichlet_xmax_values = new double[] {1.0};
        options.Dirichlet_ymin_values = new double[] {1.0};
        options.Dirichlet_ymax_values = new double[] {1.0};
        options.Dirichlet_zmin_values = new double[] {1.0};
        options.Dirichlet_zmax_values = new double[] {1.0};
        //        default_microenvironment_options.Dirichlet_all.assign( 1, true );
        //        default_microenvironment_options.Dirichlet_xmin.assign( 1, false );
        //        default_microenvironment_options.Dirichlet_xmax.assign( 1, false );
        //        default_microenvironment_options.Dirichlet_ymin.assign( 1, false );
        //        default_microenvironment_options.Dirichlet_ymax.assign( 1, false );
        //        default_microenvironment_options.Dirichlet_zmin.assign( 1, false );
        //        default_microenvironment_options.Dirichlet_zmax.assign( 1, false );

        //        default_microenvironment_options.Dirichlet_xmin_values.assign( 1, 1.0 );
        //        default_microenvironment_options.Dirichlet_xmax_values.assign( 1, 1.0 );
        //        default_microenvironment_options.Dirichlet_ymin_values.assign( 1, 1.0 );
        //        default_microenvironment_options.Dirichlet_ymax_values.assign( 1, 1.0 );
        //        default_microenvironment_options.Dirichlet_zmin_values.assign( 1, 1.0 );
        //        default_microenvironment_options.Dirichlet_zmax_values.assign( 1, 1.0 );


    }

    public double[] get(int n)
    {
        return p_density_vectors[n];
    }


    public void simulate_cell_sources_and_sinks(List<BasicAgent> basic_agent_list, double dt)
    {
        for( int i = 0; i < basic_agent_list.size(); i++ )
        {
            basic_agent_list.get( i ).simulateSecretionUptake( this, dt );
        }

        return;
    }

    public void simulate_cell_sources_and_sinks(double dt)
    {
        simulate_cell_sources_and_sinks( BasicAgent.allBasicAgents, dt );
    }

    public void simulate_diffusion_decay(double dt) throws Exception
    {
        if( solver != null )
            solver.solve( this, dt );
    }

    void apply_dirichlet_conditions()
    {
        /*
        #pragma omp parallel for 
        for( unsigned int i=0 ; i < dirichlet_indices.size() ; i++ )
        { density_vector( dirichlet_indices[i] ) = dirichlet_value_vectors[i]; }
        */

        // #pragma omp parallel for 
        for( int i = 0; i < mesh.voxels.length; i++ )
        {
            /*
            if( mesh.voxels[i].is_Dirichlet == true )
            { density_vector(i) = dirichlet_value_vectors[i]; }
            */
            if( mesh.voxels[i].isDirichlet == true )
            {
                for( int j = 0; j < dirichlet_value_vectors[i].length; j++ )
                {
                    // if( dirichlet_activation_vector[j] == true )
                    if( dirichlet_activation_vectors[i][j] == true )
                    {
                        density_vector( i )[j] = dirichlet_value_vectors[i][j];
                    }
                }

            }
        }
        return;
    }

    public void write_to_matlab(String filename)
    {
        int number_of_data_entries = mesh.voxels.length;
        int size_of_each_datum = 3 + 1 + p_density_vectors[0].length;
        File f = new File( filename );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f ) ))
        {
            write_matlab4_header( size_of_each_datum, number_of_data_entries, bw, "multiscale_microenvironment" );
            // storing data as cols 
            for( int i = 0; i < number_of_data_entries; i++ )
            {
                bw.append( String.valueOf( mesh.voxels[i].center[0] ) + "\t" );
                bw.append( String.valueOf( mesh.voxels[i].center[1] ) + "\t" );
                bw.append( String.valueOf( mesh.voxels[i].center[2] ) + "\t" );
                bw.append( String.valueOf( mesh.voxels[i].volume ) + "\t" );

                // densities  
                for( int j = 0; j < p_density_vectors[i].length; j++ )
                {
                    bw.append( String.valueOf( p_density_vectors[i][j] ) + "\t" );
                }
                bw.append( "\n" );
            }
        }
        catch( Exception ex )
        {

        }
    }

    void write_matlab_header(double[][] input, BufferedWriter bw, String variable_name) throws Exception
    {
        int number_of_data_entries = input.length;
        int size_of_each_datum = input[0].length;

        int rows = size_of_each_datum; // storing data as cols
        int cols = number_of_data_entries; // storing data as cols

        write_matlab4_header( rows, cols, bw, variable_name );

        // // storing data as rows
        bw.append( "\nStoring data as rows\n" );
        for( int j = 0; j < size_of_each_datum; j++ )
        {
            for( int i = 0; i < number_of_data_entries; i++ )
            {
                //         fwrite( (char*) &(input[i])[j] , sizeof(double), 1 , fp );
                bw.append( String.valueOf( input[i][j] ) );
                bw.append( "\t" );
            }
            bw.append( "\n" );
        }

        // storing data as cols 
        bw.append( "\nStoring data as cols\n" );
        for( int i = 0; i < number_of_data_entries; i++ )
        {
            for( int j = 0; j < size_of_each_datum; j++ )
            {
                //          bw.write( input[i])[j] );
                bw.append( String.valueOf( input[i][j] ) );
                bw.append( "\t" );
                //       fwrite( (char*) &(input[i])[j] , sizeof(double), 1 , fp );
            }
            bw.append( "\n" );
        }
    }

    void write_matlab4_header(int nrows, int ncols, BufferedWriter bw, String variable_name) throws Exception
    {
        //     FILE* fp; 
        //     fp = fopen( filename.c_str() , "wb" );
        //     if( fp == NULL )
        //     {
        //      std::cout << "Error: could not open file " << filename << "!" << std::endl;
        //      return NULL;
        //     }
        //             typedef unsigned int UINT;
        //             UINT UINTs = sizeof(UINT);
        //     
        //     UINT temp;
        //     
        //     UINT type_numeric_format = 0; // little-endian assumed for now!
        //     UINT type_reserved = 0;
        //     UINT type_data_format = 0; // doubles for all entries 
        //     UINT type_matrix_type = 0; // full matrix, not sparse
        //   
        //
        //     temp = 1000*type_numeric_format + 100*type_reserved + 10*type_data_format + type_matrix_type;
        int type_numeric_format = 0; // little-endian assumed for now!
        int type_reserved = 0;
        int type_data_format = 0; // doubles for all entries 
        int type_matrix_type = 0; // full matrix, not sparse
        int temp = 1000 * type_numeric_format + 100 * type_reserved + 10 * type_data_format + type_matrix_type;
        bw.append( "\nHEADER\n" );
        bw.append( String.valueOf( temp ) );
        bw.append( "\t" );
        //     fwrite( (char*) &temp , UINTs , 1 , fp );
        //     
        //    // UINT rows = (UINT) number_of_data_entries; // storing data as rows
        //     UINT rows = (UINT) nrows; // size_of_each_datum; // storing data as cols
        int rows = nrows;
        //     fwrite( (char*) &rows , UINTs , 1, fp );
        bw.append( String.valueOf( rows ) );
        bw.append( "\t" );
        //    // UINT cols = (UINT) size_of_each_datum; // storing data as rows
        //     UINT cols = (UINT) ncols; // number_of_data_entries; // storing data as cols
        int cols = ncols;
        //     fwrite( (char*) &cols, UINTs , 1 , fp );
        bw.append( String.valueOf( cols ) );
        bw.append( "\t" );
        //     UINT imag = 0; // no complex matrices!
        int imag = 0;
        //     fwrite( (char*) &imag, UINTs, 1 , fp );
        bw.append( String.valueOf( imag ) );
        bw.append( "\t" );
        //     UINT name_length = variable_name.size(); // strlen( variable_name );
        int name_length = variable_name.length();
        //     fwrite( (char*) &name_length, UINTs, 1 , fp );
        bw.append( String.valueOf( name_length ) );
        bw.append( "\t" );
        // this is the end of the 20-byte header 

        // write the name

        //     fwrite( variable_name.c_str() , name_length , 1 , fp );
        bw.append( variable_name );
        bw.append( "\n" );
        //     return fp; 
    }

    int voxel_index(int i, int j, int k)
    {
        return mesh.voxel_index( i, j, k );
    }

    /*! access the density vector at  [ X(i),Y(j),Z(k) ] */
    double[] density_vector(int i, int j, int k)
    {
        return p_density_vectors[voxel_index( i, j, k )];
    }

    /*! access the density vector at [x,y,z](n) */
    public double[] density_vector(int voxel_index)
    {
        return p_density_vectors[voxel_index];
    }

    int nearest_voxel_index(double[] position)
    {
        return mesh.nearest_voxel_index( position );
    }

    public void set_density(int index, String name, String units)
    {
        // fix in PhysiCell preview November 2017 
        if( index == 0 )
        {
            options.use_oxygen_as_first_field = false;//TODO: check
        }

        density_names[index] = name;
        density_units[index] = units;
    }

    public int number_of_densities()
    {
        return p_density_vectors[0].length;
    }

    public int number_of_voxels()
    {
        return mesh.voxels.length;
    }

    public Voxel voxels(int voxel_index)
    {
        return mesh.voxels[voxel_index];
    }

    public void resize_space_uniform(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end,
            double dx_new)
    {
        resize_space( x_start, x_end, y_start, y_end, z_start, z_end, dx_new, dx_new, dx_new );
    }

    void resize_space(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double x_nodes,
            double y_nodes, double z_nodes)
    {
        //        temporary_density_vectors1.assign( mesh.voxels.size() , zero ); 
        //        temporary_density_vectors2.assign( mesh.voxels.size() , zero ); 
        //        gradient_vectors.resize( mesh.voxels.size() ); 
        //        for( unsigned int k=0 ; k < mesh.voxels.size() ; k++ )
        //        {
        //            gradient_vectors[k].resize( number_of_densities() ); 
        //            for( unsigned int i=0 ; i < number_of_densities() ; i++ )
        //            {
        //                (gradient_vectors[k])[i].resize( 3, 0.0 );
        //            }
        //        }
        //        gradient_vector_computed.resize( mesh.voxels.size() , false );  
        //        dirichlet_value_vectors.assign( mesh.voxels.size(), one ); 
        //        dirichlet_activation_vectors.assign( mesh.voxels.length , dirichlet_activation_vector );  
        mesh.resize( x_start, x_end, y_start, y_end, z_start, z_end, x_nodes, y_nodes, z_nodes );
        temporary_density_vectors1 = new double[mesh.voxels.length][zero.length];
        p_density_vectors = temporary_density_vectors1;
        temporary_density_vectors2 = new double[mesh.voxels.length][zero.length];
        gradient_vectors = new double[mesh.voxels.length][number_of_densities()][3];
        gradient_vector_computed = new boolean[mesh.voxels.length];
        dirichlet_value_vectors = VectorUtil.assign( mesh.voxels.length, one );
        dirichlet_activation_vectors = VectorUtil.assign( mesh.voxels.length, dirichlet_activation_vector );
    }

    void resize_space(int x_nodes, int y_nodes, int z_nodes)
    {
        //        temporary_density_vectors1.assign( mesh.voxels.length , zero ); 
        //        temporary_density_vectors2.assign( mesh.voxels.length , zero ); 
        //        gradient_vectors.resize( mesh.voxels.length ); 

        //        for( int k=0 ; k < mesh.voxels.length ; k++ )
        //        {
        //            gradient_vectors[k] = new double[number_of_densities()][];//.resize(  ); 
        //            for( int i=0 ; i < number_of_densities() ; i++ )
        //            {
        //                gradient_vectors[k][i] = new double[3];
        ////                (gradient_vectors[k])[i].resize( 3, 0.0 );
        //            }
        //        }

        //        gradient_vector_computed = new double
        //        gradient_vector_computed.resize( mesh.voxels.length , false );  
        //        dirichlet_value_vectors.assign( mesh.voxels.length, one ); 
        //        dirichlet_activation_vectors.assign( mesh.voxels.length , dirichlet_activation_vector );
        mesh.resize( x_nodes, y_nodes, z_nodes );
        temporary_density_vectors1 = new double[mesh.voxels.length][];
        p_density_vectors = temporary_density_vectors1;
        temporary_density_vectors2 = new double[mesh.voxels.length][];
        gradient_vectors = new double[mesh.voxels.length][number_of_densities()][3];
        gradient_vector_computed = new boolean[mesh.voxels.length];
        dirichlet_value_vectors = VectorUtil.assign( mesh.voxels.length, one );
        dirichlet_activation_vectors = VectorUtil.assign( mesh.voxels.length, dirichlet_activation_vector );

    }

    public String displayInformation()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nMicroenvironment summary: " + name + ": \n" );
        sb.append( mesh.display_information() );
        sb.append( "Densities: (" + number_of_densities() + " total)" + "\n" );
        for( int i = 0; i < density_names.length; i++ )
        {
            sb.append( "   " + density_names[i] + ":" + "\n" );
            sb.append( "     units: " + density_units[i] + "\n" );
            sb.append( "     diffusion coefficient: " + diffusion_coefficients[i] );
            sb.append( " " + spatial_units + "^2 / " + time_units + "\n" );
            sb.append( "     decay rate: " + decay_rates[i] );
            sb.append( " " + time_units + "^-1" + "\n" );
            sb.append( "     diffusion length scale: " + Math.sqrt( diffusion_coefficients[i] / ( 1e-12 + decay_rates[i] ) ) );
            sb.append( " " + spatial_units + "\n" );
            //            sb.append( "     initial condition: " + options.initial_condition_vector[i] );
            sb.append( " " + density_units[i] + "\n" );
            sb.append( "     boundary condition: " + options.Dirichlet_condition_vector[i] );
            sb.append( " " + density_units[i] + " (enabled: " );
            if( dirichlet_activation_vector[i] == true )
            {
                sb.append( "true" );
            }
            else
            {
                sb.append( "false" );
            }
            sb.append( ")" + "\n" );
        }
        sb.append( "\n" );

        return sb.toString();
    }
}
