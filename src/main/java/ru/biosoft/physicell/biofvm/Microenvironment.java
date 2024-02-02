package ru.biosoft.physicell.biofvm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;

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
    private Set<BasicAgent> agents = new HashSet<>();

    public <T extends BasicAgent> Set<T> getAgents(Class<T> clazz)
    {
        return (Set<T>)agents;
    }

    public Set<BasicAgent> getAgents()
    {
        return agents;
    }

    public void addAgent(BasicAgent agent)
    {
        agents.add( agent );
    }

    public void removeAgent(BasicAgent agent)
    {
        agents.remove( agent );
    }

    public int getAgentsCount()
    {
        return agents.size();
    }

    public MicroenvironmentOptions options;
    public MicroenvironmentOptions getOptions()
    {
        return options;
    }
    private DiffusionDecaySolver solver;

    public void setSolver(DiffusionDecaySolver solver)
    {
        this.solver = solver;
    }

    /*! For internal use and accelerations in solvers */
    double[][] temporary_density_vectors1;
    /*! For internal use and accelerations in solvers */
    double[][] temporary_density_vectors2;
    /*! for internal use in bulk source/sink solvers */
    double[][] bulk_source_sink_solver_temp1;
    double[][] bulk_source_sink_solver_temp2;
    double[][] bulk_source_sink_solver_temp3;
    boolean bulk_source_sink_solver_setup_done;

    /*! stores pointer to current density solutions. Access via operator() functions. */
    public double[][] density;
    double[][][] gradient_vectors;
    boolean[] gradient_vector_computed;

    /*! helpful for solvers -- resize these whenever adding/removing substrates */
    double[] one;
    double[] zero;
    double[] one_half;
    double[] one_third;

    /*! for internal use in diffusion solvers : these make the solvers safe across microenvironments ""*/
    double[][] thomas_temp1;
    double[][] thomas_temp2;
    double[] thomas_constant1x;
    double[] thomas_constant1y;
    double[] thomas_constant1z;
    double[] thomas_neg_constant1x;
    double[] thomas_neg_constant1y;
    double[] thomas_neg_constant1z;

    boolean thomas_setup_done;
    int thomas_i_jump;
    int thomas_j_jump;
    int thomas_k_jump;

    double[] thomas_constant1;
    double[] thomas_constant1a;
    double[] thomas_constant2;
    double[] thomas_constant3;
    double[] thomas_constant3a;
    double[][] thomas_denomx;
    double[][] thomas_cx;
    double[][] thomas_denomy;
    double[][] thomas_cy;
    double[][] thomas_denomz;
    double[][] thomas_cz;
    boolean diffusion_solver_setup_done;

    // on "resize density" type operations, need to extend all of these 

    /*
    std::vector<int> dirichlet_indices; 
    std::vector< std::vector<double> > dirichlet_value_vectors; 
    std::vector<bool> dirichlet_node_map; 
    */
    double[][] dirichlet_value_vectors;
    boolean[] dirichlet_activation_vector;
    /* new in Version 1.7.0 -- activation vectors can be specified on a voxel-by-voxel basis */
    boolean[][] dirichlet_activation_vectors;

    /*! The mesh for the diffusing quantities */
    public CartesianMesh mesh;
    public AgentContainer agentContainer = new AgentContainer();
    public String spatialUnits;
    public String timeUnits;
    public String name;

    // diffusing entities 
    public String[] density_names;
    String[] density_units;

    // coefficients 
    public double[] diffusion_coefficients;
    public double[] decay_rates;
    double[][] supply_target_densities_times_supply_rates;
    double[][] supply_rates;
    double[][] uptake_rates;

    public Microenvironment(String name, String timeUnits, String spatialUnits)
    {
        this();
        this.name = name;
        this.spatialUnits = spatialUnits;
        this.timeUnits = timeUnits;
        mesh.units = spatialUnits;
    }

    public Microenvironment(String name, double size, double nodeSize, String timeUnits, String spatialUnits)
    {
        this();
        this.name = name;
        resizeSpace( 0, size, 0, size, 0, size, nodeSize, nodeSize, nodeSize );
        this.spatialUnits = spatialUnits;
        this.timeUnits = timeUnits;
        mesh.units = spatialUnits;
    }

    public Microenvironment()
    {
        name = "unnamed";
        spatialUnits = "none";
        timeUnits = "none";

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
        density = temporary_density_vectors1;

        gradient_vectors = new double[mesh.voxels.length][1][3];
        gradient_vector_computed = new boolean[mesh.voxels.length];

        density_names = new String[] {"unnamed"};
        density_units = new String[] {"none"};

        diffusion_coefficients = new double[number_of_densities()];
        decay_rates = new double[number_of_densities()];
        one_half = new double[] {0.5};
        one_third = new double[] {1.0 / 3.0};

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
    }

    public void addDensity(String name, String units, double diffusion_constant, double decay_rate)
    {
        // fix in PhysiCell preview November 2017 
        // default_microenvironment_options.use_oxygen_as_first_field = false; 

        zero = VectorUtil.push_back( zero, 0 );
        one = VectorUtil.push_back( one, 1 );
        one_half = VectorUtil.push_back( one_half, 0.5 );
        one_third = VectorUtil.push_back( one_third, 1.0 / 3.0 );

        // update units
        density_names = VectorUtil.push_back( density_names, name );
        density_units = VectorUtil.push_back( density_units, units );

        // update coefficients 
        diffusion_coefficients = VectorUtil.push_back( diffusion_coefficients, diffusion_constant );
        decay_rates = VectorUtil.push_back( decay_rates, decay_rate );

        // update sources and such 
        for( int i = 0; i < temporary_density_vectors1.length; i++ )
        {
            temporary_density_vectors1[i] = VectorUtil.push_back( temporary_density_vectors1[i], 0.0 );
            temporary_density_vectors2[i] = VectorUtil.push_back( temporary_density_vectors2[i], 0.0 );
        }

        // resize the gradient data structures 
        for( int k = 0; k < mesh.voxels.length; k++ )
        {
            gradient_vectors[k] = VectorUtil.resize( gradient_vectors[k], number_of_densities() );
            for( int i = 0; i < number_of_densities(); i++ )
            {
                gradient_vectors[k][i] = VectorUtil.resize( gradient_vectors[k][i], 3 );
            }
        }
        gradient_vector_computed = VectorUtil.resize( gradient_vector_computed, mesh.voxels.length );

        dirichlet_value_vectors = VectorUtil.assign( mesh.voxels.length, one );
        //        dirichlet_value_vectors.assign( mesh.voxels.size(), one );
        dirichlet_activation_vector = VectorUtil.push_back( dirichlet_activation_vector, false );
        //        dirichlet_activation_vectors.assign( mesh.voxels.size(), dirichlet_activation_vector );
        dirichlet_activation_vectors = VectorUtil.assign( mesh.voxels.length, dirichlet_activation_vector );

        // fix in PhysiCell preview November 2017 
        options.Dirichlet_condition_vector = VectorUtil.push_back( options.Dirichlet_condition_vector, 1.0 ); // = one; 
        options.Dirichlet_activation_vector = VectorUtil.push_back( options.Dirichlet_activation_vector, false ); // assign( number_of_densities(), false ); 

        options.initial_condition_vector = VectorUtil.push_back( options.initial_condition_vector, 1.0 );

        options.Dirichlet_all = VectorUtil.push_back( options.Dirichlet_all, true );
        options.Dirichlet_xmin = VectorUtil.push_back( options.Dirichlet_xmin, false );
        options.Dirichlet_xmax = VectorUtil.push_back( options.Dirichlet_xmax, false );
        options.Dirichlet_ymin = VectorUtil.push_back( options.Dirichlet_ymin, false );
        options.Dirichlet_ymax = VectorUtil.push_back( options.Dirichlet_ymax, false );
        options.Dirichlet_zmin = VectorUtil.push_back( options.Dirichlet_zmin, false );
        options.Dirichlet_zmax = VectorUtil.push_back( options.Dirichlet_zmax, false );

        options.Dirichlet_xmin_values = VectorUtil.push_back( options.Dirichlet_xmin_values, 1.0 );
        options.Dirichlet_xmax_values = VectorUtil.push_back( options.Dirichlet_xmax_values, 1.0 );
        options.Dirichlet_ymin_values = VectorUtil.push_back( options.Dirichlet_ymin_values, 1.0 );
        options.Dirichlet_ymax_values = VectorUtil.push_back( options.Dirichlet_ymax_values, 1.0 );
        options.Dirichlet_zmin_values = VectorUtil.push_back( options.Dirichlet_zmin_values, 1.0 );       
        options.Dirichlet_zmax_values = VectorUtil.push_back( options.Dirichlet_zmax_values, 1.0 );
    }

    public double[] get(int n)
    {
        return density[n];
    }

    public void simulate_cell_sources_and_sinks(Set<BasicAgent> agents, double dt)
    {
        for( BasicAgent agent : agents )
            agent.simulateSecretionUptake( this, dt );
    }

    public void simulate_cell_sources_and_sinks(double dt)
    {
        simulate_cell_sources_and_sinks( agents, dt );
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
                        getDensity( i )[j] = dirichlet_value_vectors[i][j];
                    }
                }

            }
        }
        return;
    }

    public void write_to_matlab(String filename)
    {
        int number_of_data_entries = mesh.voxels.length;
        int size_of_each_datum = 3 + 1 + density[0].length;
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
                for( int j = 0; j < density[i].length; j++ )
                {
                    bw.append( String.valueOf( density[i][j] ) + "\t" );
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
        return density[voxel_index( i, j, k )];
    }

    /*! access the density vector at [x,y,z](n) */
    public double[] getDensity(int voxel_index)
    {
        return density[voxel_index];
    }

    int nearest_voxel_index(double[] position)
    {
        return mesh.nearest_voxel_index( position );
    }

    /**
     * Changes index density characteristics: name, units, diffusion and decay rate
     * Microenvironment always contains one density after creation. If you need more densities use addDensity method
     */
    public void setDensity(int index, String name, String units, double diffussion, double decay)
    {
        // fix in PhysiCell preview November 2017 
        //        if( index == 0 )
        //        {
        //            options.use_oxygen_as_first_field = false;//TODO: check
        //        }

        //        VectorUtil.resize( dirichlet_activation_vector, index )
        density_names[index] = name;
        density_units[index] = units;
        diffusion_coefficients[index] = diffussion;
        decay_rates[index] = decay;
    }

    public int number_of_densities()
    {
        return density[0].length;
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
        resizeSpace( x_start, x_end, y_start, y_end, z_start, z_end, dx_new, dx_new, dx_new );
    }

    public void resizeSpace(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double x_nodes,
            double y_nodes, double z_nodes)
    {
        mesh.resize( x_start, x_end, y_start, y_end, z_start, z_end, x_nodes, y_nodes, z_nodes );
        temporary_density_vectors1 = new double[mesh.voxels.length][zero.length];
        density = temporary_density_vectors1;
        temporary_density_vectors2 = new double[mesh.voxels.length][zero.length];
        gradient_vectors = new double[mesh.voxels.length][number_of_densities()][3];
        gradient_vector_computed = new boolean[mesh.voxels.length];
        dirichlet_value_vectors = VectorUtil.assign( mesh.voxels.length, one );
        dirichlet_activation_vectors = VectorUtil.assign( mesh.voxels.length, dirichlet_activation_vector );
    }

    void resize_space(int x_nodes, int y_nodes, int z_nodes)
    {
        mesh.resize( x_nodes, y_nodes, z_nodes );
        temporary_density_vectors1 = new double[mesh.voxels.length][];
        density = temporary_density_vectors1;
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
            sb.append( " " + spatialUnits + "^2 / " + timeUnits + "\n" );
            sb.append( "     decay rate: " + decay_rates[i] );
            sb.append( " " + timeUnits + "^-1" + "\n" );
            sb.append( "     diffusion length scale: " + Math.sqrt( diffusion_coefficients[i] / ( 1e-12 + decay_rates[i] ) ) );
            sb.append( " " + spatialUnits + "\n" );
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

    public int findDensityIndex(String name)
    {
        for( int i = 0; i < density_names.length; i++ )
        {
            if( density_names[i].equals( name ) )
            {
                return i;
            }
        }
        return -1;
    }

    public void compute_all_gradient_vectors()
    {
        double two_dx = mesh.dx;
        double two_dy = mesh.dy;
        double two_dz = mesh.dz;
        boolean gradient_constants_defined = false;
        if( gradient_constants_defined == false )
        {
            two_dx *= 2.0;
            two_dy *= 2.0;
            two_dz *= 2.0;
            gradient_constants_defined = true;
        }

        //        #pragma omp parallel for 
        for( int k = 0; k < mesh.z_coordinates.length; k++ )
        {
            for( int j = 0; j < mesh.y_coordinates.length; j++ )
            {
                // endcaps 
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int i = 0;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][0] = density[n + thomas_i_jump][q];
                    gradient_vectors[n][q][0] -= density[n][q];
                    gradient_vectors[n][q][0] /= mesh.dx;

                    gradient_vector_computed[n] = true;
                }
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int i = mesh.x_coordinates.length - 1;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][0] = ( density )[n][q];
                    gradient_vectors[n][q][0] -= ( density )[n - thomas_i_jump][q];
                    gradient_vectors[n][q][0] /= mesh.dx;

                    gradient_vector_computed[n] = true;
                }

                for( int i = 1; i < mesh.x_coordinates.length - 1; i++ )
                {
                    for( int q = 0; q < number_of_densities(); q++ )
                    {
                        int n = voxel_index( i, j, k );
                        // x-derivative of qth substrate at voxel n
                        gradient_vectors[n][q][0] = ( density )[n + thomas_i_jump][q];
                        gradient_vectors[n][q][0] -= ( density )[n - thomas_i_jump][q];
                        gradient_vectors[n][q][0] /= two_dx;

                        gradient_vector_computed[n] = true;
                    }
                }

            }
        }

        //        #pragma omp parallel for 
        for( int k = 0; k < mesh.z_coordinates.length; k++ )
        {
            for( int i = 0; i < mesh.x_coordinates.length; i++ )
            {
                // endcaps 
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int j = 0;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][1] = ( density )[n + thomas_j_jump][q];
                    gradient_vectors[n][q][1] -= ( density )[n][q];
                    gradient_vectors[n][q][1] /= mesh.dy;

                    gradient_vector_computed[n] = true;
                }
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int j = mesh.y_coordinates.length - 1;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][1] = ( density )[n][q];
                    gradient_vectors[n][q][1] -= ( density )[n - thomas_j_jump][q];
                    gradient_vectors[n][q][1] /= mesh.dy;

                    gradient_vector_computed[n] = true;
                }

                for( int j = 1; j < mesh.y_coordinates.length - 1; j++ )
                {
                    for( int q = 0; q < number_of_densities(); q++ )
                    {
                        int n = voxel_index( i, j, k );
                        // y-derivative of qth substrate at voxel n
                        gradient_vectors[n][q][1] = ( density )[n + thomas_j_jump][q];
                        gradient_vectors[n][q][1] -= ( density )[n - thomas_j_jump][q];
                        gradient_vectors[n][q][1] /= two_dy;
                        gradient_vector_computed[n] = true;
                    }
                }

            }
        }

        // don't bother computing z component if there is no z-directoin 
        if( mesh.z_coordinates.length == 1 )
        {
            return;
        }

        //        #pragma omp parallel for 
        for( int j = 0; j < mesh.y_coordinates.length; j++ )
        {
            for( int i = 0; i < mesh.x_coordinates.length; i++ )
            {
                // endcaps 
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int k = 0;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][2] = ( density )[n + thomas_k_jump][q];
                    gradient_vectors[n][q][2] -= ( density )[n][q];
                    gradient_vectors[n][q][2] /= mesh.dz;

                    gradient_vector_computed[n] = true;
                }
                for( int q = 0; q < number_of_densities(); q++ )
                {
                    int k = mesh.z_coordinates.length - 1;
                    int n = voxel_index( i, j, k );
                    // x-derivative of qth substrate at voxel n
                    gradient_vectors[n][q][2] = ( density )[n][q];
                    gradient_vectors[n][q][2] -= ( density )[n - thomas_k_jump][q];
                    gradient_vectors[n][q][2] /= mesh.dz;

                    gradient_vector_computed[n] = true;
                }

                for( int k = 1; k < mesh.z_coordinates.length - 1; k++ )
                {
                    for( int q = 0; q < number_of_densities(); q++ )
                    {
                        int n = voxel_index( i, j, k );
                        // y-derivative of qth substrate at voxel n
                        gradient_vectors[n][q][2] = ( density )[n + thomas_k_jump][q];
                        gradient_vectors[n][q][2] -= ( density )[n - thomas_k_jump][q];
                        gradient_vectors[n][q][2] /= two_dz;
                        gradient_vector_computed[n] = true;
                    }
                }

            }
        }
    }

    public double[] nearest_density_vector(double[] position)
    {
        return ( density )[mesh.nearest_voxel_index( position )];
    }

    public double[] nearest_density_vector(int voxel_index)
    {
        return ( density )[voxel_index];
    }

    public void add_dirichlet_node(int voxel_index, double[] value)
    {
        mesh.voxels[voxel_index].isDirichlet = true;
        /*
        dirichlet_indices.push_back( voxel_index );
        dirichlet_value_vectors.push_back( value ); 
        */

        dirichlet_value_vectors[voxel_index] = value; // .assign( mesh.voxels.size(), one ); 
    }

    public double[][] gradient_vector(int n)
    {
        // if the gradient has not yet been computed, then do it!
        if( gradient_vector_computed[n] == false )
        {
            compute_gradient_vector( n );
        }
        
        
        return gradient_vectors[n];
    }

    void compute_gradient_vector( int n )
    {
        double two_dx = mesh.dx; 
        double two_dy = mesh.dy; 
        double two_dz = mesh.dz; 
        boolean gradient_constants_defined = false; 
        int[] indices = new int[3];//(3,0);
        
        if( gradient_constants_defined == false )
        {
            two_dx *= 2.0; 
            two_dy *= 2.0; 
            two_dz *= 2.0;
            gradient_constants_defined = true; 
        }   
        
        indices = cartesian_indices( n );
        
        // d/dx 
        if( indices[0] > 0 && indices[0] < mesh.x_coordinates.length-1 )
        {
            for( int q=0; q < number_of_densities() ; q++ )
            {
                gradient_vectors[n][q][0] = ( density )[n + thomas_i_jump][q];
                gradient_vectors[n][q][0] -= ( density )[n - thomas_i_jump][q];
                gradient_vectors[n][q][0] /= two_dx; 
                gradient_vector_computed[n] = true; 
            }
        }
        
        // don't bother computing y and z component if there is no y-direction. (1D)
        if( mesh.y_coordinates.length == 1 )
        return; 
        
        // d/dy 
        if( indices[1] > 0 && indices[1] < mesh.y_coordinates.length-1 )
        {
            for( int q = 0; q < number_of_densities(); q++ )
            {
                gradient_vectors[n][q][1] = ( density )[n + thomas_j_jump][q];
                gradient_vectors[n][q][1] -= ( density )[n - thomas_j_jump][q];
                gradient_vectors[n][q][1] /= two_dy; 
                gradient_vector_computed[n] = true; 
            }
        }
        
        // don't bother computing z component if there is no z-direction (2D) 
        if( mesh.z_coordinates.length == 1 )
         return;
        
        // d/dz 
        if( indices[2] > 0 && indices[2] < mesh.z_coordinates.length-1 )
        {
            for( int q=0; q < number_of_densities() ; q++ )
            {
                gradient_vectors[n][q][2] = ( density )[n + thomas_k_jump][q];
                gradient_vectors[n][q][2] -= ( density )[n - thomas_k_jump][q];
                gradient_vectors[n][q][2] /= two_dz; 
                gradient_vector_computed[n] = true; 
            }
        }
        
    }

    int[] cartesian_indices(int n)
    {
        return mesh.cartesian_indices( n );
    }

    public static void initialize_microenvironment(Microenvironment microenvironment)
    {
        MicroenvironmentOptions default_microenvironment_options = microenvironment.options;
        //        Microenvironment microenvironment = new Microenvironment(m.name, m.timeUnits, m.spatialUnits);
        // register the diffusion solver 
        if( microenvironment.options.simulate_2D == true )
        {
//            microenvironment.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_2D; 
        }
        else
        {
            microenvironment.solver = new ConstantCoefficientsLOD3D();// diffusion_decay_solver__constant_coefficients_LOD_3D; 
        }
        
        // set the default substrate to oxygen (with typical units of mmHg)
        if( microenvironment.options.use_oxygen_as_first_field == true )
        {
            microenvironment.setDensity(0, "oxygen" , "mmHg", 1E5, 0.1 );
//            result.diffusion_coefficients[0] = 1e5; 
//            result.decay_rates[0] = 0.1; 
        }
        
        // resize the microenvironment  
        if( default_microenvironment_options.simulate_2D == true )
        {
            default_microenvironment_options.Z_range[0] = -default_microenvironment_options.dz/2.0; 
            default_microenvironment_options.Z_range[1] = default_microenvironment_options.dz/2.0;
        }
        microenvironment.resizeSpace( default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1] , 
            default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1], 
            default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1], 
            default_microenvironment_options.dx,default_microenvironment_options.dy,default_microenvironment_options.dz );
            
        // set units
//        microenvironment.spatial_units = default_microenvironment_options.spatial_units;
//        microenvironment.time_units = default_microenvironment_options.time_units;
        microenvironment.mesh.units = default_microenvironment_options.spatial_units;

        // set the initial densities to the values set in the initial_condition_vector
        
        // if the initial condition vector has not been set, use the Dirichlet condition vector 
        if( default_microenvironment_options.initial_condition_vector.length != 
            microenvironment.number_of_densities() )
        {
           System.out.println( "BioFVM Warning: Initial conditions not set. "
                      + "                Using Dirichlet condition vector to set initial substrate values!" 
                      + "                In the future, set default_microenvironment_options.initial_condition_vector.");  
            default_microenvironment_options.initial_condition_vector = default_microenvironment_options.Dirichlet_condition_vector; 
        }

        // set the initial condition 
        for(  int n=0; n < microenvironment.number_of_voxels() ; n++ )
        {
            microenvironment.density[n] = default_microenvironment_options.initial_condition_vector.clone();
        }

        // now, figure out which sides have BCs (for at least one substrate): 

        boolean xmin = false; 
        boolean xmax = false; 
        boolean ymin = false; 
        boolean ymax = false; 
        boolean zmin = false; 
        boolean zmax = false; 
        
        if( default_microenvironment_options.outer_Dirichlet_conditions == true )
        {
            for( int n=0 ; n < microenvironment.number_of_densities() ; n++ )
            {
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_xmin[n] )
                    { xmin = true; }
                
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_xmax[n] )
                    { xmax = true; }
                
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_ymin[n] )
                    { ymin = true; }
                
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_ymax[n] )
                    { ymax = true; }
                    
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_zmin[n] )
                    { zmin = true; }
                
                if( default_microenvironment_options.Dirichlet_all[n] || 
                    default_microenvironment_options.Dirichlet_zmax[n] )
                    { zmax = true; }
            }
            
            // add the Dirichlet nodes in the right places 
            
        }
        System.out.println( "which boundaries?"); 
        System.out.println( xmin + " " + xmax + " " + ymin + " " + ymax + " " + zmin + " " + zmax);// << std::endl; 

        // add the Dirichlet nodes in the right places 
        // now, go in and set the values 
        if( default_microenvironment_options.outer_Dirichlet_conditions == true ) 
        {
            // set xmin if xmin = true or all = true 
            if( xmin == true )
            {
                for( int k = 0; k < microenvironment.mesh.z_coordinates.length; k++ )
                {
                    int I = 0; 
                    // set Dirichlet conditions along the xmin outer edges 
                    for( int j = 0; j < microenvironment.mesh.y_coordinates.length; j++ )
                    {
                        // set the value 
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(I,j,k) , default_microenvironment_options.Dirichlet_xmin_values );
                        
                        // set the activation 
                        microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(I,j,k) , 
                        default_microenvironment_options.Dirichlet_xmin ); 
                        
                    }
                }
            }           
            
            // set xmax if xmax = true or all = true 
            if( xmax == true )
            {
                for( int k = 0; k < microenvironment.mesh.z_coordinates.length; k++ )
                {
                    int I = microenvironment.mesh.x_coordinates.length - 1;
                    ;
                    // set Dirichlet conditions along the xmax outer edges 
                    for( int j = 0; j < microenvironment.mesh.y_coordinates.length; j++ )
                    {
                        // set the values 
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(I,j,k) , default_microenvironment_options.Dirichlet_xmax_values );
                        
                        // set the activation 
                        microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(I,j,k) , 
                        default_microenvironment_options.Dirichlet_xmax ); 
                    }
                }
            }           
            
            // set ymin if ymin = true or all = true 
            if( ymin == true )
            {
                for( int k = 0; k < microenvironment.mesh.z_coordinates.length; k++ )
                {
                    int J = 0; // microenvironment.mesh.x_coordinates.size()-1;; 
                    // set Dirichlet conditions along the ymin outer edges 
                    for( int i = 0; i < microenvironment.mesh.x_coordinates.length; i++ )
                    {
                        // set the values 
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,J,k) , default_microenvironment_options.Dirichlet_ymin_values );
                        
                        // set the activation 
                        microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,J,k) , 
                        default_microenvironment_options.Dirichlet_ymin ); 
                    }
                }
            }   
            
            // set ymzx if ymax = true or all = true; 
            if( ymax == true )
            {
                for( int k = 0; k < microenvironment.mesh.z_coordinates.length; k++ )
                {
                    int J = microenvironment.mesh.y_coordinates.length - 1;
                    ;
                    // set Dirichlet conditions along the ymin outer edges 
                    for( int i = 0; i < microenvironment.mesh.x_coordinates.length; i++ )
                    {
                        // set the value 
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,J,k) , default_microenvironment_options.Dirichlet_ymax_values );
                        
                        // set the activation 
                        microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,J,k) , 
                        default_microenvironment_options.Dirichlet_ymax ); 
                    }
                }
            }   
            
            // if not 2D:
            if( default_microenvironment_options.simulate_2D == false )
            {
                // set zmin if zmin = true or all = true 
                if( zmin == true )
                {
                    for( int j = 0; j < microenvironment.mesh.y_coordinates.length; j++ )
                    {
                        int K = 0; // microenvironment.mesh.z_coordinates.size()-1;; 
                        // set Dirichlet conditions along the ymin outer edges 
                        for( int i = 0; i < microenvironment.mesh.x_coordinates.length; i++ )
                        {
                            // set the value 
                            microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,K) , default_microenvironment_options.Dirichlet_zmin_values );
                        
                            // set the activation 
                            microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,K) , 
                            default_microenvironment_options.Dirichlet_zmin ); 
                        }
                    }
                }               
                
                // set zmax if zmax = true or all = true 
                if( zmax == true )
                {
                    for( int j = 0; j < microenvironment.mesh.y_coordinates.length; j++ )
                    {
                        int K = microenvironment.mesh.z_coordinates.length - 1;
                        ;
                        // set Dirichlet conditions along the ymin outer edges 
                        for( int i = 0; i < microenvironment.mesh.x_coordinates.length; i++ )
                        {
                            // set the value 
                            microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,K) , default_microenvironment_options.Dirichlet_zmax_values );
                            
                            // set the activation 
                            microenvironment.set_substrate_dirichlet_activation( microenvironment.voxel_index(i,j,K) , 
                            default_microenvironment_options.Dirichlet_zmax );                      
                        }
                    }
                }               
            }
            
        }
        
    /*
        if( default_microenvironment_options.outer_Dirichlet_conditions == true ) 
        {
            
            for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ )
            {
                // set Dirichlet conditions along the 4 outer edges 
                for( unsigned int i=0 ; i < microenvironment.mesh.x_coordinates.size() ; i++ )
                {
                    int J = microenvironment.mesh.y_coordinates.size()-1;
                    microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,0,k) , default_microenvironment_options.Dirichlet_condition_vector );
                    microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,J,k) , default_microenvironment_options.Dirichlet_condition_vector );
                }
                int I = microenvironment.mesh.x_coordinates.size()-1;
                for( unsigned int j=1; j < microenvironment.mesh.y_coordinates.size()-1 ; j++ )
                {
                    microenvironment.add_dirichlet_node( microenvironment.voxel_index(0,j,k) , default_microenvironment_options.Dirichlet_condition_vector );
                    microenvironment.add_dirichlet_node( microenvironment.voxel_index(I,j,k) , default_microenvironment_options.Dirichlet_condition_vector );
                }       
            }
            // if 3-D, also along the corresponding additional faces 
            if( default_microenvironment_options.simulate_2D == false )
            {
                int K = microenvironment.mesh.z_coordinates.size()-1; 
                for( unsigned int j=1 ; j < microenvironment.mesh.y_coordinates.size()-1 ; j++ )
                {
                    for( unsigned int i=1; i < microenvironment.mesh.x_coordinates.size()-1 ; i++ )
                    {
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,0) , default_microenvironment_options.Dirichlet_condition_vector );
                        microenvironment.add_dirichlet_node( microenvironment.voxel_index(i,j,K) , default_microenvironment_options.Dirichlet_condition_vector );
                    }   
                }   
            }
            
        }
    */
        
        // April 2023: no longer necessary after flipping our approach and doing an "additive" instead of "subtractive" DCs handling. I.e., we assume DC activation is false by default; make true on-demand.

        // // set the Dirichlet condition activation vector to match the microenvironment options 
        // for( int i=0 ; i < default_microenvironment_options.Dirichlet_activation_vector.size(); i++ )
        // {
        //  microenvironment.set_substrate_dirichlet_activation( i , default_microenvironment_options.Dirichlet_activation_vector[i] ); 
        // }
        
        microenvironment.displayInformation();
    }

    void set_substrate_dirichlet_activation( int substrate_index , boolean new_value )
    {
        dirichlet_activation_vector[substrate_index] = new_value; 
        
        for( int n = 0 ; n < mesh.voxels.length ; n++ )
        {
            dirichlet_activation_vectors[n][substrate_index] = new_value;
        }
    }

    void set_substrate_dirichlet_activation(int index, boolean[] new_value)
    {
        dirichlet_activation_vectors[index] = new_value.clone();
    }
}
