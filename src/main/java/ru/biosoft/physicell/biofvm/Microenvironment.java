package ru.biosoft.physicell.biofvm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

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
    private Map<String, Integer> substrateToIndex;
    public MicroenvironmentOptions options;

    public <T extends BasicAgent> Set<T> getAgents(Class<T> clazz)
    {
        return (Set<T>)agents;
    }

    public int getSubstrateIndex(String name)
    {
        return substrateToIndex.get( name );
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

    public MicroenvironmentOptions getOptions()
    {
        return options;
    }
    private DiffusionDecaySolver solver;

    public void setSolver(DiffusionDecaySolver solver)
    {
        this.solver = solver;
    }

    public DiffusionDecaySolver getSolver()
    {
        return solver;
    }

    /*! For internal use and accelerations in solvers */
    double[][] tempDensity1;
    /*! For internal use and accelerations in solvers */
    double[][] tempDensity2;
    /*! for internal use in bulk source/sink solvers */
    double[][] bulk_source_sink_solver_temp1;
    double[][] bulk_source_sink_solver_temp2;
    double[][] bulk_source_sink_solver_temp3;
    boolean bulk_source_sink_solver_setup_done;

    /*! stores pointer to current density solutions. Access via operator() functions. */
    public double[][] density;
    double[][][] gradients;
    boolean[] gradientComputed;

    /*! helpful for solvers -- resize these whenever adding/removing substrates */
    double[] one;
    double[] zero;
    double[] oneHalf;
    double[] oneThird;

    /*! for internal use in diffusion solvers : these make the solvers safe across microenvironments ""*/
    double[][] thomasTemp1;
    double[][] thomasTemp2;
    double[] thomas_constant1x;
    double[] thomas_constant1y;
    double[] thomas_constant1z;
    double[] thomas_neg_constant1x;
    double[] thomas_neg_constant1y;
    double[] thomas_neg_constant1z;

    boolean thomasSetup;
    int iJump;
    int jJump;
    int kJump;

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
    boolean solverSetup;

    // on "resize density" type operations, need to extend all of these 

    /*
    std::vector<int> dirichlet_indices; 
    std::vector< std::vector<double> > dirichlet_value_vectors; 
    std::vector<bool> dirichlet_node_map; 
    */
    double[][] dirichletValue;
    boolean[] dirichletActivation;
    /* new in Version 1.7.0 -- activation vectors can be specified on a voxel-by-voxel basis */
    boolean[][] dirichletActivations;

    /*! The mesh for the diffusing quantities */
    public CartesianMesh mesh;
    public AgentContainer agentContainer = new AgentContainer();
    public String spatialUnits;
    public String timeUnits;
    public String name;

    // diffusing entities 
    public String[] densityNames;
    String[] densityUnits;

    // coefficients 
    public double[] diffusionCoefficients;
    public double[] decayRates;
    double[][] supply_target_densities_times_supply_rates;
    double[][] supplyRates;
    double[][] uptakeRates;

    public double time;

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
        thomasSetup = false;
        solverSetup = false;

        //        solver = new TornadoSolverParallel();
        //        solver = new TornadoSolverParallel2();
        //        solver = new TornadoSolver();
        solver = new ConstantCoefficientsLOD3D();

        mesh = new CartesianMesh();
        mesh.resize( 1, 1, 1 );

        one = new double[] {1};
        zero = new double[] {0};

        tempDensity1 = new double[mesh.voxels.length][1];
        tempDensity2 = new double[mesh.voxels.length][1];
        density = tempDensity1;

        gradients = new double[mesh.voxels.length][1][3];
        gradientComputed = new boolean[mesh.voxels.length];

        densityNames = new String[] {"unnamed"};
        densityUnits = new String[] {"none"};

        diffusionCoefficients = new double[numberDensities()];
        decayRates = new double[numberDensities()];
        oneHalf = new double[] {0.5};
        oneThird = new double[] {1.0 / 3.0};

        dirichletValue = VectorUtil.assign( mesh.voxels.length, one );
        dirichletActivation = new boolean[] {false};
        dirichletActivations = VectorUtil.assign( 1, dirichletActivation );

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

    public void addDensity(String name, String units, double diffusionConstant, double decayRate)
    {
        // fix in PhysiCell preview November 2017 
        // default_microenvironment_options.use_oxygen_as_first_field = false; 

        zero = VectorUtil.push_back( zero, 0 );
        one = VectorUtil.push_back( one, 1 );
        oneHalf = VectorUtil.push_back( oneHalf, 0.5 );
        oneThird = VectorUtil.push_back( oneThird, 1.0 / 3.0 );

        // update units
        densityNames = VectorUtil.push_back( densityNames, name );
        densityUnits = VectorUtil.push_back( densityUnits, units );

        // update coefficients 
        diffusionCoefficients = VectorUtil.push_back( diffusionCoefficients, diffusionConstant );
        decayRates = VectorUtil.push_back( decayRates, decayRate );

        // update sources and such 
        for( int i = 0; i < tempDensity1.length; i++ )
        {
            tempDensity1[i] = VectorUtil.push_back( tempDensity1[i], 0.0 );
            tempDensity2[i] = VectorUtil.push_back( tempDensity2[i], 0.0 );
        }

        // resize the gradient data structures 
        for( int k = 0; k < mesh.voxels.length; k++ )
        {
            gradients[k] = VectorUtil.resize( gradients[k], numberDensities() );
            for( int i = 0; i < numberDensities(); i++ )
            {
                gradients[k][i] = VectorUtil.resize( gradients[k][i], 3 );
            }
        }
        gradientComputed = VectorUtil.resize( gradientComputed, mesh.voxels.length );

        dirichletValue = VectorUtil.assign( mesh.voxels.length, one );
        //        dirichlet_value_vectors.assign( mesh.voxels.size(), one );
        dirichletActivation = VectorUtil.push_back( dirichletActivation, false );
        //        dirichlet_activation_vectors.assign( mesh.voxels.size(), dirichlet_activation_vector );
        dirichletActivations = VectorUtil.assign( mesh.voxels.length, dirichletActivation );

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

    public void simulateSourcesSinks(Set<BasicAgent> agents, double dt)
    {
        for( BasicAgent agent : agents )
            agent.simulateSecretionUptake( this, dt );
    }

    public void simulateSourcesSinks(double dt)
    {
        simulateSourcesSinks( agents, dt );
    }

    public void simulateDiffusionDecay(double dt) throws Exception
    {
        if( solver != null )
            solver.solve( this, dt );
    }

    private List<Integer>[] dirichletIndices;

    void initDirichlet()
    {
        dirichletIndices = new List[densityNames.length];
        for (int i=0; i< densityNames.length; i++)
            dirichletIndices[i] = new ArrayList<>();

        for( int i = 0; i < mesh.voxels.length; i++ )
        {
            if( mesh.voxels[i].isDirichlet )
            {
                for( int j = 0; j < dirichletValue[i].length; j++ )
                {
                    if( dirichletActivations[i][j] )
                    {
                        dirichletIndices[j].add( i );
                    }
                }
            }
        }
    }

    void applyDirichletConditions()
    {
        for( int i = 0; i < dirichletIndices.length; i++ )
            for( int j = 0; j < dirichletIndices[i].size(); j++ )
            {
                int k = dirichletIndices[i].get( j );
                getDensity( k )[i] = dirichletValue[k][i];
            }
    }

    void applyDirichletConditions2()
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
                for( int j = 0; j < dirichletValue[i].length; j++ )
                {
                    // if( dirichlet_activation_vector[j] == true )
                    if( dirichletActivations[i][j] == true )
                    {
                        getDensity( i )[j] = dirichletValue[i][j];
                    }
                }
            }
        }
    }

    public void writeDensity(String filename)
    {
        int dataEntries = mesh.voxels.length;
        File f = new File( filename );
        try (BufferedWriter bw = new BufferedWriter( new FileWriter( f ) ))
        {
            bw.append( "X\tY\tZ\tVolume" );
            for( int i = 0; i < this.densityNames.length; i++ )
                bw.append( "\t" + densityNames[i] );

            bw.append( "\n" );
            for( int i = 0; i < dataEntries; i++ )
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

    int getVoxelIndex(int i, int j, int k)
    {
        return mesh.voxel_index( i, j, k );
    }

    /*! access the density vector at  [ X(i),Y(j),Z(k) ] */
    double[] getDensity(int i, int j, int k)
    {
        return density[getVoxelIndex( i, j, k )];
    }

    /*! access the density vector at [x,y,z](n) */
    public double[] getDensity(int voxel_index)
    {
        return density[voxel_index];
    }

    int nearestVoxelIndex(double[] position)
    {
        return mesh.nearestVoxelIndex( position );
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
        densityNames[index] = name;
        densityUnits[index] = units;
        diffusionCoefficients[index] = diffussion;
        decayRates[index] = decay;
    }

    public int numberDensities()
    {
        return density[0].length;
    }

    public int numberVoxels()
    {
        return mesh.voxels.length;
    }

    public Voxel voxels(int voxelIndex)
    {
        return mesh.voxels[voxelIndex];
    }

    public void resizeSpaceUniform(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double dx_new)
    {
        resizeSpace( x_start, x_end, y_start, y_end, z_start, z_end, dx_new, dx_new, dx_new );
    }

    public void resizeSpace(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end, double x_nodes,
            double y_nodes, double z_nodes)
    {
        mesh.resize( x_start, x_end, y_start, y_end, z_start, z_end, x_nodes, y_nodes, z_nodes );
        tempDensity1 = new double[mesh.voxels.length][zero.length];
        density = tempDensity1;
        tempDensity2 = new double[mesh.voxels.length][zero.length];
        gradients = new double[mesh.voxels.length][numberDensities()][3];
        gradientComputed = new boolean[mesh.voxels.length];
        dirichletValue = VectorUtil.assign( mesh.voxels.length, one );
        dirichletActivations = VectorUtil.assign( mesh.voxels.length, dirichletActivation );
    }

    void resizeSpace(int xNodes, int yNodes, int zNodes)
    {
        mesh.resize( xNodes, yNodes, zNodes );
        tempDensity1 = new double[mesh.voxels.length][];
        density = tempDensity1;
        tempDensity2 = new double[mesh.voxels.length][];
        gradients = new double[mesh.voxels.length][numberDensities()][3];
        gradientComputed = new boolean[mesh.voxels.length];
        dirichletValue = VectorUtil.assign( mesh.voxels.length, one );
        dirichletActivations = VectorUtil.assign( mesh.voxels.length, dirichletActivation );

    }

    public String displayInformation()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "\nMicroenvironment summary: " + name + ": \n" );
        sb.append( mesh.display() );
        sb.append( "Densities: (" + numberDensities() + " total)" + "\n" );
        for( int i = 0; i < densityNames.length; i++ )
        {
            sb.append( "   " + densityNames[i] + ":" + "\n" );
            sb.append( "     units: " + densityUnits[i] + "\n" );
            sb.append( "     diffusion coefficient: " + diffusionCoefficients[i] );
            sb.append( " " + spatialUnits + "^2 / " + timeUnits + "\n" );
            sb.append( "     decay rate: " + decayRates[i] );
            sb.append( " " + timeUnits + "^-1" + "\n" );
            sb.append( "     diffusion length scale: " + Math.sqrt( diffusionCoefficients[i] / ( 1e-12 + decayRates[i] ) ) );
            sb.append( " " + spatialUnits + "\n" );
            //            sb.append( "     initial condition: " + options.initial_condition_vector[i] );
            sb.append( " " + densityUnits[i] + "\n" );
            sb.append( "     boundary condition: " + options.Dirichlet_condition_vector[i] );
            sb.append( " " + densityUnits[i] + " (enabled: " );
            if( dirichletActivation[i] == true )
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
        for( int i = 0; i < densityNames.length; i++ )
        {
            if( densityNames[i].equals( name ) )
            {
                return i;
            }
        }
        return -1;
    }

    public void computeAllGradientVectors()
    {
        double two_dx = mesh.dx * 2;
        double two_dy = mesh.dy * 2;
        double two_dz = mesh.dz * 2;
        int xLength = mesh.x_coordinates.length;
        int yLength = mesh.y_coordinates.length;
        int zLength = mesh.z_coordinates.length;
        int numDens = numberDensities();

        IntStream.range( 0, zLength ).parallel().forEach( k ->
        {
            for( int j = 0; j < yLength; j++ )
            {
                for( int q = 0; q < numDens; q++ )
                {
                    int n = getVoxelIndex( 0, j, k );
                    gradients[n][q][0] = ( density[n + iJump][q] - density[n][q] ) / mesh.dx;
                    gradientComputed[n] = true;
                }
                for( int q = 0; q < numDens; q++ )
                {
                    int n = getVoxelIndex( xLength - 1, j, k );
                    gradients[n][q][0] = ( density[n][q] - density[n - iJump][q] ) / mesh.dx;
                    gradientComputed[n] = true;
                }

                for( int i = 1; i < xLength - 1; i++ )
                {
                    for( int q = 0; q < numDens; q++ )
                    {
                        int n = getVoxelIndex( i, j, k );
                        gradients[n][q][0] = ( density[n + iJump][q] - density[n - iJump][q] ) / two_dx;
                        gradientComputed[n] = true;
                    }
                }

            }
        } );

        IntStream.range( 0, zLength ).parallel().forEach( k ->
        {
            for( int i = 0; i < xLength; i++ )
            {
                for( int q = 0; q < numDens; q++ )
                {
                    int n = getVoxelIndex( i, 0, k );
                    gradients[n][q][1] = ( density[n + jJump][q] - density[n][q] ) / mesh.dy;
                    gradientComputed[n] = true;
                }
                for( int q = 0; q < numDens; q++ )
                {
                    int n = getVoxelIndex( i, yLength - 1, k );
                    gradients[n][q][1] = ( density[n][q] - density[n - jJump][q] ) / mesh.dy;
                    gradientComputed[n] = true;
                }

                for( int j = 1; j < yLength - 1; j++ )
                {
                    for( int q = 0; q < numDens; q++ )
                    {
                        int n = getVoxelIndex( i, j, k );
                        gradients[n][q][1] = ( density[n + jJump][q] - density[n - jJump][q] ) / two_dy;
                        gradientComputed[n] = true;
                    }
                }

            }
        } );
        if( zLength == 1 ) // don't bother computing z component if there is no z-directoin 
            return;

        IntStream.range( 0, yLength ).parallel().forEach( j ->
        {
            for( int i = 0; i < xLength; i++ )
            {
                for( int q = 0; q < numDens; q++ )
                {
                    int k = 0;
                    int n = getVoxelIndex( i, j, k );
                    gradients[n][q][2] = ( density[n + kJump][q] - density[n][q] ) / mesh.dz;
                    gradientComputed[n] = true;
                }
                for( int q = 0; q < numDens; q++ )
                {
                    int k = zLength - 1;
                    int n = getVoxelIndex( i, j, k );
                    gradients[n][q][2] = ( density[n][q] - density[n - kJump][q] ) / mesh.dz;
                    gradientComputed[n] = true;
                }

                for( int k = 1; k < zLength - 1; k++ )
                {
                    for( int q = 0; q < numDens; q++ )
                    {
                        int n = getVoxelIndex( i, j, k );
                        gradients[n][q][2] = ( density[n + kJump][q] - density[n - kJump][q] ) / two_dz;
                        gradientComputed[n] = true;
                    }
                }

            }
        } );
    }

    public double[] nearestDensity(double[] position)
    {
        return density[mesh.nearestVoxelIndex( position )];
    }

    public double[] nearestDensity(int voxelIndex)
    {
        return density[voxelIndex];
    }

    public void addDirichletNode(int voxelIndex, double[] value)
    {
        mesh.voxels[voxelIndex].isDirichlet = true;
        dirichletValue[voxelIndex] = value.clone();
    }

    public double[][] getGradient(int n)
    {
        if( !gradientComputed[n] )
            computeGradient( n );

        return gradients[n];
    }

    void computeGradient(int n)
    {
        double two_dx = mesh.dx * 2;
        double two_dy = mesh.dy * 2;
        double two_dz = mesh.dz * 2;
        int[] indices = mesh.cartesian_indices( n );

        if( indices[0] > 0 && indices[0] < mesh.x_coordinates.length - 1 )
        {
            for( int q = 0; q < numberDensities(); q++ )
            {
                gradients[n][q][0] = ( density[n + iJump][q] - density[n - iJump][q] ) / two_dx;
                gradientComputed[n] = true;
            }
        }

        // don't bother computing y and z component if there is no y-direction. (1D)
        if( mesh.y_coordinates.length == 1 )
            return;

        if( indices[1] > 0 && indices[1] < mesh.y_coordinates.length - 1 )
        {
            for( int q = 0; q < numberDensities(); q++ )
            {
                gradients[n][q][1] = ( density[n + jJump][q] - density[n - jJump][q] ) / two_dy;
                gradientComputed[n] = true;
            }
        }

        // don't bother computing z component if there is no z-direction (2D) 
        if( mesh.z_coordinates.length == 1 )
            return;

        if( indices[2] > 0 && indices[2] < mesh.z_coordinates.length - 1 )
        {
            for( int q = 0; q < numberDensities(); q++ )
            {
                gradients[n][q][2] = ( density[n + kJump][q] - density[n - kJump][q] ) / two_dz;
                gradientComputed[n] = true;
            }
        }

    }

    public static void initialize(Microenvironment m)
    {
        MicroenvironmentOptions options = m.options;

        if( options.simulate2D )
        {
            options.Z_range[0] = -options.dz / 2.0;
            options.Z_range[1] = options.dz / 2.0;
        }
        m.resizeSpace( options.X_range[0], options.X_range[1], options.Y_range[0], options.Y_range[1], options.Z_range[0],
                options.Z_range[1], options.dx, options.dy, options.dz );

        m.mesh.units = options.spatial_units;

        if( options.initial_condition_vector.length != m.numberDensities() )
        {
            System.out.println( "BioFVM Warning: Initial conditions not set. "
                    + "                Using Dirichlet condition vector to set initial substrate values!"
                    + "                In the future, set default_microenvironment_options.initial_condition_vector." );
            options.initial_condition_vector = options.Dirichlet_condition_vector;
        }

        // set the initial condition 
        for( int n = 0; n < m.numberVoxels(); n++ )
            m.density[n] = options.initial_condition_vector.clone();

        // now, figure out which sides have BCs (for at least one substrate): 
        boolean xmin = false;
        boolean xmax = false;
        boolean ymin = false;
        boolean ymax = false;
        boolean zmin = false;
        boolean zmax = false;

        if( options.outer_Dirichlet_conditions )
        {
            for( int n = 0; n < m.numberDensities(); n++ )
            {
                if( options.Dirichlet_all[n] || options.Dirichlet_xmin[n] )
                    xmin = true;

                if( options.Dirichlet_all[n] || options.Dirichlet_xmax[n] )
                    xmax = true;

                if( options.Dirichlet_all[n] || options.Dirichlet_ymin[n] )
                    ymin = true;

                if( options.Dirichlet_all[n] || options.Dirichlet_ymax[n] )
                    ymax = true;

                if( options.Dirichlet_all[n] || options.Dirichlet_zmin[n] )
                    zmin = true;

                if( options.Dirichlet_all[n] || options.Dirichlet_zmax[n] )
                    zmax = true;
            }
        }

        int xLength = m.mesh.x_coordinates.length;
        int yLength = m.mesh.y_coordinates.length;
        int zLength = m.mesh.z_coordinates.length;

        if( options.outer_Dirichlet_conditions )
        {
            if( xmin )
            {
                for( int k = 0; k < zLength; k++ )
                {
                    for( int j = 0; j < yLength; j++ )
                    {
                        m.addDirichletNode( m.getVoxelIndex( 0, j, k ), options.Dirichlet_xmin_values );
                        m.setSubstrateActivation( m.getVoxelIndex( 0, j, k ), options.Dirichlet_xmin );

                    }
                }
            }

            if( xmax )
            {
                for( int k = 0; k < zLength; k++ )
                {
                    for( int j = 0; j < yLength; j++ )
                    {
                        m.addDirichletNode( m.getVoxelIndex( xLength - 1, j, k ), options.Dirichlet_xmax_values );
                        m.setSubstrateActivation( m.getVoxelIndex( xLength - 1, j, k ), options.Dirichlet_xmax );
                    }
                }
            }

            if( ymin )
            {
                for( int k = 0; k < zLength; k++ )
                {
                    for( int i = 0; i < xLength; i++ )
                    {
                        m.addDirichletNode( m.getVoxelIndex( i, 0, k ), options.Dirichlet_ymin_values );
                        m.setSubstrateActivation( m.getVoxelIndex( i, 0, k ), options.Dirichlet_ymin );
                    }
                }
            }

            if( ymax )
            {
                for( int k = 0; k < zLength; k++ )
                {
                    for( int i = 0; i < xLength; i++ )
                    {
                        m.addDirichletNode( m.getVoxelIndex( i, yLength - 1, k ), options.Dirichlet_ymax_values );
                        m.setSubstrateActivation( m.getVoxelIndex( i, yLength - 1, k ), options.Dirichlet_ymax );
                    }
                }
            }

            if( !options.simulate2D )
            {
                if( zmin )
                {
                    for( int j = 0; j < yLength; j++ )
                    {
                        for( int i = 0; i < xLength; i++ )
                        {
                            m.addDirichletNode( m.getVoxelIndex( i, j, 0 ), options.Dirichlet_zmin_values );
                            m.setSubstrateActivation( m.getVoxelIndex( i, j, 0 ), options.Dirichlet_zmin );
                        }
                    }
                }

                if( zmax )
                {
                    for( int j = 0; j < yLength; j++ )
                    {
                        for( int i = 0; i < xLength; i++ )
                        {
                            m.addDirichletNode( m.getVoxelIndex( i, j, zLength - 1 ), options.Dirichlet_zmax_values );
                            m.setSubstrateActivation( m.getVoxelIndex( i, j, zLength - 1 ), options.Dirichlet_zmax );
                        }
                    }
                }
            }
        }

        m.initDirichlet();
        m.displayInformation();

        m.substrateToIndex = new HashMap<>();
        for( int i = 0; i < m.densityNames.length; i++ )
            m.substrateToIndex.put( m.densityNames[i], i );
    }

    void setDirichletActivation(int substrateIndex, boolean value)
    {
        dirichletActivation[substrateIndex] = value;
        for( int n = 0; n < mesh.voxels.length; n++ )
            dirichletActivations[n][substrateIndex] = value;
    }

    void setSubstrateActivation(int index, boolean[] value)
    {
        dirichletActivations[index] = value.clone();
    }

    void resetGradients()
    {
        for( int k = 0; k < mesh.voxels.length; k++ )
        {
            for( int i = 0; i < numberDensities(); i++ )
            {
                gradients[k][i] = new double[3];
            }
        }
        gradientComputed = new boolean[mesh.voxels.length];
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "================================" );
        sb.append( "\nMicroenvironment summary: " + name + ":" );
        sb.append( "\n================================" );
        sb.append( "\n" + mesh.display() + "\n" );
        sb.append( "\nDensities: (" + numberDensities() + " total)" );
        sb.append( "\n--------------------------------" );
        for( int i = 0; i < densityNames.length; i++ )
        {
            sb.append( "\n\t" + i + ". " + densityNames[i] + ":" );
            //            if( !density_units[i].equals( "dimensionless" ) )
            //                sb.append( "\n\tunits: " + density_units[i] );
            sb.append( "\tinitial: " + options.initial_condition_vector[i] );
            if( dirichletActivation[i] )
                sb.append( "\tboundary: " + options.Dirichlet_condition_vector[i] );
            sb.append( "\tdiffusion: " + diffusionCoefficients[i] );
            sb.append( "\tdecay: " + decayRates[i] );
            //            sb.append( "\tdiffusion length scale: "
            //                    + PhysiCellUtilities.print( Math.sqrt( diffusion_coefficients[i] / ( 1e-12 + decay_rates[i] ) ) ) );
        }
        return sb.toString();
    }
}
