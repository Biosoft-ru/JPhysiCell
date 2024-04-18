package ru.biosoft.physicell.biofvm;

import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.enums.ProfilerMode;
import uk.ac.manchester.tornado.api.types.arrays.DoubleArray;
import uk.ac.manchester.tornado.api.types.matrix.Matrix2DDouble;
import uk.ac.manchester.tornado.api.types.matrix.Matrix2DInt;

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
public class TornadoSolverParallel extends DiffusionDecaySolver
{
    int xLength;
    int yLength;
    int zLength;

    Matrix2DDouble density;
    Matrix2DDouble denomx;
    Matrix2DDouble cx;
    Matrix2DDouble denomy;
    Matrix2DDouble cy;
    Matrix2DDouble denomz;
    Matrix2DDouble cz;
    DoubleArray constant1;
    int densities;
    Matrix2DInt indices;
    Matrix2DDouble values;

    public void setup(Microenvironment m, double dt)
    {
        System.out.println( "Using implicit 3-D LOD with Thomas Algorithm ... " );

        xLength = m.mesh.x_coordinates.length;
        yLength = m.mesh.y_coordinates.length;
        zLength = m.mesh.z_coordinates.length;
        density = new Matrix2DDouble( m.density );
        densities = m.densityNames.length;

        m.thomas_denomx = new double[xLength][];
        m.thomas_cx = new double[xLength][];
        m.thomas_denomy = new double[yLength][];
        m.thomas_cy = new double[yLength][];
        m.thomas_denomz = new double[zLength][];
        m.thomas_cz = new double[zLength][];

        m.iJump = 1; //offset to get next x coordinate
        m.jJump = xLength; //offset to get next y coordinate
        m.kJump = m.jJump * yLength; //offset to get next z coordinate

        /*We solve tridiagonal matrix equation with each matrix element as a vector:
          
           [ b_1  c_1   0   ...   0  ] [ x_1 ]   [ p_1 ]
           [ a_2  b_2  c_2  ...   0  ] [ x_2 ]   [ p_2 ]
           [  0   a_3  b_3  c_3   0  ] [ x_3 ] = [ p_3 ] 
           [           ...           ] [ ... ]   [ ... ]
           [  0   ...       a_n  b_n ] [ x_n ]   [ p_n ]
        
           where a_i = c_i = dt*D/dx^2
                 b_1 = b_n = 1 + 1/3*dt*lambda + dt*D/dx^2
                 b_i = 1 + 1/3*dt*lambda + 2*dt*D/dx^2 ,  i=2,..,n-1
                 p_i - density on previous step
                 x_i - density on current step
                 
           Using Thomas algorithm
        */
        m.thomas_constant1 = VectorUtil.newProd( m.diffusionCoefficients, dt ); // c1 = dt*D/dx^2  This is -a_i = -c_i, i=0,..,n
        VectorUtil.div( m.thomas_constant1, m.mesh.dx );
        VectorUtil.div( m.thomas_constant1, m.mesh.dx );
        m.thomas_constant1a = VectorUtil.newProd( m.thomas_constant1, -1.0 ); // c1a = -dt*D/dx^2;  This is a_i = c_i, i=0,..,n
        m.thomas_constant2 = VectorUtil.newProd( m.decayRates, dt / 3.0 ); // c2 = (1/3)* dt*lambda 
        m.thomas_constant3 = VectorUtil.newSum( m.one, m.thomas_constant1 ); // c3 = 1 + 2*c1 + c2; //this is b_i, i=0,..,n-1
        VectorUtil.sum( m.thomas_constant3, m.thomas_constant1 );
        VectorUtil.sum( m.thomas_constant3, m.thomas_constant2 );
        m.thomas_constant3a = VectorUtil.newSum( m.one, m.thomas_constant1 ); //c3a = 1 + c1 + c2;  //this is b_0 = b_n
        VectorUtil.sum( m.thomas_constant3a, m.thomas_constant2 );

        // First part of forward sweep, calculate Thomas coeffecients c' and denom
        //x coefficients        
        //First step: cx (c') = c_i , denom = b_i 
        m.thomas_cx = VectorUtil.assign( xLength, m.thomas_constant1a );
        m.thomas_denomx = VectorUtil.assign( xLength, m.thomas_constant3 );
        m.thomas_denomx[0] = m.thomas_constant3a.clone();
        m.thomas_denomx[xLength - 1] = m.thomas_constant3a.clone();

        if( xLength == 1 ) //degenerate case
            m.thomas_denomx[0] = VectorUtil.newSum( m.one, m.thomas_constant2 );

        //Second step:
        // cx_i = c_i  / (b_i - a_i * cx_(i-1))
        // denom_i = (b_i - a_i * cx_(i-1)) - stored for future use
        VectorUtil.div( m.thomas_cx[0], m.thomas_denomx[0] ); // cx_1 = 
        for( int i = 1; i <= xLength - 1; i++ )
        {
            VectorUtil.axpy( m.thomas_denomx[i], m.thomas_constant1, m.thomas_cx[i - 1] );
            VectorUtil.div( m.thomas_cx[i], m.thomas_denomx[i] ); // the value at  size-1 is not actually used  
        }

        //y coefficients the same as for x
        m.thomas_cy = VectorUtil.assign( yLength, m.thomas_constant1a );
        m.thomas_denomy = VectorUtil.assign( yLength, m.thomas_constant3 );
        m.thomas_denomy[0] = m.thomas_constant3a.clone();
        m.thomas_denomy[yLength - 1] = m.thomas_constant3a.clone();
        if( yLength == 1 )
            m.thomas_denomy[0] = VectorUtil.newSum( m.one, m.thomas_constant2 );
        VectorUtil.div( m.thomas_cy[0], m.thomas_denomy[0] );

        for( int i = 1; i <= yLength - 1; i++ )
        {
            VectorUtil.axpy( m.thomas_denomy[i], m.thomas_constant1, m.thomas_cy[i - 1] );
            VectorUtil.div( m.thomas_cy[i], m.thomas_denomy[i] ); // the value at  size-1 is not actually used  
        }

        //z coefficients the same as for x
        m.thomas_cz = VectorUtil.assign( zLength, m.thomas_constant1a );
        m.thomas_denomz = VectorUtil.assign( zLength, m.thomas_constant3 );
        m.thomas_denomz[0] = m.thomas_constant3a.clone();
        m.thomas_denomz[zLength - 1] = m.thomas_constant3a.clone();
        if( zLength == 1 )
            m.thomas_denomz[0] = VectorUtil.newSum( m.one, m.thomas_constant2 );

        VectorUtil.div( m.thomas_cz[0], m.thomas_denomz[0] );

        for( int i = 1; i <= zLength - 1; i++ )
        {
            VectorUtil.axpy( m.thomas_denomz[i], m.thomas_constant1, m.thomas_cz[i - 1] );
            VectorUtil.div( m.thomas_cz[i], m.thomas_denomz[i] ); // the value at  size-1 is not actually used  
        }
        m.solverSetup = true;

        this.constant1 = DoubleArray.fromArray( m.thomas_constant1 );
        this.cx = new Matrix2DDouble( m.thomas_cx );
        this.cy = new Matrix2DDouble( m.thomas_cy );
        this.cz = new Matrix2DDouble( m.thomas_cz );
        this.denomx = new Matrix2DDouble( m.thomas_denomx );
        this.denomy = new Matrix2DDouble( m.thomas_denomy );
        this.denomz = new Matrix2DDouble( m.thomas_denomz );
        this.indices = new Matrix2DInt( m.density.length, densities );
        this.values = new Matrix2DDouble( m.density.length, densities );
        for( int i = 0; i < m.density.length; i++ )
        {
            for( int j = 0; j < densities; j++ )
            {
                if( m.dirichletActivations[i][j] )
                {
                    indices.set( i, j, 1 );
                    values.set( i, j, m.dirichletValue[i][j] );
                }
            }
        }
    }

    static double time = 0;
    public void solve(Microenvironment m, double dt) throws Exception
    {
        if( !m.solverSetup )
            setup( m, dt );

        time++;
        for( int i = 0; i < density.getNumRows(); i++ )
            for( int j = 0; j < density.getNumColumns(); j++ )
                density.set( i, j, m.density[i][j] );

        //        if( time > 1000 )
        //        {
        //            StringBuffer sb = new StringBuffer();
        //            for( int i = 0; i < xLength; i++ )
        //            {
        //                for( int j = 0; j < yLength; j++ )
        //                    sb.append( density.get( voxelIndex( i, j, 0, xLength, yLength ), 0 ) + "\t" );
        //                sb.append( "\n" );
        //            }
        //            System.out.println( sb.toString() );
        //        }
        run( density, denomx, cx, denomy, cy, denomz, cz, constant1, indices, values, xLength, yLength, zLength, densities );

        //        if( time > 1000 )
        //        {
        //            StringBuffer sb2 = new StringBuffer();
        //            for( int i = 0; i < xLength; i++ )
        //            {
        //                for( int j = 0; j < yLength; j++ )
        //                    sb2.append( density.get( voxelIndex( i, j, 0, xLength, yLength ), 0 ) + "\t" );
        //                sb2.append( "\n" );
        //            }
        //            System.out.println( sb2.toString() );
        //        }
        for( int i = 0; i < density.getNumRows(); i++ )
            for( int j = 0; j < density.getNumColumns(); j++ )
                m.density[i][j] = density.get( i, j );
    }

    //    public static void solveX(Matrix2DDouble density, Matrix2DDouble denomx, Matrix2DDouble cx, DoubleArray constant1, int xLength,
    //            int yLength, int zLength, int jump, int densities)
    //    {
    //        for( @Parallel
    //        int k = 0; k < zLength; k++ ) // #pragma omp parallel for
    //        {
    //            for( @Parallel
    //            int j = 0; j < yLength; j++ )
    //            {
    //                int n = voxelIndex( 0, j, k, xLength, xLength );
    //                for( int s = 0; s < densities; s++ )
    //                    density.set( n, s, density.get( n, s ) / denomx.get( 0, s ) );
    //
    //                for( int i = 1; i < xLength; i++ )
    //                {
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int l = 0; l < densities; l++ )
    //                    {
    //                        density.set( n, l, ( density.get( n, l ) + constant1.get( l ) * density.get( n - jump, l ) ) / denomx.get( i, l ) );
    //                    }
    //                }
    //                for( int ri = 0; ri <= xLength - 2; ri++ )
    //                {
    //                    int i = xLength - 2 - ri;
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int s = 0; s < densities; s++ )
    //                    {
    //                        density.set( n, s, density.get( n, s ) - cx.get( i, s ) * density.get( n + jump, s ) );
    //                    }
    //                }
    //            }
    //        }
    //    }
    //
    //    public static void solveY(Matrix2DDouble density, Matrix2DDouble denomy, Matrix2DDouble cy, DoubleArray constant1, int xLength,
    //            int yLength, int zLength, int jump, int densities)
    //    {
    //        for( @Parallel
    //        int k = 0; k < zLength; k++ )
    //        {
    //            for( @Parallel
    //            int i = 0; i < xLength; i++ )
    //            {
    //                int n = voxelIndex( i, 0, k, xLength, xLength );
    //                for( int s = 0; s < densities; s++ )
    //                    density.set( n, s, density.get( n, s ) / denomy.get( 0, s ) );
    //
    //                for( int j = 1; j < xLength; j++ )
    //                {
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int l = 0; l < densities; l++ )
    //                    {
    //                        density.set( n, l, ( density.get( n, l ) + constant1.get( l ) * density.get( n - jump, l ) ) / denomy.get( i, l ) );
    //                    }
    //                }
    //                for( int rj = 0; rj <= xLength - 2; rj++ )
    //                {
    //                    int j = xLength - 2 - rj;
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int s = 0; s < densities; s++ )
    //                    {
    //                        density.set( n, s, density.get( n, s ) - cy.get( i, s ) * density.get( n + jump, s ) );
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    public static void solveZ(Matrix2DDouble density, Matrix2DDouble denomz, Matrix2DDouble cz, DoubleArray constant1, int xLength,
    //            int yLength, int zLength, int jump, int densities)
    //    {
    //        for( @Parallel
    //        int j = 0; j < yLength; j++ )
    //        {
    //            for( @Parallel
    //            int i = 0; i < xLength; i++ )
    //            {
    //                int n = voxelIndex( i, j, 0, xLength, xLength );
    //                for( int s = 0; s < densities; s++ )
    //                    density.set( n, s, density.get( n, s ) / denomz.get( 0, s ) );
    //
    //                for( int k = 1; k < xLength; k++ )
    //                {
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int l = 0; l < densities; l++ )
    //                    {
    //                        density.set( n, l, ( density.get( n, l ) + constant1.get( l ) * density.get( n - jump, l ) ) / denomz.get( i, l ) );
    //                    }
    //                }
    //                for( int rk = 0; rk <= xLength - 2; rk++ )
    //                {
    //                    int k = xLength - 2 - rk;
    //                    n = voxelIndex( i, j, k, xLength, yLength );
    //                    for( int s = 0; s < densities; s++ )
    //                    {
    //                        density.set( n, s, density.get( n, s ) - cz.get( i, s ) * density.get( n + jump, s ) );
    //                    }
    //                }
    //            }
    //        }
    //    }

    public static void solveX(Matrix2DDouble density, Matrix2DDouble denomx, Matrix2DDouble cx, DoubleArray constant1, int xLength,
            int yLength, int zLength, int densities)
    {
        int xJump = 1;
        for( int k = 0; k < zLength; k++ )
        {
            for( int j = 0; j < yLength; j++ )
            {
                int n = voxelIndex( 0, j, k, xLength, xLength );
                for( int s = 0; s < densities; s++ )
                    density.set( n, s, density.get( n, s ) / denomx.get( 0, s ) );

                for( int i = 1; i < xLength; i++ )
                {
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int l = 0; l < densities; l++ )
                    {
                        density.set( n, l,
                                ( density.get( n, l ) + constant1.get( l ) * density.get( n - xJump, l ) ) / denomx.get( i, l ) );
                    }
                }
                for( int ri = 0; ri <= xLength - 2; ri++ )
                {
                    int i = xLength - 2 - ri;
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int s = 0; s < densities; s++ )
                    {
                        density.set( n, s, density.get( n, s ) - cx.get( i, s ) * density.get( n + xJump, s ) );
                    }
                }
            }
        }
        //        for( int i = 0; i < xLength; i++ )
        //        {
        //            for( int j = 0; j < xLength; j++ )
        //            {
        //                density.set( voxelIndex(i,j, 0, xLength, yLength),0, -13);
        //            }
        //        }
    }

    public static void solveY(Matrix2DDouble density, Matrix2DDouble denomy, Matrix2DDouble cy, DoubleArray constant1, int xLength,
            int yLength, int zLength, int densities)
    {
        int yJump = xLength;
        for( @Parallel
        int k = 0; k < zLength; k++ )
        {
            for( @Parallel
            int i = 0; i < xLength; i++ )
            {
                int n = voxelIndex( i, 0, k, xLength, xLength );
                for( int s = 0; s < densities; s++ )
                {
                    density.set( n, s, density.get( n, s ) / denomy.get( 0, s ) );
                }
                for( int j = 1; j < xLength; j++ )
                {
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int l = 0; l < densities; l++ )
                    {
                        double r = ( density.get( n, l ) + constant1.get( l ) * density.get( n - yJump, l ) ) / denomy.get( j, l );
                        density.set( n, l, r );
                    }
                }
                for( int rj = 0; rj <= xLength - 2; rj++ )
                {
                    int j = yLength - 2 - rj;
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int s = 0; s < densities; s++ )
                    {
                        double r = density.get( n, s ) - cy.get( j, s ) * density.get( n + yJump, s );
                        density.set( n, s, r );
                    }
                }
            }
        }
    }

    public static void solveZ(Matrix2DDouble density, Matrix2DDouble denomz, Matrix2DDouble cz, DoubleArray constant1, int xLength,
            int yLength, int zLength, int densities)
    {

        int zJump = xLength * yLength;
        for( @Parallel
        int j = 0; j < yLength; j++ )
        {
            for( @Parallel
            int i = 0; i < xLength; i++ )
            {
                int n = voxelIndex( i, j, 0, xLength, xLength );
                for( int s = 0; s < densities; s++ )
                    density.set( n, s, density.get( n, s ) / denomz.get( 0, s ) );

                for( int k = 1; k < xLength; k++ )
                {
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int l = 0; l < densities; l++ )
                    {
                        density.set( n, l,
                                ( density.get( n, l ) + constant1.get( l ) * density.get( n - zJump, l ) ) / denomz.get( k, l ) );
                    }
                }
                for( int rk = 0; rk <= xLength - 2; rk++ )
                {
                    int k = xLength - 2 - rk;
                    n = voxelIndex( i, j, k, xLength, yLength );
                    for( int s = 0; s < densities; s++ )
                    {
                        density.set( n, s, density.get( n, s ) - cz.get( k, s ) * density.get( n + zJump, s ) );
                    }
                }
            }
        }
    }

    static void applyDirichlet(Matrix2DDouble density, Matrix2DInt indices, Matrix2DDouble values, int densities)
    {
        for( @Parallel
        int i = 0; i < density.size(); i++ )
        {
            for( @Parallel
            int j = 0; j < densities; j++ )
            {
                if( indices.get( i, j ) > 0 )
                    density.set( i, j, values.get( i, j ) );
            }
        }
    }

    public static void sett(Matrix2DDouble d)
    {
        d.set( 0, 0, d.get( 0, 0 ) + 1 );
    }

    public static int voxelIndex(int i, int j, int k, int xLength, int yLength)
    {
        return ( k * yLength + j ) * xLength + i;
    }

    public static void run(Matrix2DDouble density, Matrix2DDouble denomx, Matrix2DDouble cx, Matrix2DDouble denomy, Matrix2DDouble cy,
            Matrix2DDouble denomz, Matrix2DDouble cz, DoubleArray const1, Matrix2DInt indices, Matrix2DDouble values, int xLength,
            int yLength, int zLength, int densities)
    {
        TaskGraph taskGraph = new TaskGraph( "s1" )
                .transferToDevice( DataTransferMode.EVERY_EXECUTION, density, denomx, cx, denomy, cy, denomz, cz, const1, indices, values )
                                .task( "d1", TornadoSolverParallel::applyDirichlet, density, indices, values, densities )
                .task( "t0", TornadoSolverParallel::solveX, density, denomx, cx, const1, xLength, yLength, zLength, densities )
                                .task( "d2", TornadoSolverParallel::applyDirichlet, density, indices, values, densities )
                                .task( "t1", TornadoSolverParallel::solveY, density, denomy, cy, const1, xLength, yLength, zLength, densities )
                                .task( "d3", TornadoSolverParallel::applyDirichlet, density, indices, values, densities )
                                .task( "t2", TornadoSolverParallel::solveZ, density, denomz, cz, const1, xLength, yLength, zLength, densities )
                                .task( "d4", TornadoSolverParallel::applyDirichlet, density, indices, values, densities )
                .task( "d5", TornadoSolverParallel::sett, density ).transferToHost( DataTransferMode.EVERY_EXECUTION, density );
        ImmutableTaskGraph immutableTaskGraph = taskGraph.snapshot();
        TornadoExecutionPlan executionPlan = new TornadoExecutionPlan( immutableTaskGraph );
        executionPlan.resetDevice().withProfiler( ProfilerMode.SILENT ).execute();
        executionPlan.freeDeviceMemory();
    }
}