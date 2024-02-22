package ru.biosoft.physicell.biofvm;

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
public class ConstantCoefficientsLOD3D extends DiffusionDecaySolver
{
    int xLength;
    int yLength;
    int zLength;
    double[][] density;

    public void setup(Microenvironment m, double dt)
    {
        System.out.println( "Using implicit 3-D LOD with Thomas Algorithm ... " );

        xLength = m.mesh.x_coordinates.length;
        yLength = m.mesh.y_coordinates.length;
        zLength = m.mesh.z_coordinates.length;
        density = m.density;

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
    }

    @Override
    public void solve(Microenvironment m, double dt) throws Exception
    {
        if( m.mesh.regularMesh == false || m.mesh.cartesianMesh == false )
            throw new IllegalArgumentException( "Error: This algorithm is written for regular Cartesian meshes. Try: other solvers!" );

        // define constants and pre-computed quantities 
        if( !m.solverSetup )
            setup( m, dt );


        // x-diffusion 
        m.applyDirichletConditions();
        for( int k = 0; k < zLength; k++ ) //        #pragma omp parallel for 
        {
            for( int j = 0; j < yLength; j++ )
            {
                // Thomas solver, x-direction remaining part of forward sweep, using pre-computed quantities 
                int n = m.getVoxelIndex( 0, j, k );
                VectorUtil.div( density[n], m.thomas_denomx[0] );
                for( int i = 1; i < xLength; i++ )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.axpy( density[n], m.thomas_constant1, density[n - m.iJump] );
                    VectorUtil.div( density[n], m.thomas_denomx[i] );
                }
                //back substitution
                for( int i = xLength - 2; i >= 0; i-- )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.naxpy( density[n], m.thomas_cx[i], density[n + m.iJump] );
                }
            }
        }
        // y-diffusion 
        m.applyDirichletConditions();
        for( int k = 0; k < zLength; k++ ) //        #pragma omp parallel for
        {
            for( int i = 0; i < xLength; i++ )
            {
                // Thomas solver, y-direction remaining part of forward sweep, using pre-computed quantities  
                int n = m.getVoxelIndex( i, 0, k );
                VectorUtil.div( density[n], m.thomas_denomy[0] );
                for( int j = 1; j < yLength; j++ )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.axpy( density[n], m.thomas_constant1, density[n - m.jJump] );
                    VectorUtil.div( density[n], m.thomas_denomy[j] );
                }
                // back substitution 
                for( int j = yLength - 2; j >= 0; j-- )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.naxpy( density[n], m.thomas_cy[j], density[n + m.jJump] );
                }
            }
        }
        // z-diffusion 
        m.applyDirichletConditions();
        for( int j = 0; j < yLength; j++ ) //     #pragma omp parallel for 
        {
            for( int i = 0; i < xLength; i++ )
            {
                // Thomas solver, z-direction remaining part of forward sweep, using pre-computed quantities 
                int n = m.getVoxelIndex( i, j, 0 );
                VectorUtil.div( density[n], m.thomas_denomz[0] );
                // should be an empty loop if mesh.z_coordinates.length < 2  
                for( int k = 1; k < zLength; k++ )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.axpy( density[n], m.thomas_constant1, density[n - m.kJump] );
                    VectorUtil.div( density[n], m.thomas_denomz[k] );
                }
                // back substitution
                for( int k = zLength - 2; k >= 0; k-- )
                {
                    n = m.getVoxelIndex( i, j, k );
                    VectorUtil.naxpy( density[n], m.thomas_cz[k], density[n + m.kJump] );
                }
            }
        }
        m.applyDirichletConditions();
    }
}