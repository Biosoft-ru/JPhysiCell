package ru.biosoft.physicell.biofvm;

import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.Model;

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
public class BasicAgent
{
    //    public static List<BasicAgent> allBasicAgents = new ArrayList<>();

    static int maxBasicAgentId = 0;

    private Model model;
    protected Microenvironment m;
    int selected_microenvironment;

    int current_microenvironment_voxel_index;
    public double volume;
    public double radius; //we consider agents to be balls
    boolean volumeChanged;
    public int currentVoxelIndex;

    public double[] sourceSinkTemp1;
    public double[] sourceSinkTemp2;
    double[] sourceSinkExport1;
    double[] sourceSinkExport2;
    protected double[] prevVelocity;

    double[] substrateChange;
    protected boolean isActive;

    public double[] secretionRates;
    public double[] saturationDensities;
    public double[] uptakeRates;
    public double[] netExportRates;

    public int ID;
    public int index;
    public int type;
    public String tag;

    public double[] position;
    public double[] velocity;

    public double[] internalizedSubstrates;
    public double[] fractionReleasedDeath;
    public double[] fractionTransferredIngested;

    public BasicAgent(Model model)
    {
        this.model = model;
        //give the agent a unique ID  
        ID = maxBasicAgentId;
        maxBasicAgentId++;
        // initialize position and velocity
        isActive = true;

        volume = 1.0;
        position = new double[3];
        velocity = new double[3];
        prevVelocity = new double[3];
        // link into the microenvironment, if one is defined 
        secretionRates = new double[0];//new std::vector<double>(0);
        uptakeRates = new double[0];//= new std::vector<double>(0);
        saturationDensities = new double[0];//= new std::vector<double>(0);
        netExportRates = new double[0];// = new std::vector<double>(0); 

        internalizedSubstrates = new double[0];// = new std::vector<double>(0); // 
        fractionReleasedDeath = new double[0];// = new std::vector<double>(0); 
        fractionTransferredIngested = new double[0];// = new std::vector<double>(0); 
        registerMicroenvironment( model.getMicroenvironment() );

        // these are done in register_microenvironment
        // //internalized_substrates.assign( get_default_microenvironment()->number_of_densities() , 0.0 ); 
    }

    public void setRadius(double radius)
    {
        this.radius = radius;
        this.volume = ( 4.0 / 3.0 ) * Math.PI * Math.pow( radius, 3.0 );
    }

    public void setTag(String tag)
    {
        this.tag = tag;
    }

    public String getTag()
    {
        return tag;
    }

    public double getRadius()
    {
        return radius;
    }

    public void updatePosition(double dt)
    {
        //make sure to update current_voxel_index if you are implementing this function
    };

    public boolean assignPosition(double[] newPosition) throws Exception
    {
        return assignPosition( newPosition[0], newPosition[1], newPosition[2] );
    }

    public boolean assignPosition(double x, double y, double z) throws Exception
    {
        if( !getMicroenvironment().mesh.isPositionValid( x, y, z ) )
            throw new IllegalArgumentException( "Error: the new position for agent " + ID + " is invalid: " + x + "," + y + "," + "z" );

        position[0] = x;
        position[1] = y;
        position[2] = z;
        updateVoxelIndex();

        // make sure the agent is not already registered
        getMicroenvironment().agentContainer.register_agent( this );
        return true;
    }

    protected void updateVoxelIndex()
    {
        if( !getMicroenvironment().mesh.isPositionValid( position[0], position[1], position[2] ) )
        {
            currentVoxelIndex = -1;
            isActive = false;
            return;
        }
        currentVoxelIndex = m.nearestVoxelIndex( position );
    }

    public void registerMicroenvironment(Microenvironment microenvironment)
    {
        this.m = microenvironment;
        microenvironment.addAgent( this );
        double[] density = microenvironment.getDensity( 0 );
        int length = density.length;

        secretionRates = VectorUtil.resize( secretionRates, length );
        saturationDensities = VectorUtil.resize( saturationDensities, length );
        uptakeRates = VectorUtil.resize( uptakeRates, length );
        netExportRates = VectorUtil.resize( netExportRates, length );

        // some solver temporary variables 
        sourceSinkTemp1 = VectorUtil.resize( sourceSinkTemp1, length );
        sourceSinkTemp2 = VectorUtil.resize( sourceSinkTemp2, length );
        sourceSinkExport1 = VectorUtil.resize( sourceSinkExport1, length );
        sourceSinkExport2 = VectorUtil.resize( sourceSinkExport2, length );

        // new for internalized substrate tracking 
        internalizedSubstrates = VectorUtil.resize( internalizedSubstrates, length );
        substrateChange = VectorUtil.resize( substrateChange, length );
        fractionReleasedDeath = VectorUtil.resize( fractionReleasedDeath, length );
        fractionTransferredIngested = VectorUtil.resize( fractionTransferredIngested, length );
        return;
    }

    public Microenvironment getMicroenvironment()
    {
        return m;
    }

    public Model getModel()
    {
        return model;
    }

    public static BasicAgent createBasicAgent(Model model)
    {
        BasicAgent pNew = new BasicAgent( model );
        //        allBasicAgents.add( pNew );
        pNew.index = model.getMicroenvironment().getAgentsCount();
        return pNew;
    }

    public boolean init = false;
    public void setUptakeConstants(double dt)
    {
        init = true;
        /* // new for tracking internal densities
        if( use_internal_densities_as_targets == true )
        {
            *saturation_densities = *internalized_substrates;
            *saturation_densities /= ( 1e-15 + volume ); 
        }*/

        // overall form: dp/dt = S*(T-p) - U*p 
        //   p(n+1) - p(n) = dt*S(n)*T(n) - dt*( S(n) + U(n))*p(n+1)
        //   p(n+1)*temp2 =  p(n) + temp1
        //   p(n+1) = (  p(n) + temp1 )/temp2

        double v = m.voxels( currentVoxelIndex ).volume; //voxel volume
        double discretizeConstant = dt * volume / v; // needs a fix  (?)

        // temp1 = dt*(V_cell/V_voxel)*S*T, where: S - secretion rate, T- saturation density
        sourceSinkTemp1 = VectorUtil.newProd( secretionRates, saturationDensities );
        VectorUtil.prod( sourceSinkTemp1, discretizeConstant );

        // temp2 = 1 + sdt*(V_cell/V_voxel)*( S + U )
        sourceSinkTemp2 = VectorUtil.assign( secretionRates.length, 1 );
        VectorUtil.axpy( sourceSinkTemp2, discretizeConstant, secretionRates );
        VectorUtil.axpy( sourceSinkTemp2, discretizeConstant, uptakeRates );

        // temp for net export 
        sourceSinkExport1 = VectorUtil.newProd( netExportRates, dt );
        sourceSinkExport2 = VectorUtil.newProd( sourceSinkExport1, 1.0 / v );
        volumeChanged = false;
    }

    public void simulateSecretionUptake(Microenvironment m, double dt)
    {
        if( !isActive )
            return;

        if( volumeChanged )
        {
            setUptakeConstants( dt );
            volumeChanged = false;
        }

        //Ilya: Rewritten below
        //        if( m.options.track_internalized_substrates_in_each_agent )
        //        {
        //            substrateChange = VectorUtil.assign( substrateChange.length, 1.0 ); // 1
        //            VectorUtil.diff( substrateChange, sourceSinkTemp2 );// 1-c2
        //            VectorUtil.prod( substrateChange, m.get( currentVoxelIndex ) ); // (1-c2)*rho 
        //            VectorUtil.sum( substrateChange, sourceSinkTemp1 ); // (1-c2)*rho+c1 
        //            VectorUtil.div( substrateChange, sourceSinkTemp2 );// ((1-c2)*rho+c1)/c2
        //            VectorUtil.prod( substrateChange, m.voxels( currentVoxelIndex ).volume );// W*((1-c2)*rho+c1)/c2 
        //            double b = 10;
        //                    VectorUtil.diff( internalizedSubstrates, substrateChange ); // opposite of net extracellular change
        //        }

        //   p(n+1) = (  p(n) + temp1 )/temp2
        double[] density = m.get( currentVoxelIndex );
        double[] densityChange = density.clone();
        VectorUtil.sum( density, sourceSinkTemp1 );
        VectorUtil.div( density, sourceSinkTemp2 );

        // now do net export 
        VectorUtil.sum( density, sourceSinkExport2 );

        if( m.options.track_internalized_substrates_in_each_agent )
        {
            double[] diff = VectorUtil.newDiff( density, densityChange ); //calculate change in density
            VectorUtil.prod( diff, m.voxels( currentVoxelIndex ).volume ); //recalcualte to volume
            VectorUtil.diff( internalizedSubstrates, diff ); //apply source sink change
            VectorUtil.diff( internalizedSubstrates, sourceSinkExport1 ); //apply export change
        }
    }

    public CellContainer get_container()
    {
        return null;//TODO: implement
    }

    public void releaseSubstrates()
    {
        //        Microenvironment pS = Microenvironment.get_default_microenvironment(); 

        // change in total in voxel: 
        // total_ext = total_ext + fraction*total_internal 
        // total_ext / vol_voxel = total_ext / vol_voxel + fraction*total_internal / vol_voxel 
        // density_ext += fraction * total_internal / vol_volume 

        // std::cout << "\t\t\t" << (*pS)(current_voxel_index) << "\t\t\t" << std::endl; 
        //        *internalized_substrates /=  pS->voxels(current_voxel_index).volume; // turn to density 
        //        *internalized_substrates *= *fraction_released_at_death;  // what fraction is released? 

        if( currentVoxelIndex == -1 )
            return;
        VectorUtil.div( internalizedSubstrates, m.voxels( currentVoxelIndex ).volume );// turn to density 
        VectorUtil.prod( internalizedSubstrates, fractionReleasedDeath );// what fraction is released? 
        // release this amount into the environment 

        //        pS.get( currentVoxelIndex ) += internalized_substrates;
        //        (*pS)(current_voxel_index) += *internalized_substrates; 
        VectorUtil.sum( m.get( currentVoxelIndex ), internalizedSubstrates );
        // zero out the now-removed substrates 
        internalizedSubstrates = new double[internalizedSubstrates.length];
        //        internalized_substrates->assign( internalized_substrates->size() , 0.0 ); 

        return;
    }

    public void setTotalVolume(double volume)
    {
        this.volume = volume;
        this.volumeChanged = true;
    }

    public int get_current_mechanics_voxel_index()
    {
        return -1;
    }

    public double[] nearest_gradient(int substrate_index)
    {
        return null;///microenvironment.gradient_vector( currentVoxelIndex )[substrate_index];
    }

}